# Plan: Robust MS1→MS2 Join in SecondPassSearch

## Context
- Location: `process_search_results!` in SecondPassSearch.
- Current behavior: After `filter!(x -> x.best_scan, psms)`, we keep MS2 PSM rows intended to represent the apex per `precursor_idx`. There can still be multiple MS2 rows per `precursor_idx` if different `:isotopes_captured` values are retained. We want to attach a single MS1 PSM row per precursor — the one whose RT is closest to the MS2 apex RT — and replicate that same MS1 info to all MS2 rows of that precursor (including across distinct `:isotopes_captured`).
- Current implementation sketch: Select closest MS1 row by scanning `psms` inside a `combine(groupby(ms1_psms, :precursor_idx))` body; then perform a left join on `:precursor_idx` only.

## Goals
- Deterministic selection of the “best” MS1 PSM per key (precursor, and run), defined as minimal |RT_ms1 − RT_ms2_apex|.
- Correct replication: When multiple MS2 rows exist for the same precursor (e.g., different `:isotopes_captured`), the same chosen MS1 row is attached to each MS2 row.
- Robust handling across runs: Do not cross-contaminate matches between different files/runs. Use `:ms_file_idx` alongside `:precursor_idx` where applicable.
- Performance: Avoid O(N^2) scans of `psms` inside group loops; precompute and vectorize where possible.
- Data hygiene: Optionally handle non-finite MS1 metrics (`Inf`/`NaN`) before downstream scoring/ML.

## Proposed Design
1. Define join key
   - Use `join_keys = [:ms_file_idx, :precursor_idx]` if available; otherwise fall back to `[:precursor_idx]` in strictly single-run contexts.
   - Rationale: Ensures MS1/ MS2 are matched within the same run.

2. Precompute canonical MS2 apex RT per key
   - After `best_scan` filtering, derive a single MS2 reference RT per `join_keys`:
     - If `psms` has one row per key: `ms2_apex_rt = psms.rt` directly.
     - If multiple rows per key remain (e.g., differing `:isotopes_captured`): choose a deterministic reference, e.g., the row with max `:prob` (or simply the first after sorting by `:prob` desc, then `:rt`).
   - Implementation idea:
     - `ms2_apex = combine(groupby(psms, join_keys), [:prob, :rt] => ((p, r) -> r[argmax(p)]) => :ms2_apex_rt)`
     - If `:prob` not suitable, choose `first(r)` after sorting `group_df` by `:prob` desc.

3. Join MS2 apex RT into MS1 candidates
   - Ensure `ms1_psms` has `:rt`; if missing, compute via `getRetentionTime(spectra, scan_idx)` (already present in the code).
   - Left-join: `ms1_psms_with_rt = leftjoin(ms1_psms, ms2_apex, on=join_keys)` to add `:ms2_apex_rt`.
   - Compute absolute RT difference: `:rt_diff = abs.(ms1_psms_with_rt.rt .- ms1_psms_with_rt.ms2_apex_rt)`; be careful with missing `:ms2_apex_rt` (skip such groups or keep as missing for later prune).

4. Select one MS1 row per key (closest RT)
   - Group `ms1_psms_with_rt` by `join_keys` and choose the row with the minimum `:rt_diff` per group. Deterministic tie-breaker: if multiple rows share the same minimum `:rt_diff`, pick the one with higher `:prob` (if present) or the first by `:scan_idx`.
   - Implementation idea:
     - Precompute per-group boolean mask of `is_min` using `argmin` on `:rt_diff`.
     - Or use `combine(..., :rt_diff => argmin => :idx)` then `groupby`-indexed lookup.
   - Result: `ms1_psms_best` with a single row per `join_keys`.

5. Final join to replicate MS1 across MS2 rows
   - Left-join original `psms` with `ms1_psms_best` on `join_keys` only (do NOT include `:isotopes_captured` so the same MS1 row is replicated across all MS2 isotopic rows of the precursor):
     - `psms = leftjoin(psms, ms1_psms_best, on=join_keys, makeunique=true, renamecols = "" => "_ms1")`
   - This satisfies the requirement: “give the MS1 data to both rows” when multiple MS2 rows exist for a precursor.

6. Optional: sanitize non-finite MS1 columns
   - Replace `Inf`/`-Inf`/`NaN` within the appended `_ms1` columns to `missing` (or clamp to numeric bounds) before ML/scoring stages that assume finiteness.
   - Implementation idea: find columns ending with `_ms1` and apply `replace!(col, x -> isfinite(x) ? x : missing)` for Float columns.

## Complexity and Performance
- The plan avoids nested lookups of `psms` from inside a group loop. Instead, it computes a compact `ms2_apex` table once, joins it into MS1 candidates, and uses groupwise `argmin`. This scales linearly with the number of MS1 and MS2 records and leverages DataFrames’ optimized grouping/joins.

## Edge Cases
- No MS1 candidates for a key: The left join produces `missing` MS1 columns; acceptable.
- No MS2 apex record for a key: `ms2_apex_rt` will be missing; either drop such groups from `ms1_psms_best` or allow the final join to propagate missing MS1 values.
- Multiple runs: Ensure `join_keys` includes `:ms_file_idx` (or whichever run identifier is canonical) to avoid cross-run matches.
- Multiple MS2 rows per key: Intentional — the same MS1 row is replicated to each MS2 row after the final join.

## Testing Plan
1. Unit tests (DataFrames-level)
   - Single-run, single-precursor: multiple MS1 candidates; verify nearest-by-RT selection.
   - Single-run, multiple MS2 rows (different `:isotopes_captured`): confirm the same MS1 row attaches to all.
   - Multi-run: include two runs with overlapping `:precursor_idx`; verify join respects `:ms_file_idx`.
   - Ties on RT difference: ensure deterministic tie-breaker (e.g., highest `:prob` or smallest `:scan_idx`).
   - Missing `:rt` in MS1: verify it is populated via `getRetentionTime` path.
   - Non-finite values: introduce `Inf` in an MS1 metric and verify sanitization logic (if enabled) results in `missing`.

2. Integration tests
   - Craft small synthetic MS2/ MS1 sets flowing through the relevant part of SecondPassSearch and assert the joined columns and row counts.

## Pseudocode (for orientation only)
```
# After best_scan filtering in psms
join_keys = hasproperty(psms, :ms_file_idx) ? [:ms_file_idx, :precursor_idx] : [:precursor_idx]

# 1) Compute ms2 apex RT per key (deterministic)
ms2_apex = combine(groupby(psms, join_keys)) do g
    g_sorted = sort(g, [:prob, :rt], rev=[true, false])
    (; ms2_apex_rt = g_sorted.rt[1])
end

# 2) Ensure ms1_psms has :rt
if !in(:rt, names(ms1_psms))
    ms1_psms.rt = [getRetentionTime(spectra, i) for i in ms1_psms.scan_idx]
end

# 3) Join ms2 apex rt into ms1 candidates and compute |Δrt|
ms1p = leftjoin(ms1_psms, ms2_apex, on=join_keys)
ms1p.rt_diff = abs.(ms1p.rt .- ms1p.ms2_apex_rt)

# 4) Pick closest per key
ms1p_sorted = sort(ms1p, [:rt_diff, :prob, :scan_idx], rev=[false, true, false])
ms1_best = combine(groupby(ms1p_sorted, join_keys)) do g
    g[1:1, :]
end

# 5) Final join: replicate MS1 across MS2 rows per key
psms = leftjoin(psms, ms1_best, on=join_keys, makeunique=true, renamecols => ("" => "_ms1"))

# 6) Optional: sanitize appended _ms1 float columns
```

## Rollout
- Implement as a focused refactor inside `process_search_results!`:
  - Introduce a small helper (e.g., `select_closest_ms1_per_precursor(ms1_psms, psms, spectra)`) for clarity and unit-testability.
  - Keep the join keys configurable based on the presence of run identifiers.
  - Add tests under `test/Routines/SearchDIA/...` mirroring this path.

## Notes
- This plan does not change the conceptual behavior (still “closest RT” MS1 attached), but makes it deterministic, scalable, and aligned with multi-run data. It also explicitly supports the case where multiple MS2 rows per precursor should receive the same MS1 data.

