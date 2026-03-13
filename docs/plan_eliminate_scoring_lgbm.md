# Plan: Eliminate ScoringSearch LightGBM, Use FirstPass Scores with PEP-Based Experiment-Wide FDR

## Motivation

The ScoringSearch LightGBM (20 features) adds ~6% more unique precursors over FirstPass LightGBM
(15 features) but at the cost of an entire additional ML training + scoring + global aggregation
pipeline. This plan eliminates ScoringSearch's LightGBM entirely. Instead:

1. **FirstPassSearch** performs global FDR filtering and writes global_prob / global_qval columns
   directly into the fold-split output tables.
2. **ScoringSearch** skips LightGBM training and global q-value computation entirely. It starts
   from **per-file PEP calculation → experiment-wide q-values** using the already-written
   `prec_prob` column (which comes from the FirstPass LightGBM trace aggregation).

### Why PEP for experiment-wide q-values

Scores from different files are not directly comparable (different target/decoy distributions,
different data quality). PEP values are calibrated probabilities ("probability this PSM is
incorrect"), so they **are** comparable across files. The experiment-wide q-value step must:

1. Compute PEP per-file from the `prec_prob` column (per-file calibration)
2. Pool all PSMs across files with their PEP values
3. Sort by PEP, compute experiment-wide q-values from the pooled PEPs

This is what `build_qvalue_spline_from_refs` already does when passed `compute_pep=true`.

---

## Current Pipeline (what exists today)

```
FirstPassSearch.summarize_results!()
├── aggregate_prescore_globally!()          # PEP-calibrate per file, log-odds combine top-√n, global q-values
├── Filter first_pass_psms → passing_precs  # Remove precursors failing global FDR (default 5%)
├── Add Phase 2 columns (irt_diff, prec_mz, pair_id, entrapment_group_id)
├── initialize_prob_group_features!()        # Add trace_prob=0, q_value=0, MBR columns
└── Write fold-split files to second_pass_psms/{file}_fold{0,1}.arrow

ScoringSearch.summarize_results!()
├── Step 1: Train LightGBM (20 features) on fold-split files           ← REMOVE
├── Step 1b: Merge fold files back into single files per MS run        ← KEEP (but move fold merge to FirstPassSearch)
├── Step 2: MBR filtering + trace→precursor aggregation                ← KEEP
├── Step 3: Best isotope trace selection                               ← KEEP
├── Step 4: Quantification processing (best trace filter)              ← KEEP
├── Steps 5-10: Build global_prob dicts, global q-values, per-file     ← REPLACE (global already done;
│   q-values via spline, filter by both thresholds                        only need experiment-wide q-values from PEP)
├── Step 11: Recalculate q-values on filtered data                     ← KEEP (but on experiment-wide q-values)
├── Steps 12-24: Protein inference + scoring + protein FDR             ← KEEP (unchanged)
```

## Proposed Pipeline

```
FirstPassSearch.summarize_results!()
├── aggregate_prescore_globally!()           # Same: PEP-calibrate, log-odds, global q-values
├── Filter first_pass_psms → passing_precs   # CHANGE: use 1% global FDR threshold (was 5%)
├── Write global_prob and global_qval columns into fold-split tables   ← NEW
├── Add Phase 2 columns
├── initialize_prob_group_features!()
├── Write NON-fold-split files to second_pass_psms/{file}.arrow        ← CHANGE (no more fold split)
└── Filter fragment index (unchanged)

ScoringSearch.summarize_results!()
├── Step 1: SKIP LightGBM training entirely                           ← CHANGED
├── Step 1b: SKIP fold merging (files are already merged)             ← CHANGED
├── Step 2: MBR filtering + trace→precursor aggregation               ← KEEP (but simpler: no LightGBM trace_prob)
├── Step 3: Best isotope trace selection                              ← KEEP
├── Step 4: Quantification processing (best trace filter)             ← KEEP
├── Steps 5-10: SIMPLIFIED                                            ← CHANGED
│   ├── Skip global_prob dict building (already in table as column)
│   ├── Skip global q-value computation (already in table as column)
│   ├── Build experiment-wide q-values from per-file PEP              ← KEY CHANGE
│   │   └── Use build_qvalue_spline_from_refs(refs, :prec_prob, ...; compute_pep=true)
│   │       This computes PEP per-file, pools, and builds a q-value spline
│   └── Filter by experiment-wide q-value threshold only
│       (global q-value filter already applied in FirstPassSearch)
├── Step 11: Recalculate experiment-wide q-values on filtered data    ← KEEP
├── Steps 12-24: Protein inference + scoring + protein FDR            ← KEEP (unchanged)
```

---

## Detailed Changes

### File 1: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl`

#### Change 1a: `summarize_results!()` — Write global_prob and global_qval into PSM tables

**Location**: Lines 292-444

**Current behavior**: `aggregate_prescore_globally!()` returns `passing_precs` (Set{UInt32}) and
`prec_best_scan` (Dict). The fold-split files contain no global_prob or global_qval columns.

**New behavior**:
1. `aggregate_prescore_globally!()` should also return the `global_prob_dict` (Dict{UInt32, Float32})
   and `global_qval_dict` (Dict{UInt32, Float32}) so we can write them as columns.
2. When writing the output tables, add two new columns:
   - `global_prob::Float32` — looked up from `global_prob_dict[precursor_idx]`
   - `global_qval::Float32` — looked up from `global_qval_dict[precursor_idx]`
3. Change global FDR threshold from 5% to 1% (since this is now the final global filter).

**Pseudocode for the column addition (inside the per-file loop, after Phase 2 columns)**:
```julia
# Add global prescore columns
tbl[!, :global_prob] = Float32[global_prob_dict[pid] for pid in tbl[!, :precursor_idx]]
tbl[!, :global_qval] = Float32[global_qval_dict[pid] for pid in tbl[!, :precursor_idx]]
```

#### Change 1b: Write single files instead of fold-split

**Current behavior**: Writes `{file}_fold0.arrow` and `{file}_fold1.arrow`.

**New behavior**: Write a single `{file}.arrow` file. ScoringSearch no longer needs fold-split
files since it's not training a LightGBM. The `cv_fold` column should still be present (protein
inference uses it).

```julia
# Replace the fold-split loop with:
file_path = "$(base_path).arrow"
writeArrow(file_path, tbl)
```

**Impact**: ScoringSearch's `get_valid_fold_file_paths()` will no longer find fold files, so we
need to update ScoringSearch to handle single files.

### File 2: `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl`

#### Change 2a: `aggregate_prescore_globally!()` — Return additional dicts

**Location**: Lines 1788-2039

**Current signature returns**: `(passing_precs::Set{UInt32}, prec_best_scan::Dict)`

**New signature returns**: `(passing_precs, prec_best_scan, global_prob_dict, global_qval_dict)`

The function already computes `global_probs_dict` and `global_qvals` internally. We just need to:
1. Build `global_prob_dict::Dict{UInt32, Float32}` mapping precursor_idx → aggregated probability
2. Build `global_qval_dict::Dict{UInt32, Float32}` mapping precursor_idx → global q-value
3. Return them alongside the existing return values

These dicts are already computed internally (the function builds vectors `global_pids`,
`global_probs`, `global_qvals`) — we just need to package them into dicts before returning.

### File 3: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`

#### Change 3a: Remove Step 1 (LightGBM training)

**Location**: Lines 296-335

Replace the entire `score_precursor_isotope_traces(...)` call with a simple read of the
already-merged files. Since FirstPassSearch now writes single files (not fold-split), we just
load the file references directly.

```julia
step1_time = @elapsed begin
    # No LightGBM training — FirstPass scores are used directly
    # Files are already single-file (not fold-split)
    @user_info "Skipping LightGBM training — using FirstPass scores directly"
end
```

#### Change 3b: Remove Step 1b (fold merging)

**Location**: Lines 337-378

Remove entirely — files are already single `.arrow` files from FirstPassSearch.

The existing code that builds `second_pass_refs` from merged paths needs to be replaced with
direct construction from the already-single files:

```julia
second_pass_paths = [getSecondPassPsms(getMSData(search_context), idx) * ".arrow"
                     for (idx, _) in valid_file_data]
second_pass_refs = [PSMFileReference(path) for path in second_pass_paths]
```

Wait — `setSecondPassPsms!` stores the base path (without `.arrow`). Need to check how
`get_valid_file_paths` resolves this. It may need the `.arrow` suffix.

**TODO**: Verify how `getSecondPassPsms` returns paths and ensure consistency with the new
single-file approach. Currently FirstPassSearch does `setSecondPassPsms!(ms_data, ms_file_idx, base_path)`
where `base_path = joinpath(second_pass_dir, file_name)` (no extension). ScoringSearch then
constructs fold paths as `"$(base_path)_fold0.arrow"`. With single files, we should either:
- Store the full path with `.arrow` extension, OR
- Add `.arrow` when constructing refs in ScoringSearch

#### Change 3c: Simplify Steps 5-10 (remove global q-value computation, add PEP-based experiment-wide q-values)

**Location**: Lines 414-468

**Remove**:
- `build_precursor_global_prob_dicts()` — global_prob is already a column in the table
- `build_global_qval_dict_from_scores()` — global_qval is already a column in the table
- The global_qval filtering threshold (already applied in FirstPassSearch)

**Replace with**:
1. Read `global_prob` and `global_qval` from the existing table columns (for downstream use
   by protein inference, etc.)
2. Build experiment-wide q-values from PEP:
   ```julia
   # Build experiment-wide q-value spline from per-file prec_prob via PEP
   score_col = has_mbr ? :MBR_boosted_prec_prob : :prec_prob
   spline_result = build_qvalue_spline_from_refs(filtered_refs, score_col, results.merged_quant_path;
       compute_pep=true, min_pep_points_per_bin=params.precursor_q_value_interpolation_points_per_bin,
       fdr_scale_factor=fdr_scale, temp_prefix="qval_sidecar")
   qval_spline = spline_result.qval_spline
   results.precursor_qval_interp[] = qval_spline
   results.precursor_pep_interp[] = spline_result.pep_interp
   ```
3. Filter by experiment-wide q-value only (global q-value filter already applied):
   ```julia
   combined_pipeline = TransformPipeline() |>
       add_interpolated_column(:qval, score_col, qval_spline) |>
       add_interpolated_column(:pep, score_col, results.precursor_pep_interp[]) |>
       filter_by_multiple_thresholds([
           (:qval, params.q_value_threshold)
       ])
   ```

**Note on `global_prob_dict` for downstream**: Steps 12-24 (protein inference) currently
use `results.precursor_global_qval_dict[]`. We need to either:
- Build this dict by streaming the `global_qval` column from the tables (lightweight), OR
- Store it in SearchContext from FirstPassSearch results

The simplest approach: stream-build the dict once in ScoringSearch from the existing column:
```julia
global_qval_dict = Dict{UInt32, Float32}()
for ref in filtered_refs
    tbl = Arrow.Table(file_path(ref))
    for i in 1:length(tbl.precursor_idx)
        global_qval_dict[tbl.precursor_idx[i]] = tbl.global_qval[i]
    end
end
results.precursor_global_qval_dict[] = global_qval_dict
```

#### Change 3d: Update `trace_prob` / `prec_prob` handling

**Current behavior**: ScoringSearch's LightGBM writes `trace_prob` into the fold files. After
merging, `prec_prob` is computed from `trace_prob` via Bayesian aggregation in Step 2.

**New behavior**: Since we skip LightGBM, the `trace_prob` column is initialized to 0 by
FirstPassSearch. We need to set `trace_prob = lgbm_prob` (the FirstPass LightGBM probability)
so that the Bayesian aggregation in Step 2 works correctly.

**Option A**: In FirstPassSearch, set `trace_prob = lgbm_prob` instead of 0.
**Option B**: In ScoringSearch Step 2, use `lgbm_prob` instead of `trace_prob` for aggregation.

**Recommendation**: Option A is cleaner — set `trace_prob = lgbm_prob` in FirstPassSearch when
writing the output tables. This way ScoringSearch's Step 2 (`_aggregate_trace_to_precursor_probs!`)
works without modification.

Actually, wait — the current pipeline has one PSM per precursor per file at this stage (best
scan selected in `process_search_results!`). So "trace" aggregation is trivial (1 trace per
precursor). The `prec_prob` would just equal `trace_prob`. So we can:

**Simplest approach**: In FirstPassSearch, set `trace_prob = lgbm_prob` so that:
- Step 2's `_aggregate_trace_to_precursor_probs!` produces `prec_prob = lgbm_prob`
- Everything downstream works unchanged

### File 4: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

No changes needed — the function simply won't be called.

### File 5: Parameter changes

#### `global_prescore_qvalue_threshold`

**Current**: defaults to 0.05 (5% FDR)
**New**: Keep as tunable parameter, default stays at 0.05. Users can set to 0.01 in their
config to test the tighter global filter. The experiment-wide q-value step (at
`q_value_threshold`, default 0.01) provides the final filter regardless.

**Location**: `parseParams.jl` where `global_prescore_qvalue_threshold` is parsed. No change needed.

---

## Risk Assessment

### Low Risk
- Writing `global_prob` / `global_qval` columns in FirstPassSearch — additive change
- Removing LightGBM training call — just deleting code
- Removing fold-split — simplification

### Medium Risk
- **`trace_prob` initialization**: If we set `trace_prob = lgbm_prob`, the MBR logic
  (which uses `trace_prob` and `MBR_boosted_trace_prob`) might behave differently since
  it expects trace-level rather than precursor-level probabilities. Need to verify that
  MBR's `MBR_transfer_candidate` detection works correctly with precursor-level scores.

  **Mitigation**: At this stage, there's 1 row per precursor (best scan already selected),
  so trace-level = precursor-level. The MBR detection uses `trace_prob` to identify
  missing values (trace_prob == 0 means missing in a run). With `trace_prob = lgbm_prob`,
  only true misses (precursor absent from file) have trace_prob = 0. This should be correct.

- **Experiment-wide q-value quality**: PEP-based q-values may be more or less conservative
  than the current LightGBM-rescored q-values. Need to validate with the OlsenAstral dataset.

### High Risk
- **Global FDR threshold change (5% → 1%)**: This is a significant tightening. Currently
  ScoringSearch re-filters at 1% after its own LightGBM. With the tighter FirstPass filter,
  we lose the chance for ScoringSearch to "rescue" borderline precursors with better features.

  **Mitigation**: The two_stage_fdr analysis showed FirstPass-only at 1% gives 188K unique
  targets vs full pipeline's 199K — a 6% loss. But the full pipeline uses 20 additional
  features. By using PEP-based experiment-wide q-values (which is what we're adding), we
  may recover some of that gap.

  **Alternative**: Keep global threshold at 5% and let the experiment-wide q-value filter
  (at 1%) do the final filtering. This preserves more candidates for the experiment-wide
  step. **This is the safer option and should be the default.**

---

## Recommended Implementation Order

1. **Modify `aggregate_prescore_globally!`** to return `global_prob_dict` and `global_qval_dict`
2. **Modify `FirstPassSearch.summarize_results!`**:
   a. Receive the new dicts
   b. Set `trace_prob = lgbm_prob` instead of 0
   c. Add `global_prob` and `global_qval` columns
   d. Write single files instead of fold-split
3. **Modify `ScoringSearch.summarize_results!`**:
   a. Remove Step 1 (LightGBM training)
   b. Remove Step 1b (fold merging) — replace with direct file ref construction
   c. Replace Steps 5-10 with PEP-based experiment-wide q-values
   d. Build `global_qval_dict` from table column for downstream protein steps
4. **Update parameter defaults** (`global_prescore_qvalue_threshold`)
5. **Test with ecoli integration test** (`SearchDIA("./data/ecoli_test/ecoli_test_params.json")`)
6. **Validate with OlsenAstral dataset** — compare precursor counts to current pipeline

## Open Questions

1. **~~Keep 5% or move to 1% global threshold?~~** RESOLVED: Keep as tunable parameter,
   default 5%. Users can set to 1% in config. The experiment-wide q-value step at
   `q_value_threshold` (default 1%) provides the final filter.

2. **MBR handling**: With no ScoringSearch LightGBM, MBR match quality assessment relies on
   FirstPass scores only. Is this sufficient? The MBR features (pair_prob, irt_diff,
   weight_ratio, etc.) are computed from the `trace_prob` column. With `trace_prob = lgbm_prob`,
   MBR should still work, but the transfer quality estimates may be less accurate without the
   additional 5 features that ScoringSearch's LightGBM used.

3. **Protein inference scoring**: Steps 12-24 use `prec_prob` (from ScoringSearch's LightGBM)
   for protein group scoring. With FirstPass's `lgbm_prob` flowing through as `prec_prob`,
   protein group scores will be based on 15-feature scores instead of 20-feature scores.
   This is acceptable — the protein probit regression (Step 15) adds its own features.

4. **`prec_prob` vs `lgbm_prob` naming**: After this change, `prec_prob` is derived from
   `lgbm_prob` (FirstPass). Should we rename for clarity, or keep the existing column names
   for downstream compatibility? **Recommendation**: Keep existing names for compatibility.
