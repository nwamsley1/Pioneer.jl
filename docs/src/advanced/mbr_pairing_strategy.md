**One-to-One Target–Decoy Pairing With Decoy Cloning**

Goal

- Guarantee that every pair_id identifies exactly one target and one decoy across PSMs, even when multiple targets map to a single decoy. Avoid unintended “target A vs target B” grouping in later stages (e.g., summarize_precursors!, apply_mbr_filter!). Persist pair_id to PSM files so MBR filtering operates on correct 1:1 pairs.

Problems To Solve

- Many-to-one: Multiple targets (A, B, …) may need the same decoy C. If all rows A, B, C share the same pair_id, the “pair” can accidentally include A–B together.
- Downstream grouping: summarize_precursors! and apply_mbr_filter! group by pair_id. If pair_id is not 1:1, we get incorrect within-run comparisons and false transfers.
- Training leakage: Cloned decoy rows should not distort training distributions, but they must exist later for correct pairing and MBR grouping.

Design Overview

- Decoy cloning: For each additional target that is assigned to decoy C beyond the first, create a cloned decoy row (same chromatogram features) and assign a new unique pair_id for the target–clone pair. The original decoy keeps a distinct pair_id with its first target.
- 1:1 invariant: After cloning, each pair_id appears in exactly two roles: one target and one decoy (canonical or clone).
- Exclude clones from training: Cloned decoy rows are not used when fitting models but must be present for inference and written out to per-run PSM files after scoring.
- Stratified pairing: Create target–decoy assignments within 10×10 bins in (prec_mz × iRT) to preserve similarity of pairs. Use global fallback bins when strata are sparse.
- Deterministic: Use a fixed RNG seed for reproducibility.

Data Model Additions

- New columns on the PSM DataFrame (in-memory only; all exported to files except those marked ephemeral):
  - `:pair_id::UInt32` – unique 1:1 identifier for target–decoy pairs (persisted to Arrow files).
  - `:pair_role::UInt8` – 0=target, 1=decoy (canonical), 2=decoy_clone (optional; can derive from `:target` and clone flag).
  - `:pair_clone_of::Union{Missing,UInt32}` – precursor_idx of the canonical decoy for clones; `missing` for targets and canonical decoys.
  - `:pair_training_mask::Bool` – false for clones to exclude them from training sets (ephemeral).

Where To Implement

- File: `src/utils/ML/percolatorSortOf.jl`
- Function entry: `sort_of_percolator_in_memory!`
- Timing: Very start, before `sort!(psms, [:pair_id, :isotopes_captured])`

Algorithm (Pair Generation)

1) Build precursor-level table
- Collapse PSMs to unique `precursor_idx`, keeping first `target::Bool`, `prec_mz::Float32`, and iRT column. iRT precedence: `:irt_pred` > `:irt_obs`; if both missing, use a single iRT bin.

2) Create 10×10 bins (prec_mz × iRT)
- Compute quantile-based edges for `prec_mz` and chosen iRT (0:0.1:1.0). Deduplicate and enforce strictly increasing edges; otherwise fallback to even-width `LinRange(min,max,11)`.
- Assign each precursor to a (bin_mz, bin_irt) stratum via `searchsortedlast` and clamp to [1,10]. If iRT missing, use 10×1 bins.

3) Stratified target–decoy assignment per stratum
- For each stratum s:
  - `targets_s = shuffle(target precursors in s)`
  - `decoys_s = shuffle(decoy precursors in s)`
  - If one side is empty locally, use a globally-shuffled pool from the full dataset for that side.
  - Determine which side is smaller:
    - If `length(decoys_s) < length(targets_s)`: assign each decoy to multiple targets (decoy-reuse). Maintain mapping `assignments_decoy[d] => Vector{targets}`.
    - If `length(targets_s) < length(decoys_s)`: assign each target to multiple decoys (target-reuse). Maintain mapping `assignments_target[t] => Vector{decoys}`.

4) Construct pair_id and clones (no-unpaired invariant)
- We ensure every precursor participates in a 1:1 pair_id unless one side is truly absent (no targets or no decoys even after global fallback). Choose cloning based on which side is larger in the stratum:
  - Decoys fewer than targets (decoy-reuse): For each decoy `d` with assigned targets `[t1, t2, …]`:
    - First target `t1`: `pair_id = next_id()`; assign to `t1` and canonical `d`.
    - Each additional target `ti` (i ≥ 2): clone `d`’s PSM rows per run; set `pair_id = next_id()`, `pair_clone_of = d`, `pair_training_mask = false`.
  - Targets fewer than decoys (target-reuse): For each target `t` with assigned decoys `[d1, d2, …]`:
    - First decoy `d1`: `pair_id = next_id()`; assign to canonical `t` and `d1`.
    - Each additional decoy `dj` (j ≥ 2): clone `t`’s PSM rows per run; set `pair_id = next_id()`, `pair_clone_of = t`, `pair_training_mask = false`.
  - Canonical (non-clone) rows always have `pair_training_mask = true`.

5) Update the full PSM DataFrame
- Overwrite/create `:pair_id` with the newly assigned values for all rows (canonical + clones + targets).
- Add `:pair_clone_of` and `:pair_training_mask` as above.
- Preserve original ordering, or re-sort by `[:pair_id, :isotopes_captured]` as today.

Training, Prediction, and Write-Back

6) Exclude clones when selecting training rows
- In `get_training_data_for_iteration!`, filter by `pair_training_mask` (keep `true` rows only). This avoids inflating decoy counts.

7) Inference on clones
- Keep clones in the test/inference set (fold assignment unchanged). They receive model probabilities directly. Alternatively (optimization), copy predictions from their canonical decoy; initially, compute directly for simplicity and correctness.

8) Persist to Arrow files
- `write_scored_psms_to_files!` writes out the PSMs including the cloned decoy rows and the regenerated `pair_id`. Keep `:pair_id` and drop only vector columns as before.

Downstream Changes

9) apply_mbr_filter! requires `:pair_id`
- Remove the fallback that re-derives `pair_id` from the library.
- If `:pair_id` is missing in merged_df, throw an error with a clear message: the file must be produced by the new pipeline that regenerates pair_id.

10) summarize_precursors! grouping remains the same
- With 1:1 `pair_id`, within-run target–decoy logic is correct. No further changes required aside from previously added robustness to use `:mbr_prob` or `:prob`.

apply_mbr_filter! adjustments (with 1:1 pair_id)

- Preconditions:
  - merged_df contains `:pair_id::UInt32`, `:prob::Float32`, `:target::Bool`, `:decoy::Bool`, `:MBR_transfer_candidate::Bool`, and `:ms_file_idx`.
  - Each `:pair_id` denotes exactly one target and one decoy per run (thanks to cloning during pairing).
- Remove library-derived `pair_id` code:
  - Delete any attempt to recompute `:pair_id` from the spectral library (e.g., via `getPairId`).
  - Insert a hard check: if `:pair_id` ∉ propertynames(merged_df) → `error("pair_id missing; regenerate pairs before ScoringSearch")`.
- Candidate set and non-candidates:
  - `candidate_mask = merged_df.MBR_transfer_candidate`.
  - Compute trace q-values only on non-candidates with library FDR scaling for diagnostics: `get_qvalues!(merged_df.prob[.!candidate_mask], merged_df.target[.!candidate_mask], trace_qval[.!candidate_mask])`.
- Within-run target–decoy dominance:
  - Group candidates by `[:ms_file_idx, :pair_id]`.
  - For each group, identify best target and best decoy (there should be at most one of each). If either is missing, skip the group with a warning counter.
  - If `best_decoy_prob ≥ best_target_prob`, mark the target row in this group as a transfer decoy (bad). Always mark decoy rows as transfer decoys, since they are the reference negatives.
- Build bad-mask and threshold:
  - `bad_mask = candidate_mask .& (merged_df.decoy .| target_marked_as_bad)`.
  - Compute τ with `get_ftr_threshold(merged_df.prob, merged_df.target, bad_mask, α; mask=candidate_mask)`, where α = `params.max_MBR_false_transfer_rate` (or α' if alpha-scaling experiment is enabled).
- Clamp candidate probabilities:
  - `merged_df._filtered_prob = ifelse.(candidate_mask .& (merged_df.prob .< τ), 0f0, merged_df.prob)` and return `:_filtered_prob`.
- Logging and diagnostics:
  - Count candidate groups lacking a target or decoy; log as potential data issues.
  - Log `(α, τ, #candidates, #transfer_decoys)` and a small sample of group-level decisions for QA.
- Interaction with stratified FTR experiments:
  - If later stratifying FTR calibration, reuse `:pair_id` structure but compute τ within strata (e.g., by `MBR_num_runs` or donor-count), applying the above dominance rule inside each stratum.

Validation & Diagnostics

- Before/after stats:
  - Unique `pair_id` count; distribution of group sizes (expect mostly 2, some 1).
  - Fraction of decoys cloned (and average clone count).
  - Per-stratum pairing coverage; number of strata using global fallback.
- Correctness checks:
  - For each `pair_id`, ensure at most one target and at most one decoy row per run.
  - On a small dataset, assert the 1:1 property holds globally.
- Logging:
  - Seed, bins used, stratum sizes, counts of clones, any fallbacks triggered.

Performance Considerations

- Cloning increases row count by up to the number of extra target assignments per decoy. In typical DIA datasets, this should be a moderate multiplier. Monitor memory and optionally gate on dataset size.
- If needed later, add a config toggle to disable cloning (fallback to original behavior) for memory-constrained runs.

Rollout Steps & Commits

1) Commit checkpoint (baseline):
- Commit current state before changes as `chore(MBR): checkpoint before 1:1 pairing work`.

2) Implement pairing + cloning + training mask + persistence
- Add `regenerate_pair_ids!` helper and call at top of `sort_of_percolator_in_memory!`.
- Modify `get_training_data_for_iteration!` to respect `pair_training_mask`.
- Ensure `write_scored_psms_to_files!` persists `:pair_id` (and optionally `:pair_clone_of` for debugging).

3) Update apply_mbr_filter!
- Remove library-derived `pair_id` fallback.
- Error if `:pair_id` missing in merged_df.

4) Commit implementation:
- `feat(MBR): 1:1 target–decoy pairing with decoy cloning and persisted pair_id`

5) Verification run
- Run a representative dataset; check logs and pairing stats; verify apply_mbr_filter! sees `pair_id` and completes.
