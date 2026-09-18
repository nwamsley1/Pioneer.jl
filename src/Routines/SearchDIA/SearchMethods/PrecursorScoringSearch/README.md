# Precursor scoring

`score_precursor_isotope_traces` uses the same Pass-1 training procedure with
match-between-runs (MBR) enabled or disabled. Implementation: `pass1_oom.jl`.

## Training pool

- Sample uniformly within the original CV folds, up to
  `SCORING_LGBM_MAX_TRAIN` rows per fold (currently 2.5 million), without class
  balancing. Select final reservoir row locations before loading features.
- Keep this pool fixed across semi-supervised iterations. The first iteration
  trains on all pool rows; later iterations train on all pool decoys and pool
  targets with q ≤ 0.03. Excluded targets remain available for evaluation and
  can qualify again.
- Score every pool row with the model trained on the opposite fold. Select
  models by pool targets at q ≤ 0.01; stop when the gain is below 1% or after
  eight iterations. These counts describe the pool, not the full experiment.

## Final scoring and memory

After model selection, release the pool and score every source file once.
Sidecars preserve row order and contain `precursor_idx`, `scan_idx`,
`trace_prob_prepass` (OOF), and `trace_prob_infold` (NaN when MBR is disabled).
The scores are attached while merging each run's folds into a single Arrow
file. Experiment-wide FDR and integrated MBR processing run downstream.

Retained LightGBM models contain trees without training datasets. Training
memory includes the pool, a temporary filtered subset, and native LightGBM
workspace; final prediction holds one file's feature matrices at a time.
The pool cap controls training size; `max_psms_in_memory` is currently unused.

Debug logs report sampling, per-fold fitting, pool prediction, q-value
calculation, and final prediction timings, plus file progress every 100 files.
Post-scoring logs time fold merging, cleanup, probability aggregation, score
sorting/merging, and q-value/PEP calculation.
MBR preparation also logs annotation, initial filtering, donor-threshold
calculation, donor indexing, and staging writes. Its eligibility index stores
at most two distinct donor run IDs per precursor, enough to determine whether
a donor exists outside any receiver run.
