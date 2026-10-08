# LightGBM thread scaling experiments

These experiments measure an independent native LightGBM thread budget while Julia stays at 24 workers. They contain **benchmark drivers and measured results**, not a production threading-policy implementation. Pioneer launchers, GUI, classifier defaults and search-buffer allocation are unchanged by this branch.

Source snapshot: develop `cfb45759dccbbe209c9b83f46f4005d96117763b`. Tests ran on RIS/Linux on October 8, 2026 using Julia 1.12.6 and LightGBM.jl 2.2.2. Each allocation was exclusive and separate from existing user jobs; native calls ran sequentially on startup Julia thread 1. Both hosts had 64 physical cores across two Intel Xeon Gold 6548Y+ sockets.

## Results

- [Kernel scaling and complete SWATH Pass-1 results](lgbm_threads_20261008/lgbm_thread_results.md): native thread counts 16/24/32/48/64. Large scoring prediction batches were about 1.9 times faster at 48 than 24. Training on the largest, 2.5-million-row stress matrix was 14% faster, but smaller fits improved less or slowed down. Larger stress matrices repeated source observations within the same CV fold; they are not independent large cohorts. All 270 fixed-model prediction timing trials matched Float64 outputs exactly; retraining sometimes differed at roughly 1e-15. Nine complete two-file Pass-1 trials matched stored Float32 scores and model-selection metrics.
- [Small-search results](lgbm_small_searches_20261008T185752Z/small_search_results.md): actual SCP 250 pg / 3 ms and Exploris PSMs freshly generated from current develop, tested as one-file and three-file cohorts. Five paired repeats per setting; all 40 measured trials matched stored sidecar values and model-selection metrics. Forty-eight native threads added 0.26–0.94 seconds to complete precursor Pass-1. Process CPU time roughly doubled despite the small elapsed-time difference.

| Complete precursor Pass-1 | 24 native threads | 48 native threads |
|---|---:|---:|
| SCP 250 pg / 3 ms, one file | 9.60 s | 9.86 s |
| SCP 250 pg / 3 ms, three files | 22.44 s | 23.38 s |
| Exploris 500 ng, one file | 15.55 s | 15.93 s |
| Exploris 500 ng, three files | 22.97 s | 23.59 s |

These are medians of warmed stage trials, not 48-thread whole-search timings. The small-search seed runs disabled MBR, QC plots and TSV writing while retaining PSM files. Neither these experiments nor their output parity checks validate Windows. Peak RSS is a high-water over each entire job and cannot isolate memory overhead by thread setting.

## Recheck the recorded measurements

The following lightweight Python commands regenerate summaries from the committed JSONL. They do not load Julia or run a search:

```sh
python3 scripts/benchmarks/lightgbm_threads/lgbm_threads_20261008/summarize.py
python3 scripts/benchmarks/lightgbm_threads/lgbm_small_searches_20261008T185752Z/summarize.py
```

The small-search summarizer independently checks complete sidecar signatures, model-selection metrics, actual native call thread counts and the number of measured trials. Recorded Julia driver hashes can be checked against the corresponding environment JSON. Logs have `.txt` extensions in this bundle, with trailing whitespace and line endings normalized. Measured JSONL records and the Julia drivers are unchanged.

## Rerun on RIS

The original layout and RIS paths are retained for provenance. Use a **new private scratch run directory**, then stage the Julia driver and sbatch file from the chosen experiment. Update the scratch paths in the small-search configurations and the `INPUT_DIR` in `pass1_experiment.jl` for a new run. `benchmark.sbatch` accepts the run directory and a frozen Pioneer source directory as its two arguments. The experiment directories' `experiment.md` files describe input files, libraries, resource settings and timing boundaries.

For the kernel driver, provide the two SWATH PSM files under `<run>/inputs/`. The standalone Pass-1 driver must be staged as `<run>/benchmark.jl` to match the sbatch entry point. The small-search driver expects `<run>/configs/` and three input Arrow links for each dataset; its configurations select read-only regression inputs/libraries and write all outputs under the private run directory. Raw MS, PSM and spectral-library data are intentionally excluded from this repository bundle; input paths and hashes are recorded in the results.

Julia must stay at 24 for these comparisons even though SLURM reserves 64 CPUs. The scripts use a private first depot with `--compiled-modules=existing`, BLAS=1 and one marking/one sweep GC thread. They do not perform package setup, precompile shared caches or modify other jobs. Submit to a node separate from existing searches. Review the local RIS operating instructions before submission.

## Proposed policy, pending implementation

Preserve the user's requested count independently, start search Julia with at most 24 workers, and limit eligible serial native LightGBM operations to at most 48 threads, respecting smaller user requests and the actual CPU allocation. Both CLI launchers and the GUI backend need to apply this before Julia starts. Search buffers must follow the capped Julia count, never the native count. Direct executable and Julia API entry points also need a defined policy.

Forty-eight is a proposed ceiling, not a universal optimum. Existing restrictions on native calls from Julia workers, especially the documented Windows failing path, must be retained. See the detailed report for the implementation requirements and validation still needed.
