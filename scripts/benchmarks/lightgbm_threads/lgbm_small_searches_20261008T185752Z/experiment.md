# Small-search LightGBM thread regression experiment

RIS job **3370284**, submitted October 8, 2026. Compare native LightGBM training and prediction at **24 versus 48** threads with Julia fixed at **24**. Same frozen develop source as the preceding benchmark: `cfb45759dccbbe209c9b83f46f4005d96117763b`. No production implementation or PR changes.

Datasets:

- `SCP_Astral_250pg_3ms`: all three original 250 pg / 3 ms injection-time replicates; human May 7, 2026 standard bitvec-10 library.
- `Olsen_Exploris_3P`: three E10H50Y40 500 ng / 30-SPD replicates; corresponding three-proteome May 7, 2026 standard bitvec-10 library.

Fresh 24-thread SearchDIA searches produce current-develop feature tables from these raw inputs. Only the new private scratch directories receive outputs. Configuration templates come from the local regression-config repository, with paths updated to read-only storage3 input/library sources and private scratch outputs. Temp files are retained for timing. MBR, QC plots and TSV output are disabled while generating source PSMs; these seed runs are not 48-thread whole-search comparisons.

Each dataset is then tested as a **one-file** and **three-file** cohort. Each of the four cases gets a 24-thread reference/warmup, a 48-thread warmup, and **five alternating-order paired repeats**: 24/48, 48/24, 24/48, 48/24, 24/48. Both CV folds, production sampling cap and semi-supervised stopping logic are retained. No source observations are replicated, no smaller artificial training cap is imposed, and full out-of-fold and in-fold Float32 scores are written.

Every trial uses fresh private PSM copies. Copying and verification happen outside stage timing. Timed stage work includes metadata/pool preparation, training, garbage collection, prediction, and sidecar writing. Process-local instrumentation measures each fit and prediction; the only prediction behavior change is using the classifier's native thread count instead of explicitly resetting it to Julia's count. All native calls assert they run on startup Julia thread 1. Signatures compare every typed sidecar value and column order/length; selected iterations and target/decoy counts are also compared. Parity differences are recorded rather than suppressed.

Exclusive 64-CPU node, 64 GB allocation, Julia 1.12.6 container, one marking/one sweep GC thread, BLAS=1. Initially requested `c2-node-008` to match the previous benchmark; our own pending allocation was updated to idle node `c2-node-009` to avoid waiting for that node. Actual topology is retained in `host_cpu.txt`. A private first depot reads installed caches using `--compiled-modules=existing`; no setup or precompile changes. At submission, running user jobs occupied `c2-node-077`, which was excluded. No existing allocation was interrupted or modified.

Record elapsed time, process CPU time, cumulative Julia allocations, GC time, fit/prediction call sizes, model feature lists and parity. Peak RSS is a whole-process high-water across configurations and cannot isolate memory overhead at 48. The experiment evaluates Linux Pass-1 on actual small-cohort PSMs; it does not validate Windows or time every other LightGBM stage at 48.

Remote experiment: `/scratch2/fs1/d.goldfarb/n.t.wamsley/runs/lgbm_small_searches_20261008T185752Z`.

Frozen source: `/scratch2/fs1/d.goldfarb/n.t.wamsley/src/Pioneer_lgbm_threads_20261008T174117Z`.

Files: `benchmark.jl`, `benchmark.sbatch`, `configs/*.json`, and fetched `results/` (only remote `stats/`, excluding the depot and bulk outputs).

**Completed:** job 3370284 exited successfully in 23m41s. All 40 measured trials matched the reference metrics and complete stored sidecar values. Independent local summarization verified those signatures and the actual fit/prediction thread counts. Results and limitations are in [small_search_results.md](small_search_results.md), with per-case medians/ranges in `summary.json`. At 48, median Pass-1 time increased by 0.26–0.94 s across the four cases. No production files were changed.
