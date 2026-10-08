# Independent Julia and LightGBM threads: measured results

October 8, 2026. Pioneer develop `cfb45759dccbbe209c9b83f46f4005d96117763b`.

Keeping Julia at 24 threads while allowing native LightGBM to use more cores works on the tested Linux host. Prediction benefits substantially on large batches. Training improves modestly on the largest matrix and can become slower on smaller matrices. Giving both training and prediction 48 threads did not improve the complete two-file Pass-1 stage; 64 threads made it slower. These results support separate training and prediction budgets, rather than automatically using every available core for both.

## What the current code does

There is no universal 24-thread cap in this develop snapshot. Both CLI launchers accept the requested Julia thread count, the GUI limits its selector to the machine's CPU count, and SearchDIA constructs its search context using `Threads.nthreads()`. The classifier builder defaults native LightGBM threads to that same count. Specific concurrent per-file native operations have their own restrictions; they are not evidence of a search-wide 24-thread cap.

LightGBM's native OpenMP pool can use a different count from Julia's worker pool. Increasing its count does not require creating additional Julia search workers or their search buffers. It can still increase native resource use, and this experiment does not measure the RSS difference between independent settings.

Pass-1 prediction currently explicitly supplies `Threads.nthreads()` even if the fitted classifier uses a different budget. For the complete-stage experiment, a process-local override changed only that prediction argument to `cls.num_threads`. The normal conversion to Float32 and probability clamping were preserved. Training received its count through the existing hyperparameter argument. No repository source or dependency cache was modified.

Source references in the clean experiment checkout:

- `src/build/CLI/pioneer:47`, `src/build/CLI/pioneer.bat:4`: Julia thread selection.
- `gui/src/App.tsx:418`: hardware-based selector ceiling.
- `src/Routines/SearchDIA.jl:221`: search context thread count.
- `src/utils/ML/PSMScoring/lightgbm_utils.jl:239`: default classifier thread count.
- `src/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/pass1_oom.jl:352`: explicit prediction override.

## Experiment conditions

Two separate, exclusive SLURM jobs ran on RIS `c2-node-008`: kernel benchmark **3369734** and complete Pass-1 benchmark **3369949**. Both completed with exit code zero. The host has 64 physical cores across two Intel Xeon Gold 6548Y+ sockets, one hardware thread per core, and four NUMA nodes. Each job reserved 64 CPUs and 32 GB. Julia 1.12.6 remained at 24 threads; GC used one marking and one sweep thread; BLAS used one thread. LightGBM.jl was 2.2.2. Native model calls ran sequentially from the startup Julia thread.

The source and outputs were isolated from existing searches. A private first depot read installed caches using `--compiled-modules=existing`; no package setup or shared precompile writes occurred. Existing user jobs were not stopped or modified. The clean source checkout remained unchanged.

Inputs were two saved real PSM files, `SWATH_r01.arrow` and `SWATH_r02.arrow`, containing 270,169 and 318,888 rows. Their combined 589,057 rows contained 304,630 fold-1 rows and 284,427 fold-0 rows. Each precursor's fold assignment was checked for consistency across files. Models used current production feature selectors and hyperparameters, including deterministic mode and seed 1776. These saved inputs lack some features introduced in newer develop; the measured available main and advanced feature sets contained **49 and 72 columns**, respectively. See `results/inputs.json` for exact included/missing columns and input hashes. The complete Pass-1 stage used production's feature-resolution behavior on these same files.

Each reported setting has three measured repeats following warmup, with settings rotated/reversed between repeats. The kernel benchmark compared 16, 24, 32, 48 and 64 native threads. The complete stage compared 24, 48 and 64. Times below are medians, not single best runs. Only lightweight result processing ran locally.

## Training and prediction kernels

Training times include dataset construction/binning and labels, model fitting, and production-style dataset detachment. They exclude feature-matrix preparation. Prediction times include Pioneer's normal matrix fill, native prediction and output copy, using a fixed 24-thread-trained model. Normal prediction batch size is 500,000; separate 50,000- and 1,000,000-row probes are retained in the raw results.

| Work, seconds | 16 threads | 24 threads | 32 threads | 48 threads | 64 threads |
|---|---:|---:|---:|---:|---:|
| Main model, train 250k rows / 50 trees | 0.446 | 0.404 | 0.377 | 0.411 | 0.476 |
| Scoring model, train 1M rows / 200 trees | 4.386 | 3.887 | 3.587 | 3.686 | 4.416 |
| Scoring model, train 2.5M rows / 200 trees | 9.008 | 7.445 | 6.834 | 6.377 | 6.618 |
| Main model, predict 500k rows | 0.0749 | 0.0545 | 0.0420 | 0.0370 | 0.0400 |
| 1M-trained scoring model, predict 500k rows | 0.3546 | 0.2411 | 0.1878 | 0.1275 | 0.1459 |
| 2.5M-trained scoring model, predict 500k rows | 0.3119 | 0.2215 | 0.1693 | 0.1185 | 0.1313 |

At 48 versus 24 threads, the scoring prediction batches were about **1.87–1.89 times faster**, reducing elapsed time by roughly 47%. Training the 2.5M-row matrix took **14.3% less elapsed time**, while the 1M-row matrix gained only 5.2%. The smaller Main model gained nothing at 48. Sixty-four threads were generally worse than 48 on this host.

The 250k training matrix contains distinct source rows. The 1M and 2.5M training matrices repeat observations from the same original training fold; prediction matrices larger than the held-out fold repeat held-out observations. Those are throughput stress cases, not independent million-row cohorts. They preserve the training/held-out separation but cannot establish identification behavior or performance on a heterogeneous large search.

Higher training counts consumed more total process CPU time: the 2.5M-row training case rose from **163 CPU-seconds at 24 threads to 270 at 48 and 372 at 64**. The large prediction cases had much smaller CPU-time increases. For example, the 1M-trained model's 500k-row prediction used 5.58, 5.71 and 6.25 CPU-seconds at 24, 48 and 64. This makes prediction a better candidate when the goal includes limiting total computation. Process CPU measurements include native runtime overhead and waiting.

Small prediction batches showed less reliable scaling: for the 1M-trained model at 50k rows, medians were 0.0252 s at 24, 0.0322 s at 48 and 0.0294 s at 64. An unconditional high thread count can hurt small calls.

A separate column-wise histogram probe on the 2.5M-row matrix took 8.110 s at 24 and 8.195 s at 64. It was slower than the production row-wise mode here. Histogram mode was varied separately from the primary thread comparison. This result does not support changing production's histogram mode for these inputs.

## Actual complete Pass-1 stage

The follow-up called `train_and_predict_pass1_oom!` on fresh private copies of the two PSM files, using both original CV folds, normal semi-supervised iteration logic, and both out-of-fold and in-fold sidecar scores. Training rows were not repeated. Input copying and parity checks were outside the timed stage. Timed work included pool loading, training, iteration selection, garbage collection, prediction and sidecar writing. Training and prediction both used the specified native count.

| Native threads | Median elapsed | Range over 3 trials | Median CPU-seconds | Median Julia GC |
|---|---:|---:|---:|---:|
| 24 | 17.166 s | 16.832–17.199 s | 274.99 | 3.235 s |
| 48 | 17.233 s | 17.158–17.293 s | 565.47 | 3.238 s |
| 64 | 19.287 s | 18.383–20.284 s | 837.90 | 3.591 s |

Forty-eight threads produced **no measured full-stage improvement**; the median difference was about 0.4%, within the observed spread. Sixty-four increased median time by 12.3%. Each fold has about 300k rows, so this is a real small-cohort stage test, not the 2.5M-row stress case. It demonstrates why isolated prediction speedups cannot be presented as complete-search speedups.

The larger training pool is bounded by the existing sampling cap, while the total number of rows to predict grows with a search. It is plausible that prediction's benefit will matter more for a large cohort. That is an inference, not a measured whole-search result. Neither a large complete cohort nor a train-at-24/predict-at-48 hybrid was benchmarked in this experiment.

## Result parity and memory

- All **270 fixed-model prediction timing trials** matched their reference Float64 predictions bit for bit.
- Retraining at different thread counts sometimes changed held-out Float64 probabilities by extremely small amounts. Across the primary row-wise training cases, the largest absolute difference was **7.94e-15**. Thus native retraining is not universally bitwise identical in this experiment. The separately changed histogram mode differed by up to 2.15e-14 from the row-wise reference.
- All **nine complete Pass-1 trials** matched the 24-thread warmup reference in selected iteration and counts. The selected iteration was 3; target/decoy counts at q=0.01 were 60,551/605.
- Typed-value SHA256 comparisons verified every sidecar column in both files: precursor IDs, scan IDs, out-of-fold Float32 probabilities and in-fold Float32 probabilities. Column order, types, lengths and all values matched for all 589,057 rows. This checks logical output values, not Arrow serialization bytes or downstream whole-search IDs.
- Julia cumulative allocations were approximately 1.047 GB per complete-stage trial at all three settings. Cumulative allocation is not peak live memory. SLURM reported step MaxRSS of 3,117,548 KiB for the kernel job and 2,258,932 KiB for the stage job. Each is a high-water over the whole job, so neither can establish variant-specific native memory overhead or savings.

## Implications for an implementation

The mechanism is feasible: keep Julia/search workers at the chosen search budget (24 in these experiments), and pass an independent native budget to LightGBM. The experiments do not establish 24 as an optimal search-wide cap, since Julia's search-thread scaling was deliberately held fixed.

A useful implementation would distinguish search workers, model training threads and model prediction threads. It would carry explicit budgets through classifier construction and prediction, replacing Pass-1's hard-coded prediction budget. Limits should respect the user's requested resources and the available CPU allocation. Increasing native threads should apply to serial model operations, not multiply native pools inside already-parallel file processing.

For this host and these inputs, 48 is a promising large-batch prediction budget; training is better around 24–32 for the smaller models and around 48 for the largest stress matrix. These are measurements on one host, not universal defaults. The next validation before choosing a production policy would compare a train-at-24/predict-at-48 hybrid with the baseline on a complete large cohort, with the current full feature set and final output parity checks.

Windows requires its own validation. `src/build/windows/LIGHTGBM_MSVC.md:214` documents silent identification loss when multithreaded native calls were entered from Julia worker threads. Our Linux experiments deliberately used serial native entry from the startup thread. That documented failing path should not be expanded merely because these Linux timings improve.

## Reproduction and artifacts

### Required launch policy clarified after the measurements

The user requires a requested count above 24 to increase only the eligible serial LightGBM budget. Search workers and each family of preallocated per-worker search structures must remain at most 24. The proposed LightGBM ceiling is 48. For requested counts 16, 24, 32, 48 and 64, the corresponding search/LightGBM limits would be 16/16, 24/24, 24/32, 24/48 and 24/48, further constrained by the CPU resources actually available to the process.

Enforce this before launching Julia: preserve the original requested resource count separately, launch with `JULIA_NUM_THREADS=min(requested,24)`, and provide an independent LightGBM budget of `min(requested,48)`. Do this in both CLI launchers and the GUI backend, which can launch `bin/SearchDIA` directly and therefore bypass the CLI wrapper. Resolve `auto` within the available CPU allocation before applying either limit. Derive other non-LightGBM worker settings from the capped search budget. Native LightGBM threads must never become an argument to search-context allocation.

SearchDIA currently passes Julia's thread count to `initSearchContext`, which passes it to `initSimpleSearchContexts` to create one search structure per worker. Starting Julia with 24 therefore makes the existing allocation path create 24 structures. Merely passing 24 to this constructor while leaving Julia at 48 is insufficient: ordinary `Threads.@threads` loops could still run on Julia worker IDs above 24. Direct Julia API calls and direct executable launches also need a defined entry-point policy; a runtime cannot silently reduce an already-created Julia worker pool. A guard that explains the supported launch configuration is one possible policy.

Acceptance checks for an implementation must cover both CLI syntaxes, the Windows batch launcher, GUI direct-executable and wrapper paths, inherited environment values, and `auto`. At a requested count of 48, verify the launched Julia count is 24, the search context and every per-worker buffer collection contain 24 elements, active search workers stay within that budget, and eligible native training and prediction are limited to 48. Requests below 24 must retain their smaller count. A native budget above 24 must not increase the number of search buffers. Existing single-thread native calls inside concurrent per-file operations must retain their safe execution policy, especially on Windows.

This policy is a documented implementation requirement, not a production change in the benchmark checkout. Forty-eight is a ceiling rather than an experimentally established optimum for every model or batch.

A subsequent [small-search benchmark](../lgbm_small_searches_20261008T185752Z/small_search_results.md) tested the exact 250 pg / 3 ms SCP dataset and an Exploris cohort, each as one-file and three-file cases. Five paired repeats at 24/48 showed an additional 0.26–0.94 s in complete precursor Pass-1 at 48; all 40 measured trials matched stored scores and model-selection metrics. Total process CPU time roughly doubled. This follow-up also held Julia at 24 and made no production changes.

LightGBM.jl's [published constructor documentation](https://iqvia-ml.github.io/LightGBM.jl/dev/functions/) exposes `num_threads`. Its [current prediction source and docstring](https://github.com/IQVIA-ML/LightGBM.jl/blob/master/src/predict.jl) explicitly expose a separate prediction `num_threads`, defaulting to the estimator's count. Microsoft's [thread parameter documentation](https://lightgbm.readthedocs.io/en/stable/Parameters.html#num_threads) recommends physical cores and cautions against excessive threads for small inputs; it does not prescribe a cap of 48. The [tuning guide](https://lightgbm.readthedocs.io/en/stable/Parameters-Tuning.html#add-more-computational-resources) explains that these operations use OpenMP. The numerical ceiling here comes from the Pioneer measurements and the requested resource policy, not an upstream recommendation.

All files below are under `/Users/nathanwamsley/Projects/pioneer-large-search-analysis/lgbm_threads_20261008/`:

- `benchmark.jl`, `benchmark.sbatch`: kernel experiment and allocation settings.
- `pass1_experiment.jl`: full-stage driver and process-local prediction override.
- `results/`: kernel environment, inputs, CPU topology, per-trial JSONL, log and completion marker.
- `pass1_results/`: stage environment, per-trial JSONL with complete sidecar signatures, log and completion marker.
- `summary.csv`, `summary.json`: grouped kernel medians, ranges, CPU times and parity statistics; generated by `summarize.py`.
- `pass1_summary.json`: grouped complete-stage medians and ranges, excluding the warmup.

Frozen source: `/scratch2/fs1/d.goldfarb/n.t.wamsley/src/Pioneer_lgbm_threads_20261008T174117Z`.

Remote kernel job: `/scratch2/fs1/d.goldfarb/n.t.wamsley/runs/lgbm_threads_20261008T174117Z`.

Remote stage job: `/scratch2/fs1/d.goldfarb/n.t.wamsley/runs/lgbm_pass1_threads_20261008T174117Z`.

Primary background documentation: [LightGBM parameters](https://lightgbm.readthedocs.io/en/stable/Parameters.html#num_threads), [deterministic mode](https://lightgbm.readthedocs.io/en/stable/Parameters.html#deterministic), and [histogram mode](https://lightgbm.readthedocs.io/en/stable/Parameters.html#force_col_wise). The official documentation describes independent native threading and recommends physical-core counts; our measured results determine the conclusions here.
