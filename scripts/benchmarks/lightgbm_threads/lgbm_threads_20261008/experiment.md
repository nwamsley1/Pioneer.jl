# Independent Julia and LightGBM threading experiment

RIS job **3369734**, submitted October 8, 2026. Clean develop snapshot `cfb45759dccbbe209c9b83f46f4005d96117763b`, separate from PR #535 and deferred S01/S02 branches. No production changes.

Current source inspection found no universal hard cap of 24 search threads: CLI launchers honor the selected `JULIA_NUM_THREADS`, GUI ceiling is the machine's available core count, SearchDIA sizes its context from `Threads.nthreads()`, and `build_lightgbm_classifier` defaults native `num_threads` to that count. LightGBM uses an independent OpenMP pool. Its Julia wrapper accepts an explicit `num_threads` for prediction, defaulting to the estimator's count.

The experiment holds Julia at 24 default worker threads and GC at one marking/one sweep thread, while native LightGBM uses 16, 24, 32, 48 or 64 threads in sequence. One exclusive 64-CPU SLURM allocation, 32 GB hard memory limit, pinned Julia 1.12.6 container. BLAS stays at one thread. A private first depot reads existing dependency caches with `--compiled-modules=existing`; there is no package setup or shared precompile-cache mutation. Submission excludes all current user-job nodes plus `c2-node-005` and `c2-node-056` observed during preparation.

Inputs are saved `SWATH_r01.arrow` and `SWATH_r02.arrow` from the earlier local memory audit, copied into private RIS scratch input storage. Existing precursor CV assignments are validated across files. Training uses fold 1 and check/prediction rows use fold 0. Features come from current Pioneer main/advanced sets, intersected with available saved columns, excluding non-ion-mobility features via the normal selector and removing constant enzymatic termini as production does. Missing features and complete source hashes are recorded.

Training scenarios: 250,000 rows with current 50-tree small MainSearch parameters; 1,000,000 and 2,500,000 rows with current 200-tree experiment-wide scoring parameters. Larger matrices repeat source training rows within their original fold. They are throughput stress cases, not independent large cohorts or full CV/identification validation. Primary runs retain production `deterministic=true`, `seed=1776` and forced row-wise histograms. Separate exploratory column-wise runs compare 24 and 64 native threads at 2.5 million rows.

Record dataset creation/binning+labels, tree training and production-style dataset detachment separately. Native operations execute sequentially on the startup Julia thread. A reference 24-thread model checks all retrained models' predictions on up to 100,000 held-out rows; differences are recorded, not hidden by timing. Fixed-model prediction tests score 50,000/500,000/1,000,000 rows using both prebuilt matrices and Pioneer's complete fill/predict/copy path. Every output is compared to its 24-thread reference. The 1-million-row test is a batch-size scaling probe; ordinary prediction batches currently default to 500,000.

Three measured trials per primary setting, rotated/reversed thread-count order after reference training and prediction warmups, all on the same node. Feature preparation and explicit GC are outside timed kernels. CPU time is whole-process CPU time across native threads. Allocated bytes are Julia cumulative allocations and do not include all native memory; peak RSS is a process high-water across settings and cannot establish variant-specific peak savings. SLURM step MaxRSS and actual CPU topology are retained separately.

Windows compatibility is a separate condition: the repository documents corrupted results when multithreaded native calls were entered from Julia worker threads on Windows. This experiment tests serial native entry from the startup thread on RIS/Linux. It does not establish Windows correctness or justify increasing native threads in concurrent per-file Julia tasks.

Relevant primary documentation: [LightGBM parameters](https://lightgbm.readthedocs.io/en/stable/Parameters.html#num_threads), [determinism](https://lightgbm.readthedocs.io/en/stable/Parameters.html#deterministic), [column-wise histograms](https://lightgbm.readthedocs.io/en/stable/Parameters.html#force_col_wise). The official documentation recommends real CPU cores for thread counts and suggests column-wise histograms at high thread counts. The experiment compares that mode separately rather than attributing its effect to threads alone.

Remote experiment: `/scratch2/fs1/d.goldfarb/n.t.wamsley/runs/lgbm_threads_20261008T174117Z`.

Remote source: `/scratch2/fs1/d.goldfarb/n.t.wamsley/src/Pioneer_lgbm_threads_20261008T174117Z`.

Both kernel job 3369734 and follow-up complete Pass-1 job 3369949 completed successfully. Results, parity checks, limitations, and implementation implications are recorded in [lgbm_thread_results.md](lgbm_thread_results.md). The follow-up used a process-local prediction override because Pass-1 otherwise explicitly resets native prediction threads to Julia's count. It ran on the same host with Julia at 24, compared 24/48/64 native threads, and verified all stored sidecar values across nine measured trials.

Fetched artifacts are in `results/` and `pass1_results/`; kernel summaries are in `summary.csv` and `summary.json`; stage medians are in `pass1_summary.json`. Only stats were fetched. No production source changes were made.
