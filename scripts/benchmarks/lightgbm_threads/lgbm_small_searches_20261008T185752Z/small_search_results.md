# 24 versus 48 LightGBM threads on smaller searches

October 8, 2026. RIS job **3370284**. Pioneer develop `cfb45759dccbbe209c9b83f46f4005d96117763b`. Julia stayed at **24** workers throughout; native LightGBM training and prediction used either **24 or 48**. No production changes.

## Complete precursor Pass-1 timings

Five measured repeats per configuration, paired in alternating order, following warmups of both configurations. Times are medians. Positive change means 48 was slower.

| Actual workload | PSM rows | 24 threads, s | 48 threads, s | Difference of medians, s | Median paired change |
|---|---:|---:|---:|---:|---:|
| SCP 250 pg / 3 ms, 1 file | 55,038 | 9.600 | 9.859 | +0.259 | +3.1% |
| SCP 250 pg / 3 ms, 3 files | 168,223 | 22.441 | 23.380 | +0.939 | +3.9% |
| Exploris 500 ng, 1 file | 209,427 | 15.545 | 15.931 | +0.386 | +2.3% |
| Exploris 500 ng, 3 files | 625,152 | 22.969 | 23.590 | +0.621 | +0.9% |

The percentage column is the median of the five individual paired ratios; it is not calculated by dividing the two median times. Per-trial ranges and paired changes are preserved in `summary.json`.

## Where the difference comes from

Fit timings include preparation of labels, native dataset construction/training and production dataset detachment. Prediction timings include native prediction, Float32 conversion and clamping for every pool-validation and final-file prediction call. Matrix gathering, q-value calculations, garbage collection and sidecar writing also contribute to the complete stage.

| Workload | Sum of fit calls, s (24 → 48) | Sum of predict calls, s (24 → 48) | Process CPU-seconds (24 → 48) | Julia GC, s (24 → 48) |
|---|---:|---:|---:|---:|
| SCP 250 pg / 3 ms, 1 file | 2.134 → 2.457 | 0.129 → 0.074 | 59.9 → 124.9 | 7.192 → 7.188 |
| SCP 250 pg / 3 ms, 3 files | 8.132 → 9.222 | 0.708 → 0.405 | 212.3 → 445.9 | 13.145 → 13.195 |
| Exploris 500 ng, 1 file | 7.846 → 8.490 | 0.632 → 0.351 | 195.8 → 402.9 | 6.768 → 6.689 |
| Exploris 500 ng, 3 files | 13.234 → 14.939 | 1.884 → 1.025 | 331.1 → 689.1 | 6.707 → 6.714 |

These columns are separately aggregated medians and need not sum to the median stage duration. The much larger CPU-time cost is material if scalability includes total computation. GC time was substantial after generating the fresh searches; the complete-stage figures include that cost. The per-call fit/prediction figures expose native behavior independently of its contribution to overall time.

## Correctness checks

All 40 measured stage trials matched their 24-thread reference in every recorded output and model-selection metric: **True**. Checks covered selected semi-supervised iteration, pool target/decoy counts at q=0.01, and typed-value SHA256 signatures of every sidecar column for every row: precursor IDs, scan IDs, out-of-fold Float32 scores and in-fold Float32 scores. Column order, types and lengths were included. The trial data and actual native thread counts were also independently checked during local summarization.

This validates stored precursor-stage results on these inputs. Pool target counts are not final pipeline identification counts. It does not establish universal bitwise equality of native retraining, final pipeline parity at 48, or Windows correctness.

## Inputs and method

The exact low-input dataset was **SCP_Astral_250pg_3ms**, with three original 250 pg / 3 ms injection-time replicates. The conventional-instrument dataset was **Olsen_Exploris_3P**, using three E10H50Y40 500 ng / 30-SPD replicates. Each was first searched from raw Arrow MS inputs with current develop, then tested as both a one-file subset and a three-file cohort. No observations were repeated and the training sampling cap was not reduced. Actual available advanced feature sets contained 72 columns; full lists are in the JSONL.

The seed searches used available May 7, 2026 bitvec-10 human and three-proteome standard libraries, respectively. They retained temp files and disabled MBR, QC plots and TSV writing while producing the source PSMs. They were not timed 24-versus-48 whole-search comparisons. Their elapsed times (292.9 s SCP and 143.0 s Exploris) include initial compilation and other pipeline work and should not be treated as warmed baseline search times.

Baseline final identification counts confirm the low-ID workload: the three SCP runs identified 3,767 / 3,652 / 3,763 precursors and 1,375 / 1,390 / 1,391 protein groups. The three Exploris runs identified 53,639 / 53,254 / 54,122 precursors and 7,171 / 7,109 / 7,179 protein groups. These are from the seed search reports, not 48-thread final-output validation.

Fresh private copies of those PSMs were used for every timing trial. Copying and output verification were outside timing. Actual production `train_and_predict_pass1_oom!` used both original precursor CV folds, its normal maximum training pool, production semi-supervised stopping logic and both out-of-fold/in-fold sidecars. A process-local override changed Pass-1 prediction from explicitly using Julia's count to using the classifier's native count. A second override instrumented the original fit behavior. Neither changed source files or caches. All native calls asserted startup Julia thread 1.

Both settings ran on the same exclusive RIS host, `c2-node-009`: Intel Xeon Gold 6548Y+, 64 physical cores, two sockets, one hardware thread per core, four NUMA nodes. Allocation: 64 CPUs, 64 GB. Runtime: Julia 1.12.6 container, Julia 24, GC one marking/one sweep thread, BLAS one thread, OpenMP dynamic adjustment disabled. A private first depot read existing caches with `--compiled-modules=existing`. Other user jobs were not interrupted or changed.

## Interpretation and limits

The low-input SCP cases show a modest elapsed-time penalty at 48: approximately 0.26 s for one file and 0.94 s for three files in this measured stage. Prediction is faster at 48, but slower training offsets that saving. CPU time rises substantially despite the small wall-time change. This supports treating 48 as an optional maximum, with a modest small-cohort penalty here; it does not show that 48 is faster or cheaper for these searches.

The measured workload is the complete precursor Pass-1 stage, including all its model fits and predictions. It does not benchmark MainSearch models, protein models, MBR models, or total SearchDIA runtime at 48. A 48-thread ceiling remains distinct from a fixed count of 48 for every model and batch. Results on this host do not establish an optimum on other hardware.

Process RSS is a high-water across both datasets and settings, so it cannot isolate native memory overhead at 48. Cumulative Julia allocations also do not measure peak live memory. The Julia worker count was fixed at 24, preserving the existing search-context allocation path's 24 sets of per-worker structures during seed searches.

Job 3370284 completed successfully in 23m41s; the srun step reported MaxRSS of 11,155,204 KiB (about 10.64 GiB). The accounting record is retained in `results/job_accounting.txt`.

## Artifacts

All files are under `/Users/nathanwamsley/Projects/pioneer-large-search-analysis/lgbm_small_searches_20261008T185752Z`:

- `benchmark.jl`, `benchmark.sbatch`, `configs/`, `local_manifest.json`: exact driver, allocation, configs and hashes.
- `results/samples.jsonl`: all seed, warmup and measured records, call timings, feature names, parity signatures and input hashes.
- `results/environment.json`, `results/host_cpu.txt`: runtime, source commit and hardware.
- `results/SCP_Astral_250pg_3ms/` and `results/Olsen_Exploris_3P/`: baseline search logs and reports, including final per-run counts.
- `summary.json`: five-repeat medians, ranges and paired differences, generated by `summarize.py`.
- `results/slurm.3370284.txt`, `results/COMPLETED`: job log and completion marker.

Remote scratch directory: `/scratch2/fs1/d.goldfarb/n.t.wamsley/runs/lgbm_small_searches_20261008T185752Z`. Only stats were fetched; bulk inputs, outputs and the private depot remain on RIS.
