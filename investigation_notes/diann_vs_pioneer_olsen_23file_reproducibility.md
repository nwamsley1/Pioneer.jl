# DIA-NN vs Pioneer — Olsen 3-Proteome 23-file reproducibility (2026-05-15)

## Datasets & configuration

- **MS data**: `RegressionTestsLite/Olsen_3P_Exploris` — 23 raw files,
  20220909 EXPL8 Evo5 ZY MixedSpecies 500 ng, mixed-species (Human + Yeast + E. coli).
- **DIA-NN**: dataset-specific library
  (`OlsenExplorisThreeProteome500ng-11-24-2025-speclib.parquet`), report
  `pr_matrix.tsv` (75,570 precursor rows × 23 file columns).
- **Pioneer**: branch `feature/mbr-batch-f` HEAD `9273f250` (n_correlated_fragments
  >0.7 + MS1 100 ppm default + maxthreadid fix). Library
  `altimeter_3P_len7o40_ch2o3_mc1_OlsenExploris_mzsorted.poin`, MBR on, q ≤ 0.01.
- **Pioneer launch**: `julia --project=. -t 10 --gcthreads=5,1 ...`. Earlier
  single-thread runs crashed inside `compute_mbr_features_dual!` (heap
  corruption from sizing per-thread buffers with `Threads.nthreads()` instead
  of `Threads.maxthreadid()`; fixed in commit `f8273334`).

Three Pioneer outputs are referenced below:

| Run | Path | Notes |
|---|---|---|
| `100p0` | `bench_ms1ppm_100p0/`            | MS1 100 ppm, `_70` not yet added. 1,207,248 prec / 145,915 PG. |
| `100p0_keep` | `bench_ms1ppm_100p0_keep/`  | Same config + `delete_temp=false` (for temp_data inspection) + `_70` added. |
| `100p0_dump` | `bench_ms1ppm_100p0_dump/`  | Same as keep + `PIONEER_DUMP_PRE_LGBM=1` (pre-LGBM PSM dump). |

Precursor identity matched across tools by canonical key
`(stripped_sequence, sorted [(position, UniMod-id)], charge)`. Modifications
parsed from DIA-NN `Modified.Sequence` (UniMod) and Pioneer `structural_mods`.

## 1. Headline overlap (human-only target precursors, q ≤ 0.01)

| | DIA-NN | Pioneer (100p0) |
|---|---:|---:|
| Unique human-only precursors | 62,014 | 50,968 |
| In all 23 files | **37,998** | **28,092** |
| In ≥ 18 files | 49,039 (79%) | 40,034 (79%) |
| In ≥ 12 files | 55,217 (89%) | 46,412 (91%) |
| In ≥ 6  files | 59,588 (96%) | 49,616 (97%) |
| In ≥ 3  files | 61,081 (99%) | 50,209 (99%) |

DIA-NN finds ~22% more human precursors overall and ~35% more in every
file. The gap narrows at the low-N thresholds (≥3 / ≥6 files) — fraction
of identifications reaching deeper reproducibility is similar.

### 23-file set Venn

| | Count |
|---|---:|
| DIA-NN ∩ Pioneer (both in all 23 files) | 26,682 |
| DIA-NN only in 23 (not in Pioneer's 23-set) | **11,316** |
| Pioneer only in 23 (not in DIA-NN's 23-set) | 1,410 |

Library coverage: of the 11,316 DIA-NN-only-in-23 keys, **11,198 (99.0%)
are in the Pioneer library**. The gap is *not* a library coverage problem.

## 2. Where do Pioneer's "extras" land?

Of the 11,316 DIA-NN-only-in-23:

- **9,402 (83.1%) are seen in Pioneer at q ≤ 0.01 in at least one file.**
- 1,914 (16.9%) are never identified by Pioneer at q ≤ 0.01 in any file
  (93.8% of those are still in the library).

Distribution of the 9,402 by # Pioneer files at q ≤ 0.01:

| Pioneer files at q ≤ 0.01 | Count | % |
|---:|---:|---:|
| 22/23 (missing exactly 1) | **3,193** | **34.0%** |
| 21/23 | 1,645 | 17.5% |
| 20/23 | 1,088 | 11.6% |
| 19/23 | 807 | 8.6% |
| 18/23 | 633 | 6.7% |
| 12-17/23 ("often") | 1,592 | 16.9% |
| 6-11/23 ("sometimes") | 334 | 3.6% |
| 1-5/23 ("rare") | 110 | 1.2% |

**78.3% miss 1-5 files; 34% miss exactly one file.** This is a
marginal-call problem on a small per-precursor handful of "bad" files,
not a discovery problem.

## 3. Per-file PEP investigation — why are those slots missing?

Using `temp_data/main_search_psms/{file}_fold{0,1}.arrow` (PSMs that
passed Pioneer's per-file PEP ≤ 0.9 filter) and `temp_data/pre_lgbm_psms_diag/{file_idx}.arrow`
(every candidate PSM before LGBM scoring — i.e., what survived BitVec +
fragment-index + deconv).

A "missing slot" = a (precursor_idx, file_name) pair where the precursor
is in the 9k-group **and** the file is one in which Pioneer did *not*
find it at q ≤ 0.01.

### Three-way per-pair split (35,120 missing slots from the `_dump` run)

| Cause | Pairs | % |
|---|---:|---:|
| **BitVec / fragment-index filtered** (no candidate at all in `pre_lgbm_psms_diag`) | 9,252 | **26.3%** |
| **PEP-failed** (had a candidate, per-file LGBM PEP > 0.9) | **15,646** | **44.6%** |
| **PEP-passed but global qval rejected** (in `main_search_psms`, experiment-wide qval > 0.01) | 10,222 | 29.1% |

Each mode is large enough to be worth attacking independently.

### Per-precursor — pure-mode buckets (n = 9,360)

| Mode | Precursors | % |
|---|---:|---:|
| ALL missing files = BitVec-fail (precursor never considered in those files) | 862 | 9.2% |
| ALL missing files = PEP-fail (had candidates, always scored too low) | 1,879 | 20.1% |
| ALL missing files = pure global-qval rejection | 2,519 | 26.9% |
| Mixed | rest (~44%) | |

**≥1 missing file by category**:
- ≥1 BitVec-fail: 3,047 precursors (32.6%)
- ≥1 PEP-fail:    5,708 precursors (61.0%)
- ≥1 global-qval-rejected (PEP-passed): 5,422 (57.9%)

## 4. PEP distribution in the global-qval-rejected subset (9,187 pairs)

For PSMs in *present* files (Pioneer at q ≤ 0.01) — n = 175,516:

| | p5 | p25 | median | p75 | p95 |
|---|---:|---:|---:|---:|---:|
| per-file PEP | 0.0008 | 0.0073 | **0.036** | 0.132 | 0.339 |

For PSMs in *missing* files where Pioneer passed the per-file PEP ≤ 0.9 filter:

| | p5 | p25 | median | p75 | p95 |
|---|---:|---:|---:|---:|---:|
| per-file PEP | 0.025 | 0.116 | **0.235** | 0.342 | 0.404 |

Bucketed:

| per-file PEP | Count | % of 9,187 |
|---|---:|---:|
| ≤ 0.01 | **170** | 1.8% |
| 0.01 – 0.05 | 813 | 8.8% |
| 0.05 – 0.10 | 995 | 10.8% |
| 0.10 – 0.20 | 1,960 | 21.3% |
| 0.20 – 0.50 | 5,249 | 57.1% |
| 0.50 – 0.80 | 0 | 0 |
| 0.80 – 0.90 | 0 | 0 |

**None of these PSMs are anywhere near the 0.9 PEP cutoff** — the
worst is PEP ≈ 0.45. The 0.9 filter is not the proximate threshold for
this subgroup. They're rejected by the experiment-wide qval calculation
even though per-file scoring is decent.

**Anomaly: 170 cases have per-file PEP ≤ 0.01 yet global qval > 0.01.**
Per-file LGBM ranked them as "extremely confident" yet the
experiment-wide ScoringSearch Pass-1 LGBM (different feature set,
trained across all files) demoted them. Worth a separate look — could
indicate the experiment-wide model is dropping signal that the per-file
model captures, or that the global-qval aggregation has a bug.

## 5. Quantification CV (raw `peak_area` vs DIA-NN raw quants)

For 9,374 precursors that DIA-NN finds in all 23 files *and* Pioneer
sees at q ≤ 0.01 in ≥ 2 files:

| | n | median CV | p25 | p75 | p90 |
|---|---:|---:|---:|---:|---:|
| DIA-NN raw quant CV | 9,374 | **22.4%** | 18.1% | 28.9% | 41.1% |
| Pioneer raw peak_area CV | 9,374 | **32.1%** | 23.9% | 45.4% | 67.4% |

Pioneer's raw peak_area CV is ~10 pp higher than DIA-NN's raw quant CV
on the same precursors. Dennis flagged a bug in chromatogram integration
(`fast_df_sort!(chromatograms, [:precursor_idx, :isotopes_captured, :rt])`
makes rt non-monotonic within each precursor block; `integrate_chrom`
sees a non-monotonic time vector) that likely accounts for some of this
gap. Not fixed yet.

## 6. Bug surfaced & fixed during this investigation

`compute_mbr_features_dual!` (`mbr_features.jl`) sized per-thread counter
buffers with `Threads.nthreads()` but indexed them inside `@inbounds`
with `Threads.threadid()`. Default Julia 1.12 launches have
`nthreads() == 1` but `maxthreadid() == 2` (interactive thread), and
`@threads :static` schedules iterations onto threadid 2. Heap corruption
manifests later as bus errors / abort traps / TypeErrors. Fixed in commit
`f8273334` by switching to `Threads.maxthreadid()`. No other
threadid-indexed buffers in the codebase had the same pattern.

Without this fix, fresh-Julia runs reliably crashed inside
`compute_mbr_features_dual!` regardless of MS1 ppm setting.

## 7. Related parameter / feature decisions

- **MS1 ppm default raised 10 → 100 ppm** (commit `9273f250`). Sweep at
  3–100 ppm on SCP_250pg_3/24ms (≤1% ID variation across the range) and
  10 vs 100 ppm on 23-file Olsen (1,214,243 → 1,207,248 prec, −0.58%;
  PGs −0.77%) — MS1 ppm has near-zero effect on IDs. The
  `ms1_m0_mass_err_ppm` feature carries the real signal; the lookup
  window just needs to be wide enough to find the peak.

- **`n_correlated_fragments` threshold 0.9 → 0.7** (commit `9273f250`).
  Direct head-to-head: 0.7 carries ~11× the LGBM gain in ScoringSearch
  Pass-1 (1,838 vs 168) and 8-15× in per-file MainSearch. IDs:
  precursor +0.12%, PG +1.09%.

## 8. Follow-ups (not yet done)

1. **862 "all BitVec-fail" precursors** — which specific files lose them,
   is it mass-error/iRT drift in those files, or a deterministic BitVec
   LUT effect.
2. **170 "per-file PEP ≤ 0.01 but global qval > 0.01" anomaly** — does
   the experiment-wide ScoringSearch LGBM systematically demote a subset
   of per-file-confident precursors? Feature-importance comparison.
3. **Feature-distribution proxy for the 15,646 PEP-failed candidates**
   (gof / weight / matched_ratio etc. from pre_lgbm dump) to characterize
   whether they're genuinely weak or borderline-LGBM calls.
4. **Chromatogram integration bug** (Dennis): `fast_df_sort!` produces
   non-monotonic rt within a precursor block when `SeperateTraces()` is
   active (always, per `IntegrateChromatogramsSearch.jl:73`). Plausibly
   the main driver of the Pioneer/DIA-NN CV gap.

## 9. Scripts (under /tmp; not checked in)

- `/tmp/dianno_pioneer_overlap.jl` — initial 23-file overlap + CV.
- `/tmp/diann_extra_histogram.jl` — histogram of the 9,402 group.
- `/tmp/per_file_pep_investigation.jl` — pair-level PEP-pass/fail
  categorization using `main_search_psms`.
- `/tmp/per_prec_categorize.jl` — per-precursor PEP-fail categorization.
- `/tmp/pep_distribution_subset.jl` — PEP distributions present vs
  missing files.
- `/tmp/bitvec_vs_pep_split.jl` — three-way per-pair split using
  `pre_lgbm_psms_diag`.
- `/tmp/bench_olsen_ms1ppm.jl`, `/tmp/bench_olsen_ncorr70.jl`,
  `/tmp/bench_olsen_keep_temp.jl`, `/tmp/bench_olsen_dump_pre_lgbm.jl` —
  Pioneer benchmark launch scripts.
