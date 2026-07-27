# A q-value-calibrated cross-run count feature — how the current ones work, and the plan

Context: branch `feat/shadow-global-max-improvements`, 2026-07-27. Nothing implemented yet.

The proposal: *"how many runs was that precursor seen in with a first-round q-value better than
1%?"* — a calibrated version of a feature that already exists in proxy form.

---

## 1. What exists today

### 1.1 The feature family

They are the **GROUP cross-run features**, `GROUP_COLS` in
`PrecursorScoringSearch/two_round_scoring.jl`:

```julia
const GROUP_COLS = [:global_max, :global_top3_logodds, :cluster_max, :cluster_top3_logodds,
                    :global_n_present, :global_n_confident, :cluster_n_present, :cluster_n_confident]

two_round_features() = vcat(GROUP_COLS, [:delta_irt])
```

Four of the eight are already count-of-runs features:

| feature | meaning |
|---|---|
| `global_n_present` | # **other** runs where the precursor was observed at all |
| `global_n_confident` | # **other** runs where round-1 OOF score `s1 ≥ GROUP_CONF_S1` |
| `cluster_n_present` | same as `global_n_present`, restricted to the run's cluster |
| `cluster_n_confident` | same as `global_n_confident`, restricted to the run's cluster |

And the threshold is, verbatim:

```julia
const GROUP_CONF_S1 = 0.99f0  # round-1 OOF score threshold for the n_confident count (q<0.01 proxy)
```

So `global_n_confident` *is* the proposed feature, with a hardcoded score of 0.99 standing in for
q<0.01. The comment already concedes it is a proxy. The same pattern exists in the global model
(`n_runs_observed`, `n_prob_gt_0_5`, `n_prob_gt_0_9`, `n_prob_gt_0_99`) — all score-thresholded.

### 1.2 How they are computed — `_write_group_features!(file_paths, max_pid)`

Called by `write_two_round_feature_columns!` after round 1 has written each file's
`.pass1_sidecar.arrow` (containing OOF `s1 = trace_prob_prepass`). Streaming, one file in memory
at a time, accumulators sized by library precursor count (`max_pid + 1`, indexed `pid + 1`).

**Per CV fold** (0 and 1 separately — `fold_of[f]` is read from a constant `cv_fold` column):

- `files = [f for f in 1:nfiles if fold_of[f] == fold]`, and **R = length(files)** are the "runs"
  for that fold. `file_paths` here are **fold-split** files — one Arrow per (run × fold) — so
  **R = the full run count**, not half of it (see §1.3).
- **Pass 1** — stream each file, accumulating per precursor:
  - `gcnt_p[pid]` += 1 for every run where present
  - `gcnt_c[pid]` += 1 when `s1 ≥ GROUP_CONF_S1`
  - `grec[pid]` — a `TopK4` record holding the top-4 `(score, run)` pairs
  - `best_s1[pid]` / `ref_irt[pid]` — iRT of the highest-`s1` instance, for `delta_irt`
  - `run_present[ri]` — presence lists, used for clustering
- **Clustering** — `k = _group_size(R)` where `_group_size(R) = R < 9 ? 0 : min(9, R ÷ 2)`.
  If `k > 0`, runs are clustered and a parallel set of accumulators (`crec`, `ccnt_p`, `ccnt_c`)
  is built per cluster. **If `k == 0` the cluster columns simply mirror the global ones.**
- **Pass 2** — re-stream each file and emit, with **leave-one-out** arithmetic so a row never sees
  its own contribution:

```julia
ownc    = s >= GROUP_CONF_S1 ? Int32(1) : Int32(0)
gnp[i]  = Float32(gcnt_p[p + 1] - Int32(1))    # exclude own row
gnc[i]  = Float32(gcnt_c[p + 1] - ownc)        # exclude own row only if it was confident
gmax[i] = _tk_loo_max(grec[p + 1], ru)         # leave-one-out max over the top-4 record
```

- **Output** — the 9 columns go to a **row-aligned `.group.sidecar.arrow`**
  (`add_columns_via_sidecar!(ref, ...; tag = "group")`), *not* back into the ~80-column fold file.
  Round-2 read sites resolve them via `_feature_columns`, and Step 1b's `PSMFileReference`
  auto-discovers the sidecar.

### 1.3 Resolution: R is the FULL run count

**Corrected 2026-07-27.** An earlier draft of this plan claimed the fold split halves the run
count. It does not. **CV folds partition PRECURSORS, not runs:** every PSM of a given
`precursor_idx` is assigned to the same fold, verified on both dumps — **0 of 378,007** precursors
(KEAP1) and **0 of 192,108** (Olsen) appear in more than one fold. MainSearch therefore writes each
run's PSMs as two files (`*_fold0.arrow`, `*_fold1.arrow`), and a fold's file set still covers
*every run*. What is halved is the precursor population per fold, not the run count — which is the
right design, since halving the runs would cripple cross-run features by construction.

So `R = n_runs`, and the leave-one-out counts span `{0 … R-1}`:

| dataset | runs | R per fold | LOO count range | cluster on? (`_group_size(R) = R<9 ? 0 : min(9, R÷2)`) |
|---|---|---|---|---|
| SCP 250pg | 3 | 3 | {0, 1, 2} | no (R < 9) |
| KEAP1 | 6 | 6 | {0 … 5} | no (R < 9) |
| Olsen Exploris 3P | 23 | 23 | {0 … 22} | **yes, k = min(9, 11) = 9** |
| APMS20 | 60 | 60 | {0 … 59} | yes, k = 9 |

This is better for the proposal than the earlier draft implied. KEAP1 at R = 6 has real
resolution (6 distinct count values), so it is a usable testbed rather than a blind one. **Olsen
3P (R = 23) remains the preferred test** — more resolution in the count, and it is the smallest
dataset where clustering is actually active, so it is the only one of the three where
`cluster_n_confident` is not a literal duplicate of `global_n_confident`. SCP at R = 3 is still
close to blind and should not be used to judge this feature.

## 2. How shadow targets/decoys interact — and why they are the reason to expect this to work

The concern with any count-of-runs feature is the MBR trap: "seen in more runs ⇒ more likely real"
is exactly how false transfers get manufactured. The shadow-decoy design is what defuses it, and
it is already wired to cover new GROUP columns automatically.

### 2.1 The graft mechanism

In `_sample_both_folds(...; inject_shadows = true)` (`pass1_oom.jl`), for every sampled real row
`i` at slot `pos`, a **shadow** is placed at slot `k + pos`:

- `s = partner[i]` — the **nearest-iRT opposite-class row** (`_shadow_partners`)
- `y[sp] = target_c[s]` — the shadow therefore carries the **opposite label**
- `shadow_src[file] ← (s, sp)` — all **non-graft** columns are copied from the partner `s`
- `shadow_graft[file] ← (i, sp)` — the **graft** columns are copied from the **original row `i`**
- `is_graft = Bool[f in graft_cols for f in features]`, with
  `graft_cols = _mbr_graft_cols() = Tuple(c for c in GROUP_COLS)`

So a shadow is *an opposite-class row wearing the original row's cross-run features*. Every GROUP
column therefore has a **1:1 target/decoy marginal by construction**, which makes it useless as a
standalone signal and forces the model to use it only jointly with within-run evidence.

### 2.2 Why this makes the plan cheap

**Adding a column to `GROUP_COLS` wires everything at once:** it is grafted (because
`_mbr_graft_cols()` is derived from `GROUP_COLS`) *and* it enters the round-2 feature set (because
`two_round_features()` is derived from `GROUP_COLS`). One-line change, both behaviours.

The counter-example is instructive: **`delta_irt` is deliberately *not* in `GROUP_COLS`** — it is
appended separately in `two_round_features()` precisely so it is **not** grafted, since each shadow
should keep its source row's own iRT self-consistency. So "grafted or not" is a real design choice.
For a count feature, grafted is correct: the whole point is to deny it standalone use.

### 2.3 What the prior evidence says about expectations

- **Supporting:** `shadow + global_max` was the winning configuration in earlier work — so the
  graft demonstrably suppresses standalone abuse *while preserving usable signal*.
- **Bounding:** the counterfactual twin control (`offline_score_twin.jl`) — a deliberately
  *maximal* adversary that symmetrically scalar-copies every MBR feature — neutralised **all** of
  them back to baseline (+0.6–1.4%). That is the ceiling on what a symmetric-graft scheme can
  preserve, and it says the graft's strength is a dial, not a free lunch.
- **Cautioning:** in the global-model monotone sweep, `n_runs_observed` was the single feature
  responsible for essentially the entire loss — freeing it recovered everything. So a count
  feature's *conditional* relationship to truth is **not** the simple monotone increase its
  marginal statistics suggest. Better calibration could sharpen the trap as easily as fix it.

---

## 3. Why `GROUP_CONF_S1 = 0.99` is the weak point

A fixed cutoff on a raw round-1 OOF score is **not calibrated**. The `s1` distribution shifts with
library, instrument, run count, sample load, and data volume, so `s1 ≥ 0.99` corresponds to a
*different true FDR* in every experiment — and even between the two CV folds of one experiment.
The feature's meaning therefore drifts across datasets, which is precisely the property that makes
a cross-run feature generalise badly.

A q-value threshold is calibrated by construction: 1% FDR means 1% FDR everywhere.

---

## 4. Implementation plan

### 4.1 Compute the calibrated threshold

Add one cheap streaming pass at the top of `_write_group_features!`, before the per-fold loop:

```
for each file:  read (s1 from .pass1_sidecar.arrow, target from the fold file)
pool into two vectors, then:
    get_qvalues!(s1_pooled, target_pooled, q; doSort = true)
    s1_star = the smallest s1 whose q <= 0.01          # calibrated threshold
```

Everything needed is already at hand — `read_file` already loads `s1`; `target` is a column of the
fold file. Cost is one extra O(N) pass plus one sort, negligible against the two passes already
there.

**Use an experiment-wide threshold, not per-fold.** Pioneer's q-value contexts are experiment-wide
and global; there is no per-run context, and a per-fold threshold would make the feature mean
different things in fold 0 and fold 1 — reintroducing the drift we are trying to remove.

### 4.2 Emit the new columns

Accumulate `gcnt_q[pid]` (and `ccnt_q[pid]`) exactly as `gcnt_c`/`ccnt_c` are, but gated on
`s1 ≥ s1_star` instead of `≥ GROUP_CONF_S1`, with the same leave-one-out subtraction:

```julia
ownq   = s >= s1_star ? Int32(1) : Int32(0)
gnq[i] = Float32(gcnt_q[p + 1] - ownq)
```

Then add `:global_n_confident_q` and `:cluster_n_confident_q` to `add_columns_via_sidecar!` and to
`GROUP_COLS`.

### 4.3 Add alongside, do not replace

Keep `global_n_confident` (the 0.99 version) in place for the first experiments. Two reasons: the
comparison is only interpretable if the baseline arm is unchanged, and if the two features turn out
near-identical (likely at small R, see §1.3) we want to *see* that rather than infer it. Drop the
proxy only if the calibrated version demonstrably dominates.

### 4.4 Sizing the change

- `two_round_scoring.jl`: +1 pooled-q pass, +2 accumulators, +2 emitted columns, +2 `GROUP_COLS`
  entries. No new inputs, no new files, no signature changes.
- Grafting and round-2 feature inclusion come free via `_mbr_graft_cols()` / `two_round_features()`.
- Round-2 feature count goes 9 → 11 cross-run columns (80 → 82 total with ADVANCED_FEATURE_SET).

---

## 5. Evaluation

- **Testbed: Olsen Exploris 3P (23 runs) primarily, KEAP1 (6 runs) as a secondary.** Per §1.3,
  Olsen gives counts 0–22 and is the smallest dataset with clustering active. APMS20 (60 runs)
  would be better still. Avoid judging on SCP (R = 3).
- **⚠ The preserved dumps are MERGED, and the round-2 path needs FOLD-SPLIT files.** Verified: each
  file in the KEAP1 and Olsen dumps contains both `cv_fold` values, because Step 1b merges
  `*_fold{0,1}.arrow` into one Arrow per run after scoring. But `_write_group_features!` calls
  `_file_fold(tbl.cv_fold)`, which reads **row 1** as the whole file's fold — so pointing the
  round-2 path at a merged dump silently mislabels every file's fold and produces garbage GROUP
  features (no error). This is also a latent trap in the stale `offline_score.jl` harness.
  **Fix: reconstruct fold-split files by filtering each merged file on `cv_fold` and writing two
  Arrows per run.** Exact and cheap — `cv_fold` is a row column — and needs no new search.
- **Harness:** the offline round-1/round-2 path on a dumped `main_search_psms`, ~40–100 s/arm.
  Baseline must be the production path with `match_between_runs = true` so round 2 actually runs.
- **Arms:** (A) current, (B) + `*_n_confident_q`, (C) `*_n_confident_q` replacing the 0.99 version,
  (D) shadows off, as a check that the graft is doing the work we think it is.
- **Seeds:** vary `bagging_seed`; match seed count to the dataset's CV (Olsen CV ≈ 0.02–0.05%, so
  5 seeds resolve fractions of a percent). Pair every arm against the baseline at the same seed.
- **Metrics — raw IDs are NOT sufficient here.** A count feature can raise IDs by transferring
  false positives. Required: the **KEAP1-into-KO leak** count locally as a false-transfer proxy,
  and **entrapment EFDR on RIS** for the real answer (the `ent-s1` libraries are not local).

---

## 6. Risks and honest expectations

1. **Effect size scales with R.** At R = 3 (SCP) the 0.99 and q≤0.01 thresholds will usually
   select the same runs, so the two features are near-duplicates there. Resolution grows with run
   count, so the gain — if any — lives on Olsen/APMS20-scale experiments. This is a weaker
   objection than the earlier draft claimed, since R is the full run count, not half.
2. **Count features are the MBR-trap lever.** See §2.3. The realistic failure mode is +IDs, −EFDR.
3. **Feature-count dilution.** Adding 2 correlated columns to ~80 may simply add noise; the
   Tier-2/Tier-4 ablations showed this model is not short of features.
4. **The graft is a dial, not a guarantee** (§2.3, counterfactual result). If the calibrated count
   *does* help, arm (D) is what tells us whether the shadows were load-bearing.

**Decision gate:** implement §4 (it is small), test on Olsen with the leak check, and only pursue
further if arm B beats A on IDs **without** worsening the leak. If IDs rise and the leak rises
too, that is the MBR trap and the feature should be dropped regardless of the ID number.
