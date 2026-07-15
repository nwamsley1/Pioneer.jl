# Two-Round Experiment-Wide Scoring — Cross-Run Consistency Features

Design & CV protocol for a second round of PSM scoring whose only new inputs are two
cross-run features derived from the first round. Offline on EWZ human/yeast (32-feature
base, **precursor-keyed folds matching production**): single-twin `twin_score` recovered
**+102,532 experiment-wide IDs at +304 false yeast** (fy rate flat: 0.130% -> 0.140%),
whereas the naive "score − best-across-all-runs" variant *tripled* false transfers. The
distinction is condition-matching: reference the score in a *comparable* run, not the
global best.

NOTE ON LEAKAGE: earlier prototype numbers (e.g. "+124,865") used random per-PSM folds,
which leak through the cross-run lookups and inflate the gain ~6-10%; the figures here are
the leakage-safe (precursor-keyed) re-measurements. A **cluster-consensus** variant (mean
`s1` over the run's whole condition-cluster instead of one twin) roughly *doubles* the
single-twin gain (**+214,631 IDs / +422 fy**, fy rate 0.138%) and is the recommended future
form; a one-hot cluster indicator is inert (run-constant). This file documents the shipped
single-twin version; see the conversation log for the cluster-consensus follow-up.

## 1. Defining "most similar" run

**Profile.** Let `V` be the precursor vocabulary (precursors observed in ≥2 runs). For each
run `R`, build a sparse profile over `V`:

```
v_R[pi] = best-per-precursor round-1 score s1 of precursor pi in run R
        = 0                                            if pi absent in R
```

The score profile encodes which precursors are confident where — i.e. sample composition.

**Similarity = cosine of L2-normalized profiles** (absent = 0 is treated as "not present"
without mean-centering the structural zeros):

```
sim(R,S) = (v_R . v_S) / (||v_R|| ||v_S||)
twin(R)  = argmax over S != R of sim(R,S)
```

**Efficient computation.** Stack L2-normalized profiles as sparse `P (n_runs x |V|)`:

```
G      = P * P^T          # n_runs x n_runs cosine matrix (one sparse Gram product)
G[R,R] = -Inf             # forbid self
twin   = argmax(G, dims=cols)
```

O(nnz·n_runs) time, O(n_runs²) memory — milliseconds at 40 runs / ~6.5M nnz, trivial to
~1000 runs. Restrict `V` to ≥2-run precursors (singletons contribute nothing).

**EWZ validation:** all 40/40 runs' most-similar neighbor is the same condition;
within-condition cosine 0.822 vs across 0.748.

**Edge cases:** precursor absent in twin → `twin_score = 0` (decoy-like, the intended
signal); single-run experiment → `twin_score = 0` (inert); ties → lowest run index; twin
need not be mutual.

## 2. Cross-run features (from round-1 scores)

```
twin_score(p) = s1( instance of precursor pi(p) in run twin(R(p)) )   # 0 if absent
```

A false transfer's twin is a comparable run where the precursor is absent, so
`twin_score ≈ 0` (decoy-like); a genuine precursor is present in its comparable run
(target-like). The target/decoy round-2 model therefore *down-weights* false transfers —
the opposite of "best-across-all-runs", whose max reaches into the yeast runs.

```
ref_irt(pi)  = irt_obs of the highest-s1 instance of precursor pi   # most confident elution
Delta_iRT(p) = | irt_obs(p) - ref_irt(pi(p)) |
```

`Delta_iRT` uses the best instance's **iRT only, never its score**, so it carries none of
the "confident-somewhere" pathology. Safe in isolation (a one-directional RT signal); adds
a small increment on top of twin_score / cluster-consensus. Left as the validated
best-instance version; a twin/cluster-referenced variant is a future A/B.

## 3. Two-round CV protocol

The "model" is two model sets (`{M1_f}`, `{M2_f}`) plus the twin map. Scores must be
out-of-fold through **both**. Two invariants make this correct:

- **Invariant A — precursor-consistent folds.** CV fold is keyed to `precursor_idx`: every
  PSM across all runs for a given precursor shares one fold. (Pioneer already does this.)
  This is *required*: the cross-run features of PSM `p` look up other instances of the
  **same precursor**, and consistent folds put those instances in the same held-out fold as
  `p`, so their `s1` is OOF-consistent with `p`.
- **Invariant B — identical folds across rounds.** A PSM's fold never changes between
  rounds 1 and 2.

```
# ---- Round 1 (standard k-fold OOF) ----
for f in folds:
    M1[f] = train_round1( PSMs where fold != f )
    s1[fold==f] = M1[f].predict( PSMs where fold == f )   # OOF

# ---- Twin map + derived features (from OOF s1) ----
P    = run x precursor profile of best-per-prec s1
twin = argmax cosine(P) off-diagonal
for every PSM p:
    twin_score[p] = s1 of pi(p) in run twin(R(p))   (0 if absent)
    Delta_iRT[p]  = | irt_obs[p] - irt_obs(argmax_s1 instance of pi(p)) |

# ---- Round 2 (SAME folds) ----
for f in folds:
    X = [ base_features , twin_score , Delta_iRT ]
    M2[f] = train_round2( X where fold != f )
    s2[fold==f] = M2[f].predict( X where fold == f )      # OOF -> experiment-wide q
```

**No leakage:** `M2[f]` excluded fold `f`; the derived features are built from `s1` of
same-precursor instances (all in fold `f`), whose `s1` came from `M1[f]` which also excluded
fold `f`. `twin(R)` is a run-level aggregate over the whole vocabulary (negligible per-PSM
influence) and only *selects* which neighbor to read. With random per-PSM folds this would
leak — the twin instance could land in a fold whose `M1` had trained on `p`; precursor-keyed
folds close that path.

## 4. Implementation checklist

- Reuse the existing `precursor_idx`-keyed CV fold assignment; assert identical across rounds.
- Persist `s1` as a per-PSM OOF column before computing any cross-run feature.
- Build the twin map once from OOF `s1` via the sparse Gram; store `twin(R)` per run.
- Compute `twin_score` and `Delta_iRT` as plain columns (twin_score=0 / single-run inert).
- Round 2 trains on base + the two columns with the SAME folds; output `s2` → exp-wide q.
- Validation gate: confirm `s2` target/decoy calibration holds and the false-yeast **rate**
  (not just count) does not rise vs the round-1 / prec_prob baseline.

## 5. Production result (EWZ, full search) — the key finding

Two-round's benefit depends **entirely on the filter ordering**. Measured end-to-end (40 files:
20 GO113 human-only + 20 GO114 human+yeast; `yeast_human.poin`; MBR off; log-odds global; ~12 min).

| config | target inst | unique | false-yeast inst | fy rate |
|---|---:|---:|---:|---:|
| global-first, vanilla | 2,039,698 | 62,946 | 10,547 | 0.517% |
| global-first, two-round | 2,045,219 | 62,295 | 10,956 | 0.536% |
| **AND, vanilla** | 1,797,588 | 62,946 | **1,494** | 0.083% |
| **AND, two-round** | **1,891,551** | 62,295 | **1,821** | 0.096% |
| DIA-NN (reference) | 2,047,659 | 63,560 | 1,648 | 0.080% |

- **Global-first** (default: filter `global_q≤1%` FIRST, then Step-11 recalc+filter `ew_q` on the
  survivors) — the ew filter is **inert** (survivors are already ~0.5% FDR), so two-round's ew-score
  improvement filters nothing → **neutral / slightly negative** (+5,521 inst, −651 uniq, +409 fy).
- **AND** (`EWFULL_AND_THRESH=0.01`: Step-10 filters `global_q≤1%` **AND** *pre-global* `ew_q≤1%`
  simultaneously, then Step-11 recalc) — the ew filter now **bites**:
  - AND alone cuts false yeast 10,547 → **1,494** (to DIA-NN's level), costing **−242k inst** (the
    "recovery band": global-passing but ew-failing).
  - **Two-round recovers +93,963 of those instances** (1,797,588 → 1,891,551) at only **+327 false
    yeast** — marginal fy rate 0.35% (99.65% clean). Closes ~40% of the ID gap to DIA-NN while
    holding false yeast near DIA-NN's level.

**Feature gains (round-2 LGBM, from `log_pass_importance`):** `twin_score` **69.4%** of total gain
(dominant, 5× the next feature — round-2 is essentially a re-rank on the nearest-run score);
`delta_irt` only **0.7%** (redundant with `twin_score` + base `irt_error` — dispensable).

**Takeaway:** the **AND filter is the false-yeast lever** (→ DIA-NN level); **two-round is the
ID-recovery lever that only expresses under the AND**. Ship them together.

### Reproduction

Base config = `EWZ_local/ewz_combined_run/config.json` (MBR off — forced on branch
`feat/two-round-twin-scoring`; log-odds global; `q_value_threshold` 0.01; no iRT/protein filters),
with `paths.results` pointed at a fresh dir. Then:

```bash
# AND + two-round (the winning config)
TWO_ROUND=1 EWFULL_AND_THRESH=0.01 \
  julia --project=. --threads 10 --gcthreads 5,1 \
  -e 'using Pioneer; SearchDIA("<config>.json")'
# baselines: drop TWO_ROUND for vanilla+AND; drop EWFULL_AND_THRESH for global-first.
```

Metrics from `<results>/precursors_long.arrow`: passing = `target & qval≤0.01 & global_qval≤0.01`;
false yeast = `species=="YEAST" & file_name startswith "GO113"`. DIA-NN reference from
`.../DIANN/both/report.parquet` (`Q.Value≤0.01`, `Decoy==0`; yeast = `Protein.Names` contains
`YEAST` & not `HUMAN`; GO113 runs). Round-2 gains print at `@user_info` via `log_pass_importance`.
