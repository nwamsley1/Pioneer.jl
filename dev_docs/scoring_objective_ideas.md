# Three ideas for the precursor scorer, ranked by value-per-risk

Status: **proposal.** Idea 1 is partly shipped (`71d80fa48`); ideas 2 and 3 are unbuilt.
Context: branch `feat/shadow-global-max-improvements`, 2026-07-27.

---

## 0. The diagnosis these ideas share

We train the experiment-wide precursor model on **binary logloss**, but we select and use it on
**targets at q ≤ 0.01** — a metric that depends only on the extreme tail of the score
distribution. Logloss spends model capacity uniformly, including on the vast bulk of PSMs whose
classification is never in doubt and never affects the reported ID count.

Today's HP result is direct evidence of that mismatch. Raising `learning_rate` 0.10 → 0.20 gained
**+1.06% to +1.56%** OOF targets across three datasets — i.e. the model was *under-fit for the
tail* under a loss that presumably looked converged in the bulk. If the loss and the metric were
aligned, a 2× step size would not be free money.

All three ideas below are attempts to close that gap. They differ enormously in cost and risk.

### Why the ranking is ordered by risk, not by elegance

The monotone-constraint experiment is the cautionary tale. It was principled, cheap to reason
about, and *wrong*: constraints cost **−0.49% to −2.74%** with a clean dose-response of about
−0.05% per constrained feature. The minimal arm — six directions that follow straight from the
formulas (`gof`, `fitted_hellinger`, `fitted_manhattan_distance` at +1; `poisson`, `err_norm`,
`irt_error` at −1) — already lost 0.49%. The lesson is that these features are genuinely
non-monotone at the extremes and that non-monotonicity is real, OOF-validated signal.

So: **cheap probes first, and let each one earn the next.**

---

## 1. Finish HP tuning — proven, zero risk, not exhausted

**Status: partly done.** `learning_rate` 0.10 → 0.20 shipped in `71d80fa48`; end-to-end KEAP1
went 705,113 → 711,404 precursors (+0.89%).

The decomposition matters and is worth not forgetting: the learning rate was the *free* half.

| lever | Olsen 3P | SCP 250pg | train time | stability (CV) |
|---|---|---|---|---|
| `learning_rate` 0.10 → 0.20 | +1.19% | +0.87% | 66 s vs 70 s (free) | 0.04% (SCP: *better* than baseline) |
| `num_leaves` 63→127 + `max_depth` 8→12 | +0.91% | +0.68% | 95 s (+36%) | worst of any arm |

The two are sub-additive (0.91 + 1.19 = 2.10% vs 1.56% together), so they are partly redundant
routes to the same effective complexity — which is why only the cheap one shipped.

**What is still untested**, and why each is plausible:

- **`min_gain_to_split = 1.0`** — the strongest candidate. This is *not a tuned value*: it comes
  from `build_lightgbm_classifier`'s kwarg default and is silently inherited, since
  `SCORING_LGBM_HP` never sets it. (`GLOBAL_PRECURSOR_LGBM_HP` explicitly sets `0.0`, which shows
  the default was never a deliberate choice for the scorer.) A gain floor of 1.0 on a model we
  just demonstrated is under-fit is a lot of pruning pressure.
- **`min_data_in_leaf = 300`** — with a 2.5M-row sample this forbids leaves finer than ~0.01% of
  the data. Tail-relevant structure may live in smaller pockets.
- **`num_iterations = 200`** — never re-swept after the learning-rate change. Doubling the step
  size halves the effective number of useful boosting rounds, so the optimum has certainly moved.
- **`feature_fraction = 0.8`**, `lambda_l1/l2 = 1.0`.

**Cost:** hours of compute, no code. **Risk:** none — HP changes cannot break the pipeline, only
the IDs move. **Expected value:** unknown but positive; the one HP we probed returned ~1%.

This is the only item on the list where I would predict a gain with confidence.

---

## 2. A tail-focused ranking objective — TESTED, AND IT FAILS

> ### ⛔ RESULT (2026-07-27): do not pursue. Measured on two dumps; every configuration is
> worse than logloss, and **the more tail-focused the objective, the worse it gets** — the
> opposite of the hypothesis. Details in §2a. The design write-up below is kept because the
> mechanism it got wrong is the useful part.

### 2a. What the probe found

Standalone probe (`rank_objective_probe.jl`): same dump, same 71 features, same `cv_fold`
split, same metric — only the objective differs. 2-fold OOF, metric = targets at q ≤ 0.01 via
Pioneer's own `get_qvalues!`. Note no calibration was needed after all: q-values depend only on
score ORDER, so a raw ranking score is directly comparable here.

**SCP 250pg** (229,550 rows; baseline 10,270):

| arm | k/B | OOF targets | vs logloss |
|---|---|---|---|
| A logloss (production) | — | 10,270 | — |
| lambdarank B1000 k300 | 0.30 | 10,043 | −2.2% (best) |
| lambdarank B100 k10 | 0.10 | 9,728 | −5.3% |
| lambdarank B50 k5 | 0.10 | 9,519 | −7.3% |
| lambdarank B200 k10 | 0.05 | 8,928 | −13.1% |
| lambdarank B1000 k30 | 0.03 | 7,710 | −24.9% |
| rank_xendcg B200 k10 | 0.05 | 5,336 | −48.0% |
| rank_xendcg B1000 k30 | 0.03 | 2 | −100% |
| rank_xendcg B5000 k100 | 0.02 | 1 | −100% |

**KEAP1** (2,213,582 rows; baseline 666,364) reproduces it:

| arm | OOF targets | vs logloss | train time |
|---|---|---|---|
| A logloss | 666,364 | — | 13 s |
| lambdarank B1000 k300 | 664,646 | −0.26% | 31 s |
| lambdarank B1000 k30 | 532,168 | −20.1% | 3 s |

Two clean signals, consistent across a 10× data range:

1. **Nothing beats logloss.** The best ranking arm merely *approaches* it (−0.26% on KEAP1) — and
   it does so by becoming the *least* tail-focused, i.e. by degenerating toward plain pairwise
   AUC — while costing 2.4× the training time.
2. **Tail concentration monotonically hurts.** Lower `k/B` is uniformly worse, all the way to
   −25% (lambdarank) and −100% (rank_xendcg). `rank_xendcg` is unusable here at any setting.

### 2b. Why it fails — the mechanism I got wrong

**Top-k within a random block is not the tail of the global score distribution.** These are
different things and §2's design silently conflated them.

Our target:decoy ratio is roughly 2:1, so a random block of 1000 rows holds ~660 targets. The
top 30 of that block are *easy* targets — they say nothing about where the global 1% FDR
threshold sits. The rows that actually determine that threshold are the **marginal** ones, and in
a random block they sit in the middle. So top-k weighting concentrates the objective on precisely
the least informative rows, and the tighter the truncation the more capacity is wasted there.
That is exactly the monotone degradation the sweep shows.

IR ranking objectives assume *sparse* relevance — a handful of relevant documents per query among
many irrelevant ones. At ~66% positives per group, NDCG@k is saturated on arrival.

**The useful corollary:** to focus on the global tail you would need blocks **stratified by
score** (a block containing only marginal rows), not random blocks. Which is essentially what our
**semi-supervised loop already does** — it trains on q ≤ 0.03 targets plus all decoys, selecting
by *global* score position. So Pioneer already has a better-targeted tail-focusing device than
blocked NDCG could provide, and the "mismatch" is smaller than §0 assumed.

### 2c. Original design write-up (for the record)

The idea: stop optimizing average-case classification and optimize *"get targets into the top of
the ranking"* directly. LightGBM ships two objectives that do this, and **LightGBM.jl already
supports them** — `LGBMRanking` exposes `lambdarank_truncation_level` and `label_gain`, and
`fit!` accepts `group::Vector{Int}`, which it passes through as
`LGBM_DatasetSetField(ds, "group", group)`. No custom loss, no `ccall`.

### ⚠ Correction to my first sketch of this

I originally proposed *"a single query group with `lambdarank_truncation_level` at the tail
size."* **That is wrong and would not have worked.** The docs specify the truncation level should
be "slightly higher than the desirable cutoff `k` in NDCG@k (e.g. `k + 3`)", and lambdarank's
pairwise cost scales with truncation level × group size. Our tail is ~700k targets out of a
2.5M-row sample, so a single group implies k ≈ 700,000 and a cost that is not remotely tractable.

### The design that does work: many modest blocks

Partition the training sample into **random blocks** of a few hundred to a few thousand rows,
each block a query group, with `label = target ? 1 : 0` and a small truncation level:

```
groups        = random blocks, |block| = B  (few hundred .. few thousand)
label         = 1 for target, 0 for decoy
label_gain    = [0, 1]                       # binary relevance, not the 2^i - 1 default
truncation    = k, small (e.g. 10 .. 30)
objective     = "rank_xendcg"                # see below
```

Optimising NDCG@k inside random blocks is a **top-weighted pairwise objective**: it pays for
ranking targets above decoys, weighted toward the top of each block. That is a practical
surrogate for *partial* AUC, which is the statistically correct surrogate for "how many targets
clear a high-specificity threshold". The `k/B` ratio is the tuning dial for how hard the objective
concentrates on the tail — a small `k` with a large `B` is more tail-focused.

**Prefer `rank_xendcg` over `lambdarank`.** The LightGBM docs describe it as "faster than and
achieves the similar performance as lambdarank", and it has no truncation-level cliff to tune.
Run `lambdarank` as the comparison arm, not the default.

### ⚠ The real cost: the output is no longer a probability

A ranking objective emits **raw scores**, not calibrated probabilities. That is fine for parts of
our pipeline and fatal for others:

| consumer | needs | ranking score OK? |
|---|---|---|
| `get_qvalues!` | ordering only (it sorts) | **yes** |
| `get_PEP!` | empirical binning of scores | mostly — bin edges shift, estimate is empirical |
| `MAIN_PEP_FILTER_THR`, `clamp(s, 1e-6, 1-1e-4)` | values in [0,1] | **no** |
| `_pos_logit` / `_logodds_from_topk` (global + cross-run) | `log(p/(1-p))` | **no — breaks hard** |

The cross-run log-odds combination is the blocker: it takes logits of probabilities, and feeding
it an unbounded raw score is meaningless rather than merely miscalibrated.

**Fix:** a monotone calibration step — logistic (Platt) or isotonic regression from OOF ranking
score to probability. This preserves the ranking exactly (so q-values are untouched) and restores
probability semantics for everything downstream. It **must be fit out-of-fold**, or it leaks.
This is standard and well understood, but it is real work and it means this idea is *not* the
"zero code risk" probe I first called it.

**Cost:** moderate — group construction, a calibration step, and the offline harness already
exists. **Risk:** moderate and well-understood. **Expected value:** genuinely unknown. This is a
probe of a hypothesis, not a known win.

---

## 3. nnPU — most principled, biggest lift

### The framing correction that motivates it

Our setup is usually described as semi-supervised, but more precisely it is the **mirror image of
positive-unlabeled learning**. Classic PU is *labeled positives + unlabeled*. Ours is:

- **labeled negatives** — decoys, generated, genuinely known-negative
- **a contaminated positive set** — targets, a mixture of true and false at unknown proportion

That is negative-unlabeled learning, symmetric to PU, so the PU machinery transfers with the roles
swapped. Percolator's iterate-and-relabel loop — which is what our semi-supervised loop is — is a
**heuristic** approximation of the PU risk estimator: it repeatedly guesses which "positives" are
actually negatives and hard-relabels them.

**nnPU** (Kiryo et al. 2017) replaces the heuristic with a principled non-negative risk estimator.
Its whole motivation is that the *unbiased* PU risk overfits badly with flexible models — which
describes a GBDT exactly. There is also a boosting-specific PU framework (AdaPU, Stat & Computing
2024) that is closer still to what we would want.

### Why this is more attractive here than in a typical PU problem

nnPU needs the class prior π (the true-positive fraction of the unlabeled set), and **estimating
π is normally the hard part of PU learning** — there is a whole literature on it. We get it
nearly free: the decoy-based null already tells us what fraction of targets are false at any score
threshold. That is our entire FDR apparatus. The usual blocker does not apply to us.

I could not find published work applying PU theory to DIA rescoring — the proteomics literature is
still largely Percolator-lineage. So this is open ground, with the corresponding risk that it is
open for a reason.

### What it costs to build

- **`ccall`.** LightGBM.jl exposes no custom-objective hook (no `fobj`/`feval` anywhere in the
  wrapper), but the bundled binary **does** export `LGBM_BoosterUpdateOneIterCustom` (confirmed
  via `nm` on the artifact dylib). The docs confirm the contract: custom objectives mean
  "gradients and hessians not computed directly by LightGBM" and "must be passed through
  parameters explicitly in the C API". So we hand-roll the boosting loop: build the booster, then
  per iteration compute grad/hess in Julia and call `UpdateOneIterCustom`.
- **Speed is not the obstacle.** Grad/hess is O(N) at ~10–20 flops/row: roughly 10 ms/iteration at
  2.5M rows, so ~2 s across 200 iterations against a ~70 s fit. 200 `ccall`s are free. Tree
  construction is unchanged and still dominates.
- **The same probability problem as idea 2, plus more.** With no built-in objective there is no
  sigmoid transform, so `predict` returns raw margins (expected, but note the docs do *not* state
  this — verify empirically). Additionally `is_unbalance`, `scale_pos_weight`, `sigmoid`, and
  `boost_from_average` all stop applying, and we currently rely on `is_unbalance = true`.
- We lose `fit!`'s scaffolding (validation sets, early stopping).

**Cost:** high. **Risk:** high — a hand-rolled training loop is a silent-wrong failure mode, and
it interacts with the semi-supervised loop it is meant to replace. **Expected value:** highest
ceiling of the three, and the only one that fixes a *conceptual* error rather than a tuning one.

---

## 4. Ranking and recommended sequence

| # | idea | cost | risk | evidence it helps | verdict |
|---|---|---|---|---|---|
| 1 | Finish HP tuning (`min_gain_to_split`, `num_iterations`, `min_data_in_leaf`) | compute only | none | ~1% from the one HP already probed | **do next** |
| 2 | Tail-focused ranking objective (blocked NDCG) | moderate | moderate | **TESTED — fails, −0.26% to −100%** | **abandon** (§2a) |
| 3 | nnPU custom loss via `ccall` | high | high | theory only | only if 2 shows the mismatch is real |

**Sequence (updated 2026-07-27 after running idea 2):** idea 2 has been tested and it fails, so
the sequence collapses to: **finish idea 1**, and treat idea 3 as a separate decision.

The probe did its job — it cost a couple of hours and killed a moderate-risk build. That is the
same logic that made the global model the right first monotone test: spend the cheap experiment to
decide whether the expensive one is worth running.

**What idea 2's failure does and does not tell us about idea 3.** It refutes *tail-weighting the
loss* as a route, and it weakens §0's framing: since the semi-supervised loop already focuses
training on the marginal region by global score, the objective/metric gap is smaller than assumed.
But nnPU's motivation is a **different** defect — bias from treating a contaminated positive set
as clean — not tail focus. So idea 3 is not refuted by this result; it is merely no longer
supported by it, and now has to stand on its own theory. Given it is the high-cost, high-risk item
on the list, I would want an independent reason to believe the contamination bias is material
before building it.

---

## 5. Evaluation protocol (learned the hard way this week)

Any of these must be judged the same way, or the result is noise:

- **Offline harness, Pioneer's real trainer.** `train_and_predict_pass1_oom!` on a dumped
  `main_search_psms` folder; ~40–100 s per arm versus a full search. The baseline arm must be the
  production path, not a reimplementation.
- **OOF metric.** `target_q01` from `_predict_sample_oof` (each fold scored by the opposite fold's
  booster). Never in-sample — a more complex model wins in-sample by memorising.
- **Vary `bagging_seed`, not `seed`.** The master `seed` is **inert** in our stack: LightGBM only
  derives unset sub-seeds from it, and LightGBM.jl always emits `bagging_seed`,
  `feature_fraction_seed`, and `data_random_seed` explicitly with its own defaults.
- **Match seed count to the dataset's CV, and pair on seed.** SCP 250pg read +0.13% at 3 seeds
  and +1.32% at 12, because its CV is 1.2–1.7%. KEAP1 and Olsen were fine at 3–5 seeds because
  their CVs are 0.08% and 0.02%. Report mean ± SEM and a sign test; never `max − min`.
- **More than one dataset, spanning the data-volume range.** SCP 250pg (3 files) → KEAP1 (6) →
  Olsen Exploris 3P (23). Extra capacity is most likely to overfit in the low-data regime, so SCP
  is the sensitive arm.
- **EFDR is still owed.** Every number in this document is decoy-based and therefore **blind to
  false transfer**. The `ent-s1` entrapment libraries are not on the local machine
  (`Pioneer_Human_canon_ent-s1.poin`, `Pioneer_HYE_canon_ent-s1.poin` are both missing), so
  entrapment validation has to happen on RIS. A ~1% ID gain that costs EFDR is not a gain.

---

## Sources

- [LightGBM Parameters](https://lightgbm.readthedocs.io/en/latest/Parameters.html) —
  `lambdarank_truncation_level` (default 30, "slightly higher than k"), `rank_xendcg`
  ("faster than and achieves the similar performance as lambdarank"), `label_gain`, and the
  custom-objective contract
- [LightGBM Advanced Topics](https://lightgbm.readthedocs.io/en/latest/Advanced-Topics.html)
- [Kiryo et al. 2017 — PU Learning with a Non-Negative Risk Estimator (nnPU)](https://arxiv.org/pdf/1703.00593)
- [A boosting framework for positive-unlabeled learning (AdaPU), Stat & Computing 2024](https://link.springer.com/article/10.1007/s11222-024-10529-y)
- [Positive Unlabeled Gradient Boosting](https://www.researchgate.net/publication/357907652_Positive_Unlabeled_Gradient_Boosting)
- [Class-prior estimation for PU data](https://arxiv.org/pdf/1611.01586) and
  [prior-shift variants](https://link.springer.com/article/10.1007/s10994-022-06190-z)
- [Käll et al. 2007 — Percolator (semi-supervised peptide identification)](https://www.nature.com/articles/nmeth1113)
- [ProteoTorch — deep semi-supervised rescoring](https://www.biorxiv.org/content/10.1101/2020.11.12.380881v2.full)
- [awesome-ml-pu-learning](https://github.com/JointEntropy/awesome-ml-pu-learning)
