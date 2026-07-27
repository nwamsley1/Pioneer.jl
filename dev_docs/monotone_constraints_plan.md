# Monotone constraints for Pioneer's LightGBM models — plan + proposed directions

Status: **proposal, nothing implemented.** Feasibility verified 2026-07-24 on branch
`feat/shadow-global-max-improvements`. Directions below are measured on KEAP1 (6 files) and
still need review — especially the ones flagged ⚠.

---

## 1. What the three LightGBM knobs actually do

### `monotone_constraints` (alias `mc`)
A vector of `-1 / 0 / +1`, **one entry per feature, in column order**.
`+1` = "holding all other features constant, increasing this feature can never *decrease* the
predicted probability". `-1` = never *increase* it. `0` = unconstrained.

### `monotone_constraints_method` — `basic` | `intermediate` | `advanced`
All three enforce the *same* constraint; they differ in how much **legitimate** model capacity
they destroy while doing it. When LightGBM splits a node under a constraint it must guarantee no
future descendant split can violate monotonicity, so it propagates min/max bounds down the tree:

- **`basic`** — applies one conservative global bound per constrained split. Free (no speed cost)
  but *over-constrains*: it rejects many splits that would not actually have violated anything.
- **`intermediate`** — tracks tighter per-node bounds. Docs: "may slow down training very
  slightly… much less constraining… should significantly improve the results."
- **`advanced`** — tightest bookkeeping, least over-constraining, "may slow down training."

**Measured here** (250k × 71, `SCORING_LGBM_HP`): none 0.82 s / basic 1.00 s / intermediate
0.81 s / advanced 0.85 s. `intermediate` came out *faster* than the unconstrained baseline, so
these differences are run-to-run noise at our size — **the method choice is effectively free, so
use `advanced` (or `intermediate`) and never `basic`.** Using `basic` would confound "monotonicity
hurt the model" with "the enforcement method threw away good splits".

### `monotone_penalty` (alias `mc_penalty`), default `0.0`
Forbids monotone splits on the first `floor(X)` levels of each tree, with a continuous increasing
penalty as a function of X. Effect: constrained features cannot form the top of the tree, so the
root structure is built from unconstrained features and constraints only shape the leaves.

Why we might want it: a constrained split near the root imposes its bounds on the entire subtree
below, which is where over-constraining does the most damage. A small penalty (1–2) is a way to
get the regularization benefit lower in the tree while leaving the global structure free.
Recommendation: leave at `0.0` for the first experiments, add it only as a knob if constrained
models underperform — otherwise it's a second variable confounding the first result.

### `interaction_constraints`
A list of feature-index groups; the features used within a tree must all come from a single group.
A feature may belong to several groups. Not a monotonicity tool — a separate regularizer that
composes well with one. See §6 for the one use that looks genuinely valuable here.

⚠ **Implementation trap unique to this parameter:** `interaction_constraints` **is** in
LightGBM.jl's `INDEXPARAMS` (`utils.jl:2`), so the wrapper converts Julia 1-based indices to
0-based automatically — pass 1-based indices and do *not* pre-subtract. `monotone_constraints` is
NOT in that list, so its values pass through verbatim (correct — they're values, not indices).
If anyone ever adds `monotone_constraints` to `INDEXPARAMS`, every constraint silently corrupts
(1→0, 0→−1).

---

## 2. Traps that must be handled before any of this is trusted

1. **Column-order alignment.** The vector is positional. A wrong-*length* vector fails loudly
   (LightGBM asserts `num_total_features() == monotone_constraints.size()`), but a
   **right-length, wrong-order** vector is silently accepted and mis-assigns every constraint.
   Our column order is `available_features = filter(f -> hasproperty(psms, f), features)` — a
   *dynamic* filter. So the vector MUST be built by mapping over `available_features`
   (`[get(DIRECTIONS, f, 0) for f in available_features]`), never from the static feature list.

2. **A wrong sign does not degrade gracefully — it annihilates the feature.** Verified: on data
   with a genuinely negative feature→label relationship, unconstrained predictions ran
   0.949→0.053, while `mc=[+1]` produced a **flat 0.502** everywhere. There is no partial credit
   and no error message.

3. **`missing → 0.0f0` imputation.** `_fill_column!` (lightgbm_utils.jl) imputes missing as
   `0.0f0`. Where 0 means "not measured" rather than "worst", a monotone constraint asserts
   something false about the entire missing subpopulation. This is the likely explanation for
   `ms1_m0_mass_err_ppm` measuring AUC 0.574 (i.e. "higher mass error ⇒ more target-like",
   physically backwards): 0 = no MS1 peak found, so the marginal signal is really
   *detected vs not detected*, not error magnitude. Such features must be `0` until the
   missingness is encoded separately.

4. **Out-of-range sentinels.** `OTHER_WINDOW_CORR_SENTINEL = -1.0` and
   `OTHER_WINDOW_APEX_DELTA_SENTINEL = 100.0` place "absent" at an extreme of the numeric range.
   A constraint then encodes a claim about the sentinel, not the feature. Force `0`.

5. **Marginal ≠ conditional.** Everything measured below is a *marginal* AUC. Constraints act
   *conditionally* (all else equal). These can disagree (Simpson's paradox) — the suspected case
   is `n_runs_observed`, marginally +1 but plausibly conditionally negative at equal top-k
   probability, which is exactly the false-transfer mechanism.

6. **Redundant features.** With 6 files `top_run_count = isqrt(6) = 2`, so
   `empirical_global_score` **is** `top2_logodds_score` (identical AUC 0.6387 confirms it). Any
   4–8-file experiment has 13 distinct global features, not 15.

---

## 3. Proposed directions — global model (15 features)

All 15 measure "more/better cross-run evidence", all 15 measure AUC > 0.5, and physics agrees
with data for every one. **So the first arm needs no per-feature judgement: all `+1`.** That
sidesteps sign-error risk entirely, which makes the global model the right place to start.

| # | feature | AUC | dir | reasoning |
|---|---|---|---|---|
| 1 | `empirical_global_score` | 0.639 | +1 | log-odds over top-k run probs; more evidence ⇒ higher |
| 2 | `top1_prec_prob` | 0.643 | +1 | best single-run probability |
| 3 | `top2_prec_prob` | 0.625 | +1 | 2nd-best run probability |
| 4 | `top3_prec_prob` | 0.615 | +1 | 3rd-best run probability |
| 5 | `top2_logodds_score` | 0.639 | +1 | = #1 when top_run_count=2 (redundant) |
| 6 | `top3_logodds_score` | 0.639 | +1 | log-odds over top 3 |
| 7 | `mean_prec_prob` | 0.643 | +1 | average evidence across runs |
| 8 | `std_prec_prob` | 0.577 | +1 ⚠ | spread. Marginally +1, but "inconsistent across runs" is arguably *bad*. Weakest physics case of the 15 — candidate for 0 |
| 9 | `min_prec_prob` | 0.630 | +1 | worst-run floor; higher floor ⇒ better |
| 10 | `top1_top2_gap` | 0.542 | +1 ⚠ | a large gap = one run much better than the rest = possible single-run fluke. Physics ambiguous; weak effect. Candidate for 0 |
| 11 | `top2_top3_gap` | 0.566 | +1 ⚠ | same as #10 |
| 12 | `n_runs_observed` | 0.615 | +1 ⚠⚠ | **the false-transfer lever.** +1 hard-wires "seen in more runs ⇒ more likely real", the exact mechanism shadow decoys exist to neutralize. Must be judged on EFDR, not IDs |
| 13 | `n_prob_gt_0_5` | 0.633 | +1 | count of runs with confident observation |
| 14 | `n_prob_gt_0_9` | 0.625 | +1 | stricter count |
| 15 | `n_prob_gt_0_99` | 0.610 | +1 | strictest count |

Arms to test: **all-`+1`** (simplest), then **`+1` except {8, 10, 11} = 0** (drop the
spread/gap features whose physics is ambiguous), then **`+1` except #12 = 0** (isolates the
false-transfer lever).

---

## 4. Proposed directions — scoring features (71)

Measured by AUC vs target/decoy on the KEAP1 dump (378,007 rows, 251,891 tgt / 126,116 dec),
then overridden by physics where the marginal is confounded (traps 3–4 above). Grouped by
category, since the *reason* differs by category.

### 4a. Match quality — direction from the definition (high confidence)
| feature | AUC | dir | reasoning |
|---|---|---|---|
| `gof` | 0.702 | +1 | `-log2(residuals/fitted)` — negated log, so higher = better fit |
| `fitted_manhattan_distance` | 0.706 | +1 | also a negated log despite the name "distance" |
| `fitted_hellinger` | 0.706 | +1 | negated log of Hellinger² |
| `smoothed_2d_shadow_hellinger` | 0.721 | +1 | strongest single feature; same negated-log form |
| `log2_intensity_explained` | 0.631 | +1 | more of the spectrum explained |
| `total_ions` | 0.682 | +1 | more matched M0 fragments |
| `y_count` | 0.652 | +1 | more matched y-ions |
| `longest_y` | 0.543 | +1 | longer contiguous y-series |
| `n_frags_detected_union` | 0.610 | +1 | more fragments detected |
| `n_frags_detected_intersection` | 0.528 | +1 | detected in every scan |
| `poisson` | 0.318 | **−1** | log-Poisson tail; range [−71, −1.5]. AUC<0.5 ⇒ higher is decoy-like |
| `err_norm` | 0.425 | **−1** | `2^error / total_ions` — a normalized error |
| `irt_error` | 0.359 | **−1** | deviation from predicted iRT |
| `top3_ms2_mass_error_mean` | 0.379 | **−1** | mean fragment mass error (magnitude, min 0) |
| `weight_rank_at_scan` | 0.357 | **−1** | rank, 1 = highest weight, so lower is better |
| `weight_ratio_at_scan` | 0.629 | +1 | ratio to the scan's max weight, 1.0 = best |
| `weight` | 0.612 | +1 | deconvolved abundance |
| `missed_cleavage` | 0.460 | **−1** | fewer missed cleavages = more tryptic = more real |

### 4b. Chromatographic / elution shape
| feature | AUC | dir | reasoning |
|---|---|---|---|
| `n_contiguous_scans` | 0.679 | +1 | a real peak spans consecutive scans |
| `n_correlated_fragments` | 0.628 | +1 | co-eluting fragments |
| `frag_corr_strength` | 0.644 | +1 | correlation strength |
| `frag_corr_effective_n` | 0.641 | +1 | effective count of correlated frags |
| `frag_corr_best_m0` | 0.584 | +1 | best fragment–M0 correlation |
| `n_scans` | 0.562 | +1 | more scans observed |
| `frag_apex_dispersion_irt` | 0.425 | **−1** | fragments should apex together; dispersion is bad |
| `smoothness` | 0.539 | +1 ⚠ | huge range [0.002, 8850]; check the definition's polarity |
| `irt_fwhm` | 0.507 | **0** | peak width — too narrow *and* too wide are both bad |
| `delta_frame_peak_center` | 0.486 | **0** | signed deviation, symmetric ⇒ non-monotone |
| `flanking_*` (8 features) | 0.52–0.70 | +1, except ⚠ | flanking support. `flanking_signal_support` measures 0.441 ⇒ **−1**, which contradicts its name — resolve before use. `flanking_frag_corr_strength` (0.512) and `..._effective_n` (0.521) are weak ⇒ 0 |

### 4c. MS1 isotope — mostly confounded by 0-as-missing ⚠
| feature | AUC | dir | reasoning |
|---|---|---|---|
| `ms1_isotope_dotp_m0_m1_m2` | 0.631 | +1 | isotope-pattern dot product |
| `ms1_m0_m1_m2_window_fraction` | 0.627 | +1 | fraction of the isotope envelope in-window |
| `ms1_m0_m1_m2_window_fraction_pc` | 0.623 | +1 | same, precursor-corrected |
| `ms1_m0_intensity`, `ms1_m1_intensity` | 0.611–0.613 | +1 | 0 = absent is legitimately the worst case |
| `frag1..8_smoothed_intensity` | 0.588–0.610 | +1 | same: 0 = absent = worst |
| `ms1_ms2_explained_delta`(`_pc`) | 0.605 / 0.603 | +1 ⚠ | a *delta*; verify it isn't symmetric |
| `ms1_m0_mass_err_ppm` | 0.574 | **0** ⚠⚠ | AUC says "higher error ⇒ more target", physically backwards. 0 = no MS1 peak, so the marginal is detected-vs-not (trap 3). Needs a missing indicator before it can be constrained |
| `ms1_m1_to_m0_ratio` | 0.607 | **0** ⚠ | the *deviation from predicted* is informative, not the level |
| `ms1_m1_to_m0_pred` | 0.606 | **0** | a library prediction = a coordinate, not evidence |
| `ms1_weight_apex_to_m0_apex_irt` | 0.456 | **−1** | MS1/MS2 apex disagreement |
| `ms1_m0_peak_frag_intensity_fraction` | 0.612 | +1 | fraction of the MS1 peak explained |
| `ms1_m0_peak_n_precursors` | 0.588 | **0** ⚠ | more competing precursors should be *worse*; AUC disagrees ⇒ confounded, leave free |

### 4d. Bitvec-rank features — need the rank convention confirmed
| feature | AUC | dir |
|---|---|---|
| `frag_apex_gt2x_flank_bitvec_rank` | 0.283 | **−1** (2nd-strongest feature overall) |
| `n_frags_detected_union_bitvec_rank` | 0.386 | **−1** |
| `n_correlated_fragments_bitvec_rank` | 0.396 | **−1** |
| `n_frags_detected_intersection_bitvec_rank` | 0.463 | **−1** |

All four agree: lower rank = better. Consistent with a rank where 1 is best. Confirm against
`_bitvec_pattern_rank` before committing.

### 4e. Coordinate / identity — force 0 regardless of AUC
| feature | AUC | dir | reasoning |
|---|---|---|---|
| `prec_mz` | 0.530 | **0** | m/z carries no directional evidence; useful only via interactions |
| `sequence_length` | 0.494 | **0** | already ~0.5 |
| `spectrum_peak_count` | 0.482 | **0** | spectrum property, not evidence |
| `precursor_fraction_transmitted` | 0.525 | **0** | quad-transmission normalization |
| `scan_prec_mz_n_precursors` | 0.497 | **0** | already ~0.5 |
| `log_by_ratio_m0` | 0.482 | **0** | b/y ratio — both extremes are bad |
| `Mox` | 0.467 | **−1** | **not** a coordinate: a prior on which peptides are real. Oxidised Met is more common among decoys/false targets. This is the case where naive "sequence property ⇒ neutral" gets the sign wrong |

### 4f. Sentinel-encoded — force 0 (trap 4)
`other_window_weight_corr` (sentinel −1.0), `other_window_apex_delta_irt` (sentinel 100.0),
`n_scans_other_windows`. All three measure AUC ≈ 0.5000 on this file because they're almost
entirely sentinel here.

**Tally:** ~41 `+1`, ~11 `−1`, ~19 `0`.

---

## 5. Experiment plan

Judge on **entrapment EFDR**, never raw IDs — constraints on cross-run/count features can raise
raw IDs while increasing false transfers.

**Step 0 — code.** Add `monotone_constraints` / `_method` / `_penalty` passthrough to
`build_lightgbm_classifier` (its kwarg list is currently closed). Build the vector from
`available_features`. Default off, so every existing call is unchanged and bit-identical.

**Step 1 — global model, offline, no search needed.** Rebuild the 15-feature table from the
preserved dump via `_accumulate_global_inputs!` + `_compute_global_feature_columns`
(1,040,374 precursors in seconds; `entrapment_group_id` is present for EFDR). Arms: unconstrained
baseline, all-`+1`, and the two ablations from §3. Report OOF targets at q≤0.01 + EFDR.

**Step 2 — the capacity test (the actual prize).** 2×2: {unconstrained, constrained} ×
{current HP, aggressive HP (higher `learning_rate`, deeper `max_depth`, more `num_leaves`)}.
The hypothesis worth testing is *not* "constraints help at fixed HP" but "constraints let a more
aggressive model avoid overfitting":
  - unconstrained + current  = baseline
  - unconstrained + aggressive = expected to overfit (EFDR worsens)
  - constrained + current = small effect either way
  - constrained + aggressive = **the win, if there is one**

**Step 3 — precursor scorer (71 features).** Only if Step 2 shows a gain. Use the existing
offline harness (~2 min/experiment). Stage it: 4a match-quality only (highest confidence) →
add 4b → add 4c/4d. Keep every ⚠ at 0 until its definition is confirmed.

**Step 4 — `interaction_constraints`** as a separate follow-up, not mixed into the above. See §6:
the within-run × cross-run split is arguably a *better* fit to the false-transfer problem than the
monotone work, and is testable on the same harness.

Full-search KEAP1 confirmation only after an offline arm wins.

---

## 6. Is `interaction_constraints` actually useful here? — yes, for one specific thing

Most framings of this parameter ("forbid spurious interactions") are too vague to act on. But
there is one use that maps exactly onto a problem this branch already fights.

### 6a. The compelling use: forbid within-run × cross-run interactions in the round-2 model

Round 2 trains on `[ADVANCED_FEATURE_SET ; GROUP_COLS ; delta_irt]`. The false-transfer failure
mode is not "the model uses cross-run evidence" — it's specifically:

> *when the within-run evidence is weak **but** the cross-run evidence is strong, call it a target*

That is an **interaction** between weak within-run features and strong cross-run features, and it
is exactly how MBR-style transfer manufactures false positives. Splitting the features into two
groups:

- group 1 = all within-run features (ADVANCED_FEATURE_SET)
- group 2 = the cross-run features (`GROUP_COLS` + `delta_irt`)

forbids any single tree from using both. The model then becomes **additive** across the two
domains: `score ≈ f(within-run) + g(cross-run)`. `g` can still raise a row's score, but only
*uniformly* — it can no longer apply a boost **conditional on the row looking bad within-run**.
The conditional-rescue pattern becomes structurally unrepresentable rather than merely penalised.

Why this is attractive relative to what we already do:
- It attacks a **different axis** than the symmetric shadow decoys. Shadows neutralise *standalone*
  use of cross-run features (they make the feature label-uninformative by construction, which is
  why the counterfactual twin experiment flattened every MBR feature to baseline). Interaction
  constraints instead permit the feature to carry real signal additively while removing the
  conditional interaction. The two are complementary, and this one does not require synthesising
  rows.
- It composes naturally with monotone `+1` on the cross-run features: together they give a
  GAM-like "cross-run evidence contributes a monotone additive bump", which is close to the
  behaviour we actually want from MBR.
- It is cheap to test — the same offline harness, no new features, no data synthesis.

Caveat to check: with two groups the model may simply allocate most trees to group 1 and
underuse group 2, in which case the result looks like "cross-run features did nothing." Compare
feature gains per group, not just the final metric.

### 6b. A weaker, secondary use: let coordinate features calibrate without bridging domains

Coordinate/context features (`prec_mz`, `spectrum_peak_count`, `sequence_length`,
`scan_prec_mz_n_precursors`) have a legitimate calibration role — a `gof` of 3 does not mean the
same thing at 400 *m/z* as at 1000 — but they also let a tree carve out a narrow
m/z × length × peak-count region and fit noise in it. Because a feature may appear in several
groups, they can be added to each domain group (so they may calibrate within a domain) while
still being unable to bridge domains. Lower expected value than 6a and more groups to get wrong;
do it only if 6a pays off.

### 6c. What NOT to constrain

- **MS1 × MS2.** Their agreement is genuinely informative and the interaction is real signal —
  MS1 support matters *more* when MS2 is ambiguous. Forbidding this would remove signal, not noise.
- **The 16 `frag*_int` / `frag*_smoothed_intensity` features.** These jointly encode the fragment
  intensity pattern; interactions among them are the point.
- **Anything, at fine granularity.** Many small groups will starve the split search
  (`feature_fraction = 0.8` already samples features per tree). Two groups is the experiment;
  five is a different, riskier one.

### 6d. Practical notes

Group membership is by **column index**, so the groups must be derived from `available_features`
exactly like the monotone vector — `[[findfirst(==(f), available_features) for f in group] ...]`
— using 1-based indices (LightGBM.jl converts them; see §1). Because interaction constraints act
per tree and we train 200 trees, the resulting model is still a sum over both groups, which is
precisely the additive structure intended.
