# Sciex ZT — Shape-Feature Experiments (a/b/c)

*Branch `feat/sciex-zt-scanning` · draft for review · 2026-07-10*

## Goal

Improve the within-metascan **shape** signal — how well the ±k-bin profiles match the
transmission shape — by testing better similarity *measures*, a better *aggregate*, and
better *templates*. We add candidates as **extra columns** and let the two LightGBM
models' split gains adjudicate head-to-head, then prune the losers. Breadth is cheap:
every candidate is a handful of length-13 ops inside the collapse loop we already run
(the profiles `w` and `Fbuf[b]` are already in hand), so this costs ~nothing vs. the
irreducible fragment gather.

## Design principles (non-negotiable)

- **Additive & env-gated.** All new columns behind one flag `PIONEER_ZT_SHAPE_EXP=1`.
  Default off ⇒ collapse output byte-identical to today; non-ZT datasets never collapse.
- **Keep the incumbent alongside challengers.** Don't replace `zt_tri_cosine` yet — emit
  it *and* the new measures as separate columns so gains are directly comparable in one
  run. Prune only after the data speaks.
- **Let gains adjudicate, then A/B the survivors.** 1-file run ranks all candidates by
  gain in both models; 3-file run confirms IDs/PGs for the kept subset. No feature stays
  on gain alone — it must not regress the 3-file precursor/PG numbers (75,870 / 12,398).
- **Reuse existing helpers.** `_frag_pcor` (centered Pearson), `_positive_corr_summary`,
  the `tri`/`tnorm` already built in `collapse_to_metascans`.

---

## (a) Template-match *measure* — is uncentered cosine the right metric?

`zt_tri_cosine` today (`metascan_collapse.jl:72`) is **uncentered** cosine of the weight
profile `w` vs a linear triangle `tri`. Weakness: a flat (decoy) profile still scores
cosine ≈ **0.89**, squeezing target-vs-decoy into `[0.89, 1.0]`.

Add, over the same `w` and `tri`:

| New column | Definition | Why |
|---|---|---|
| `zt_tri_pcor` | `_frag_pcor(w, tri)` — **mean-centered** Pearson | Flat → ~0, triangle → 1: full `[0,1]` spread. Consistent with the fragment axis (already Pearson). |
| `zt_tri_spectral_angle` | `1 − (2/π)·acos(clamp(cosine,0,1))` | Linearizes near 1 (standard in spectral matching); better resolution among near-perfect matches. Cheap byproduct of the cosine we already compute. |

*Prediction:* `zt_tri_pcor` ≥ `zt_tri_cosine` in gain (esp. main-search, where the flat
decoys live). If so, `zt_tri_cosine` becomes a prune candidate.

---

## (b) Summed-fragment aggregate — one high-SNR trace instead of 8 noisy ones

Let `Sfrag = Σ_b Fbuf[b]` (element-wise sum of the 8 fragment profiles across the 13
bins) — a single high-SNR "total fragment" transmission profile. Optionally a
rank-weighted sum `Sfrag_rw = Σ_b (1/√rank_b)·Fbuf[b]` (down-weight low-rank fragments).

| New column | Definition | Question it answers |
|---|---|---|
| `frag_sum_corr_weight_shape` | `_frag_pcor(Sfrag, w)` | Do the fragments **collectively** ride the precursor's own transmission profile? |
| `frag_sum_corr_tri_shape` | `_frag_pcor(Sfrag, tri)` | Do they collectively ride the **ideal** triangle? (interference-robust — independent of a possibly-contaminated `w`.) |
| `frag_sum_corr_weight_rw_shape` *(opt.)* | `_frag_pcor(Sfrag_rw, w)` | Same, rank-weighted — test whether weighting helps vs. the plain sum. |

*Rationale:* the current `frag_corr_strength_shape` rank-sums **8 separate per-fragment
corrs**, each estimated from a 13-point vector (noisy). Correlating the **summed** profile
once is higher-SNR and cheaper (1 corr vs 8). It may capture the same signal more robustly
— or it may lose the "how many fragments agree" information that `effective_n` carries.
Testing both tells us which. *Prediction:* `frag_sum_corr_weight_shape` is competitive
with `frag_corr_strength_shape` and possibly complementary (keep both if both gain).

---

## (c) Template *shape* — linear triangle vs. trapezoid vs. empirical

The template feeds `zt_tri_*` **and** the `frag_sum_corr_tri_shape` (b). A scanning-quad
transmission ⊗ isolation window is closer to a **trapezoid** than a sharp triangle; the
memory flags an "empirical-triangle template" TODO. Make the template selectable and add
template-specific match columns.

| Template | Definition | Notes |
|---|---|---|
| **Linear triangle** (incumbent) | `tri[j] = max(0, 1 − \|j\|/(k+1))` | baseline |
| **Trapezoid** `zt_trap_pcor` | flat top for `\|j\|·S ≤ W/2`, linear decay to 0 at the transmission edge; parameterize by `W` (effective ≈7 m/z, *not* the filename "5") and `S=1.022` | matches a scanning quad better; W,S via env `PIONEER_ZT_W`, default from the measured ≈7 |
| **Empirical** `zt_emp_pcor` | mean of `metascan_w_*` over confident targets (q<0.01) from the main-search's own high-prob set, normalized; Pearson of `w` vs this template | data-driven, no shape assumption; second wave (needs a within-file high-confidence pool — reuse the RT-spline target set already collected in main-search) |

All template-match columns use **centered Pearson** (from (a)'s conclusion) so we're
comparing templates, not metrics. *Prediction:* trapezoid ≥ triangle if the true profile
has a flat top; empirical ≥ both if the shape is non-parametric — but empirical risks
overfitting the calibration pool, so gate it on the 3-file (cross-run) numbers, not gain.

---

## Wiring (one place, one flag)

1. **`metascan_collapse.jl`** — inside the existing per-meta-PSM block (where `w`, `Fbuf`,
   `tri`, `tnorm` already exist), when `PIONEER_ZT_SHAPE_EXP` is set, also compute:
   `zt_tri_pcor`, `zt_tri_spectral_angle`, `zt_trap_pcor`, `frag_sum_corr_weight_shape`,
   `frag_sum_corr_tri_shape` (+ optional `_rw`, `zt_emp_pcor` in wave 2). Push to new
   accumulators; write as columns. Build the trapezoid template once (like `tri`).
2. **`MainSearch.jl`** — append the experimental symbols to the main-search feature vcat
   only when the flag is on.
3. **`model_config.jl`** — add the same symbols to `ADVANCED_FEATURE_SET`; `hasproperty`
   filters them out when the flag is off / non-ZT (existing mechanism, no risk).
4. Define a `const ZT_SHAPE_EXP_FEATURES` list so the flag toggles one vector in all three
   sites (mirrors `ZT_SHAPE_FEATURES`).

## Experiment protocol — **one candidate at a time** (decided 2026-07-10)

Test each candidate *in isolation* so gains and 3-file numbers attribute cleanly to that
one feature (the A1 lesson: isolated changes, judged on cross-run PGs). Infrastructure is
wired once — a single env flag `PIONEER_ZT_SHAPE_EXP` and a `ZT_SHAPE_EXP_FEATURES` vector
that holds **only the candidate currently under test**. Incumbents stay in for head-to-head.

Per candidate:
1. Set `ZT_SHAPE_EXP_FEATURES = [<this candidate>]`; compute its column in the collapse.
2. **Byte-identity check, flag OFF:** `zt_cmp_dumps.jl` center_only vs the current committed
   dump — every existing column bit-identical (guards "additive").
3. **1-file REP1, flag ON:** capture both LGBM gain tables; read the candidate's gain and
   rank vs its incumbent (`zt_tri_cosine` for a/c, `frag_corr_strength_shape` for b).
4. **3-file, flag ON:** IDs/PGs vs 75,870 / 12,398. Accept only if it holds-or-improves.
5. **Verdict:** keep (fold into `ZT_SHAPE_FEATURES`, possibly retiring the incumbent it
   beats) or drop. Record in the feature panel + memory. Then move to the next candidate.

Order: (a) `zt_tri_pcor` → (a) `zt_tri_spectral_angle` → (b) `frag_sum_corr_weight_shape`
→ (b) `frag_sum_corr_tri_shape` → (c) `zt_trap_pcor` → [wave 2] `zt_emp_pcor`.
`zt_tri_pcor` goes first because its verdict (centered vs uncentered) sets the metric for
the (b)/(c) template-match columns.

## Guardrails

- **Perf:** all candidates are O(13) per meta-PSM; expect collapse time flat (±1 s). Keep
  the split-timing log (`f641960d2`) on to confirm we didn't reintroduce cost.
- **Non-ZT frozen:** `_zt_k == 0` path untouched; verify one Astral/Exploris config
  byte-identical.
- **No silent template magic:** if the empirical template is used, `log()` the pool size
  and that it's per-file, so a small-pool file doesn't quietly get a garbage template.
- **Overfit watch:** a challenger that gains high in 1-file but regresses 3-file PGs is
  overfitting the single file — reject (this is exactly how A1's `frag_corr_best_shape`
  lesson generalizes: judge on cross-run PGs, not gain).

## Candidate summary (wave 1)

```
(a)  zt_tri_pcor                    Pearson(w, triangle)          — vs uncentered cosine
     zt_tri_spectral_angle          SA(w, triangle)               — linearized cosine
(b)  frag_sum_corr_weight_shape     Pearson(Σfrag, w)             — aggregate vs weight
     frag_sum_corr_tri_shape        Pearson(Σfrag, triangle)      — aggregate vs template
     [opt] frag_sum_corr_weight_rw_shape  Pearson(Σ(1/√r)·frag, w)
(c)  zt_trap_pcor                   Pearson(w, trapezoid(W,S))    — better template
     [wave 2] zt_emp_pcor           Pearson(w, empirical target profile)
```

7 new columns, one env flag, one collapse pass. Cheap to try, easy to prune.
