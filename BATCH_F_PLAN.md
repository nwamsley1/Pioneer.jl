# Batch F — Paired-counterfactual MBR + q-value FTR

Author note (2026-05-12): replaces the Phase-6 FTR controller (currently
live as of commit `f459df89`) with a paired-counterfactual training
scheme that lets the FTR LGBM learn the discriminative MBR feature
distribution directly from a 50/50 real-vs-fake training set per
candidate.

## Result on 8-file Olsen (commit 09d3916a, 2026-05-12)

| Metric | Live MBR | Batch F | Δ |
|---|---:|---:|---:|
| Precursor rows (long) | 421,419 | **426,566** | **+5,147** |
| PG rows (long) | 51,906 | **52,047** | **+141** |
| MaxLFQ protein groups | 7,390 | **7,405** | **+15** |
| Candidates | 150,248 (14.06%) | 72,806 (6.81%) | -77k |
| Recovered | 59,319 (39.48%) | 68,920 (94.66%) | +9,601 |
| Recovered targets | ~55,606 | 63,174 | +7,568 |
| Recovered decoys | ~3,713 | 5,746 | +2,033 |
| τ | 0.9464 | 0.9166 |  |

Batch F wins despite half the candidate pool — paired training drives a
much tighter threshold (94.66% recovery rate vs 39.48% live). FTR-model
importance concentrates sharply on the dominant MBR feature
(`MBR_max_pair_prob_true` gain 360,090 vs live FTR's 6,538 → 55× — the
"Siamese isolation" effect we hoped for).

Decoy recovery rose ~55%; global qval ≤ 0.01 cut absorbs most. Bears
watching but the +15 MaxLFQ groups suggests the global FDR is still
healthy.

2-file Olsen is hostile to Batch F (53% true-donor coverage on 2 files
vs ~94% on 8 files), so 2-file regressed −1,170 IDs. The proper test is
8-file (multi-condition) and the proper test wins.

## Entrapment FDR validation (8-file Olsen, entr1 library)

Re-ran the same 8-file pipeline against the entrapment library
`altimeter_3P_len7o40_ch2o3_mc1_OlsenExploris_mzsorted_entr1.poin`
(13.6M precursors, half target / half entrapment). Output:
`mbr_entr_run_01/efdr_out/`.

Search totals: 386,017 IDs / 49,376 PGs / 7,083 MaxLFQ groups.
FTR: 56,989 candidates (6.43%) → 54,983 recovered (96.48%).

Calibration errors (mean abs |EFDR − decoyFDR| over q ≤ 0.01 range):

| Method | prec_prob (per-file) | global_prob (global) |
|---|---:|---:|
| Combined EFDR | 0.0005 | 0.0004 |
| Paired EFDR | 0.0005 | 0.0004 |

At decoy FDR = 1%:
- per-file EFDR ≈ 0.70–0.73% (slightly conservative)
- global EFDR ≈ 1.0% (essentially on y=x)

Plots: `efdr_comparison_prec_prob.{png,pdf}`, `efdr_comparison_global_prob.{png,pdf}`.

Conclusion: Batch F's MBR + FTR is FDR-calibrated. The reported q-values
are honest; the α=0.01 FTR budget does not inflate the empirical FDR.

## YeastMBR Dennis FTR + EFDR (Batch F iRT-NN, entr1 lib, 2026-05-12)

20-file YeastMBR subset: 10 HelaOnly (GO113 Inj1-10) + 10 Hela+Yeast
(GO114 Inj1-10). Dennis method = compare yeast IDs in HelaOnly files
between Run A (HelaOnly alone) and Run B (Combined). Difference is the
MBR-caused false transfers.

| | Run A (HelaOnly alone) | Run B (Combined 20-file) |
|---|---:|---:|
| Total q≤.01 target rows | 402,796 | 854,261 |
| Yeast in HelaOnly | 874 | 1,280 |
| Yeast in HelaYeast | — | 40,106 |
| MBR candidates | 66,416 (7.18%) | 189,235 (9.48%) |
| MBR recovered | 64,171 (96.62%) | 181,973 (96.16%) |
| MaxLFQ groups | 5,682 | 6,464 |

Dennis FTR:
- Baseline yeast in HelaOnly (Run A): 874
- Combined yeast in HelaOnly (Run B): 1,280
- Net false transfers = 406
- Dennis FTR (vs HelaOnly total in B): 406 / 409,709 = **0.10%**
- Simple FTR (B-only):                 1,280 / 854,261 = 0.15%

Both an order of magnitude below the α=0.01 budget.

EFDR calibration (mean abs |EFDR − decoyFDR| at q ≤ 0.01):

| Method | Run A prec / global | Run B prec / global |
|---|---|---|
| Combined EFDR | 0.0004 / 0.0003 | 0.0005 / 0.0003 |
| Paired EFDR | 0.0004 / 0.0002 | 0.0004 / 0.0002 |

All within 0.0005 of perfect.

Of the 6,913 added IDs in HelaOnly between Runs A and B, only 406 are
yeast false transfers (~6%) — the rest are real human MBR rescues
within the HelaOnly pool.

## Baseline to beat

8-file Olsen Exploris MBR run (commit `f459df89`, 2026-05-12):

| | Value |
|---|---:|
| Total precursor rows (long) | 421,419 |
| Total protein-group rows (long) | 51,906 |
| MaxLFQ protein groups | 7,390 |
| Candidates | 150,248 (14.06% of PSMs) |
| τ | 0.9464 |
| Recovered | 59,319 (39.48% of candidates) |
| Empirical FTR (recovered) | T←D 0.14% + D←T 0.86% = 1.00% (= α=0.01 dead on) |
| FTR model top feature | `MBR_best_irt_diff` (gain 44,921) |

## Conceptual change

**Today** (live, mbr_ftr.jl): the FTR LGBM trains on the candidate cohort
with label `is_bad_transfer = (T←D) | (D←T)`. Most candidates are T←T
(legit) so the bad-class is a minority; the model has to find the bad
signal in a high-noise environment.

**Batch F**: for every candidate, manufacture a *counterfactual*
training row using a per-precursor decoy partner. Each candidate now
contributes TWO rows to the FTR training set:
- one row labeled "real MBR" with `MBR_*_true` features (donor = same precursor in another file)
- one row labeled "fake MBR" with `MBR_*_false` features (donor = the counterfactual partner — a decoy if candidate is target, vice versa)

All other features are identical across the pair. LGBM gain is concentrated
on the MBR feature axis by construction. Threshold is set q-value style on
the *doubled* training frame.

## Confirmed design choices

1. **Threshold direction** — reject `is_false_mbr=true` rows at α=0.01,
   computed q-value-style with monotonization. Reuse the existing
   `get_qvalues!` in `src/utils/ML/fdrUtilities.jl`: swap
   `is_decoy → is_false_mbr` (so "decoys" in q-value parlance are the
   fake-MBR rows). The threshold τ_q is the highest score for which
   the q-value (= empirical fraction of `is_false_mbr=true` rows
   among rows scoring ≥ that score, monotonized) is ≤ α. Recovery:
   candidates whose `_true`-row q-value ≤ α.

2. **Partner stability — fully experiment-global.** Build the partner
   map once over **every unique precursor_idx that appears in at least
   one PSM in the experiment** (= unique values in `best_psms.precursor_idx`).
   No per-file partnering; a precursor has one global partner used by
   every row of that precursor across every file. Use `get!` to write
   each partner exactly once for determinism.

3. **`_false` lookup window** — exclude file=k for symmetry with `_true`.

4. **No lone precursors** — every precursor must have a partner.
   Algorithm: stratify on (cv_fold × mz_decile × 0.5-iRT bin); pair within
   stratum using M:N wrap (`mod1`). If a stratum has only targets or only
   decoys, fall back outward:
     a. nearest non-empty neighbor stratum in the same cv_fold + mz_decile
        (search by widening the iRT window until a partner is found),
     b. then the same cv_fold globally (any iRT, any mz),
     c. then experiment-globally (any cv_fold).
   Log how many precursors fell back to each tier — if many fall to (b)/(c)
   the 0.5-iRT bin is too fine and should be widened. Hard requirement:
   `partner_map` covers every unique `precursor_idx` in the PSM table.

5. **All 5 MBR features get the `_true`/`_false` split:**
   - `MBR_max_pair_prob_true` / `_false`
   - `MBR_log2_weight_ratio_true` / `_false`
   - `MBR_log2_explained_ratio_true` / `_false`
   - `MBR_best_irt_diff_true` / `_false`
   - `MBR_is_missing_true` / `_false`

6. **Non-MBR features in FTR LGBM** — keep them all (~48 features). In
   the doubled training frame their values are identical between paired
   rows, so trees splitting on them gain zero IG within a pair.
   Cross-candidate variation may still inform the model (e.g., "marginal
   MBR is OK iff the candidate has high spectral_contrast"). Drop later
   if importance is zero.

7. **Drop Pass 2 entirely.** `:trace_prob_mbr` exists today only to feed
   the FTR LGBM as one input feature. In Batch F the FTR LGBM learns the
   MBR signal directly from the `_true`/`_false` split, so Pass 2 is
   redundant. Cuts one full LightGBM training (~30s of the 173s
   precursor-scoring step on 8-file Olsen). Removes `:trace_prob_mbr`
   from `FTR_FEATURES` too.

## File-level changes

| File | Action | LOC |
|---|---|---:|
| `mbr_pairing.jl` | Add `regenerate_counterfactual_partners!` (M:N partner map). Keep `regenerate_pair_ids!` — it's used by `compute_mbr_features!` and by downstream Probit features. | +50 |
| `mbr_features.jl` | New `compute_mbr_features_dual!` — builds two parallel donor lookups (per `precursor_idx` and per `partner_pid`) and writes 10 `MBR_*_{true,false}` columns. | +120 (rewrite) |
| `mbr_ftr.jl` | New `apply_mbr_filter_paired!` — row-doubled FTR training + q-value threshold. Old `apply_mbr_filter!` stays for A/B comparison. | +180 |
| `model_config.jl::ADVANCED_FEATURE_SET` | Replace 5 MBR features with the 5 `_true` variants (do NOT add `_false` — they're FTR-only). | trivial |
| `score_psms.jl` | Add `regenerate_counterfactual_partners!` call after `regenerate_pair_ids!`. Drop Pass 2 call + `:trace_prob_mbr` assignment. Switch from `apply_mbr_filter!` → `apply_mbr_filter_paired!`. | -20 +10 |
| `PrecursorScoringSearch.jl` | Unchanged (qval bypass step keys on `:mbr_recovered` regardless of which FTR variant produced it). | 0 |
| `ftrUtilities.jl` | Possibly add `get_ftr_qvalues!` (q-value style with monotonization). Or reuse `get_qvalues!` directly with `is_false_mbr` in place of `is_decoy`. | +20 |

## Pair regeneration with finer iRT bins

Keep the **strict 1:1 pair_id assignment** as today (`regenerate_pair_ids!`)
since the qval-bypass step and diagnostic logging still use `pair_id`.
Add a **separate experiment-global partner-map** for the counterfactual:

```julia
function regenerate_counterfactual_partners!(
    psms::DataFrame,
    precursors::LibraryPrecursors;
    rng_seed::Int = 1846,           # different seed from pair_id
    irt_bin_width::Float32 = 0.5f0,
)
    rng = MersenneTwister(rng_seed)
    # ── 1. Build precursor-level table over UNIQUE precursors in psms ──
    # (every unique precursor_idx in psms.precursor_idx is eligible).
    seen = Set{UInt32}(); plist = (pid, target, mz, irt, cv_fold)[]  # one row per unique pid

    # ── 2. Stratify (cv_fold × mz_decile × 0.5-iRT bin) ──
    # mz_decile from 10-quantile binning of plist.mz
    # irt_bin   from floor((plist.irt - irt_min) / 0.5)
    stratum_key = 10_000_000 * cv_fold + 10_000 * mz_decile + irt_bin

    # ── 3. Build target/decoy lists per stratum ──
    strata_t[s], strata_d[s] = targets, decoys in stratum s

    # ── 4. Partner assignment with hard "no-lone" guarantee ──
    pid_to_partner = Dict{UInt32, UInt32}()
    sizehint!(pid_to_partner, length(plist))
    n_in_stratum = 0; n_fb_neighbor = 0; n_fb_fold = 0; n_fb_global = 0

    for s in keys(merge(strata_t, strata_d))
        ts = copy(get(strata_t, s, UInt32[]))
        ds = copy(get(strata_d, s, UInt32[]))
        shuffle!(rng, ts); shuffle!(rng, ds)
        if !isempty(ts) && !isempty(ds)
            # M:N wrap within stratum
            n = max(length(ts), length(ds))
            for k in 1:n
                t = ts[mod1(k, length(ts))]
                d = ds[mod1(k, length(ds))]
                get!(pid_to_partner, t, d)   # write-on-first
                get!(pid_to_partner, d, t)
            end
            n_in_stratum += length(ts) + length(ds)
        end
        # else: one-sided stratum — handled in the fallback pass below
    end

    # ── 5. Fallback for any precursor still missing a partner ──
    # (one-sided strata leave their occupants with no partner.)
    # Tiers: (a) nearest non-empty neighbor stratum in same cv_fold+mz_decile,
    #        (b) same cv_fold (any mz/irt),
    #        (c) experiment-globally.
    # For each unpartnered precursor, sample one OPPOSITE-class precursor
    # from the smallest-tier non-empty pool that has one. Write with get!.

    for pid in unpartnered:
        partner = find_partner_widening(pid, target_pools, decoy_pools)
        get!(pid_to_partner, pid, partner)
        # bump n_fb_* counter

    # ── 6. Apply to PSM rows ──
    psms[!, :counterfactual_partner_pid] = UInt32[
        pid_to_partner[p] for p in psms.precursor_idx
    ]

    @user_info "MBR Batch F — counterfactual partner map:"
    @user_info "  unique precursors: $(length(plist))"
    @user_info "  partnered within stratum:    $n_in_stratum"
    @user_info "  fallback to neighbor:         $n_fb_neighbor"
    @user_info "  fallback to fold-global:      $n_fb_fold"
    @user_info "  fallback to experiment-global: $n_fb_global"
    @assert all(haskey(pid_to_partner, p) for p in plist.pid) \
        "Counterfactual partner map missing some precursors"
end
```

The hard `@assert` at the bottom enforces the "no lone precursors"
requirement. If any precursor falls through all three fallback tiers
(which can only happen if the experiment has *zero* opposite-class
precursors anywhere — degenerate), the run errors out loudly.

## Dual donor lookup

Replace `compute_mbr_features!` with a function that builds two
precursor-keyed top-2 dicts (not pair_id-keyed):

```julia
struct _MBRDonorEntry  # same as today
    trace_prob::Float32
    weight::Float32
    log2_intensity_explained::Float32
    irt_residual::Float32
    ms_file_idx::UInt32
end
const _SENTINEL = _MBRDonorEntry(typemin(Float32), 0, 0, 0, 0)

# Build dict: precursor_idx → (top-1, top-2) by trace_prob_prepass
top2_per_pid = Dict{UInt32, NTuple{2, _MBRDonorEntry}}()
for i in 1:n
    pid = precursor_idx[i]
    e = _MBRDonorEntry(score[i], weight[i], log2ie[i], irt_res[i], file[i])
    cur = get(top2_per_pid, pid, (_SENTINEL, _SENTINEL))
    if e.trace_prob > cur[1].trace_prob
        top2_per_pid[pid] = (e, cur[1])
    elseif e.trace_prob > cur[2].trace_prob
        top2_per_pid[pid] = (cur[1], e)
    end
end

# For each row at (precursor_idx=p, file=k, partner=q=psms.counterfactual_partner_pid[i]):
#   true_donor  = top entry of top2_per_pid[p] with file != k
#   false_donor = top entry of top2_per_pid[q] with file != k
# Compute MBR_*_true from true_donor, MBR_*_false from false_donor.
# Sentinels (no donor): -1.0, MBR_is_missing_{true,false} = true.
```

Output columns (10 numeric + 2 bool):
- `MBR_max_pair_prob_true`, `MBR_max_pair_prob_false`
- `MBR_log2_weight_ratio_true`, `MBR_log2_weight_ratio_false`
- `MBR_log2_explained_ratio_true`, `MBR_log2_explained_ratio_false`
- `MBR_best_irt_diff_true`, `MBR_best_irt_diff_false`
- `MBR_is_missing_true`, `MBR_is_missing_false`

No `MBR_is_best_decoy` needed — that label exists today only for the
`is_bad_transfer` definition which Batch F replaces.

## Row-doubled FTR training

```julia
function apply_mbr_filter_paired!(
    psms::DataFrame;
    alpha::Float32 = 0.01f0,
    q_thresh::Float32 = 0.01f0,
)
    # 1. Pre-MBR q-values (unchanged: from :trace_prob_prepass)
    pre_qvals = get_qvalues!(:trace_prob_prepass, :target)

    # Donor floor (same as today)
    prob_thresh = min(:trace_prob_prepass[pre_qvals ≤ 0.01 & target])

    # 2. Candidate gate — now requires BOTH donors:
    cand = (pre_qvals > q_thresh) &
           (MBR_max_pair_prob_true  ≥ prob_thresh) &
           (.!MBR_is_missing_true) &
           (.!MBR_is_missing_false)         # ← new: needs a counterfactual too
    n_cand = count(cand)

    # 3. Build doubled training frame
    sub = psms[cand, [non_MBR_features ∪ MBR_features ∪ :cv_fold]]
    # Top half: real-MBR rows (label is_false_mbr=false)
    # Bottom half: fake-MBR rows (label is_false_mbr=true) — copy the row,
    #              replace MBR_*_true cols with MBR_*_false values.
    X_real  = feature_matrix(sub, ftr_features_using_true_cols)
    X_fake  = feature_matrix(sub, ftr_features_using_false_cols)
    X       = vcat(X_real, X_fake)
    y_real  = .!is_false_mbr_top_half       # all true (label is_real)
    y       = vcat(falses(n_cand), trues(n_cand))   # is_false_mbr labels
    cv_double = vcat(sub.cv_fold, sub.cv_fold)       # paired-row fold stays

    # 4. 2-fold CV LGBM. Train to predict y (is_false_mbr).
    #    Output: ftr_score on every doubled row; ftr_score on row i_real
    #    should be LOW (high "is_real" prob); on i_fake should be HIGH.
    # We'll predict score = P(is_false_mbr=true). Lower = more "real".
    ftr_score = oof_predict(LGBM, X, y, cv_double)

    # 5. Q-value-style threshold on doubled frame
    # Apply get_qvalues! treating is_false_mbr=true as the "decoy" indicator
    # and using the negated score (so q-vals are computed walking
    # ftr_score ascending = most-real-first).
    qvals_double = similar(ftr_score, Float32)
    get_qvalues!(.-ftr_score, .!y, qvals_double)
    # qvals_double[i] = empirical fraction of doubled rows scoring ≤
    # ftr_score[i] that have is_false_mbr=true (monotonized).
    # τ_q = max ftr_score for which q ≤ α.

    # 6. Recovery — only the REAL-MBR rows are recovered (the top half).
    # A candidate is recovered iff its real-MBR row has q ≤ α.
    qvals_real = qvals_double[1:n_cand]      # top half corresponds to _true rows
    recovered_in_cand = qvals_real .≤ alpha
    psms[!, :mbr_recovered] = falses(n)
    psms[cand, :mbr_recovered] = recovered_in_cand
    psms[!, :MBR_transfer_candidate] = cand
    psms[!, :ftr_qval_true] = fill(NaN32, n)
    psms[cand, :ftr_qval_true] = qvals_real
end
```

## Q-value computation detail

We already have `get_qvalues!(scores, targets, qvals; doSort=true)` in
`fdrUtilities.jl`. The semantics:
- Walk scores descending
- For each row, `q(score) = (# negatives with score ≥ this) / (# total with score ≥ this)`, monotonized (running min)

For our FTR application, swap the polarity:
- "negatives" = `is_false_mbr` (= true means fake donor)
- "positives" = `.!is_false_mbr` (= true means real donor)
- Higher `ftr_score` ⇒ model thinks "more fake"; lower ⇒ "more real"
- We want a q-value where "real-and-clear" rows have low q, monotonized

So pass `scores = .-ftr_score` (negate so lower-score → higher-rank) and
`targets = .!is_false_mbr` (real = target). Then `get_qvalues!` gives us
the empirical fraction of fakes among top-ranked rows, monotonized.

Recovery: rows with `q ≤ α=0.01`.

## Implementation order

1. **Add `:counterfactual_partner_pid` column** via `regenerate_counterfactual_partners!`. Run pipeline; verify column appears in `second_pass_psms`.
2. **Add `compute_mbr_features_dual!`** with `_true`/`_false` split. Run pipeline; verify columns appear and donor-hit rates are sane.
3. **Add `apply_mbr_filter_paired!`** behind a flag. Run pipeline both ways and compare:
   - 2-file Olsen for fast iteration
   - 8-file Olsen for the final number
4. **Drop Pass 2 + `:trace_prob_mbr`** if Batch F wins. Otherwise keep both code paths for further experimentation.

## Diagnostics to log per run

- Counterfactual partner coverage by tier: stratum / fold-neighbor /
  fold-global / experiment-global (hard requirement: 100%, asserted)
- Per-stratum target/decoy sizes — distribution after fine binning
- Candidate count (note: `_false` requirement always satisfied since
  every precursor now has a partner)
- Doubled-frame size + LGBM training time
- Q-value distribution among `_true` rows; histogram among `_false` rows
- Recovery composition: T←T / T←D / D←T / D←D
- Empirical FTR among recovered rows (≤ α=0.01 by construction of
  q-value monotonization)

## Risks / things to watch

1. **Stratum starvation at 0.5-iRT bins.** 0.5-iRT bins on Olsen
   (iRT range ~100 units after normalization) ⇒ ~200 bins per (mz_decile,
   cv_fold) ⇒ many bins with one-sided occupancy. The fallback tiers
   guarantee 100% partner coverage, but if too many fall to tier (b) or
   (c) the partners are no longer "comparable" to the candidate's true
   donor (mz/irt mismatch). If fallback-fraction > ~30%, widen the iRT
   bin (1.0 or 2.0 iRT) to push the failure rate down.
2. **Doubled-frame size** for 150k candidates → 300k training rows.
   Should be fine; current FTR LGBM trains on the 150k candidates in
   ~1s; doubling stays sub-minute.
3. **Within-pair feature identity.** Most LGBM gain comes from the MBR
   features (today: MBR_best_irt_diff dominates with gain 45k/4k for the
   next-best feature). The 1:1 paired training should sharpen this
   ratio further. If LGBM somehow finds split signal in the
   "identical-between-pair" features, that's a bug — investigate.
4. **Q-value floor.** If `n_false_donors` in a top-K slice is 0 (highly
   confident reals only), the FTR q-value is 0 and we recover everyone.
   That's correct behavior but verify on small-n edge cases.
