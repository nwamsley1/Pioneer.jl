# Batch F — Paired-counterfactual MBR + q-value FTR

Author note (2026-05-12): replaces the Phase-6 FTR controller (currently
live as of commit `f459df89`) with a paired-counterfactual training
scheme that lets the FTR LGBM learn the discriminative MBR feature
distribution directly from a 50/50 real-vs-fake training set per
candidate.

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

## Open questions — resolved defaults

1. **Threshold direction** — reject `is_false_mbr=true` rows at α=0.01,
   computed q-value-style with monotonization (use the same machinery
   as `get_qvalues!` in `src/utils/ML/fdrUtilities.jl`, just swap
   `is_decoy → is_false_mbr` and `target ↔ is_real`). Recovery:
   candidates whose `_true`-row score has q ≤ 0.01.
2. **Partner stability** — assign once at pair-regen time, stable across
   all files and all FTR runs. Re-shuffling per candidate would make the
   FTR LGBM's training set non-deterministic and would defeat the
   paired-counterfactual logic.
3. **`_false` lookup window** — exclude file=k for symmetry with
   `_true`. Otherwise a same-file partner could "win" by chance and
   leak into the donor.
4. **Lone precursors** — skip from FTR (no counterfactual = no training
   row pair). Candidate set tightens to "has a true donor AND has a
   counterfactual partner". Diagnostic: log how many candidates we lose.
5. **All 5 MBR features get the `_true`/`_false` split:**
   - `MBR_max_pair_prob_true` / `_false`
   - `MBR_log2_weight_ratio_true` / `_false`
   - `MBR_log2_explained_ratio_true` / `_false`
   - `MBR_best_irt_diff_true` / `_false`
   - `MBR_is_missing_true` / `_false`
6. **Non-MBR features in FTR LGBM** — keep them all (today's
   `FTR_FEATURES` ≈ 48). In the doubled training frame their values are
   identical between paired rows, so LGBM trees splitting on them gain
   zero IG within a pair. Across-candidate variation may still inform the
   model (e.g., "marginal MBR is OK iff the candidate has high
   spectral_contrast"). Cost is small; if importance is zero we'll drop
   them.
7. **Drop Pass 2 entirely** — `:trace_prob_mbr` exists today only to
   feed the FTR LGBM as one input feature. In Batch F the FTR LGBM
   learns MBR signal directly from the `_true`/`_false` split, so Pass 2
   is redundant. Cuts one full LGBM training (~30s of the 173s
   precursor-scoring step on 8-file Olsen).

## File-level changes

| File | Action | LOC |
|---|---|---:|
| `mbr_pairing.jl` | Add `regenerate_counterfactual_partners!` (M:N partner map). Keep `regenerate_pair_ids!` — it's used by `compute_mbr_features!` and by downstream Probit features. | +50 |
| `mbr_features.jl` | New `compute_mbr_features_dual!` — builds two parallel donor lookups (per `precursor_idx` and per `partner_pid`) and writes 10 `MBR_*_{true,false}` columns. | +120 (rewrite) |
| `mbr_ftr.jl` | New `apply_mbr_filter_paired!` — row-doubled FTR training + q-value threshold. Old `apply_mbr_filter!` stays for A/B comparison. | +180 |
| `model_config.jl::ADVANCED_FEATURE_SET` | Replace 5 MBR features with the 5 `_true` variants (do NOT add `_false` — they're FTR-only). | trivial |
| `score_psms.jl` | Drop Pass 2 call + `:trace_prob_mbr` assignment. Switch from `apply_mbr_filter!` → `apply_mbr_filter_paired!`. | -20 +10 |
| `PrecursorScoringSearch.jl` | Unchanged (qval bypass step keys on `:mbr_recovered` regardless of which FTR variant produced it). | 0 |
| `ftrUtilities.jl` | Possibly add `get_ftr_qvalues!` (q-value style with monotonization). Or reuse `get_qvalues!` directly with `is_false_mbr` in place of `is_decoy`. | +20 |

## Pair regeneration with finer iRT bins

Modify the stratum key in `regenerate_pair_ids!` and add a parallel
partner-map computation:

```julia
# Replace 10-decile iRT binning with fixed 0.5-iRT-wide bins:
irt_min, irt_max = extrema(plist_irt)
n_irt_bins = max(1, ceil(Int, (irt_max - irt_min) / 0.5f0))
bin_irt[i] = clamp(1 + floor(Int, (plist_irt[i] - irt_min) / 0.5f0), 1, n_irt_bins)

# Stratum key (unchanged structure, finer iRT):
stratum[i] = 1_000_000 * Int(cv_fold) +
             10_000   * mz_decile     +   # 1..10
                        bin_irt[i]        # 1..n_irt_bins
```

Keep the **strict 1:1 pair_id assignment** as today (each pid gets one
pair_id; pairs have one target + one decoy where possible; lone pids
get standalone pair_ids). This stays useful for:
- the qval-bypass step (which keys on rows, not on partners)
- diagnostic logging
- future re-introduction of pair-based MBR features if the
  counterfactual approach turns out worse

Then add a **separate M:N partner map** for the counterfactual:

```julia
function regenerate_counterfactual_partners!(
    psms::DataFrame,
    precursors::LibraryPrecursors;
    rng_seed::Int = 1846,           # different seed from pair_id
    irt_bin_width::Float32 = 0.5f0,
)
    # Build precursor-level table same as regenerate_pair_ids! but stratum
    # uses irt_bin_width-wide bins.
    # For each stratum:
    ts, ds = targets, decoys
    isempty(ts) || isempty(ds) && continue
    n = max(length(ts), length(ds))    # cover both sides
    shuffle!(rng, ts); shuffle!(rng, ds)
    for k in 1:n
        t = ts[mod1(k, length(ts))]
        d = ds[mod1(k, length(ds))]
        # First wins: don't overwrite an existing partner for stability.
        get!(pid_to_partner, t, d)
        get!(pid_to_partner, d, t)
    end
    # Apply to PSM rows
    psms[!, :counterfactual_partner_pid] = UInt32[
        get(pid_to_partner, p, UInt32(0)) for p in prec_idx_col
    ]
    # Log: % of pids that got a partner, distribution of stratum sizes
end
```

Use `get!` (write only on first encounter) so the partner is stable —
later targets in the same stratum that rotate around `ds` will overwrite
the FIRST occurrence; using `get!` prevents that. Result: each target's
partner is fixed at the first pass through its stratum.

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

- Counterfactual partner coverage: % of precursors with a partner assigned
- Per-stratum sizes — distribution of `min(|T|, |D|)` after fine binning
- Candidate count with vs without the `_false` requirement
- Doubled-frame size + LGBM training time
- Q-value distribution of `_true` rows; histogram of `_false` rows
- Recovery composition: T←T / T←D / D←T / D←D (same as today's log)
- Empirical FTR among recovered rows (should be ≤ α=0.01 by construction
  of q-value monotonization)

## Risks / things to watch

1. **Partner-coverage cliff at finer iRT bins.** 0.5-iRT bins on Olsen
   (iRT range ~100 units after normalization) ⇒ ~200 bins per (mz_decile,
   cv_fold) ⇒ many bins with <2 targets or <2 decoys. If counterfactual
   partner coverage drops below ~80%, widen the bin (e.g., 1.0 iRT).
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
