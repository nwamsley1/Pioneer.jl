# LightGBM training (with low-data probit fallback) and best-per-precursor
# selection for MainSearch. The training core (`train_psm_classifier_with_fallback`)
# is shared with PrecursorScoringSearch — both stages use the same hyperparameters,
# the same 2-fold CV, and the same low-data fallback so the classifier behaves
# identically wherever it's applied.

# Shared LightGBM hyperparameters — used by both MainSearch (per-file) and
# PrecursorScoringSearch (experiment-wide).
# Tuned 2026-05-10. Final config picked to balance ID count vs wall time on
# Olsen Exploris one-file with the full MS1 + frag-chromatogram feature set.
#
# Iteration / learning-rate sweep results (paired EFDR, ~46 features incl.
# MS1Corr + frag-chrom features, 250k rows × 10 threads):
#
#   Config             wall    q<=.001  q<=.01   EFDR@.01  Notes
#   3000 iter, lr=.007 130 s   52,600   75,236   0.0106    Reference (heavy)
#   200 iter, lr=.10    9.6s   51,729   74,400   0.0100    Older default
#   100 iter, lr=.20    ~5s    -        -        -         CURRENT (post 2026-05-18 sweep)
#   100 iter, lr=.15    5.5s   51,099   73,235   0.0099    -2.7% IDs, 24x faster
#   50  iter, lr=.20    3.5s   51,023   73,226   0.0101    Cheap mode, -2.7% IDs, 37x faster
#   50  iter, lr=.10    3.4s   45,607   71,327   0.0100    Too-few iter for the regularization
#
# Updated 2026-05-18 (round 2): 100/0.20 → 50/0.10 after the cheap-MainSearch
# sweep on Astral / Exploris / MTAC (6 files each).
# Per-file MainSearch LGBM at 50/0.10 (lr_total = 5) vs 100/0.20 (lr_total = 20):
#   Astral:   -28s wall, -0.00% IDs, +0.18% PGs
#   Exploris: -9s  wall, +0.41% IDs, +0.23% PGs
#   MTAC:     -7s  wall, +1.28% IDs, +0.02% PGs
# IDs neutral-to-positive on all three datasets; ScoringSearch's 200/0.10
# recovers any quality lost in MainSearch (the per-file LGBM's job is just
# best-scan selection + the PEP > 0.9 gate, not final scoring).
const _SHARED_LGBM_HP_BASE = (num_iterations=50, learning_rate=0.10, max_depth=8,
                              num_leaves=63, min_data_in_leaf=300, feature_fraction=0.8,
                              bagging_fraction=0.8, bagging_freq=1, is_unbalance=false,
                              max_bin=255, lambda_l1=1.0, lambda_l2=1.0)

# Env overrides for HP sweep without recompile:
#   PIONEER_LGBM_ITERS, PIONEER_LGBM_LR, PIONEER_LGBM_MIN_DATA, PIONEER_LGBM_LEAVES
# Default (no env set) preserves the baked-in values above. Evaluated at
# call time (not module-load) so env changes between Julia sessions pick up
# without invalidating the precompile cache.
function current_shared_lgbm_hp()
    hp = _SHARED_LGBM_HP_BASE
    haskey(ENV, "PIONEER_LGBM_ITERS")    && (hp = merge(hp, (num_iterations    = parse(Int, ENV["PIONEER_LGBM_ITERS"]),)))
    haskey(ENV, "PIONEER_LGBM_LR")       && (hp = merge(hp, (learning_rate     = parse(Float64, ENV["PIONEER_LGBM_LR"]),)))
    haskey(ENV, "PIONEER_LGBM_MIN_DATA") && (hp = merge(hp, (min_data_in_leaf  = parse(Int, ENV["PIONEER_LGBM_MIN_DATA"]),)))
    haskey(ENV, "PIONEER_LGBM_LEAVES")   && (hp = merge(hp, (num_leaves        = parse(Int, ENV["PIONEER_LGBM_LEAVES"]),)))
    return hp
end
# Kept for back-compat — anything still referencing the const sees the base value.
const SHARED_LGBM_HP = _SHARED_LGBM_HP_BASE

# Per-experiment scoring LGBM hyperparams (used by PrecursorScoringSearch).
# Same shape as SHARED_LGBM_HP but lower learning rate × more iterations:
# lr_total stays in the lr_total≈20 sweet-spot identified for the per-file HP,
# while the slower lr lets boosting refine more decision boundaries on the
# larger, mixed-file PSM pool. Tunable independently of per-file via this
# constant; per-file HP is intentionally untouched.
const SCORING_LGBM_HP = (num_iterations=200, learning_rate=0.10, max_depth=8,
                         num_leaves=63, min_data_in_leaf=300, feature_fraction=0.8,
                         bagging_fraction=0.8, bagging_freq=1, is_unbalance=false,
                         max_bin=255, lambda_l1=1.0, lambda_l2=1.0)

# Per-fold training cap; folds larger than this are random-subsampled.
const SHARED_LGBM_MAX_TRAIN = 250_000

# Adaptive HP selection trigger. After running the default LGBM, if its OOF
# target count at q≤0.01 (summed across folds) is below this threshold, the
# dataset is small enough that competing HP variants of differing complexity
# is cheap and worth trying. Above the threshold the default already has
# plenty of signal to dominate.
#
# Per-file PSM counts are typically 500k-6M even for SCP, so a row-count
# trigger doesn't discriminate small from large. The OOF target count is
# the real "did the default find enough signal" signal — and it's already
# computed for the probit-vs-LGBM low-data fallback.
const ADAPTIVE_HP_TARGET_THRESHOLD = 5_000

# Adaptive HP candidates (used when n_psms < ADAPTIVE_HP_THRESHOLD).
# All inherit env-overridable defaults from current_shared_lgbm_hp() and
# only override min_data_in_leaf and num_leaves.
const ADAPTIVE_HP_OVERRIDES = (
    medium = (min_data_in_leaf = 50, num_leaves = 31),
    simple = (min_data_in_leaf = 20, num_leaves = 15),
)
# Low-data threshold per fold: below this we CV-select between LightGBM and probit.
const SHARED_LGBM_LOW_DATA_THRESHOLD = 10_000

# Diagnostic dump: when nonzero, dump pre-reduction psms (with lgbm_score) to
# DIAG_DUMP_DIR/<file_idx>.arrow inside train_lgbm_and_select_best.
const DIAG_DUMP_FILE_IDX = Ref{Int64}(0)
const DIAG_DUMP_DIR      = Ref{String}("")

"""
    train_psm_classifier_with_fallback(psms; features) -> (scores, last_classifier, info)

Two-fold CV LightGBM training with the shared hyperparameters
(`SHARED_LGBM_HP`). The function:

1. Builds a feature matrix from `features` (filtered to columns present in `psms`).
2. Splits PSMs into the existing `cv_fold` column's 0/1 folds.
3. Trains LightGBM per fold (sub-sampling each fold to at most
   `SHARED_LGBM_MAX_TRAIN`).
4. If `min(|fold0|, |fold1|) < SHARED_LGBM_LOW_DATA_THRESHOLD`, also trains a
   probit regression and picks whichever scores more targets at 1% OOF FDR.
5. Returns:
   - `all_scores :: Vector{Float64}` — fold-assigned per-PSM scores (psms order)
   - `last_classifier` — the LightGBM classifier object from the last fold,
     or `nothing` if a fold was unanimous-class
   - `info :: NamedTuple` — `(slice, fit, predict, low_data, lgbm_psms_at_1pct,
     probit_psms_at_1pct, winner)` for diagnostics

Used by both `train_lgbm_and_select_best` (MainSearch) and
`score_precursor_isotope_traces` (PrecursorScoringSearch).
"""
function train_psm_classifier_with_fallback(
    psms::DataFrame;
    features::Vector{Symbol},
    lgbm_hp = current_shared_lgbm_hp(),
    compute_infold::Bool = false,
)
    targets_col = psms[!, :target]
    n_total = nrow(psms)

    # Filter to available features
    apply_feature_blacklist!(features)
    available_features = filter(f -> hasproperty(psms, f), features)

    # Build feature matrix
    X_all = feature_matrix(psms, available_features)

    # Two-fold cross-validation using existing cv_fold column
    cv_fold = psms[!, :cv_fold]
    idx0 = findall(cv_fold .== 0)
    idx1 = findall(cv_fold .== 1)
    all_scores = Vector{Float64}(undef, n_total)
    infold_all = compute_infold ? Vector{Float64}(undef, n_total) : nothing

    MAX_TRAIN = SHARED_LGBM_MAX_TRAIN
    LOW_DATA_THRESHOLD = SHARED_LGBM_LOW_DATA_THRESHOLD

    # Fold order: (train_idx, test_idx)
    #   fold_pairs[1] = (idx1, idx0): trained on cv_fold==1, OOF for idx0,
    #                                 in-fold for idx1.
    #   fold_pairs[2] = (idx0, idx1): trained on cv_fold==0, OOF for idx1,
    #                                 in-fold for idx0.
    fold_pairs = [(idx1, idx0), (idx0, idx1)]

    _sample_pos(n_avail) = n_avail > MAX_TRAIN ? randperm(n_avail)[1:MAX_TRAIN] : collect(1:n_avail)
    sub_positions = [_sample_pos(length(tr)) for (tr, _) in fold_pairs]
    min_fold_size = min(length(idx0), length(idx1))
    low_data = min_fold_size < LOW_DATA_THRESHOLD

    # LightGBM CV. Slices train/test matrices on demand (transient peak) to avoid
    # retaining two ~400MB per-fold matrices across the whole function.
    # When compute_infold=true, also predicts on the FULL train fold per
    # iteration (rtv3: "memorization gap" — pairing OOF + in-fold scores lets
    # the downstream MBR-FTR LGBM see how overconfident the Pass-1 model is on
    # rows it was trained on, which is independent of the OOF probability).
    function _lgbm_cv(hp_overrides::NamedTuple = NamedTuple())
        hp_eff = isempty(hp_overrides) ? lgbm_hp : merge(lgbm_hp, hp_overrides)
        fold_scores  = Vector{Vector{Float64}}(undef, 2)
        infold_scores = compute_infold ? Vector{Vector{Float64}}(undef, 2) : nothing
        last_cls = nothing
        t_slice = 0.0; t_fit = 0.0; t_predict = 0.0
        for (fi, (train_idx, test_idx)) in enumerate(fold_pairs)
            sub_pos = train_idx[sub_positions[fi]]
            ts = time()
            X_tr = X_all[sub_pos, :]
            y_lbl = _prepare_labels(targets_col[sub_pos])
            t_slice += time() - ts
            if length(unique(y_lbl)) == 1
                fold_scores[fi] = fill(y_lbl[1] == 0 ? 0.0 : 1.0, length(test_idx))
                if compute_infold
                    infold_scores[fi] = fill(y_lbl[1] == 0 ? 0.0 : 1.0, length(train_idx))
                end
            else
                cls = build_lightgbm_classifier(; hp_eff...)
                tf = time(); LightGBM.fit!(cls, X_tr, y_lbl; verbosity = -1); t_fit += time() - tf
                ts2 = time(); X_te = X_all[test_idx, :]; t_slice += time() - ts2
                tp = time(); raw = LightGBM.predict(cls, X_te); t_predict += time() - tp
                fold_scores[fi] = ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
                if compute_infold
                    ts3 = time(); X_tr_full = X_all[train_idx, :]; t_slice += time() - ts3
                    tp2 = time(); raw_in = LightGBM.predict(cls, X_tr_full); t_predict += time() - tp2
                    infold_scores[fi] = ndims(raw_in) == 2 ? dropdims(raw_in; dims=2) : raw_in
                end
                last_cls = cls
            end
        end
        fold_scores, infold_scores, last_cls, (slice=t_slice, fit=t_fit, predict=t_predict)
    end

    # Probit CV — only runs in low-data branch, so build fold DataFrames lazily here.
    function _probit_cv()
        # Build Float64 fold DataFrames with intercept (ProbitRegression's z_score_bounds is typed Float64).
        function _mk_df(idx)
            df = DataFrame(Float64.(X_all[idx, :]), :auto)
            df[!, :intercept] = ones(Float64, nrow(df))
            df
        end
        df_folds = [_mk_df(idx0), _mk_df(idx1)]  # indexed by fold number (0/1) → +1
        fold_scores  = Vector{Vector{Float64}}(undef, 2)
        infold_scores = compute_infold ? Vector{Vector{Float64}}(undef, 2) : nothing
        for (fi, (train_idx, test_idx)) in enumerate(fold_pairs)
            train_fold = cv_fold[train_idx[1]] + 1  # 1 or 2
            test_fold  = cv_fold[test_idx[1]] + 1
            df_tr_full = df_folds[train_fold]
            df_te = df_folds[test_fold]
            sub_pos = sub_positions[fi]
            tr = df_tr_full[sub_pos, :]
            y_bool = Vector{Bool}(targets_col[train_idx[sub_pos]])
            β = zeros(Float64, ncol(tr))
            tr_chunks = Iterators.partition(1:length(y_bool), max(1, cld(length(y_bool), Threads.nthreads())))
            ProbitRegression(β, tr, y_bool, tr_chunks; max_iter=15)
            s = zeros(Float64, nrow(df_te))
            te_chunks = Iterators.partition(1:nrow(df_te), max(1, cld(nrow(df_te), Threads.nthreads())))
            ModelPredictProbs!(s, df_te, β, te_chunks)
            fold_scores[fi] = s
            if compute_infold
                s_in = zeros(Float64, nrow(df_tr_full))
                tr_full_chunks = Iterators.partition(1:nrow(df_tr_full),
                                                     max(1, cld(nrow(df_tr_full), Threads.nthreads())))
                ModelPredictProbs!(s_in, df_tr_full, β, tr_full_chunks)
                infold_scores[fi] = s_in
            end
        end
        fold_scores, infold_scores
    end

    function _fold_oof_psms_at_1pct(y_fold::AbstractVector, scores_on_test::AbstractVector)
        probs = Float32.(scores_on_test)
        qvals = zeros(Float16, length(probs))
        get_qvalues!(probs, y_fold, qvals)
        count(i -> qvals[i] <= Float16(0.01) && y_fold[i], eachindex(qvals))
    end

    @debug_l1 "  LightGBM CV: fold0=$(length(idx0)) fold1=$(length(idx1)) PSMs; train $(length.(sub_positions))"

    # Always run the default-HP LGBM (single source of truth for big datasets
    # and the safest fallback for small ones).
    lgbm_scores, lgbm_infold_scores, last_classifier, lgbm_timings = _lgbm_cv()
    @debug_l1 "  LightGBM timings: slice=$(round(lgbm_timings.slice, digits=2))s fit=$(round(lgbm_timings.fit, digits=2))s predict=$(round(lgbm_timings.predict, digits=2))s"

    # Helper to score a candidate by OOF targets at q≤0.01 across both folds.
    _oof_count(fs) = _fold_oof_psms_at_1pct(targets_col[idx0], fs[1]) +
                     _fold_oof_psms_at_1pct(targets_col[idx1], fs[2])

    # Adaptive model selection: if the default LGBM's OOF target count at q≤0.01
    # is below ADAPTIVE_HP_TARGET_THRESHOLD, the dataset is small enough that
    # competing alternative HP configs (medium=md50/leaves31, simple=md20/leaves15)
    # might learn better despite the default's complexity advantage on bigger
    # data. Below LOW_DATA_THRESHOLD per fold, probit also joins the competition.
    # Winner = highest OOF target count at 1% FDR.
    n_lgbm_default = _oof_count(lgbm_scores)
    adaptive = n_lgbm_default < ADAPTIVE_HP_TARGET_THRESHOLD
    candidates = [(name="lgbm_default", scores=lgbm_scores, infold=lgbm_infold_scores,
                   last=last_classifier, oof=n_lgbm_default)]

    if adaptive
        for (variant_name, hp_over) in pairs(ADAPTIVE_HP_OVERRIDES)
            sc, infold, lc, _ = _lgbm_cv(hp_over)
            push!(candidates, (name="lgbm_$(variant_name)", scores=sc, infold=infold,
                               last=lc, oof=_oof_count(sc)))
        end
    end

    info_n_probit = -1
    if low_data
        probit_scores, probit_infold_scores = _probit_cv()
        n_probit = _oof_count(probit_scores)
        info_n_probit = n_probit
        push!(candidates, (name="probit", scores=probit_scores, infold=probit_infold_scores,
                           last=nothing, oof=n_probit))
    end

    # Pick winner by OOF target count.
    best_idx = argmax([c.oof for c in candidates])
    winner = candidates[best_idx]
    info_winner = winner.name
    info_n_lgbm = n_lgbm_default

    if length(candidates) > 1
        compete_str = join(["$(c.name)=$(c.oof)" for c in candidates], " ")
        @debug_l1 "  Model selection (n=$n_total adaptive=$adaptive low_data=$low_data): $compete_str → $info_winner"
    end

    all_scores[idx0] .= winner.scores[1]
    all_scores[idx1] .= winner.scores[2]
    if compute_infold && winner.infold !== nothing
        infold_all[idx1] .= winner.infold[1]
        infold_all[idx0] .= winner.infold[2]
    end
    # Replace last_classifier with the winner's so importance() reflects the
    # chosen model (mostly cosmetic for diagnostics — only LGBM winners populate
    # this; probit winners keep nothing).
    last_classifier = winner.last

    info = (
        slice = lgbm_timings.slice,
        fit = lgbm_timings.fit,
        predict = lgbm_timings.predict,
        low_data = low_data,
        adaptive = adaptive,
        lgbm_psms_at_1pct = info_n_lgbm,
        probit_psms_at_1pct = info_n_probit,
        winner = info_winner,
        available_features = available_features,
        candidate_oof = Dict(c.name => c.oof for c in candidates),
    )
    return all_scores, infold_all, last_classifier, info
end

"""
    train_lgbm_and_select_best(psms; features) -> (best_psms, scores, timings)

Train LightGBM (with low-data probit fallback) on ALL PSMs (all scans) using
the shared `train_psm_classifier_with_fallback` helper, then select the best
scan per precursor by LightGBM score and log feature importances.

Returns:
- best_psms: DataFrame with one row per precursor (best by LightGBM score)
- scores: Vector{Float32} of LightGBM probabilities for best_psms
- timings: NamedTuple with timing breakdowns
"""
function train_lgbm_and_select_best(
    psms::DataFrame;
    features::Vector{Symbol} = collect(PRESCORE_FEATURES),
)
    t0 = time()
    # Per-precursor PSM count, broadcast to every row so MainSearch's per-file
    # LightGBM can use :n_scans as a feature.
    if nrow(psms) > 0 && !hasproperty(psms, :n_scans)
        counts = Dict{UInt32, UInt32}()
        @inbounds for pid in psms[!, :precursor_idx]::Vector{UInt32}
            counts[pid] = get(counts, pid, UInt32(0)) + UInt32(1)
        end
        psms[!, :n_scans] = UInt32[counts[pid] for pid in psms[!, :precursor_idx]]
    end
    all_scores, _, last_classifier, info = train_psm_classifier_with_fallback(psms; features=features)
    t_train_cv = time()

    # Surface the adaptive model selection result so we can see per-file picks
    # in the log without grepping debug_l1.
    if hasproperty(info, :adaptive) && info.adaptive
        oof_str = join(["$k=$v" for (k,v) in sort(collect(info.candidate_oof))], " ")
        @user_info "  Adaptive model selection (n=$(nrow(psms))): $oof_str → $(info.winner)"
    end

    model = if last_classifier !== nothing
        LightGBMModel(last_classifier, info.available_features, nothing)
    else
        LightGBMModel(nothing, info.available_features, 0.0f0)
    end

    # Add scores to psms for best-per-precursor selection
    psms[!, :lgbm_score] = Float32.(all_scores)

    # Diagnostic dump: pre-reduction PSMs WITH lgbm_score
    if DIAG_DUMP_FILE_IDX[] != 0
        diag_dir = DIAG_DUMP_DIR[]
        isdir(diag_dir) || mkpath(diag_dir)
        cols_to_keep = [:precursor_idx, :scan_idx, :target, :weight, :gof, :fitted_manhattan_distance, :lgbm_score]
        keep = intersect(cols_to_keep, propertynames(psms))
        Arrow.write(joinpath(diag_dir, "$(DIAG_DUMP_FILE_IDX[]).arrow"), psms[!, keep])
    end

    # Select best scan per precursor by LightGBM score
    psms = select_best_per_precursor!(psms, :lgbm_score)

    # Extract scores for best PSMs
    scores = psms[!, :lgbm_score]
    t_best = time()

    # Feature importances — top 15 by gain. Promoted from debug_l2 to user_info
    # for the tuning phase; revert to @debug_l2 once feature set stabilizes.
    let imp = importance(model)
        if imp !== nothing
            sorted_imp = sort(imp, by = x -> -x[2])
            lines = ["MainSearch per-file LGBM feature gains (all $(length(sorted_imp))):"]
            for (fname, gain) in sorted_imp
                push!(lines, "    $(rpad(string(fname), 40)) $(round(Int, gain))")
            end
            @user_info join(lines, "\n")
        end
    end

    timings = (
        matrix = 0.0,            # subsumed into train_cv now
        train_cv = t_train_cv - t0,
        best = t_best - t_train_cv,
    )

    return psms, Vector{Float32}(scores), timings
end

"""
    select_best_per_precursor!(psms::DataFrame, score_col::Symbol) -> DataFrame

Keeps one row per precursor_idx. Uses sortperm for contiguous group processing:
per group, selects the highest-weight PSM among those with score ≥ p75 (if ≥4 PSMs),
otherwise the highest-score PSM. Computes `irt_fwhm` (iRT span of scans with
weight ≥ 50% of peak weight), `n_above_hm`, `rt_fwhm`, and `best_rt` per precursor.
"""
function select_best_per_precursor!(psms::DataFrame, score_col::Symbol)
    scores = psms[!, score_col]::Vector{Float32}
    prec_ids = psms[!, :precursor_idx]::Vector{UInt32}
    has_irt = hasproperty(psms, :irt_obs)
    has_weight = hasproperty(psms, :weight)
    has_rt = hasproperty(psms, :rt)
    irt_obs = has_irt ? psms[!, :irt_obs]::Vector{Float32} : nothing
    rt_vals = has_rt ? psms[!, :rt]::Vector{Float32} : nothing
    weights = (has_irt && has_weight) ? psms[!, :weight]::Vector{Float32} : nothing
    # Per-scan columns used for max-across-scans aggregation features.
    gof_vec  = hasproperty(psms, :gof) ? psms[!, :gof] : nothing
    fmd_vec  = hasproperty(psms, :fitted_manhattan_distance) ? psms[!, :fitted_manhattan_distance] : nothing
    n = nrow(psms)

    # sortperm groups PSMs by precursor_idx for contiguous processing
    perm = sortperm(prec_ids)

    # Output arrays — one element per unique precursor
    keep_rows = Vector{Int}()
    sizehint!(keep_rows, n ÷ 10)
    compute_fwhm = weights !== nothing  # implies has_irt
    compute_rt = has_rt

    out_irt_fwhm = compute_fwhm ? Vector{Float32}() : nothing
    out_n_above_hm = compute_fwhm ? Vector{UInt16}() : nothing
    out_rt_fwhm = (compute_fwhm && compute_rt) ? Vector{Float32}() : nothing
    out_best_rt = compute_rt ? Vector{Float32}() : nothing
    # Re-added 2026-05-11 (orphaned in 2025-03 consolidation):
    # `smoothness` = Σ(((Δw_left + Δw_right)/w_apex)²) — squared second-derivative
    # of weight chromatogram; real peaks → smooth → low value, noise → jagged → high.
    # `num_scans` = group_len (count of PSMs / MS2 scans for this precursor).
    out_smoothness = (compute_fwhm && compute_rt) ? Vector{Float32}() : nothing
    out_num_scans  = Vector{UInt16}()
    # Per-precursor max-across-scans aggregations (added 2026-05-11 so the
    # experiment-wide LightGBM gets the same max_* signal it had pre-consolidation).
    out_max_weight = weights !== nothing ? Vector{Float32}() : nothing
    out_max_gof    = gof_vec  !== nothing ? Vector{eltype(gof_vec)}()  : nothing
    out_max_fmd    = fmd_vec  !== nothing ? Vector{eltype(fmd_vec)}()  : nothing

    if compute_fwhm
        sizehint!(out_irt_fwhm, n ÷ 10)
        sizehint!(out_n_above_hm, n ÷ 10)
    end
    if out_rt_fwhm !== nothing
        sizehint!(out_rt_fwhm, n ÷ 10)
    end
    if out_best_rt !== nothing
        sizehint!(out_best_rt, n ÷ 10)
    end
    if out_smoothness !== nothing
        sizehint!(out_smoothness, n ÷ 10)
    end
    sizehint!(out_num_scans, n ÷ 10)
    out_max_weight !== nothing && sizehint!(out_max_weight, n ÷ 10)
    out_max_gof    !== nothing && sizehint!(out_max_gof,    n ÷ 10)
    out_max_fmd    !== nothing && sizehint!(out_max_fmd,    n ÷ 10)

    # Reusable buffers
    p75_buf = Vector{Float32}(undef, 128)
    smooth_w_buf  = Vector{Float32}(undef, 128)  # weights sorted by rt
    smooth_rt_buf = Vector{Float32}(undef, 128)  # rt sorted ascending
    smooth_ord_buf = Vector{Int}(undef, 128)     # per-group sort permutation

    # Single pass: process each precursor group contiguously
    gi = 1
    while gi <= n
        pid = prec_ids[perm[gi]]
        group_start = gi
        while gi <= n && prec_ids[perm[gi]] == pid
            gi += 1
        end
        group_len = gi - group_start

        # --- Sub-pass 1: find max_weight, max_gof, max_fmd, and best-score row ---
        best_s = typemin(Float32)
        best_row = perm[group_start]
        mw = 0f0
        max_gof_val = out_max_gof !== nothing ? typemin(eltype(out_max_gof)) : nothing
        max_fmd_val = out_max_fmd !== nothing ? typemin(eltype(out_max_fmd)) : nothing

        @inbounds for k in 0:(group_len - 1)
            row = perm[group_start + k]
            s = scores[row]
            if s > best_s
                best_s = s
                best_row = row
            end
            if weights !== nothing
                w = weights[row]
                if w > mw
                    mw = w
                end
            end
            if max_gof_val !== nothing
                g = gof_vec[row]
                g > max_gof_val && (max_gof_val = g)
            end
            if max_fmd_val !== nothing
                f = fmd_vec[row]
                f > max_fmd_val && (max_fmd_val = f)
            end
        end

        # --- Sub-pass 1b: p75 re-selection by weight ---
        if weights !== nothing && group_len >= 4
            if length(p75_buf) < group_len
                resize!(p75_buf, group_len)
            end
            @inbounds for k in 0:(group_len - 1)
                p75_buf[k + 1] = scores[perm[group_start + k]]
            end
            p75_idx = ceil(Int, 0.75 * group_len)
            partialsort!(view(p75_buf, 1:group_len), p75_idx)
            score_threshold = p75_buf[p75_idx]

            best_w = typemin(Float32)
            @inbounds for k in 0:(group_len - 1)
                row = perm[group_start + k]
                scores[row] >= score_threshold || continue
                w = weights[row]
                if w > best_w
                    best_w = w
                    best_row = row
                end
            end
        end

        push!(keep_rows, best_row)

        # --- Sub-pass 2: FWHM bounds ---
        if compute_fwhm
            half_max = 0.5f0 * mw
            irt_lo = typemax(Float32)
            irt_hi = typemin(Float32)
            rt_lo = typemax(Float32)
            rt_hi = typemin(Float32)
            n_hm = UInt16(0)

            @inbounds for k in 0:(group_len - 1)
                row = perm[group_start + k]
                weights[row] >= half_max || continue
                n_hm += UInt8(1)
                irt = irt_obs[row]
                irt_lo = min(irt_lo, irt)
                irt_hi = max(irt_hi, irt)
                if compute_rt
                    rt = rt_vals[row]
                    rt_lo = min(rt_lo, rt)
                    rt_hi = max(rt_hi, rt)
                end
            end

            push!(out_irt_fwhm, n_hm > 0 ? irt_hi - irt_lo : 0f0)
            push!(out_n_above_hm, n_hm)
            if out_rt_fwhm !== nothing
                push!(out_rt_fwhm, n_hm > 0 ? rt_hi - rt_lo : 0f0)
            end
        end

        # --- Sub-pass 3: smoothness (squared second-derivative of weight chrom) ---
        if out_smoothness !== nothing
            if length(smooth_w_buf) < group_len
                resize!(smooth_w_buf,  group_len)
                resize!(smooth_rt_buf, group_len)
                resize!(smooth_ord_buf, group_len)
            end
            # Populate buffers and sort by RT
            @inbounds for k in 0:(group_len - 1)
                row = perm[group_start + k]
                smooth_w_buf[k+1]  = weights[row]
                smooth_rt_buf[k+1] = rt_vals[row]
                smooth_ord_buf[k+1] = k + 1
            end
            sort!(view(smooth_ord_buf, 1:group_len), by = ki -> smooth_rt_buf[ki])
            # Compute the roughness sum
            rough = 0f0
            if mw > 0f0
                if group_len == 1
                    # Single point — original code: rough = (−2w/w_apex)² = 4
                    rough = 4f0
                else
                    @inbounds for k in 1:group_len
                        ki = smooth_ord_buf[k]
                        w_i = smooth_w_buf[ki]
                        if k == 1
                            ki_r = smooth_ord_buf[k+1]
                            dt_r = smooth_rt_buf[ki_r] - smooth_rt_buf[ki]
                            d_r = dt_r > 0 ? (smooth_w_buf[ki_r] - w_i) / dt_r : 0f0
                            d_l = dt_r > 0 ? (-w_i) / dt_r : 0f0
                            rough += ((d_l + d_r) / mw)^2
                        elseif k == group_len
                            ki_l = smooth_ord_buf[k-1]
                            dt_l = smooth_rt_buf[ki] - smooth_rt_buf[ki_l]
                            d_l = dt_l > 0 ? (smooth_w_buf[ki_l] - w_i) / dt_l : 0f0
                            d_r = dt_l > 0 ? (-w_i) / dt_l : 0f0
                            rough += ((d_l + d_r) / mw)^2
                        else
                            ki_l = smooth_ord_buf[k-1]
                            ki_r = smooth_ord_buf[k+1]
                            dt_l = smooth_rt_buf[ki] - smooth_rt_buf[ki_l]
                            dt_r = smooth_rt_buf[ki_r] - smooth_rt_buf[ki]
                            d_l = dt_l > 0 ? (smooth_w_buf[ki_l] - w_i) / dt_l : 0f0
                            d_r = dt_r > 0 ? (smooth_w_buf[ki_r] - w_i) / dt_r : 0f0
                            rough += ((d_l + d_r) / mw)^2
                        end
                    end
                end
            end
            push!(out_smoothness, rough)
        end

        push!(out_num_scans, UInt16(group_len))

        if out_best_rt !== nothing
            push!(out_best_rt, rt_vals[best_row])
        end
        out_max_weight !== nothing && push!(out_max_weight, mw)
        out_max_gof    !== nothing && push!(out_max_gof,    max_gof_val)
        out_max_fmd    !== nothing && push!(out_max_fmd,    max_fmd_val)
    end  # while gi <= n

    # Build result from selected rows
    result = psms[keep_rows, :]

    # Attach computed columns
    if out_irt_fwhm !== nothing
        result[!, :irt_fwhm] = out_irt_fwhm
        result[!, :n_above_hm] = out_n_above_hm
    end
    if out_rt_fwhm !== nothing
        result[!, :rt_fwhm] = out_rt_fwhm
    end
    if out_best_rt !== nothing
        result[!, :best_rt] = out_best_rt
    end
    if out_smoothness !== nothing
        result[!, :smoothness] = out_smoothness
    end
    result[!, :num_scans] = out_num_scans
    out_max_weight !== nothing && (result[!, :max_weight] = out_max_weight)
    out_max_gof    !== nothing && (result[!, :max_gof]    = out_max_gof)
    out_max_fmd    !== nothing && (result[!, :max_fitted_manhattan_distance] = out_max_fmd)

    return result
end
