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
const SHARED_LGBM_HP = (num_iterations=50, learning_rate=0.20, max_depth=8,
                        num_leaves=63, min_data_in_leaf=300, feature_fraction=0.8,
                        bagging_fraction=0.8, bagging_freq=1, is_unbalance=true,
                        max_bin=255, lambda_l1=1.0, lambda_l2=1.0)

# Lightweight per-file MainSearch model (config C): smaller trees + coarser bins.
# Applies ONLY to the MainSearch per-file best-scan/PEP-gate model
# (train_lgbm_and_select_best) — NOT pass1_oom, the FTR model, or ScoringSearch,
# which keep SHARED_LGBM_HP / SCORING_LGBM_HP. Matches the env-sweep "C" arm.
const MAINSEARCH_LGBM_HP = merge(SHARED_LGBM_HP,
    (max_depth = 4, num_leaves = 15, max_bin = 63))

# Per-experiment scoring LGBM hyperparams (used by PrecursorScoringSearch).
# Same shape as SHARED_LGBM_HP but lower learning rate × more iterations:
# the slower lr lets boosting refine more decision boundaries on the
# larger, mixed-file PSM pool.
const SCORING_LGBM_HP = (num_iterations=200, learning_rate=0.10, max_depth=8,
                         num_leaves=63, min_data_in_leaf=300, feature_fraction=0.8,
                         bagging_fraction=0.8, bagging_freq=1, is_unbalance=true,
                         max_bin=255, lambda_l1=1.0, lambda_l2=1.0)

# Per-stage per-fold subsample caps.
#
# MainSearch per-file LGBM: its job is best-per-precursor selection + the PEP
# > 0.9 gate. Coarse classification — 250k subsample is enough.
#
# ScoringSearch experiment-wide LGBM: this is the final classifier driving
# target/decoy discrimination at q≤0.01. The 2026-05-19 sweep showed +0.88%
# IDs on Exploris and +1.47% on MTAC when bumping 250k → 1M; gains plateau
# around 1M-2M.
const MAIN_LGBM_MAX_TRAIN    = 250_000
const SCORING_LGBM_MAX_TRAIN = 2_500_000

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
# All inherit defaults from SHARED_LGBM_HP and only override
# min_data_in_leaf and num_leaves.
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
    train_psm_classifier_with_fallback(psms; features) -> (scores, infold_scores, last_classifier, info)

Two-fold CV LightGBM training with the shared hyperparameters
(`SHARED_LGBM_HP`). The function:

1. Builds a feature matrix from `features` (filtered to columns present in `psms`).
2. Splits PSMs into the existing `cv_fold` column's 0/1 folds.
3. Trains LightGBM per fold (sub-sampling each fold to at most `max_train`).
4. If `min(|fold0|, |fold1|) < SHARED_LGBM_LOW_DATA_THRESHOLD`, also trains a
   probit regression and picks whichever scores more targets at 1% OOF FDR.
5. Returns:
   - `all_scores :: Vector{Float64}` — fold-assigned per-PSM scores (psms order)
   - `last_classifier` — the LightGBM classifier object from the last fold,
     or `nothing` if a fold was unanimous-class
   - `info :: NamedTuple` — diagnostics plus `available_features` and a
     reusable fold predictor for reapplying scores after feature refreshes

Used by both `train_lgbm_and_select_best` (MainSearch) and
`score_precursor_isotope_traces` (PrecursorScoringSearch).

`buffers` supplies reusable backing stores for the feature matrices; the default
allocates a fresh set per call, so callers that don't reuse buffers are
unaffected. MainSearch passes a per-search set (see `MainSearchResults`) so the
population-scaled matrices are built once and reused across every file.
"""
function train_psm_classifier_with_fallback(
    psms::DataFrame;
    features::Vector{Symbol},
    lgbm_hp = SHARED_LGBM_HP,
    compute_infold::Bool = false,
    max_train::Int = MAIN_LGBM_MAX_TRAIN,
    training_mask::Union{Nothing, AbstractVector{Bool}} = nothing,
    buffers::LGBMMatrixBuffers = LGBMMatrixBuffers(),
)
    # The feature matrices built below are `unsafe_wrap` views into these backing
    # vectors, and the inner CV closures capture the wrapped matrices — a capture
    # boxes the variable and does NOT root the backing store. Preserve the
    # vectors here, in a frame no closure captures, for the whole call.
    b_all, b_train, b_test, b_infold =
        buffers.all, buffers.train, buffers.test, buffers.infold
    return GC.@preserve b_all b_train b_test b_infold _train_psm_classifier_with_fallback(
        psms;
        features = features,
        lgbm_hp = lgbm_hp,
        compute_infold = compute_infold,
        max_train = max_train,
        training_mask = training_mask,
        buffers = buffers,
    )
end

function _train_psm_classifier_with_fallback(
    psms::DataFrame;
    features::Vector{Symbol},
    lgbm_hp,
    compute_infold::Bool,
    max_train::Int,
    training_mask::Union{Nothing, AbstractVector{Bool}},
    buffers::LGBMMatrixBuffers,
)
    targets_col = psms[!, :target]
    n_total = nrow(psms)

    available_features = filter(f -> hasproperty(psms, f), features)
    if :num_enzymatic_termini in available_features
        enzymatic_termini = psms[!, :num_enzymatic_termini]
        first_value = isempty(enzymatic_termini) ? nothing : first(enzymatic_termini)
        if isempty(enzymatic_termini) ||
           all(value -> isequal(value, first_value), enzymatic_termini)
            # Full-specific libraries (and some small semi-specific files) have
            # no information in this column. Keep their model matrix identical
            # to the pre-specificity path and avoid a collinear probit input.
            deleteat!(
                available_features,
                findfirst(==(:num_enzymatic_termini), available_features),
            )
        end
    end
    n_features = length(available_features)

    # Build feature matrix
    X_all = feature_matrix!(buffers.all, psms, available_features)

    # Two-fold cross-validation using existing cv_fold column
    cv_fold = psms[!, :cv_fold]
    idx0 = findall(cv_fold .== 0)
    idx1 = findall(cv_fold .== 1)
    all_scores = Vector{Float64}(undef, n_total)
    infold_all = compute_infold ? Vector{Float64}(undef, n_total) : nothing

    MAX_TRAIN = max_train
    LOW_DATA_THRESHOLD = SHARED_LGBM_LOW_DATA_THRESHOLD

    function _fit_idx(full_train_idx)
        training_mask === nothing && return full_train_idx
        return full_train_idx[training_mask[full_train_idx]]
    end

    # Fold order: (fit_idx, full_train_idx, test_idx)
    #   fold_pairs[1] = trained on masked cv_fold==1 rows, OOF for idx0,
    #                   in-fold for all idx1 rows when requested.
    #   fold_pairs[2] = trained on masked cv_fold==0 rows, OOF for idx1,
    #                   in-fold for all idx0 rows when requested.
    fold_pairs = [(_fit_idx(idx1), idx1, idx0), (_fit_idx(idx0), idx0, idx1)]

    function _sample_pos(n_avail)
        n_avail > MAX_TRAIN || return collect(1:n_avail)
        return randperm(n_avail)[1:MAX_TRAIN]
    end
    sub_positions = [_sample_pos(length(fit_idx)) for (fit_idx, _, _) in fold_pairs]
    min_fit_size = minimum(length(fit_idx) for (fit_idx, _, _) in fold_pairs)
    low_data = min_fit_size < LOW_DATA_THRESHOLD

    # Size the slice buffers to the largest slice each will ever hold on this
    # call, up front and once. Every row count is already known here (both folds'
    # sub-sample sizes and both folds' full sizes), so nothing below resizes a
    # buffer — a resize can move the data and would dangle a live wrap.
    _size_matrix_buffer!(buffers.train, maximum(length, sub_positions; init = 0), n_features)
    max_fold_rows = max(length(idx0), length(idx1))
    _size_matrix_buffer!(buffers.test, max_fold_rows, n_features)
    compute_infold && _size_matrix_buffer!(buffers.infold, max_fold_rows, n_features)

    # LightGBM CV. Train/test matrices are gathered into the caller's pre-sized
    # buffers rather than sliced into fresh matrices — same values, no per-file
    # allocation. `bufs` is an explicit parameter (not a capture) so the closure
    # doesn't box it.
    # When compute_infold=true, also predicts on the FULL train fold per
    # iteration (rtv3: "memorization gap" — pairing OOF + in-fold scores lets
    # the downstream MBR-FTR LGBM see how overconfident the Pass-1 model is on
    # rows it was trained on, which is independent of the OOF probability).
    function _lgbm_cv(bufs::LGBMMatrixBuffers, hp_overrides::NamedTuple = NamedTuple())
        hp_eff = isempty(hp_overrides) ? lgbm_hp : merge(lgbm_hp, hp_overrides)
        fold_scores  = Vector{Vector{Float64}}(undef, 2)
        infold_scores = compute_infold ? Vector{Vector{Float64}}(undef, 2) : nothing
        fold_predictors = Vector{Any}(undef, 2)
        last_cls = nothing
        t_slice = 0.0; t_fit = 0.0; t_predict = 0.0
        for (fi, (fit_idx, full_train_idx, test_idx)) in enumerate(fold_pairs)
            sub_pos = fit_idx[sub_positions[fi]]
            ts = time()
            X_tr = gather_rows!(bufs.train, X_all, sub_pos)
            y_lbl = _prepare_labels(targets_col[sub_pos])
            t_slice += time() - ts
            if isempty(y_lbl) || length(unique(y_lbl)) == 1
                constant_score = isempty(y_lbl) || y_lbl[1] == 0 ? 0.0 : 1.0
                fold_scores[fi] = fill(constant_score, length(test_idx))
                fold_predictors[fi] = (
                    kind = :constant,
                    value = constant_score,
                    model = nothing,
                    beta = nothing,
                )
                if compute_infold
                    infold_scores[fi] = fill(constant_score, length(full_train_idx))
                end
            else
                cls = build_lightgbm_classifier(; hp_eff...)
                tf = time(); LightGBM.fit!(cls, X_tr, y_lbl; verbosity = -1); t_fit += time() - tf
                ts2 = time(); X_te = gather_rows!(bufs.test, X_all, test_idx); t_slice += time() - ts2
                tp = time(); raw = LightGBM.predict(cls, X_te); t_predict += time() - tp
                fold_scores[fi] = ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
                fold_predictors[fi] = (
                    kind = :lgbm,
                    value = NaN,
                    model = cls,
                    beta = nothing,
                )
                if compute_infold
                    ts3 = time(); X_tr_full = gather_rows!(bufs.infold, X_all, full_train_idx); t_slice += time() - ts3
                    tp2 = time(); raw_in = LightGBM.predict(cls, X_tr_full); t_predict += time() - tp2
                    infold_scores[fi] = ndims(raw_in) == 2 ? dropdims(raw_in; dims=2) : raw_in
                end
                last_cls = cls
            end
        end
        return (
            fold_scores,
            infold_scores,
            fold_predictors,
            last_cls,
            (slice = t_slice, fit = t_fit, predict = t_predict),
        )
    end

    # Probit CV — only runs in low-data branch, so build fold DataFrames lazily here.
    function _probit_cv()
        # Build Float64 fold DataFrames with intercept (ProbitRegression's z_score_bounds is typed Float64).
        function _mk_df(idx)
            df = DataFrame(Float64.(X_all[idx, :]), :auto)
            df[!, :intercept] = ones(Float64, nrow(df))
            df
        end
        fold_scores  = Vector{Vector{Float64}}(undef, 2)
        infold_scores = compute_infold ? Vector{Vector{Float64}}(undef, 2) : nothing
        fold_predictors = Vector{Any}(undef, 2)
        for (fi, (fit_idx, full_train_idx, test_idx)) in enumerate(fold_pairs)
            sub_pos = sub_positions[fi]
            train_sub_idx = fit_idx[sub_pos]
            tr = _mk_df(train_sub_idx)
            df_te = _mk_df(test_idx)
            y_bool = Vector{Bool}(targets_col[train_sub_idx])
            if isempty(y_bool) || length(unique(y_bool)) == 1
                constant_score = isempty(y_bool) || !y_bool[1] ? 0.0 : 1.0
                fold_scores[fi] = fill(constant_score, length(test_idx))
                fold_predictors[fi] = (
                    kind = :constant,
                    value = constant_score,
                    model = nothing,
                    beta = nothing,
                )
                if compute_infold
                    infold_scores[fi] = fill(constant_score, length(full_train_idx))
                end
                continue
            end
            β = zeros(Float64, ncol(tr))
            tr_chunks = Iterators.partition(1:length(y_bool), max(1, cld(length(y_bool), Threads.nthreads())))
            ProbitRegression(β, tr, y_bool, tr_chunks; max_iter=15)
            fold_predictors[fi] = (
                kind = :probit,
                value = NaN,
                model = nothing,
                beta = copy(β),
            )
            s = zeros(Float64, nrow(df_te))
            te_chunks = Iterators.partition(1:nrow(df_te), max(1, cld(nrow(df_te), Threads.nthreads())))
            ModelPredictProbs!(s, df_te, β, te_chunks)
            fold_scores[fi] = s
            if compute_infold
                df_tr_full = _mk_df(full_train_idx)
                s_in = zeros(Float64, nrow(df_tr_full))
                tr_full_chunks = Iterators.partition(1:nrow(df_tr_full),
                                                     max(1, cld(nrow(df_tr_full), Threads.nthreads())))
                ModelPredictProbs!(s_in, df_tr_full, β, tr_full_chunks)
                infold_scores[fi] = s_in
            end
        end
        return fold_scores, infold_scores, fold_predictors
    end

    function _fold_oof_psms_at_1pct(y_fold::AbstractVector, scores_on_test::AbstractVector)
        probs = Float32.(scores_on_test)
        qvals = zeros(Float16, length(probs))
        get_qvalues!(probs, y_fold, qvals)
        count(i -> qvals[i] <= Float16(0.01) && y_fold[i], eachindex(qvals))
    end

    @debug_l1 "  LightGBM CV: fold0=$(length(idx0)) fold1=$(length(idx1)) PSMs; " *
               "fit_folds=$([length(fit_idx) for (fit_idx, _, _) in fold_pairs]) " *
               "train=$(length.(sub_positions))  MAX_TRAIN=$MAX_TRAIN"

    # Always run the default-HP LGBM (single source of truth for big datasets
    # and the safest fallback for small ones).
    lgbm_scores, lgbm_infold_scores, lgbm_predictors, last_classifier, lgbm_timings = _lgbm_cv(buffers)
    @debug_l1 "  LightGBM timings: slice=$(round(lgbm_timings.slice, digits=2))s " *
               "fit=$(round(lgbm_timings.fit, digits=2))s " *
               "predict=$(round(lgbm_timings.predict, digits=2))s"

    # Helper to score a candidate by OOF targets at q≤0.01 across both folds.
    _oof_count(fs) = _fold_oof_psms_at_1pct(targets_col[idx0], fs[1]) +
                     _fold_oof_psms_at_1pct(targets_col[idx1], fs[2])

    # Adaptive model selection: if the default LGBM's OOF target count at q≤0.01
    # is below ADAPTIVE_HP_TARGET_THRESHOLD (5_000), the dataset is small enough
    # that competing alternative HP configs (medium=md50/leaves31,
    # simple=md20/leaves15) might learn better despite the default's complexity
    # advantage on bigger data. Below LOW_DATA_THRESHOLD per fold, probit also
    # joins the competition. Winner = highest OOF target count at 1% FDR.
    #
    # Perf shortcut: on large datasets the OOF count at q≤0.01 will be tens of
    # thousands of targets, well above the 5_000 adaptive threshold — adaptive
    # will always be false.
    skip_oof_probe = min_fit_size > 50_000
    if skip_oof_probe
        n_lgbm_default = -1   # sentinel: not computed
        adaptive = false
    else
        n_lgbm_default = _oof_count(lgbm_scores)
        adaptive = n_lgbm_default < ADAPTIVE_HP_TARGET_THRESHOLD
    end
    candidate_type = @NamedTuple{
        name::String,
        scores::Vector{Vector{Float64}},
        infold::Union{Nothing, Vector{Vector{Float64}}},
        last::Union{Nothing, LightGBM.LGBMClassification},
        predictors::Vector{Any},
        oof::Int64,
    }
    candidates = candidate_type[
        (name="lgbm_default", scores=lgbm_scores, infold=lgbm_infold_scores,
         last=last_classifier, predictors=lgbm_predictors, oof=n_lgbm_default)
    ]

    if adaptive
        for (variant_name, hp_over) in pairs(ADAPTIVE_HP_OVERRIDES)
            sc, infold, preds, lc, _ = _lgbm_cv(buffers, hp_over)
            push!(candidates, (name="lgbm_$(variant_name)", scores=sc, infold=infold,
                               last=lc, predictors=preds, oof=_oof_count(sc)))
        end
    end

    info_n_probit = -1
    if low_data
        probit_scores, probit_infold_scores, probit_predictors = _probit_cv()
        n_probit = _oof_count(probit_scores)
        info_n_probit = n_probit
        push!(candidates, (name="probit", scores=probit_scores, infold=probit_infold_scores,
                           last=nothing, predictors=probit_predictors, oof=n_probit))
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
        predictor = (
            available_features = available_features,
            fold_predictors = winner.predictors,
        ),
    )
    return all_scores, infold_all, last_classifier, info
end

function _predict_psm_classifier_fold(fold_predictor, X::AbstractMatrix)
    n = size(X, 1)
    if fold_predictor.kind == :constant
        return fill(Float64(fold_predictor.value), n)
    elseif fold_predictor.kind == :lgbm
        raw = LightGBM.predict(fold_predictor.model, X)
        return ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
    elseif fold_predictor.kind == :probit
        df = DataFrame(Float64.(X), :auto)
        df[!, :intercept] = ones(Float64, nrow(df))
        scores = zeros(Float64, nrow(df))
        chunks = Iterators.partition(1:nrow(df), max(1, cld(nrow(df), Threads.nthreads())))
        ModelPredictProbs!(scores, df, fold_predictor.beta, chunks)
        return scores
    else
        throw(ArgumentError("unknown PSM classifier predictor kind $(fold_predictor.kind)"))
    end
end

function predict_psm_classifier_scores(
    psms::DataFrame,
    predictor;
    buffers::LGBMMatrixBuffers = LGBMMatrixBuffers(),
)
    # As in train_psm_classifier_with_fallback: the matrices below are
    # `unsafe_wrap` views into these vectors, so root them for the whole call.
    b_all, b_test = buffers.all, buffers.test
    return GC.@preserve b_all b_test _predict_psm_classifier_scores(psms, predictor, buffers)
end

function _predict_psm_classifier_scores(
    psms::DataFrame,
    predictor,
    buffers::LGBMMatrixBuffers,
)
    available_features = Vector{Symbol}(predictor.available_features)
    X_all = feature_matrix!(buffers.all, psms, available_features)
    cv_fold = psms[!, :cv_fold]
    test_indices = [findall(cv_fold .== 0), findall(cv_fold .== 1)]
    scores = Vector{Float64}(undef, nrow(psms))
    # Both folds gather into the same buffer, pre-sized to the larger of the two
    # so the loop never resizes it. Each fold's matrix is fully consumed by the
    # predict below before the next iteration overwrites it.
    _size_matrix_buffer!(
        buffers.test,
        max(length(test_indices[1]), length(test_indices[2])),
        length(available_features),
    )
    for fi in 1:2
        idx = test_indices[fi]
        isempty(idx) && continue
        scores[idx] .= _predict_psm_classifier_fold(
            predictor.fold_predictors[fi],
            gather_rows!(buffers.test, X_all, idx),
        )
    end
    return scores
end

"""
    train_lgbm_and_select_best(psms; features) -> (best_psms, scores, timings, predictor)

Train LightGBM (with low-data probit fallback) on ALL PSMs (all scans) using
the shared `train_psm_classifier_with_fallback` helper, then select the best
scan per precursor by LightGBM score and log feature importances.

Returns:
- best_psms: DataFrame with one row per precursor (best by LightGBM score)
- scores: Vector{Float32} of LightGBM probabilities for best_psms
- timings: NamedTuple with timing breakdowns
- predictor: fold predictor from the fitted classifier; after iRT refinement,
  refreshed candidate features are scored with this same model before
  best-scan selection is repeated
"""
function train_lgbm_and_select_best(
    psms::DataFrame;
    features::Vector{Symbol} = collect(PRESCORE_FEATURES),
    lgbm_hp = MAINSEARCH_LGBM_HP,
    center_mzs = nothing,
    isolation_widths = nothing,
    buffers::LGBMMatrixBuffers = LGBMMatrixBuffers(),
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
    all_scores, _, last_classifier, info = train_psm_classifier_with_fallback(
        psms;
        features = features,
        lgbm_hp = lgbm_hp,
        buffers = buffers,
    )
    t_train_cv = time()

    # Surface the adaptive model selection result so we can see per-file picks
    # in the log without grepping debug_l1.
    if hasproperty(info, :adaptive) && info.adaptive
        oof_str = join(["$k=$v" for (k,v) in sort(collect(info.candidate_oof))], " ")
        @debug_l1 "  Adaptive model selection (n=$(nrow(psms))): $oof_str → $(info.winner)"
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
    psms = select_best_per_precursor!(
        psms,
        :lgbm_score;
        center_mzs = center_mzs,
        isolation_widths = isolation_widths,
    )

    # Extract scores for best PSMs
    scores = psms[!, :lgbm_score]
    t_best = time()

    # Feature importances — diagnostic only.
    let imp = importance(model)
        if imp !== nothing
            sorted_imp = sort(imp, by = x -> -x[2])
            lines = ["MainSearch per-file LGBM feature gains (all $(length(sorted_imp))):"]
            for (fname, gain) in sorted_imp
                push!(lines, "    $(rpad(string(fname), 40)) $(round(Int, gain))")
            end
            @debug_l1 join(lines, "\n")
        end
    end

    timings = (
        matrix = 0.0,            # subsumed into train_cv now
        train_cv = t_train_cv - t0,
        best = t_best - t_train_cv,
    )

    return psms, Vector{Float32}(scores), timings, info.predictor
end

function reapply_psm_classifier_and_select_best!(
    psms::DataFrame,
    predictor;
    center_mzs = nothing,
    isolation_widths = nothing,
    buffers::LGBMMatrixBuffers = LGBMMatrixBuffers(),
)
    t0 = time()
    if nrow(psms) > 0 && !hasproperty(psms, :n_scans)
        counts = Dict{UInt32, UInt32}()
        @inbounds for pid in psms[!, :precursor_idx]::Vector{UInt32}
            counts[pid] = get(counts, pid, UInt32(0)) + UInt32(1)
        end
        psms[!, :n_scans] = UInt32[counts[pid] for pid in psms[!, :precursor_idx]]
    end

    all_scores = predict_psm_classifier_scores(psms, predictor; buffers = buffers)
    t_predict = time()
    psms[!, :lgbm_score] = Float32.(all_scores)
    best_psms = select_best_per_precursor!(
        psms,
        :lgbm_score;
        center_mzs = center_mzs,
        isolation_widths = isolation_widths,
    )
    scores = best_psms[!, :lgbm_score]
    t_best = time()

    timings = (
        predict = t_predict - t0,
        best = t_best - t_predict,
    )
    return best_psms, Vector{Float32}(scores), timings
end

function _mainsearch_peps_and_pass_mask(
    scores::AbstractVector{<:Real},
    targets::AbstractVector{Bool};
    pep_threshold::Float32 = MAIN_PEP_FILTER_THR,
)
    n = length(scores)
    if pep_threshold >= 1.0f0
        return zeros(Float32, n), trues(n)
    end
    score_f32 = scores isa AbstractVector{Float32} ? scores : Float32.(scores)
    peps = Vector{Float32}(undef, n)
    get_PEP!(score_f32, targets, peps; doSort = true, fdr_scale_factor = 1.0f0)
    return peps, peps .<= pep_threshold
end

function add_precursor_fraction_transmitted!(
    best_psms::DataFrame,
    quad_transmission_model::QuadTransmissionModel,
    search_data::Vector{SimpleLibrarySearch{IsotopeSplineModel{Float32}}},
    prec_charge::AbstractArray{UInt8},
    prec_mz::AbstractArray{Float32},
    sulfur_count::AbstractArray{UInt8},
    centerMz::AbstractVector{Union{Missing, Float32}},
    isolationWidthMz::AbstractVector{Union{Missing, Float32}},
)
    prec_ids = best_psms[!, :precursor_idx]::Vector{UInt32}
    scan_idxs = best_psms[!, :scan_idx]::Vector{UInt32}
    n = nrow(best_psms)
    out_fraction = Vector{Float32}(undef, n)
    iso_splines = getIsoSplines(search_data[1])
    precursor_transmission = zeros(Float32, 5)

    @inbounds for i in 1:n
        pid = prec_ids[i]
        scan_id = scan_idxs[i]
        _, _, scan_mz, window_width = _compute_scan_window(scan_id, centerMz, isolationWidthMz)
        out_fraction[i] = _compute_fraction_transmitted(
            iso_splines,
            quad_transmission_model,
            prec_mz[pid],
            prec_charge[pid],
            sulfur_count[pid],
            scan_mz,
            window_width,
            precursor_transmission,
        )
    end

    best_psms[!, :precursor_fraction_transmitted] = out_fraction
    return best_psms
end

const MAINSEARCH_FRAGMENT_INTENSITY_COLUMNS = (
    :frag1_int, :frag2_int, :frag3_int, :frag4_int,
    :frag5_int, :frag6_int, :frag7_int, :frag8_int,
)

const OTHER_WINDOW_CORR_SENTINEL = -1.0f0
const OTHER_WINDOW_APEX_DELTA_SENTINEL = 100.0f0

mutable struct MainSearchOtherWindowScratch
    selected_cycles::Vector{UInt32}
    selected_weight_values::Vector{Float32}
    window_to_idx::Dict{Tuple{Float32, Float32}, Int}
    window_best_peps::Vector{Float32}
    window_cycles::Vector{Vector{UInt32}}
    window_weight_values::Vector{Vector{Float32}}
    window_apex_weight::Vector{Float32}
    window_apex_irt::Vector{Float32}
    window_smooth_prev_row::Vector{Int}
    window_smooth_same_row::Vector{Int}
    window_smooth_next_row::Vector{Int}
end

function MainSearchOtherWindowScratch()
    return MainSearchOtherWindowScratch(
        UInt32[],
        Float32[],
        Dict{Tuple{Float32, Float32}, Int}(),
        Float32[],
        Vector{UInt32}[],
        Vector{Float32}[],
        Float32[],
        Float32[],
        Int[],
        Int[],
        Int[],
    )
end

function _reset_other_window_scratch!(scratch::MainSearchOtherWindowScratch)
    empty!(scratch.selected_cycles)
    empty!(scratch.selected_weight_values)
    empty!(scratch.window_to_idx)
    empty!(scratch.window_best_peps)
    @inbounds for cycles in scratch.window_cycles
        empty!(cycles)
    end
    @inbounds for values in scratch.window_weight_values
        empty!(values)
    end
    empty!(scratch.window_apex_weight)
    empty!(scratch.window_apex_irt)
    empty!(scratch.window_smooth_prev_row)
    empty!(scratch.window_smooth_same_row)
    empty!(scratch.window_smooth_next_row)
    return nothing
end

function _other_window_idx!(
    scratch::MainSearchOtherWindowScratch,
    key::Tuple{Float32, Float32},
)
    idx = get(scratch.window_to_idx, key, 0)
    idx != 0 && return idx

    idx = length(scratch.window_best_peps) + 1
    scratch.window_to_idx[key] = idx
    push!(scratch.window_best_peps, Inf32)
    if idx > length(scratch.window_cycles)
        push!(scratch.window_cycles, UInt32[])
        push!(scratch.window_weight_values, Float32[])
    end
    push!(scratch.window_apex_weight, typemin(Float32))
    push!(scratch.window_apex_irt, NaN32)
    push!(scratch.window_smooth_prev_row, 0)
    push!(scratch.window_smooth_same_row, 0)
    push!(scratch.window_smooth_next_row, 0)
    return idx
end

@inline function _store_other_window_smooth_row!(
    scratch::MainSearchOtherWindowScratch,
    other_idx::Int,
    cycle::UInt32,
    center_cycle::UInt32,
    row::Int,
)
    if center_cycle > UInt32(1) && cycle == center_cycle - UInt32(1)
        scratch.window_smooth_prev_row[other_idx] = row
    elseif cycle == center_cycle
        scratch.window_smooth_same_row[other_idx] = row
    elseif cycle == center_cycle + UInt32(1)
        scratch.window_smooth_next_row[other_idx] = row
    end
    return nothing
end

@inline function _trace_push_cycle_sum!(
    cycles::Vector{UInt32},
    values::Vector{Float32},
    cycle::UInt32,
    value::Float32,
)
    len = length(cycles)
    if len == 0 || cycle > cycles[end]
        push!(cycles, cycle)
        push!(values, value)
        return nothing
    elseif cycle == cycles[end]
        values[end] += value
        return nothing
    end

    i = 1
    @inbounds while i <= len && cycles[i] < cycle
        i += 1
    end
    if i <= len && cycles[i] == cycle
        values[i] += value
    else
        insert!(cycles, i, cycle)
        insert!(values, i, value)
    end
    return nothing
end

function _trace_aligned_corr_sorted(
    cycles_a::Vector{UInt32},
    values_a::Vector{Float32},
    cycles_b::Vector{UInt32},
    values_b::Vector{Float32},
)
    n_cycles_a = length(cycles_a)
    n_cycles_b = length(cycles_b)
    n_aligned_cycles = 0
    mean_a = 0f0
    mean_b = 0f0

    index_a = 1
    index_b = 1
    @inbounds while index_a <= n_cycles_a || index_b <= n_cycles_b
        if index_b > n_cycles_b || (index_a <= n_cycles_a && cycles_a[index_a] < cycles_b[index_b])
            mean_a += values_a[index_a]
            index_a += 1
        elseif index_a > n_cycles_a || cycles_b[index_b] < cycles_a[index_a]
            mean_b += values_b[index_b]
            index_b += 1
        else
            mean_a += values_a[index_a]
            mean_b += values_b[index_b]
            index_a += 1
            index_b += 1
        end
        n_aligned_cycles += 1
    end

    n_aligned_cycles < 2 && return 0f0
    inv_aligned_cycles = 1f0 / Float32(n_aligned_cycles)
    mean_a *= inv_aligned_cycles
    mean_b *= inv_aligned_cycles

    sum_sq_a = 0f0
    sum_sq_b = 0f0
    sum_cross = 0f0
    index_a = 1
    index_b = 1
    @inbounds while index_a <= n_cycles_a || index_b <= n_cycles_b
        value_a = 0f0
        value_b = 0f0
        if index_b > n_cycles_b || (index_a <= n_cycles_a && cycles_a[index_a] < cycles_b[index_b])
            value_a = values_a[index_a]
            index_a += 1
        elseif index_a > n_cycles_a || cycles_b[index_b] < cycles_a[index_a]
            value_b = values_b[index_b]
            index_b += 1
        else
            value_a = values_a[index_a]
            value_b = values_b[index_b]
            index_a += 1
            index_b += 1
        end
        delta_a = value_a - mean_a
        delta_b = value_b - mean_b
        sum_sq_a += delta_a * delta_a
        sum_sq_b += delta_b * delta_b
        sum_cross += delta_a * delta_b
    end

    denom = sqrt(sum_sq_a * sum_sq_b)
    return denom > 0f0 ? Float32(sum_cross / denom) : 0f0
end

function _mainsearch_flanking_core_bounds!(
    order::Vector{Int},
    group_rows::Vector{Int},
    group_len::Int,
    scan_idxs::AbstractVector{UInt32},
    cycle_idxs::AbstractVector,
    pass_mask::AbstractVector{Bool},
    selected_scan::UInt32,
)
    if length(order) < group_len
        resize!(order, group_len)
    end

    best_row = group_rows[1]
    @inbounds for i in 1:group_len
        row = group_rows[i]
        order[i] = row
        if scan_idxs[row] == selected_scan
            best_row = row
        end
    end

    best_scan = scan_idxs[best_row]
    pass_mask[best_row] || return (
        scan_min = best_scan,
        scan_max = best_scan,
        n_scans = UInt16(1),
        best_row = best_row,
        left_row = 0,
        right_row = 0,
    )

    ord = view(order, 1:group_len)
    sort!(
        ord;
        lt = (a, b) -> begin
            ca = UInt32(cycle_idxs[a])
            cb = UInt32(cycle_idxs[b])
            ca == cb ? scan_idxs[a] < scan_idxs[b] : ca < cb
        end,
    )

    best_pos = 1
    @inbounds for i in 1:group_len
        if ord[i] == best_row
            best_pos = i
            break
        end
    end

    lo = best_pos
    @inbounds while lo > 1
        current = ord[lo]
        previous = ord[lo - 1]
        pass_mask[previous] || break
        UInt32(cycle_idxs[previous]) == UInt32(cycle_idxs[current]) - UInt32(1) || break
        lo -= 1
    end

    hi = best_pos
    @inbounds while hi < group_len
        current = ord[hi]
        next_row = ord[hi + 1]
        pass_mask[next_row] || break
        UInt32(cycle_idxs[next_row]) == UInt32(cycle_idxs[current]) + UInt32(1) || break
        hi += 1
    end

    scan_min = typemax(UInt32)
    scan_max = UInt32(0)
    @inbounds for i in lo:hi
        scan = scan_idxs[ord[i]]
        scan_min = min(scan_min, scan)
        scan_max = max(scan_max, scan)
    end
    n_scans = UInt16(min(hi - lo + 1, Int(typemax(UInt16))))
    return (
        scan_min = scan_min,
        scan_max = scan_max,
        n_scans = n_scans,
        best_row = best_row,
        left_row = best_pos > lo ? ord[best_pos - 1] : 0,
        right_row = best_pos < hi ? ord[best_pos + 1] : 0,
    )
end

@inline function _mainsearch_fragment_presence_mask(row::Int, frag_cols)
    mask = UInt16(0)
    bit = UInt16(1)
    @inbounds for col in frag_cols
        if Float32(col[row]) > 0f0
            mask |= bit
        end
        bit <<= 1
    end
    return mask
end

@inline function _mainsearch_fragment_intensity_sum(row::Int, frag_cols)
    total = 0f0
    @inbounds for col in frag_cols
        total += max(Float32(col[row]), 0f0)
    end
    return total
end

function _mainsearch_radius1_smooth_fragments!(
    out_cols,
    out_row::Int,
    apex_row::Int,
    left_row::Int,
    right_row::Int,
    frag_cols,
)
    denom = 2f0
    left_row > 0 && (denom += 1f0)
    right_row > 0 && (denom += 1f0)

    @inbounds for rank in 1:8
        signal = 2f0 * max(Float32(frag_cols[rank][apex_row]), 0f0)
        left_row > 0 && (signal += max(Float32(frag_cols[rank][left_row]), 0f0))
        right_row > 0 && (signal += max(Float32(frag_cols[rank][right_row]), 0f0))
        out_cols[rank][out_row] = Float32(signal / denom)
    end
    return nothing
end

function _mainsearch_radius1_2d_shadow_hellinger(
    apex_row::Int,
    left_row::Int,
    right_row::Int,
    other_prev_row::Int,
    other_same_row::Int,
    other_next_row::Int,
    fitted_cols,
    shadow_cols,
)
    smooth_denom = 2f0
    left_row > 0 && (smooth_denom += 1f0)
    right_row > 0 && (smooth_denom += 1f0)
    other_prev_row > 0 && (smooth_denom += 0.5f0)
    other_same_row > 0 && (smooth_denom += 1f0)
    other_next_row > 0 && (smooth_denom += 0.5f0)

    sum_fitted = 0f0
    sum_shadow = 0f0
    bc_sum = 0f0
    @inbounds for rank in 1:8
        fitted = max(Float32(fitted_cols[rank][apex_row]), 0f0)
        shadow = 2f0 * max(Float32(shadow_cols[rank][apex_row]), 0f0)
        left_row > 0 && (shadow += max(Float32(shadow_cols[rank][left_row]), 0f0))
        right_row > 0 && (shadow += max(Float32(shadow_cols[rank][right_row]), 0f0))
        other_prev_row > 0 && (shadow += 0.5f0 * max(Float32(shadow_cols[rank][other_prev_row]), 0f0))
        other_same_row > 0 && (shadow += max(Float32(shadow_cols[rank][other_same_row]), 0f0))
        other_next_row > 0 && (shadow += 0.5f0 * max(Float32(shadow_cols[rank][other_next_row]), 0f0))
        shadow = Float32(shadow / smooth_denom)
        sum_fitted += fitted
        sum_shadow += shadow
        bc_sum += sqrt(fitted * shadow)
    end

    denom = sqrt(sum_fitted * sum_shadow)
    hellinger_sq = denom > 0f0 ? 1f0 - bc_sum / denom : 1f0
    return Float32(-log2(max(hellinger_sq, 1f-10)))
end

function add_trace_and_fragment_features!(
    best_psms::DataFrame,
    psms::DataFrame,
    pass_mask::AbstractVector{Bool};
    bitvec_rank_table,
    center_mzs = nothing,
    isolation_widths = nothing,
    pep_values = nothing,
)
    best_prec_ids = best_psms[!, :precursor_idx]::Vector{UInt32}
    best_scan_idxs = best_psms[!, :scan_idx]::Vector{UInt32}

    prec_ids = psms[!, :precursor_idx]::Vector{UInt32}
    scan_idxs = psms[!, :scan_idx]::Vector{UInt32}
    cycle_idxs = psms[!, :cycle_idx]::Vector{UInt32}
    weights = psms[!, :weight]::Vector{Float32}
    irt_obs = psms[!, :irt_obs]::Vector{Float32}
    ms1_m0_intensities = psms[!, :ms1_m0_intensity]::Vector{Float32}
    # `Tuple(psms[!, c] for c in ...)` inferred as an abstractly-typed tuple, because `df[!, col]`
    # infers as `AbstractVector`. Unlike the gather helpers in features.jl -- whose result crosses a
    # function boundary and so specialises on the concrete runtime tuple -- these are consumed by the
    # per-row loops *in this same function*, so the abstract element type survived and every
    # `Float32(col[row])` boxed and dispatched dynamically: 8 per row, per fragment helper.
    # `ntuple` over the const 8-column tuple plus the assertion makes the type statically concrete.
    # The columns are Float32 by construction (MainUnscoredPSM{Float32}), consistent with the
    # `::Vector{Float32}` assertions already made above on this same table.
    frag_cols = ntuple(
        i -> psms[!, MAINSEARCH_FRAGMENT_INTENSITY_COLUMNS[i]]::Vector{Float32},
        length(MAINSEARCH_FRAGMENT_INTENSITY_COLUMNS),
    )
    fitted_frag_cols = ntuple(
        i -> psms[!, FITTED_FRAGMENT_INTENSITY_COLUMNS[i]]::Vector{Float32},
        length(FITTED_FRAGMENT_INTENSITY_COLUMNS),
    )
    shadow_frag_cols = ntuple(
        i -> psms[!, SHADOW_FRAGMENT_INTENSITY_COLUMNS[i]]::Vector{Float32},
        length(SHADOW_FRAGMENT_INTENSITY_COLUMNS),
    )

    n_best = nrow(best_psms)
    out_flanking_core_scan_min = Vector{UInt32}(undef, n_best)
    out_flanking_core_scan_max = Vector{UInt32}(undef, n_best)
    out_n_contiguous_scans = Vector{UInt16}(undef, n_best)
    out_flanking_core_ms1_m0_signal = Vector{Float32}(undef, n_best)
    out_flanking_core_frag_signal = Vector{Float32}(undef, n_best)
    out_union = Vector{UInt8}(undef, n_best)
    out_intersection = Vector{UInt8}(undef, n_best)
    out_union_rank = Vector{UInt16}(undef, n_best)
    out_intersection_rank = Vector{UInt16}(undef, n_best)
    out_smooth_frag_cols = ntuple(_ -> Vector{Float32}(undef, n_best), 8)
    out_2d_shadow_hellinger = Vector{Float32}(undef, n_best)
    out_n_scans_other_windows = Vector{UInt32}(undef, n_best)
    out_other_window_weight_corr = Vector{Float32}(undef, n_best)
    out_other_window_apex_delta_irt = Vector{Float32}(undef, n_best)

    flanking_core_order = Int[]
    selected_window_rows = Int[]
    use_window_groups = center_mzs !== nothing && isolation_widths !== nothing
    compute_other_windows = use_window_groups && pep_values !== nothing
    other_window_scratch = MainSearchOtherWindowScratch()

    row = 1
    n = nrow(psms)
    @inbounds for i in 1:n_best
        pid = best_prec_ids[i]
        while prec_ids[row] < pid
            row += 1
        end
        group_start = row
        while row <= n && prec_ids[row] == pid
            row += 1
        end
        group_stop = row - 1

        empty!(selected_window_rows)
        selected_window_key = use_window_groups ?
            _scan_window_key(best_scan_idxs[i], center_mzs, isolation_widths) :
            (0f0, 0f0)
        for group_row in group_start:group_stop
            if !use_window_groups ||
               _scan_window_key(scan_idxs[group_row], center_mzs, isolation_widths) == selected_window_key
                push!(selected_window_rows, group_row)
            end
        end

        flanking_core = _mainsearch_flanking_core_bounds!(
            flanking_core_order,
            selected_window_rows,
            length(selected_window_rows),
            scan_idxs,
            cycle_idxs,
            pass_mask,
            best_scan_idxs[i],
        )
        out_flanking_core_scan_min[i] = flanking_core.scan_min
        out_flanking_core_scan_max[i] = flanking_core.scan_max
        out_n_contiguous_scans[i] = flanking_core.n_scans

        _mainsearch_radius1_smooth_fragments!(
            out_smooth_frag_cols,
            i,
            flanking_core.best_row,
            flanking_core.left_row,
            flanking_core.right_row,
            frag_cols,
        )
        smoothed_2d_shadow_hellinger = _mainsearch_radius1_2d_shadow_hellinger(
            flanking_core.best_row,
            flanking_core.left_row,
            flanking_core.right_row,
            0,
            0,
            0,
            fitted_frag_cols,
            shadow_frag_cols,
        )

        core_ms1_m0_signal = 0f0
        core_frag_signal = 0f0
        union_mask = UInt16(0)
        intersection_mask = UInt16(0x00ff)

        for group_row in selected_window_rows
            mask = _mainsearch_fragment_presence_mask(group_row, frag_cols)
            union_mask |= mask
            intersection_mask &= mask
            if flanking_core.scan_min <= scan_idxs[group_row] <= flanking_core.scan_max
                core_ms1_m0_signal += max(ms1_m0_intensities[group_row], 0f0)
                core_frag_signal += _mainsearch_fragment_intensity_sum(group_row, frag_cols)
            end
        end

        out_flanking_core_ms1_m0_signal[i] = core_ms1_m0_signal
        out_flanking_core_frag_signal[i] = core_frag_signal
        out_union[i] = UInt8(count_ones(union_mask))
        out_intersection[i] = UInt8(count_ones(intersection_mask))
        out_union_rank[i] = _bitvec_pattern_rank(bitvec_rank_table, union_mask)
        out_intersection_rank[i] = _bitvec_pattern_rank(bitvec_rank_table, intersection_mask)

        n_scans_other_windows = UInt32(0)
        other_window_weight_corr = OTHER_WINDOW_CORR_SENTINEL
        other_window_apex_delta_irt = OTHER_WINDOW_APEX_DELTA_SENTINEL

        if compute_other_windows
            _reset_other_window_scratch!(other_window_scratch)
            selected_apex_weight = typemin(Float32)
            selected_apex_irt = NaN32
            center_cycle = UInt32(cycle_idxs[flanking_core.best_row])

            for group_row in group_start:group_stop
                pass_mask[group_row] || continue
                w = max(weights[group_row], 0f0)
                cycle = UInt32(cycle_idxs[group_row])
                group_window_key = _scan_window_key(scan_idxs[group_row], center_mzs, isolation_widths)

                if group_window_key == selected_window_key
                    _trace_push_cycle_sum!(
                        other_window_scratch.selected_cycles,
                        other_window_scratch.selected_weight_values,
                        cycle,
                        w,
                    )
                    if w > selected_apex_weight
                        selected_apex_weight = w
                        selected_apex_irt = irt_obs[group_row]
                    end
                else
                    n_scans_other_windows += UInt32(1)
                    other_idx = _other_window_idx!(other_window_scratch, group_window_key)
                    other_window_scratch.window_best_peps[other_idx] =
                        min(other_window_scratch.window_best_peps[other_idx], Float32(pep_values[group_row]))
                    _store_other_window_smooth_row!(
                        other_window_scratch,
                        other_idx,
                        cycle,
                        center_cycle,
                        group_row,
                    )
                    _trace_push_cycle_sum!(
                        other_window_scratch.window_cycles[other_idx],
                        other_window_scratch.window_weight_values[other_idx],
                        cycle,
                        w,
                    )
                    if w > other_window_scratch.window_apex_weight[other_idx]
                        other_window_scratch.window_apex_weight[other_idx] = w
                        other_window_scratch.window_apex_irt[other_idx] = irt_obs[group_row]
                    end
                end
            end

            if n_scans_other_windows > UInt32(0) &&
               !isempty(other_window_scratch.selected_cycles) &&
               isfinite(selected_apex_irt)
                chosen_other_idx = 0
                chosen_other_pep = Inf32
                @inbounds for other_idx in eachindex(other_window_scratch.window_best_peps)
                    pep = other_window_scratch.window_best_peps[other_idx]
                    if pep < chosen_other_pep
                        chosen_other_pep = pep
                        chosen_other_idx = other_idx
                    end
                end

                if chosen_other_idx != 0
                    other_window_weight_corr = _trace_aligned_corr_sorted(
                        other_window_scratch.selected_cycles,
                        other_window_scratch.selected_weight_values,
                        other_window_scratch.window_cycles[chosen_other_idx],
                        other_window_scratch.window_weight_values[chosen_other_idx],
                    )
                    other_prev_row = other_window_scratch.window_smooth_prev_row[chosen_other_idx]
                    other_same_row = other_window_scratch.window_smooth_same_row[chosen_other_idx]
                    other_next_row = other_window_scratch.window_smooth_next_row[chosen_other_idx]
                    if other_prev_row != 0 || other_same_row != 0 || other_next_row != 0
                        smoothed_2d_shadow_hellinger = _mainsearch_radius1_2d_shadow_hellinger(
                            flanking_core.best_row,
                            flanking_core.left_row,
                            flanking_core.right_row,
                            other_prev_row,
                            other_same_row,
                            other_next_row,
                            fitted_frag_cols,
                            shadow_frag_cols,
                        )
                    end
                    other_apex_irt = other_window_scratch.window_apex_irt[chosen_other_idx]
                    if isfinite(other_apex_irt)
                        other_window_apex_delta_irt = abs(selected_apex_irt - other_apex_irt)
                    end
                end
            end
        end

        out_2d_shadow_hellinger[i] = smoothed_2d_shadow_hellinger
        out_n_scans_other_windows[i] = n_scans_other_windows
        out_other_window_weight_corr[i] = other_window_weight_corr
        out_other_window_apex_delta_irt[i] = other_window_apex_delta_irt
    end

    best_psms[!, :flanking_core_scan_min] = out_flanking_core_scan_min
    best_psms[!, :flanking_core_scan_max] = out_flanking_core_scan_max
    best_psms[!, :n_contiguous_scans] = out_n_contiguous_scans
    best_psms[!, :flanking_core_ms1_m0_signal] = out_flanking_core_ms1_m0_signal
    best_psms[!, :flanking_core_frag_signal] = out_flanking_core_frag_signal
    best_psms[!, :n_frags_detected_union] = out_union
    best_psms[!, :n_frags_detected_intersection] = out_intersection
    best_psms[!, :n_frags_detected_union_bitvec_rank] = out_union_rank
    best_psms[!, :n_frags_detected_intersection_bitvec_rank] = out_intersection_rank
    best_psms[!, :n_scans_other_windows] = out_n_scans_other_windows
    best_psms[!, :other_window_weight_corr] = out_other_window_weight_corr
    best_psms[!, :other_window_apex_delta_irt] = out_other_window_apex_delta_irt
    @inbounds for rank in 1:8
        best_psms[!, SMOOTHED_FRAGMENT_INTENSITY_COLUMNS[rank]] = out_smooth_frag_cols[rank]
    end
    best_psms[!, :smoothed_2d_shadow_hellinger] = out_2d_shadow_hellinger
    select!(best_psms, Not([
        FITTED_FRAGMENT_INTENSITY_COLUMNS...,
        SHADOW_FRAGMENT_INTENSITY_COLUMNS...,
    ]))
    return best_psms
end

"""
    select_best_per_precursor!(psms::DataFrame, score_col::Symbol) -> DataFrame

Keeps one row per precursor_idx. Uses sortperm for contiguous group processing:
per group, selects the highest-weight PSM among those with score ≥ p75 (if ≥4 PSMs),
otherwise the highest-score PSM. Computes `irt_fwhm` (iRT span of scans with
weight ≥ 50% of peak weight), `n_above_hm`, `rt_fwhm`, and `best_rt`. When scan
window arrays are supplied, those shape features use only the selected PSM's
isolation window.
"""
function select_best_per_precursor!(
    psms::DataFrame,
    score_col::Symbol;
    center_mzs = nothing,
    isolation_widths = nothing,
)
    scores = psms[!, score_col]::Vector{Float32}
    prec_ids = psms[!, :precursor_idx]::Vector{UInt32}
    use_window_groups = center_mzs !== nothing && isolation_widths !== nothing
    scan_idxs = use_window_groups ? psms[!, :scan_idx]::Vector{UInt32} : nothing
    has_irt = hasproperty(psms, :irt_obs)
    has_weight = hasproperty(psms, :weight)
    has_rt = hasproperty(psms, :rt)
    irt_obs = has_irt ? psms[!, :irt_obs]::Vector{Float32} : nothing
    rt_vals = has_rt ? psms[!, :rt]::Vector{Float32} : nothing
    weights = (has_irt && has_weight) ? psms[!, :weight]::Vector{Float32} : nothing
    # Type-assert the Float16 scoring columns. Without `::Vector{Float16}`
    # the DataFrame accessor returns an abstract type, which forces dynamic
    # dispatch inside the sub-pass-1 hot loop and costs ~3.5s/file on
    # Astral. Identified 2026-05-19.
    n = nrow(psms)

    # sortperm groups PSMs by precursor_idx for contiguous processing
    perm = sortperm(prec_ids)

    # Output arrays — one element per unique precursor
    keep_rows = Vector{Int}()
    sizehint!(keep_rows, n ÷ 10)
    compute_fwhm = weights !== nothing  # implies has_irt
    compute_rt = has_rt

    out_irt_fwhm = compute_fwhm ? Vector{Float32}() : nothing
    # n_above_hm + rt_fwhm not LGBM features but still consumed by
    # prescore_aggregation.jl, so keep them as outputs.
    out_n_above_hm = compute_fwhm ? Vector{UInt16}() : nothing
    out_rt_fwhm = (compute_fwhm && compute_rt) ? Vector{Float32}() : nothing
    out_best_rt = compute_rt ? Vector{Float32}() : nothing
    # Re-added 2026-05-11 (orphaned in 2025-03 consolidation):
    # `smoothness` = Σ(((Δw_left + Δw_right)/w_apex)²) — squared second-derivative
    # of weight chromatogram; real peaks → smooth → low value, noise → jagged → high.
    # `num_scans` = group_len (count of PSMs / MS2 scans for this precursor).
    out_smoothness = (compute_fwhm && compute_rt) ? Vector{Float32}() : nothing
    # 2026-05-21: max_* / min_* / n_above_hm / num_scans / rt_fwhm outputs
    # removed (compute+write cost wasn't worth the ~+1% ID gain).

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

    # Reusable buffers
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

        # --- Sub-pass 1: find best-score row and max weight (max weight is
        # used by sub-pass 2 for the FWHM half-max threshold, even though
        # it's not output as a feature). ---
        best_s = typemin(Float32)
        best_row = perm[group_start]
        mw = 0f0

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
        end

        # --- Sub-pass 1b (Method 1: score-margin 0.1): among scans scoring
        # within 0.1 of the max score, pick the highest-weight one. Falls back
        # to the max-score row (best_row) when only that row qualifies. ---
        if weights !== nothing
            score_threshold = best_s - 0.1f0
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

        selected_window_key = use_window_groups ?
            _scan_window_key(scan_idxs[best_row], center_mzs, isolation_widths) :
            (0f0, 0f0)
        window_mw = mw
        if use_window_groups && weights !== nothing
            window_mw = 0f0
            @inbounds for k in 0:(group_len - 1)
                row = perm[group_start + k]
                _scan_window_key(scan_idxs[row], center_mzs, isolation_widths) == selected_window_key || continue
                w = weights[row]
                if w > window_mw
                    window_mw = w
                end
            end
        end

        push!(keep_rows, best_row)

        # --- Sub-pass 2: FWHM bounds. irt_fwhm is a LGBM feature;
        # n_above_hm + rt_fwhm are not LGBM features but are still consumed
        # by prescore_aggregation.jl. ---
        if compute_fwhm
            half_max = 0.5f0 * window_mw
            irt_lo = typemax(Float32)
            irt_hi = typemin(Float32)
            rt_lo = typemax(Float32)
            rt_hi = typemin(Float32)
            n_hm = UInt16(0)

            @inbounds for k in 0:(group_len - 1)
                row = perm[group_start + k]
                if use_window_groups
                    _scan_window_key(scan_idxs[row], center_mzs, isolation_widths) == selected_window_key || continue
                end
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
            smooth_len = group_len
            if use_window_groups
                smooth_len = 0
                @inbounds for k in 0:(group_len - 1)
                    row = perm[group_start + k]
                    _scan_window_key(scan_idxs[row], center_mzs, isolation_widths) == selected_window_key || continue
                    smooth_len += 1
                end
            end
            if length(smooth_w_buf) < smooth_len
                resize!(smooth_w_buf,  smooth_len)
                resize!(smooth_rt_buf, smooth_len)
                resize!(smooth_ord_buf, smooth_len)
            end
            @inbounds let j = 0
                for k in 0:(group_len - 1)
                    row = perm[group_start + k]
                    if use_window_groups
                        _scan_window_key(scan_idxs[row], center_mzs, isolation_widths) == selected_window_key || continue
                    end
                    j += 1
                    smooth_w_buf[j] = weights[row]
                    smooth_rt_buf[j] = rt_vals[row]
                    smooth_ord_buf[j] = j
                end
            end
            sort!(view(smooth_ord_buf, 1:smooth_len), by = ki -> smooth_rt_buf[ki])
            # Compute the roughness sum
            rough = 0f0
            if window_mw > 0f0
                if smooth_len == 1
                    # Single point — original code: rough = (−2w/w_apex)² = 4
                    rough = 4f0
                else
                    @inbounds for k in 1:smooth_len
                        ki = smooth_ord_buf[k]
                        w_i = smooth_w_buf[ki]
                        if k == 1
                            ki_r = smooth_ord_buf[k+1]
                            dt_r = smooth_rt_buf[ki_r] - smooth_rt_buf[ki]
                            d_r = dt_r > 0 ? (smooth_w_buf[ki_r] - w_i) / dt_r : 0f0
                            d_l = dt_r > 0 ? (-w_i) / dt_r : 0f0
                            rough += ((d_l + d_r) / window_mw)^2
                        elseif k == smooth_len
                            ki_l = smooth_ord_buf[k-1]
                            dt_l = smooth_rt_buf[ki] - smooth_rt_buf[ki_l]
                            d_l = dt_l > 0 ? (smooth_w_buf[ki_l] - w_i) / dt_l : 0f0
                            d_r = dt_l > 0 ? (-w_i) / dt_l : 0f0
                            rough += ((d_l + d_r) / window_mw)^2
                        else
                            ki_l = smooth_ord_buf[k-1]
                            ki_r = smooth_ord_buf[k+1]
                            dt_l = smooth_rt_buf[ki] - smooth_rt_buf[ki_l]
                            dt_r = smooth_rt_buf[ki_r] - smooth_rt_buf[ki]
                            d_l = dt_l > 0 ? (smooth_w_buf[ki_l] - w_i) / dt_l : 0f0
                            d_r = dt_r > 0 ? (smooth_w_buf[ki_r] - w_i) / dt_r : 0f0
                            rough += ((d_l + d_r) / window_mw)^2
                        end
                    end
                end
            end
            push!(out_smoothness, rough)
        end

        if out_best_rt !== nothing
            push!(out_best_rt, rt_vals[best_row])
        end
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

    return result
end
