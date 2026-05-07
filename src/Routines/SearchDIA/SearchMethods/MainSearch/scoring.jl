# LightGBM training and best-per-precursor selection. PrecursorScoringSearch
# still uses `train_psm_classifier_with_fallback`; MainSearch's per-run
# prescore model uses the semi-supervised loop below.

# Shared LightGBM hyperparameters — used by both MainSearch (per-file) and
# PrecursorScoringSearch (experiment-wide). Config "D" (depth=8, 50 iter, lr=0.10)
# was chosen via OOF 1%-FDR sweep over model complexity × training size.
const SHARED_LGBM_HP = (num_iterations=50, learning_rate=0.10, max_depth=8,
                        num_leaves=63, min_data_in_leaf=30, feature_fraction=0.8,
                        bagging_fraction=0.8, bagging_freq=1, is_unbalance=false,
                        max_bin=255)
# Per-fold training cap; folds larger than this are random-subsampled.
const SHARED_LGBM_MAX_TRAIN = 250_000
# Low-data threshold per fold: below this we CV-select between LightGBM and probit.
const SHARED_LGBM_LOW_DATA_THRESHOLD = 10_000
const MAINSEARCH_SEMISUPERVISED_ITERATIONS = 1
const MAINSEARCH_SEMISUPERVISED_POSITIVE_QVALUE = 0.10f0
const MAINSEARCH_FDR_LOG_THRESHOLDS = Float32[0.01, 0.05, 0.10]

"""
    _mainsearch_training_mask(targets, previous_scores)

Return rows used for the next MainSearch semi-supervised LightGBM iteration.
The first iteration (`previous_scores === nothing`) uses every target and decoy.
Subsequent iterations retain all decoys and only target rows passing the
configured FDR threshold under the previous iteration's out-of-fold scores.
"""
function _mainsearch_training_mask(
    targets::AbstractVector{Bool},
    previous_scores::Union{Nothing,AbstractVector},
)
    previous_scores === nothing && return trues(length(targets))

    qvals = Vector{Float32}(undef, length(targets))
    get_qvalues!(Float32.(previous_scores), targets, qvals; doSort=true)
    mask = Vector{Bool}(undef, length(targets))
    @inbounds for i in eachindex(targets)
        mask[i] = !targets[i] || qvals[i] <= MAINSEARCH_SEMISUPERVISED_POSITIVE_QVALUE
    end
    return mask
end

function _mainsearch_label_counts(
    targets::AbstractVector{Bool},
    train_mask::AbstractVector{Bool},
)
    positives = 0
    negatives = 0
    @inbounds for i in eachindex(targets, train_mask)
        train_mask[i] || continue
        if targets[i]
            positives += 1
        else
            negatives += 1
        end
    end
    return (positives = positives, negatives = negatives)
end

function _mainsearch_scale_pos_weight(
    targets::AbstractVector{Bool},
    train_mask::AbstractVector{Bool},
)
    counts = _mainsearch_label_counts(targets, train_mask)
    counts.positives == 0 && return 1.0f0
    return Float32(counts.negatives / counts.positives)
end

function _mainsearch_target_fdr_counts(
    scores::AbstractVector,
    targets::AbstractVector{Bool},
    precursor_idx::AbstractVector{UInt32};
    thresholds::Vector{Float32} = MAINSEARCH_FDR_LOG_THRESHOLDS,
)
    qvals = Vector{Float32}(undef, length(targets))
    get_qvalues!(Float32.(scores), targets, qvals; doSort=true)

    counts = Dict{Float32, NamedTuple{(:target_psms, :unique_target_precursors), Tuple{Int, Int}}}()
    for threshold in thresholds
        target_psms = 0
        target_precursors = Set{UInt32}()
        @inbounds for i in eachindex(targets, precursor_idx, qvals)
            if targets[i] && qvals[i] <= threshold
                target_psms += 1
                push!(target_precursors, precursor_idx[i])
            end
        end
        counts[threshold] = (
            target_psms = target_psms,
            unique_target_precursors = length(target_precursors),
        )
    end
    return counts
end

function _log_mainsearch_semisupervised_iteration!(
    iteration::Int,
    scores::AbstractVector,
    targets::AbstractVector{Bool},
    precursor_idx::AbstractVector{UInt32},
    train_mask::AbstractVector{Bool},
)
    label_counts = _mainsearch_label_counts(targets, train_mask)
    fdr_counts = _mainsearch_target_fdr_counts(scores, targets, precursor_idx)
    @debug_l1 "  MainSearch semi-supervised LightGBM iter $iteration training: " *
              "positives=$(label_counts.positives) negatives=$(label_counts.negatives)"
    for threshold in MAINSEARCH_FDR_LOG_THRESHOLDS
        c = fdr_counts[threshold]
        pct = round(100 * threshold, digits = 1)
        @debug_l1 "  MainSearch semi-supervised LightGBM iter $iteration q≤$(pct)%: " *
                  "target_psms=$(c.target_psms) unique_target_precursors=$(c.unique_target_precursors)"
    end
    return nothing
end

const MAINSEARCH_IRT_REFINEMENT_QVALUE_THRESHOLD = 0.01f0
const MAINSEARCH_IRT_REFINEMENT_MIN_PRECURSORS = 2
const MAINSEARCH_ENABLE_IRT_REFINEMENT = false

struct MainSearchIrtRefinementModel
    intercept::Float32
    irt_pred_coefficient::Float32
    irt_pred_squared_coefficient::Float32
    coefficients::Dict{String, Float32}
end

@inline function _mainsearch_irt_pred_basis(current_irt_pred::Float32)
    x = Float64(current_irt_pred)
    return x, x * x
end

function _mainsearch_precursor_token_counts(
    sequence::AbstractString,
    structural_mods::Union{Missing, AbstractString},
)
    counts = Dict{String, Float32}()
    position_mods = Dict{Int, Vector{String}}()
    residue_tokens = String[]

    if !ismissing(structural_mods) && !isempty(structural_mods)
        mod_regex = r"\((\d+),([A-Z]|[nc]),([^,\)]+)\)"
        for m in eachmatch(mod_regex, structural_mods)
            position = parse(Int, m.captures[1])
            site = only(m.captures[2])
            mod_name = m.captures[3]

            if site == 'n' || site == 'c'
                token = string(site, "|", mod_name)
                counts[token] = get(counts, token, 0.0f0) + 1.0f0
            elseif 1 <= position <= length(sequence)
                push!(get!(position_mods, position, String[]), mod_name)
            end
        end
    end

    for (position, aa) in enumerate(sequence)
        mods_here = get(position_mods, position, nothing)
        token = if isnothing(mods_here) || isempty(mods_here)
            string(aa)
        else
            string(aa, "|", join(sort(mods_here), "&"))
        end
        push!(residue_tokens, token)
        counts[token] = get(counts, token, 0.0f0) + 1.0f0
    end

    if !isempty(residue_tokens)
        n_term_token = "NTERM|" * first(residue_tokens)
        c_term_token = "CTERM|" * last(residue_tokens)
        counts[n_term_token] = get(counts, n_term_token, 0.0f0) + 1.0f0
        counts[c_term_token] = get(counts, c_term_token, 0.0f0) + 1.0f0
    end

    return counts
end

function _mainsearch_precursor_token_counts(
    precursors::LibraryPrecursors,
    precursor_idx::UInt32,
    token_cache::Dict{UInt32, Dict{String, Float32}},
)
    return get!(token_cache, precursor_idx) do
        sequence = String(getSequence(precursors)[precursor_idx])
        structural_mods = getStructuralMods(precursors)[precursor_idx]
        _mainsearch_precursor_token_counts(sequence, structural_mods)
    end
end

function _fit_mainsearch_irt_refinement_model(
    precursors::LibraryPrecursors,
    precursor_ids::Vector{UInt32},
    irt_pred_inputs::Vector{Float32},
    irt_corrections::Vector{Float32},
    token_cache::Dict{UInt32, Dict{String, Float32}},
)
    n = length(precursor_ids)
    if n != length(irt_pred_inputs) || n != length(irt_corrections)
        throw(DimensionMismatch("iRT refinement inputs must have equal lengths"))
    end
    n < MAINSEARCH_IRT_REFINEMENT_MIN_PRECURSORS && return nothing

    token_to_col = Dict{String, Int}()
    precursor_counts = Vector{Dict{String, Float32}}(undef, n)
    for i in eachindex(precursor_ids)
        counts = _mainsearch_precursor_token_counts(precursors, precursor_ids[i], token_cache)
        precursor_counts[i] = counts
        for token in keys(counts)
            if !haskey(token_to_col, token)
                token_to_col[token] = length(token_to_col) + 1
            end
        end
    end

    use_quadratic_basis = n >= 3
    irt_basis_cols = use_quadratic_basis ? 2 : 1
    X = zeros(Float64, n, length(token_to_col) + 1 + irt_basis_cols)
    X[:, 1] .= 1.0
    for i in eachindex(precursor_counts)
        linear_term, quadratic_term = _mainsearch_irt_pred_basis(irt_pred_inputs[i])
        X[i, 2] = linear_term
        if use_quadratic_basis
            X[i, 3] = quadratic_term
        end
        for (token, count) in precursor_counts[i]
            X[i, token_to_col[token] + 1 + irt_basis_cols] = Float64(count)
        end
    end

    β = X \ Float64.(irt_corrections)

    coefficients = Dict{String, Float32}()
    for (token, col_idx) in token_to_col
        coefficients[token] = Float32(β[col_idx + 1 + irt_basis_cols])
    end

    return MainSearchIrtRefinementModel(
        Float32(β[1]),
        Float32(β[2]),
        use_quadratic_basis ? Float32(β[3]) : 0.0f0,
        coefficients,
    )
end

function _predict_mainsearch_irt_refinement(
    model::MainSearchIrtRefinementModel,
    counts::Dict{String, Float32},
    current_irt_pred::Float32,
)
    linear_term, quadratic_term = _mainsearch_irt_pred_basis(current_irt_pred)
    prediction = model.intercept +
                 model.irt_pred_coefficient * Float32(linear_term) +
                 model.irt_pred_squared_coefficient * Float32(quadratic_term)
    for (token, count) in counts
        prediction += get(model.coefficients, token, 0.0f0) * count
    end
    return prediction
end

function _predict_mainsearch_irt_refinement(
    precursors::LibraryPrecursors,
    model::MainSearchIrtRefinementModel,
    precursor_idx::UInt32,
    current_irt_pred::Float32,
    token_cache::Dict{UInt32, Dict{String, Float32}},
)
    return _predict_mainsearch_irt_refinement(
        model,
        _mainsearch_precursor_token_counts(precursors, precursor_idx, token_cache),
        current_irt_pred,
    )
end

function _mainsearch_passing_irt_training_precursors(
    precursor_idx::AbstractVector,
    targets::AbstractVector{Bool},
    scores::AbstractVector{<:AbstractFloat},
    irt_pred::AbstractVector,
    irt_obs::AbstractVector;
    q_value_threshold::Float32 = MAINSEARCH_IRT_REFINEMENT_QVALUE_THRESHOLD,
)
    n = length(scores)
    n == 0 && return UInt32[], Float32[], Float32[]

    q_values = Vector{Float32}(undef, n)
    get_qvalues!(Float32.(scores), Bool.(targets), q_values; doSort=true)

    first_pred = Dict{UInt32, Float32}()
    sum_obs = Dict{UInt32, Float64}()
    count_obs = Dict{UInt32, Int}()
    @inbounds for i in eachindex(scores)
        if targets[i] && q_values[i] <= q_value_threshold
            pid = UInt32(precursor_idx[i])
            if !haskey(first_pred, pid)
                first_pred[pid] = Float32(irt_pred[i])
            end
            sum_obs[pid] = get(sum_obs, pid, 0.0) + Float64(irt_obs[i])
            count_obs[pid] = get(count_obs, pid, 0) + 1
        end
    end

    precursor_ids = collect(keys(sum_obs))
    sort!(precursor_ids)
    irt_pred_inputs = Float32[first_pred[pid] for pid in precursor_ids]
    irt_corrections = Float32[
        (sum_obs[pid] / count_obs[pid]) - first_pred[pid]
        for pid in precursor_ids
    ]
    return precursor_ids, irt_pred_inputs, irt_corrections
end

function _summary_quantile(values::Vector{Float32}, q::Real)
    isempty(values) && return 0.0f0
    return Float32(quantile(values, q))
end

function _format_mainsearch_irt_refinement_summary(label::AbstractString, summary)
    r(x) = round(Float64(x), digits=4)
    return "  MainSearch iRT refinement [$label]: " *
           "refined=$(summary.refined), models=$(summary.models_fit), " *
           "training_precursors=$(summary.training_precursors), " *
           "irt_error med/p90 $(r(summary.error_median_before))/$(r(summary.error_p90_before)) -> " *
           "$(r(summary.error_median_after))/$(r(summary.error_p90_after))"
end

function _refine_mainsearch_irt_after_first_iteration!(
    psms::AbstractDataFrame,
    precursors::LibraryPrecursors,
    scores::AbstractVector{<:AbstractFloat};
    label::AbstractString = "MainSearch",
)
    n = nrow(psms)
    length(scores) == n || throw(DimensionMismatch("iRT refinement scores must match PSM rows"))

    irt_error_before = Float32.(psms[!, :irt_error])
    cv_folds = sort!(unique(UInt8.(psms[!, :cv_fold])))
    token_cache = Dict{UInt32, Dict{String, Float32}}()
    models = Dict{UInt8, MainSearchIrtRefinementModel}()
    training_precursors = 0

    for fold in cv_folds
        train_idx = findall(!=(fold), psms[!, :cv_fold])
        precursor_ids, irt_pred_inputs, irt_corrections =
            _mainsearch_passing_irt_training_precursors(
                psms[train_idx, :precursor_idx],
                psms[train_idx, :target],
                scores[train_idx],
                psms[train_idx, :irt_pred],
                psms[train_idx, :irt_obs],
            )
        training_precursors += length(precursor_ids)
        model = _fit_mainsearch_irt_refinement_model(
            precursors,
            precursor_ids,
            irt_pred_inputs,
            irt_corrections,
            token_cache,
        )
        model === nothing || (models[fold] = model)
    end

    current_pred = Float32.(psms[!, :irt_pred])
    refined_pred = copy(current_pred)
    precursor_idx = psms[!, :precursor_idx]
    cv_fold = psms[!, :cv_fold]
    correction_cache = Dict{Tuple{UInt8, UInt32, Float32}, Float32}()

    @inbounds for row in 1:n
        model = get(models, UInt8(cv_fold[row]), nothing)
        model === nothing && continue
        pid = UInt32(precursor_idx[row])
        fold = UInt8(cv_fold[row])
        current = current_pred[row]
        correction = get!(correction_cache, (fold, pid, current)) do
            _predict_mainsearch_irt_refinement(
                precursors,
                model,
                pid,
                current,
                token_cache,
            )
        end
        refined_pred[row] = current + correction
    end

    refined = !isempty(models)
    if refined
        psms[!, :irt_pred] = refined_pred
        psms[!, :irt_error] = abs.(Float32.(psms[!, :irt_obs]) .- refined_pred)
    end

    irt_error_after = refined ? Float32.(psms[!, :irt_error]) : irt_error_before
    summary = (
        refined = refined,
        models_fit = length(models),
        training_precursors = training_precursors,
        error_median_before = _summary_quantile(irt_error_before, 0.5),
        error_p90_before = _summary_quantile(irt_error_before, 0.9),
        error_median_after = _summary_quantile(irt_error_after, 0.5),
        error_p90_after = _summary_quantile(irt_error_after, 0.9),
    )
    @debug_l1 _format_mainsearch_irt_refinement_summary(label, summary)
    return summary
end

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

Used by `score_precursor_isotope_traces` (PrecursorScoringSearch). MainSearch
uses `train_mainsearch_semisupervised_lgbm`.
"""
function train_psm_classifier_with_fallback(
    psms::DataFrame;
    features::Vector{Symbol},
)
    targets_col = psms[!, :target]
    n_total = nrow(psms)

    # Filter to available features
    available_features = filter(f -> hasproperty(psms, f), features)

    # Build feature matrix
    X_all = feature_matrix(psms, available_features)

    # Two-fold cross-validation using existing cv_fold column
    cv_fold = psms[!, :cv_fold]
    idx0 = findall(cv_fold .== 0)
    idx1 = findall(cv_fold .== 1)
    all_scores = Vector{Float64}(undef, n_total)

    lgbm_hp = SHARED_LGBM_HP
    MAX_TRAIN = SHARED_LGBM_MAX_TRAIN
    LOW_DATA_THRESHOLD = SHARED_LGBM_LOW_DATA_THRESHOLD

    # Fold order: (train_idx, test_idx)
    fold_pairs = [(idx1, idx0), (idx0, idx1)]

    _sample_pos(n_avail) = n_avail > MAX_TRAIN ? randperm(n_avail)[1:MAX_TRAIN] : collect(1:n_avail)
    sub_positions = [_sample_pos(length(tr)) for (tr, _) in fold_pairs]
    min_fold_size = min(length(idx0), length(idx1))
    low_data = min_fold_size < LOW_DATA_THRESHOLD

    # LightGBM CV. Slices train/test matrices on demand (transient peak) to avoid
    # retaining two ~400MB per-fold matrices across the whole function.
    function _lgbm_cv()
        fold_scores = Vector{Vector{Float64}}(undef, 2)
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
            else
                cls = build_lightgbm_classifier(; lgbm_hp...)
                tf = time(); LightGBM.fit!(cls, X_tr, y_lbl; verbosity = -1); t_fit += time() - tf
                ts2 = time(); X_te = X_all[test_idx, :]; t_slice += time() - ts2
                tp = time(); raw = LightGBM.predict(cls, X_te); t_predict += time() - tp
                fold_scores[fi] = ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
                last_cls = cls
            end
        end
        fold_scores, last_cls, (slice=t_slice, fit=t_fit, predict=t_predict)
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
        fold_scores = Vector{Vector{Float64}}(undef, 2)
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
        end
        fold_scores
    end

    function _fold_oof_psms_at_1pct(y_fold::AbstractVector, scores_on_test::AbstractVector)
        probs = Float32.(scores_on_test)
        qvals = zeros(Float16, length(probs))
        get_qvalues!(probs, y_fold, qvals)
        count(i -> qvals[i] <= Float16(0.01) && y_fold[i], eachindex(qvals))
    end

    @debug_l1 "  LightGBM CV: fold0=$(length(idx0)) fold1=$(length(idx1)) PSMs; train $(length.(sub_positions))"

    lgbm_scores, last_classifier, lgbm_timings = _lgbm_cv()
    @debug_l1 "  LightGBM timings: slice=$(round(lgbm_timings.slice, digits=2))s fit=$(round(lgbm_timings.fit, digits=2))s predict=$(round(lgbm_timings.predict, digits=2))s"

    info_winner = "lgbm"
    info_n_lgbm = -1
    info_n_probit = -1
    if low_data
        probit_scores = _probit_cv()
        n_lgbm   = _fold_oof_psms_at_1pct(targets_col[idx0], lgbm_scores[1])   + _fold_oof_psms_at_1pct(targets_col[idx1], lgbm_scores[2])
        n_probit = _fold_oof_psms_at_1pct(targets_col[idx0], probit_scores[1]) + _fold_oof_psms_at_1pct(targets_col[idx1], probit_scores[2])
        winner = n_lgbm >= n_probit ? "lgbm" : "probit"
        @debug_l1 "  Low-data CV selection (min fold=$min_fold_size < $LOW_DATA_THRESHOLD): lgbm=$n_lgbm probit=$n_probit → $winner"
        chosen = n_lgbm >= n_probit ? lgbm_scores : probit_scores
        all_scores[idx0] .= chosen[1]
        all_scores[idx1] .= chosen[2]
        info_winner = winner
        info_n_lgbm = n_lgbm
        info_n_probit = n_probit
    else
        all_scores[idx0] .= lgbm_scores[1]
        all_scores[idx1] .= lgbm_scores[2]
    end

    info = (
        slice = lgbm_timings.slice,
        fit = lgbm_timings.fit,
        predict = lgbm_timings.predict,
        low_data = low_data,
        lgbm_psms_at_1pct = info_n_lgbm,
        probit_psms_at_1pct = info_n_probit,
        winner = info_winner,
        available_features = available_features,
    )
    return all_scores, last_classifier, info
end

"""
    train_mainsearch_semisupervised_lgbm(psms; features) -> (scores, last_classifier, info)

MainSearch-only per-run LightGBM training. The first two-fold CV iteration
trains on all targets and decoys; any additional configured iterations train
on all decoys plus target rows passing the configured FDR threshold under the
previous iteration scores.
"""
function train_mainsearch_semisupervised_lgbm(
    psms::DataFrame;
    features::Vector{Symbol},
    precursors::Union{Nothing,LibraryPrecursors} = nothing,
    label::AbstractString = "MainSearch",
)
    targets_col = psms[!, :target]
    precursor_idx = psms[!, :precursor_idx]
    n_total = nrow(psms)

    available_features = filter(f -> hasproperty(psms, f), features)
    X_all = feature_matrix(psms, available_features)

    cv_fold = psms[!, :cv_fold]
    idx0 = findall(cv_fold .== 0)
    idx1 = findall(cv_fold .== 1)
    fold_pairs = [(idx1, idx0), (idx0, idx1)]

    all_scores = Vector{Float64}(undef, n_total)
    previous_scores = nothing
    last_classifier = nothing
    lgbm_hp = SHARED_LGBM_HP
    max_train = SHARED_LGBM_MAX_TRAIN

    t_slice = 0.0
    t_fit = 0.0
    t_predict = 0.0

    @debug_l1 "  MainSearch semi-supervised LightGBM CV: fold0=$(length(idx0)) fold1=$(length(idx1)) PSMs; iterations=$MAINSEARCH_SEMISUPERVISED_ITERATIONS"

    for iteration in 1:MAINSEARCH_SEMISUPERVISED_ITERATIONS
        train_mask = _mainsearch_training_mask(targets_col, previous_scores)
        label_counts = _mainsearch_label_counts(targets_col, train_mask)
        if label_counts.positives == 0 || label_counts.negatives == 0
            @user_warn "MainSearch semi-supervised LightGBM iteration has insufficient label diversity; reusing all targets and decoys" iteration=iteration positives=label_counts.positives negatives=label_counts.negatives
            train_mask .= true
        end

        for (train_idx, test_idx) in fold_pairs
            isempty(test_idx) && continue
            selected_train_idx = train_idx[train_mask[train_idx]]
            if isempty(selected_train_idx)
                all_scores[test_idx] .= 0.0
                continue
            end

            if length(selected_train_idx) > max_train
                selected_train_idx = selected_train_idx[randperm(length(selected_train_idx))[1:max_train]]
            end

            ts = time()
            X_tr = X_all[selected_train_idx, :]
            y_lbl = _prepare_labels(targets_col[selected_train_idx])
            t_slice += time() - ts

            if length(unique(y_lbl)) == 1
                all_scores[test_idx] .= y_lbl[1] == 0 ? 0.0 : 1.0
            else
                fold_targets = targets_col[selected_train_idx]
                fold_mask = trues(length(fold_targets))
                scale_pos_weight = _mainsearch_scale_pos_weight(fold_targets, fold_mask)
                @debug_l1 "  MainSearch semi-supervised LightGBM iter $iteration fold train: " *
                          "positives=$(count(identity, fold_targets)) negatives=$(count(!, fold_targets)) " *
                          "scale_pos_weight=$(round(scale_pos_weight, digits=3))"
                cls = build_lightgbm_classifier(; lgbm_hp..., scale_pos_weight=scale_pos_weight)
                tf = time()
                LightGBM.fit!(cls, X_tr, y_lbl; verbosity = -1)
                t_fit += time() - tf

                ts2 = time()
                X_te = X_all[test_idx, :]
                t_slice += time() - ts2

                tp = time()
                raw = LightGBM.predict(cls, X_te)
                t_predict += time() - tp
                scores = ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
                all_scores[test_idx] .= scores
                last_classifier = cls
            end
        end

        _log_mainsearch_semisupervised_iteration!(
            iteration, all_scores, targets_col, precursor_idx, train_mask
        )
        if MAINSEARCH_ENABLE_IRT_REFINEMENT && iteration == 1 && precursors !== nothing &&
           hasproperty(psms, :irt_pred) && hasproperty(psms, :irt_error) && hasproperty(psms, :irt_obs)
            summary = _refine_mainsearch_irt_after_first_iteration!(
                psms,
                precursors,
                all_scores;
                label = label,
            )
            if summary.refined && any(f -> f === :irt_error || f === :irt_pred, available_features)
                X_all = feature_matrix(psms, available_features)
            end
        end
        previous_scores = copy(all_scores)
    end

    @debug_l1 "  MainSearch semi-supervised LightGBM timings: slice=$(round(t_slice, digits=2))s fit=$(round(t_fit, digits=2))s predict=$(round(t_predict, digits=2))s"

    info = (
        slice = t_slice,
        fit = t_fit,
        predict = t_predict,
        low_data = false,
        lgbm_psms_at_1pct = -1,
        probit_psms_at_1pct = -1,
        winner = "lgbm_semisupervised",
        available_features = available_features,
        iterations = MAINSEARCH_SEMISUPERVISED_ITERATIONS,
    )
    return all_scores, last_classifier, info
end

"""
    _add_mainsearch_prescore_scan_features!(psms)

Broadcast per-precursor scan support features to every PSM row before the
MainSearch per-run LightGBM model is trained. `n_scans` is the number of PSM
rows for the precursor. `n_consecutive_cycles` is the length of the adjacent
cycle-index run containing this PSM, deduplicating repeated PSM rows within the
same cycle.
"""
function _add_mainsearch_prescore_scan_features!(psms::DataFrame)
    n = nrow(psms)
    n == 0 && return psms

    prec_ids = psms[!, :precursor_idx]::Vector{UInt32}
    cycle_ids = hasproperty(psms, :cycle_idx) ? psms[!, :cycle_idx] : nothing
    n_scans = Vector{UInt32}(undef, n)
    n_consecutive_cycles = cycle_ids === nothing ? nothing : Vector{UInt32}(undef, n)

    perm = sortperm(prec_ids)
    cycle_buf = UInt32[]
    run_len_buf = UInt32[]
    gi = 1
    while gi <= n
        pid = prec_ids[perm[gi]]
        group_start = gi
        while gi <= n && prec_ids[perm[gi]] == pid
            gi += 1
        end
        group_len = gi - group_start
        group_n = UInt32(group_len)

        unique_cycle_count = 0
        if cycle_ids !== nothing
            if length(cycle_buf) < group_len
                resize!(cycle_buf, group_len)
            end
            @inbounds for k in 0:(group_len - 1)
                row = perm[group_start + k]
                cycle_buf[k + 1] = UInt32(cycle_ids[row])
            end
            sort!(view(cycle_buf, 1:group_len))

            previous = UInt32(0)
            have_previous = false
            @inbounds for k in 1:group_len
                cycle = cycle_buf[k]
                if !have_previous || cycle != previous
                    unique_cycle_count += 1
                    cycle_buf[unique_cycle_count] = cycle
                    have_previous = true
                end
                previous = cycle
            end

            if length(run_len_buf) < unique_cycle_count
                resize!(run_len_buf, unique_cycle_count)
            end
            run_start = 1
            while run_start <= unique_cycle_count
                run_end = run_start
                @inbounds while run_end < unique_cycle_count &&
                                cycle_buf[run_end + 1] == cycle_buf[run_end] + UInt32(1)
                    run_end += 1
                end
                run_len = UInt32(run_end - run_start + 1)
                @inbounds for k in run_start:run_end
                    run_len_buf[k] = run_len
                end
                run_start = run_end + 1
            end
        end

        @inbounds for k in 0:(group_len - 1)
            row = perm[group_start + k]
            n_scans[row] = group_n
            if n_consecutive_cycles !== nothing
                cycle = UInt32(cycle_ids[row])
                cycle_pos = searchsortedfirst(view(cycle_buf, 1:unique_cycle_count), cycle)
                n_consecutive_cycles[row] = run_len_buf[cycle_pos]
            end
        end
    end

    psms[!, :n_scans] = n_scans
    if n_consecutive_cycles !== nothing
        psms[!, :n_consecutive_cycles] = n_consecutive_cycles
    end
    return psms
end

"""
    train_lgbm_and_select_best(psms; features) -> (best_psms, scores, timings)

Train the MainSearch semi-supervised LightGBM on all candidate PSMs, then select
the best scan per precursor by LightGBM score and log feature importances.

Returns:
- best_psms: DataFrame with one row per precursor (best by LightGBM score)
- scores: Vector{Float32} of LightGBM probabilities for best_psms
- timings: NamedTuple with timing breakdowns
"""
function train_lgbm_and_select_best(
    psms::DataFrame;
    features::Vector{Symbol} = collect(PRESCORE_FEATURES),
    precursors::Union{Nothing,LibraryPrecursors} = nothing,
    label::AbstractString = "MainSearch",
)
    t0 = time()
    _add_mainsearch_prescore_scan_features!(psms)
    all_scores, last_classifier, info = train_mainsearch_semisupervised_lgbm(
        psms;
        features = features,
        precursors = precursors,
        label = label,
    )
    t_train_cv = time()

    model = if last_classifier !== nothing
        LightGBMModel(last_classifier, info.available_features, nothing)
    else
        LightGBMModel(nothing, info.available_features, 0.0f0)
    end

    # Add scores to psms for best-per-precursor selection
    psms[!, :lgbm_score] = Float32.(all_scores)

    # Select best scan per precursor by LightGBM score
    psms = select_best_per_precursor!(psms, :lgbm_score)

    # Extract scores for best PSMs
    scores = psms[!, :lgbm_score]
    t_best = time()

    # Feature importances (debug_l2 only — `importance()` allocates a Dict per call)
    if DEBUG_CONSOLE_LEVEL[] >= 2
        imp = importance(model)
        if imp !== nothing
            sorted_imp = sort(imp, by = x -> -x[2])
            @debug_l2 "MainSearch LightGBM Feature Importances (gain):"
            for (fname, gain) in sorted_imp
                @debug_l2 "  $(rpad(fname, 40)) $(round(gain, digits=2))"
            end
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

    # Reusable buffer for p75 score computation (avoids per-group allocation)
    p75_buf = Vector{Float32}(undef, 128)

    # Single pass: process each precursor group contiguously
    gi = 1
    while gi <= n
        pid = prec_ids[perm[gi]]
        group_start = gi
        while gi <= n && prec_ids[perm[gi]] == pid
            gi += 1
        end
        group_len = gi - group_start

        # --- Sub-pass 1: find max_weight and best-score row ---
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

    return result
end
