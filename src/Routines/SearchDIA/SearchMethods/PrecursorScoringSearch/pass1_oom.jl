# Pass-1 LightGBM: sample a fixed pool, fit/select models on pool OOF scores,
# then release the pool and stream final predictions over every source file.

using Random

# Per-fold sampling cap. Resolved at call time (SCORING_LGBM_MAX_TRAIN isn't
# bound at parse time — MainSearch/scoring.jl loads AFTER this file in
# importScripts.jl).
default_pass1_oom_k_per_fold() = SCORING_LGBM_MAX_TRAIN

function _pass1_log_phase(phase::AbstractString, started::Real)
    @debug_l1 "Pass-1 $phase complete: $(round(time() - started, digits = 2))s"
end

function _pass1_log_progress(phase::AbstractString, done::Int, total::Int)
    if done % 100 == 0 || done == total
        @debug_l1 "Pass-1 $phase: files=$done/$total"
    end
end

function _count_psm_rows_by_fold(file_paths::Vector{String})
    n0 = 0; n1 = 0
    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        n = length(tbl.cv_fold)
        n == 0 && continue
        fold_c = tbl.cv_fold
        @inbounds for i in 1:n
            f = UInt8(fold_c[i])
            f == UInt8(0) && (n0 += 1; continue)
            f == UInt8(1) && (n1 += 1; continue)
        end
    end
    return n0, n1
end

# Determine which features in `requested` are actually present in the per-file
# Arrow tables. We trust the first non-empty file's schema (all per-file
# arrows produced by MainSearch share the same columns). The optional digestion
# feature gets one lightweight UInt8 scan so full-specific libraries retain the
# same model matrix they used before the feature existed.
function _resolve_available_features(file_paths::Vector{String}, requested::Vector{Symbol})
    available = Symbol[]
    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        length(tbl.precursor_idx) == 0 && continue
        available = filter(f -> hasproperty(tbl, f), requested)
        break
    end

    if :num_enzymatic_termini in available
        first_value = nothing
        found_value = false
        varies = false
        for fpath in file_paths
            tbl = Arrow.Table(fpath)
            hasproperty(tbl, :num_enzymatic_termini) || continue
            for value in tbl.num_enzymatic_termini
                if !found_value
                    first_value = value
                    found_value = true
                elseif !isequal(value, first_value)
                    varies = true
                    break
                end
            end
            varies && break
        end
        if !varies
            deleteat!(
                available,
                findfirst(==(:num_enzymatic_termini), available),
            )
        end
    end
    return available
end

# Pass 2: select final reservoir sources, then gather one file at a time.
# Selection preserves the serial RNG sequence and original reservoir slots.
# Only the final occupant of each slot is copied: replaced rows never require
# feature reads, and no Arrow columns remain referenced across files.
function _sample_both_folds(
    file_paths::Vector{String},
    features::Vector{Symbol},
    k_per_fold::Int,
    n0::Int,
    n1::Int,
    rng::AbstractRNG,
)
    nfeat = length(features)
    k0 = min(n0, k_per_fold)
    k1 = min(n1, k_per_fold)
    n_files = length(file_paths)
    max(n_files, k0, k1) <= typemax(Int32) || throw(ArgumentError(
        "Pass-1 sample file and reservoir indexes must fit in Int32",
    ))

    t_select_start = time()
    @debug_l1 "Pass-1 sample selection started: files=$n_files sample_fold0=$k0 sample_fold1=$k1"
    plan0, y0, plan1, y1 = _select_pass1_sample_plans(
        file_paths, k0, k1, n0, n1, rng,
    )
    _pass1_log_phase("sample selection", t_select_start)

    X0 = Matrix{Float32}(undef, k0, nfeat)
    X1 = Matrix{Float32}(undef, k1, nfeat)
    t_gather_start = time()
    @debug_l1 "Pass-1 sample gather started: files=$n_files rows=$(k0 + k1)"
    for file_idx in eachindex(file_paths)
        if !isempty(plan0[file_idx]) || !isempty(plan1[file_idx])
            _gather_pass1_sample_file!(
                X0, X1, file_paths[file_idx], features, plan0[file_idx], plan1[file_idx],
            )
        end
        _pass1_log_progress("sample gather", file_idx, n_files)
    end
    _pass1_log_phase("sample gather", t_gather_start)
    return X0, y0, X1, y1
end

function _select_pass1_sample_plans(
    file_paths::Vector{String}, k0::Int, k1::Int, n0::Int, n1::Int,
    rng::AbstractRNG,
)
    # A slot holds (source file, source row), replacing previous occupants.
    sources0 = Vector{Tuple{Int32, Int32}}(undef, k0)
    sources1 = Vector{Tuple{Int32, Int32}}(undef, k1)
    y0 = Vector{Bool}(undef, k0)
    y1 = Vector{Bool}(undef, k1)
    n_files = length(file_paths)
    seen0 = 0
    seen1 = 0
    for (file_idx, fpath) in enumerate(file_paths)
        seen0, seen1 = _select_pass1_sample_file!(
            sources0, y0, sources1, y1, fpath, Int32(file_idx),
            seen0, seen1, rng,
        )
        _pass1_log_progress("sample selection", file_idx, n_files)
    end
    @assert seen0 == n0 "Sampler saw $seen0 rows for fold=0 but Pass-1 metadata counted $n0"
    @assert seen1 == n1 "Sampler saw $seen1 rows for fold=1 but Pass-1 metadata counted $n1"

    return _group_pass1_sample_sources(sources0, n_files), y0,
           _group_pass1_sample_sources(sources1, n_files), y1
end

# Keep each table inside this scope; specialize the hot loop on its concrete
# column types instead of indexing Arrow.Table's heterogeneous columns there.
function _select_pass1_sample_file!(
    sources0, y0, sources1, y1, fpath::String, file_idx::Int32,
    seen0::Int, seen1::Int, rng::AbstractRNG,
)
    tbl = Arrow.Table(fpath)
    fold_c = tbl.cv_fold
    target_c = tbl.target
    n = length(fold_c)
    n <= typemax(Int32) || throw(ArgumentError(
        "Pass-1 source row indexes must fit in Int32: $fpath has $n rows",
    ))
    return _update_pass1_sample_sources!(
        sources0, y0, sources1, y1, fold_c, target_c, file_idx,
        seen0, seen1, rng,
    )
end

function _update_pass1_sample_sources!(
    sources0, y0, sources1, y1, fold_c::AbstractVector, target_c::AbstractVector,
    file_idx::Int32, seen0::Int, seen1::Int, rng::AbstractRNG,
)
    k0 = length(sources0)
    k1 = length(sources1)
    @inbounds for i in eachindex(fold_c)
        f = UInt8(fold_c[i])
        if f == UInt8(0)
            seen0 += 1
            pos = if seen0 <= k0
                seen0
            else
                draw = rand(rng, 1:seen0)
                draw > k0 && continue
                draw
            end
            sources0[pos] = (file_idx, Int32(i))
            y0[pos] = Bool(target_c[i])
        elseif f == UInt8(1)
            seen1 += 1
            pos = if seen1 <= k1
                seen1
            else
                draw = rand(rng, 1:seen1)
                draw > k1 && continue
                draw
            end
            sources1[pos] = (file_idx, Int32(i))
            y1[pos] = Bool(target_c[i])
        end
    end
    return seen0, seen1
end

function _group_pass1_sample_sources(sources::Vector{Tuple{Int32, Int32}}, n_files::Int)
    counts = zeros(Int, n_files)
    @inbounds for (file_idx, _) in sources
        counts[file_idx] += 1
    end
    plans = [Vector{Tuple{Int32, Int32}}(undef, n) for n in counts]
    fill!(counts, 0)
    @inbounds for dst_i in eachindex(sources)
        file_idx, src_i = sources[dst_i]
        counts[file_idx] += 1
        plans[file_idx][counts[file_idx]] = (src_i, Int32(dst_i))
    end
    return plans
end

function _gather_pass1_sample_file!(
    X0::Matrix{Float32}, X1::Matrix{Float32}, fpath::String,
    features::Vector{Symbol}, p0::Vector{Tuple{Int32, Int32}}, p1::Vector{Tuple{Int32, Int32}},
)
    tbl = Arrow.Table(fpath)
    feat_cols = AbstractVector[getproperty(tbl, f) for f in features]

    # One column per task keeps writes disjoint. The barrier specializes the
    # heterogeneous Arrow feature columns without per-value boxing.
    Threads.@threads for j in eachindex(features)
        _copy_sampled_column!(X0, X1, feat_cols[j], p0, p1, j)
    end
    return nothing
end

# Copy planned (src_row -> dst_row) values from one feature column into the
# sampled fold matrices. Isolated as its own method to act as a function barrier:
# `col` is specialized to its concrete type here, so `Float32(col[src_i])` does
# not box (see caller).
@inline function _copy_sampled_column!(
    X0::Matrix{Float32}, X1::Matrix{Float32}, col::AbstractVector,
    p0::Vector{Tuple{Int32, Int32}}, p1::Vector{Tuple{Int32, Int32}}, j::Int,
)
    @inbounds for k in eachindex(p0)
        src_i, dst_i = p0[k]
        X0[dst_i, j] = Float32(col[src_i])
    end
    @inbounds for k in eachindex(p1)
        src_i, dst_i = p1[k]
        X1[dst_i, j] = Float32(col[src_i])
    end
    return
end

# Fit a prediction-only model; no training dataset is retained with the trees.
function _fit_pass1_booster(
    X::Matrix{Float32}, y::Vector{Bool}, lgbm_hp::NamedTuple,
)
    if isempty(y) || length(unique(y)) == 1
        constant_score = isempty(y) || !y[1] ? 0.0f0 : 1.0f0
        return (kind = :constant, value = constant_score, model = nothing)
    end
    cls = build_lightgbm_classifier(; lgbm_hp...)
    LightGBM.fit!(cls, X, _prepare_labels(y); verbosity = -1)
    # LightGBM.jl implements Booster deepcopy by serializing/reloading the trees.
    cls.booster = deepcopy(cls.booster)
    return cls
end

_pass1_importance_classifier(cls) = cls isa LightGBM.LGBMClassification ? cls : nothing

function _pass1_fold_row_indices(fold_col::AbstractVector)
    n0 = 0
    n1 = 0
    @inbounds for i in eachindex(fold_col)
        f = UInt8(fold_col[i])
        f == UInt8(0) && (n0 += 1; continue)
        f == UInt8(1) && (n1 += 1; continue)
    end

    idx0 = Vector{Int}(undef, n0)
    idx1 = Vector{Int}(undef, n1)
    pos0 = 0
    pos1 = 0
    @inbounds for i in eachindex(fold_col)
        f = UInt8(fold_col[i])
        if f == UInt8(0)
            pos0 += 1
            idx0[pos0] = i
        elseif f == UInt8(1)
            pos1 += 1
            idx1[pos1] = i
        end
    end
    return (fold0 = idx0, fold1 = idx1)
end

function _fill_fold_column!(
    M::Matrix{Float32},
    j::Int,
    col::AbstractVector{T},
    row_idx::AbstractVector{<:Integer},
) where {T<:Real}
    @inbounds @simd for k in eachindex(row_idx)
        M[k, j] = Float32(col[Int(row_idx[k])])
    end
end

function _fill_fold_column!(
    M::Matrix{Float32},
    j::Int,
    col::AbstractVector{<:Union{Missing, T}},
    row_idx::AbstractVector{<:Integer},
) where {T<:Real}
    @inbounds for k in eachindex(row_idx)
        v = col[Int(row_idx[k])]
        M[k, j] = v === missing ? 0.0f0 : Float32(v)
    end
end

_fill_fold_column!(::Matrix{Float32}, ::Int, col::AbstractVector, ::AbstractVector) =
    throw(ArgumentError("Unsupported feature type $(eltype(col)) for LightGBM"))

function _pass1_fold_feature_matrices(
    tbl,
    features::Vector{Symbol},
    idx0::AbstractVector{<:Integer},
    idx1::AbstractVector{<:Integer},
)
    nfeat = length(features)
    X0 = Matrix{Float32}(undef, length(idx0), nfeat)
    X1 = Matrix{Float32}(undef, length(idx1), nfeat)
    feat_cols = AbstractVector[getproperty(tbl, features[j]) for j in 1:nfeat]
    Threads.@threads for j in 1:nfeat
        col = feat_cols[j]
        _fill_fold_column!(X0, j, col, idx0)
        _fill_fold_column!(X1, j, col, idx1)
    end
    return X0, X1
end

function _predict_pass1_block(cls, X::Matrix{Float32})
    n = size(X, 1)
    n == 0 && return Float32[]
    if cls isa NamedTuple && cls.kind == :constant
        return fill(clamp(Float32(cls.value), 1f-6, 1f0 - 1f-4), n)
    end
    # Call from the main thread, with files/folds processed sequentially, to
    # avoid concurrent entry into libomp. Each predict uses the Julia core count.
    raw = LightGBM.predict(cls, X; num_threads = Threads.nthreads())
    scores = ndims(raw) == 2 ? dropdims(raw; dims = 2) : raw
    return clamp.(Float32.(scores), 1f-6, 1f0 - 1f-4)
end

# Each file keeps its original row order. OOF uses the opposite-fold model;
# optional in-fold scores use the model trained on that row's fold.
function _predict_pass1_to_sidecar(
    fpath::String, features::Vector{Symbol},
    cls_trained_on::NamedTuple, compute_infold::Bool,
)
    tbl = Arrow.Table(fpath)
    n = length(tbl.precursor_idx)
    oof = Vector{Float32}(undef, n)
    infold = fill(NaN32, n)
    if n > 0
        fold_rows = _pass1_fold_row_indices(tbl.cv_fold)
        idx0, idx1 = fold_rows.fold0, fold_rows.fold1
        X0, X1 = _pass1_fold_feature_matrices(tbl, features, idx0, idx1)
        oof[idx0] .= _predict_pass1_block(cls_trained_on.fold1, X0)
        oof[idx1] .= _predict_pass1_block(cls_trained_on.fold0, X1)
        if compute_infold
            infold[idx0] .= _predict_pass1_block(cls_trained_on.fold0, X0)
            infold[idx1] .= _predict_pass1_block(cls_trained_on.fold1, X1)
        end
    end
    writeArrow(fpath * PASS1_SIDECAR_SUFFIX, DataFrame(
        precursor_idx = collect(UInt32.(tbl.precursor_idx)),
        scan_idx = collect(UInt32.(tbl.scan_idx)),
        trace_prob_prepass = oof,
        trace_prob_infold = infold,
    ))
    return nothing
end

function _predict_pass1_files(
    file_paths::Vector{String}, features::Vector{Symbol},
    cls_trained_on::NamedTuple, compute_infold::Bool,
)
    # Keep native prediction on the main thread; each call already uses all cores.
    for (file_idx, path) in enumerate(file_paths)
        _predict_pass1_to_sidecar(path, features, cls_trained_on, compute_infold)
        _pass1_log_progress("full-data prediction", file_idx, length(file_paths))
    end
    return nothing
end

# Fit one fold's currently eligible pool rows. Keep this in a function so the
# temporary filtered matrix cannot be retained by an iteration's saved state.
function _fit_pass1_pool_fold(
    X::Matrix{Float32}, y::Vector{Bool}, training_mask,
    lgbm_hp::NamedTuple,
)
    if training_mask === nothing || all(training_mask)
        X_fit, y_fit = X, y
    else
        rows = findall(training_mask)
        X_fit, y_fit = X[rows, :], y[rows]
    end
    n_targets = count(y_fit)
    n_decoys = length(y_fit) - n_targets
    cls = _fit_pass1_booster(X_fit, y_fit, lgbm_hp)
    return (classifier = cls, targets = n_targets, decoys = n_decoys)
end

# Return only models/scalars so the pool and its score/mask buffers can be
# collected before final prediction. Excluded targets still get pool OOF scores.
function _train_pass1_pool(
    X0::Matrix{Float32}, y0::Vector{Bool},
    X1::Matrix{Float32}, y1::Vector{Bool};
    lgbm_hp::NamedTuple = SHARED_LGBM_HP,
    semisupervised::Bool = false,
    semisupervised_train_q_threshold::Float32 = SCORING_SEMISUPERVISED_TRAIN_QVALUE_THRESHOLD,
    semisupervised_stop_q_threshold::Float32 = SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD,
    semisupervised_min_gain::Float32 = SCORING_SEMISUPERVISED_MIN_TARGET_GAIN,
    semisupervised_max_iterations::Int = SCORING_SEMISUPERVISED_MAX_ITERATIONS,
)
    n0, n1 = length(y0), length(y1)
    n_pool = n0 + n1
    targets = vcat(y0, y1)
    scores = Vector{Float32}(undef, n_pool)
    training_mask = nothing
    best_state = nothing
    previous_target_q01 = -1
    max_iterations = semisupervised ? semisupervised_max_iterations : 1

    for iter_idx in 1:max_iterations
        mask0 = training_mask === nothing ? nothing : @view training_mask[1:n0]
        mask1 = training_mask === nothing ? nothing : @view training_mask[(n0 + 1):n_pool]
        @debug_l1 "Pass-1 pool iter $iter_idx: fold 0 fit starting"
        started = time()
        fit0 = _fit_pass1_pool_fold(X0, y0, mask0, lgbm_hp)
        # The filtered copy leaves scope only when the helper returns. Reclaim
        # it before the other fold allocates its own training subset.
        GC.gc()
        _pass1_log_phase("pool iter $iter_idx fold 0 fit", started)

        @debug_l1 "Pass-1 pool iter $iter_idx: fold 1 fit starting"
        started = time()
        fit1 = _fit_pass1_pool_fold(X1, y1, mask1, lgbm_hp)
        GC.gc()
        _pass1_log_phase("pool iter $iter_idx fold 1 fit", started)
        cls_trained_on = (fold0 = fit0.classifier, fold1 = fit1.classifier)
        # The full pool stays available for prediction, including rows left
        # out of these fits, and for all subsequent iterations.

        @debug_l1 "Pass-1 pool iter $iter_idx: OOF prediction starting; rows=$n_pool"
        started = time()
        scores[1:n0] .= _predict_pass1_block(cls_trained_on.fold1, X0)
        scores[(n0 + 1):n_pool] .= _predict_pass1_block(cls_trained_on.fold0, X1)
        _pass1_log_phase("pool iter $iter_idx OOF prediction", started)

        @debug_l1 "Pass-1 pool iter $iter_idx: q-value sorting and masks starting; rows=$n_pool"
        started = time()
        metrics = _scoring_semisupervised_metrics_and_mask(
            scores, targets;
            train_q_threshold = semisupervised_train_q_threshold,
            stop_q_threshold = semisupervised_stop_q_threshold,
        )
        _pass1_log_phase("pool iter $iter_idx q-value sorting and masks", started)
        n_train_targets = fit0.targets + fit1.targets
        n_train_decoys = fit0.decoys + fit1.decoys
        @debug_l1 "  ScoringSearch semi-supervised iter $iter_idx (Pass-1 training pool): " *
                   "train targets=$n_train_targets decoys=$n_train_decoys; " *
                   "pool_fold0=$n0 pool_fold1=$n1; " *
                   "pool q≤$semisupervised_stop_q_threshold targets=$(metrics.target_q01) decoys=$(metrics.decoy_q01)"

        importance_classifier = _pass1_importance_classifier(fit0.classifier)
        importance_classifier === nothing &&
            (importance_classifier = _pass1_importance_classifier(fit1.classifier))
        state = (
            iter = iter_idx,
            cls_trained_on = cls_trained_on,
            last_classifier = importance_classifier,
            target_q01 = metrics.target_q01,
            decoy_q01 = metrics.decoy_q01,
        )
        best_state = _scoring_better_iteration_state(best_state, state)
        if semisupervised && iter_idx > 1 && !_scoring_target_gain_sufficient(
            previous_target_q01, state.target_q01;
            min_fraction = semisupervised_min_gain,
        )
            @debug_l1 "  ScoringSearch semi-supervised stopping (Pass-1 training pool): " *
                       "iter $iter_idx pool targets=$(state.target_q01) did not improve by " *
                       "$(round(100 * semisupervised_min_gain, digits = 2))% over $previous_target_q01; " *
                       "using iter $(best_state.iter) with pool targets=$(best_state.target_q01)"
            break
        elseif !semisupervised
            break
        elseif iter_idx == max_iterations
            @debug_l1 "  ScoringSearch semi-supervised stopping (Pass-1 training pool): " *
                       "hit max iterations $max_iterations; using iter $(best_state.iter) " *
                       "with pool targets=$(best_state.target_q01)"
            break
        end
        previous_target_q01 = state.target_q01
        training_mask = metrics.training_mask
    end
    return best_state
end

"""
    train_and_predict_pass1_oom!(file_paths; features, compute_infold,
        k_per_fold=default_pass1_oom_k_per_fold(), semisupervised=false, ...)

Sample up to `k_per_fold` rows per original CV fold once, then fit and select
models using pool OOF scores. Later iterations retain all pool decoys and
q-value-passing pool targets. Returned discovery counts describe the pool.
After selection, release the pool and write OOF/optional in-fold sidecars for
all input rows; final experiment-wide FDR remains downstream.
"""
function train_and_predict_pass1_oom!(
    file_paths::Vector{String};
    features::Vector{Symbol},
    compute_infold::Bool,
    lgbm_hp::NamedTuple = SHARED_LGBM_HP,
    k_per_fold::Int = default_pass1_oom_k_per_fold(),
    rng::AbstractRNG = MersenneTwister(1776),
    semisupervised::Bool = false,
    semisupervised_train_q_threshold::Float32 = SCORING_SEMISUPERVISED_TRAIN_QVALUE_THRESHOLD,
    semisupervised_stop_q_threshold::Float32 = SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD,
    semisupervised_min_gain::Float32 = SCORING_SEMISUPERVISED_MIN_TARGET_GAIN,
    semisupervised_max_iterations::Int = SCORING_SEMISUPERVISED_MAX_ITERATIONS,
)
    k_per_fold > 0 || throw(ArgumentError("k_per_fold must be positive"))
    semisupervised_max_iterations > 0 ||
        throw(ArgumentError("semisupervised_max_iterations must be positive"))
    t_total_start = time()
    @debug_l1 "Pass-1 metadata starting: files=$(length(file_paths))"
    available = _resolve_available_features(file_paths, features)
    if isempty(available)
        error("train_and_predict_pass1_oom!: no requested features are present " *
              "in the per-file Arrow schema")
    end
    total_n0, total_n1 = _count_psm_rows_by_fold(file_paths)
    n_total = total_n0 + total_n1
    _pass1_log_phase("metadata", t_total_start)
    if n_total == 0
        @user_warn "train_and_predict_pass1_oom!: no PSM rows found across $(length(file_paths)) files"
        return (n_total = 0, last_classifier = nothing, available_features = available)
    end

    @debug_l1 "Pass-1 OOM training:"
    @debug_l1 "  rows: total=$n_total  fold0=$total_n0  fold1=$total_n1"
    @debug_l1 "  fixed representative pool: k_per_fold=$k_per_fold"
    @debug_l1 "  features in use: $(length(available)) / $(length(features))"
    X0, y0, X1, y1 = _sample_both_folds(
        file_paths, available, k_per_fold, total_n0, total_n1, rng,
    )
    n_pool_fold0, n_pool_fold1 = length(y0), length(y1)
    best_state = _train_pass1_pool(
        X0, y0, X1, y1;
        lgbm_hp = lgbm_hp,
        semisupervised = semisupervised,
        semisupervised_train_q_threshold = semisupervised_train_q_threshold,
        semisupervised_stop_q_threshold = semisupervised_stop_q_threshold,
        semisupervised_min_gain = semisupervised_min_gain,
        semisupervised_max_iterations = semisupervised_max_iterations,
    )
    # No model retains training data. Drop the pool before allocating per-file
    # prediction matrices; the helper's scores and masks are already out of scope.
    X0 = nothing; X1 = nothing; y0 = nothing; y1 = nothing
    GC.gc()

    @debug_l1 "Pass-1 full-data prediction starting: files=$(length(file_paths)) rows=$n_total; " *
              "selected_iter=$(best_state.iter)"
    started = time()
    _predict_pass1_files(
        file_paths, available, best_state.cls_trained_on, compute_infold,
    )
    _pass1_log_phase("full-data prediction", started)
    elapsed_total = time() - t_total_start
    @debug_l1 "Pass-1 OOM elapsed: total=$(round(elapsed_total, digits = 2))s"
    return (
        n_total = n_total,
        n_fold0 = total_n0,
        n_fold1 = total_n1,
        k_per_fold = k_per_fold,
        n_pool_fold0 = n_pool_fold0,
        n_pool_fold1 = n_pool_fold1,
        metric_scope = :training_pool,
        last_classifier = best_state.last_classifier,
        cls_trained_on = best_state.cls_trained_on,
        available_features = available,
        semisupervised_iter = best_state.iter,
        target_q01 = best_state.target_q01,
        decoy_q01 = best_state.decoy_q01,
        elapsed_s = elapsed_total,
    )
end
