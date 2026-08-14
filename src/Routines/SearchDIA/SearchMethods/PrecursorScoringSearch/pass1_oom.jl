# Pass-1 LightGBM training, out-of-memory variant.
#
# Replaces the previous load-everything-into-best_psms approach with three
# streaming passes over the per-file Arrow tables:
#
#   Pass 1 (metadata): Count rows per cv_fold; pick the feature list. Reads
#                      Arrow metadata plus the optional UInt8 digestion column
#                      to remove it when it is constant.
#   Pass 2 (sample):   Reservoir-sample up to K rows per fold into a fixed-
#                      size feature matrix. Train one LightGBM booster per
#                      fold on its sampled matrix.
#   Pass 3 (predict):  For each file, build the file's feature matrix and
#                      predict OOF + in-fold trace_prob scores, writing
#                      them straight to the file's .pass1_sidecar.arrow.
#
# Memory profile: bounded by the per-fold sample buffers (K × n_features ×
# 4 B each, default ~50 MB per fold) plus two LightGBM boosters and one
# file's feature matrix at a time. Independent of total PSM count.
#
# Probit fallback is not implemented here — at OOM scale the per-fold
# sample is always >> the low-data threshold. Caller should keep the
# in-memory train_psm_classifier_with_fallback for non-MBR / tiny-data runs.

using Random

# Per-fold sampling cap. Resolved at call time (SCORING_LGBM_MAX_TRAIN isn't
# bound at parse time — MainSearch/scoring.jl loads AFTER this file in
# importScripts.jl).
default_pass1_oom_k_per_fold() = SCORING_LGBM_MAX_TRAIN

# Pass 1: count rows per cv_fold across all files. Returns (n_fold0, n_fold1).
# When `training_masks` is provided, only rows with a true mask entry are
# counted for the training sample; prediction still scores every row.
function _count_psm_rows_by_fold(file_paths::Vector{String}, training_masks = nothing)
    n0 = 0; n1 = 0
    for (file_idx, fpath) in enumerate(file_paths)
        tbl = Arrow.Table(fpath)
        n = length(tbl.cv_fold)
        n == 0 && continue
        mask = training_masks === nothing ? nothing : training_masks[file_idx]
        fold_c = tbl.cv_fold
        @inbounds for i in 1:n
            mask !== nothing && !mask[i] && continue
            f = UInt8(fold_c[i])
            f == UInt8(0) && (n0 += 1; continue)
            f == UInt8(1) && (n1 += 1; continue)
        end
    end
    return n0, n1
end

function _pass1_file_row_counts(file_paths::Vector{String})
    row_counts = Vector{Int}(undef, length(file_paths))
    total_rows = 0
    for (file_idx, fpath) in enumerate(file_paths)
        n = length(Arrow.Table(fpath).precursor_idx)
        row_counts[file_idx] = n
        total_rows += n
    end
    return row_counts, total_rows
end

function _pass1_file_offsets(row_counts::AbstractVector{<:Integer})
    offsets = Vector{Int}(undef, length(row_counts) + 1)
    offsets[1] = 0
    @inbounds for i in eachindex(row_counts)
        offsets[i + 1] = offsets[i] + Int(row_counts[i])
    end
    return offsets
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

# Pass 2: single-pass reservoir-sample both folds simultaneously.
# Walks each file ONCE (vs the previous two-pass version that walked every
# file twice — once per fold), routing each row to its fold-specific
# reservoir based on `cv_fold`. When `k_per_fold >= n_avail_in_fold` (the
# default regime — k_per_fold = 2.5M, typical n_avail = 0.5-2M), every
# row is accepted, no random draws happen, and the function is purely an
# I/O + feature-gather pass. When `k_per_fold < n_avail_in_fold`, standard
# reservoir replace-with-prob k/seen kicks in for that fold.
function _sample_both_folds(
    file_paths::Vector{String},
    features::Vector{Symbol},
    k_per_fold::Int,
    n0::Int,
    n1::Int,
    rng::AbstractRNG,
    training_masks = nothing,
)
    nfeat = length(features)
    k0 = min(n0, k_per_fold)
    k1 = min(n1, k_per_fold)
    X0 = Matrix{Float32}(undef, k0, nfeat)
    X1 = Matrix{Float32}(undef, k1, nfeat)
    y0 = Vector{Bool}(undef, k0)
    y1 = Vector{Bool}(undef, k1)
    n_files = length(file_paths)

    # Phase A (serial): run the reservoir-sampling decision logic to build a
    # per-file plan of (src_row_in_file => dst_row_in_matrix) pairs. We also
    # fill y0/y1 here. Keeping this serial preserves the existing RNG
    # sequence so IDs are bit-stable vs the prior column-major loop.
    plan0 = Vector{Vector{Tuple{Int32, Int32}}}(undef, n_files)
    plan1 = Vector{Vector{Tuple{Int32, Int32}}}(undef, n_files)
    all_feat_cols = Vector{Vector{AbstractVector}}(undef, n_files)
    seen0 = 0
    seen1 = 0
    for (file_idx, fpath) in enumerate(file_paths)
        tbl = Arrow.Table(fpath)
        n = length(tbl.cv_fold)
        plan0[file_idx] = Tuple{Int32, Int32}[]
        plan1[file_idx] = Tuple{Int32, Int32}[]
        all_feat_cols[file_idx] = AbstractVector[getproperty(tbl, f) for f in features]
        n == 0 && continue
        mask = training_masks === nothing ? nothing : training_masks[file_idx]
        fold_c = tbl.cv_fold
        target_c = tbl.target
        @inbounds for i in 1:n
            mask !== nothing && !mask[i] && continue
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
                push!(plan0[file_idx], (Int32(i), Int32(pos)))
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
                push!(plan1[file_idx], (Int32(i), Int32(pos)))
                y1[pos] = Bool(target_c[i])
            end
        end
    end
    @assert seen0 == n0 "Sampler saw $seen0 rows for fold=0 but Pass-1 metadata counted $n0"
    @assert seen1 == n1 "Sampler saw $seen1 rows for fold=1 but Pass-1 metadata counted $n1"

    # Phase B (column-parallel): execute the plan. Each thread takes one
    # output column j, walks every file's plan, and copies the source values
    # at planned indices into X0[:, j] and X1[:, j]. Plans are independent
    # of thread, so all threads write disjoint matrix locations.
    Threads.@threads for j in 1:nfeat
        for file_idx in 1:n_files
            # Function barrier: `all_feat_cols[file_idx]` is a Vector{AbstractVector}
            # (heterogeneous feature column types), so `col` has an abstract static
            # type and indexing it boxes a value per element. Passing it to a method
            # specializes the copy loop on the concrete column type, eliminating the
            # per-element boxing in this hot loop.
            _copy_sampled_column!(X0, X1, all_feat_cols[file_idx][j],
                                  plan0[file_idx], plan1[file_idx], j)
        end
    end

    return X0, y0, X1, y1
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

# Fit a LightGBM classifier on the sampled (X, y). Returns the booster.
function _fit_pass1_booster(
    X::Matrix{Float32}, y::Vector{Bool}, lgbm_hp::NamedTuple,
)
    if isempty(y) || length(unique(y)) == 1
        constant_score = isempty(y) || !y[1] ? 0.0f0 : 1.0f0
        return (kind = :constant, value = constant_score, model = nothing)
    end
    cls = build_lightgbm_classifier(; lgbm_hp...)
    LightGBM.fit!(cls, X, _prepare_labels(y); verbosity = -1)
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

# Pass 3: predict each file's rows and write Pass-1 sidecar.
#
# `cls_trained_on` is a NamedTuple keyed by the fold the booster was trained
# on, e.g. `(fold0 = booster_fit_on_fold0_rows, fold1 = booster_fit_on_fold1_rows)`.
#   - OOF for fold0 rows uses cls_trained_on.fold1 (the model that didn't see them)
#   - In-fold for fold0 rows uses cls_trained_on.fold0 (only if compute_infold)
# and symmetrically for fold1.
#
# Sidecar layout matches the post-integration MBR staging contract:
#   (:precursor_idx, :scan_idx, :trace_prob_prepass, :trace_prob_infold)
#
# Scores are clamped to [1e-6, 1 - 1e-4] (matches the in-memory pipeline's
# clamp at score_psms.jl).
function _predict_pass1_to_sidecar(
    fpath::String, features::Vector{Symbol},
    cls_trained_on::NamedTuple,
    compute_infold::Bool,
    write_sidecar::Bool = true;
    scores_out = nothing,
    targets_out = nothing,
    return_predictions::Bool = true,
)
    tbl = Arrow.Table(fpath)
    n = length(tbl.precursor_idx)
    if n == 0
        # Empty file — write an empty sidecar to keep downstream invariants.
        if write_sidecar
            empty_df = DataFrame(
                precursor_idx       = UInt32[],
                scan_idx            = UInt32[],
                trace_prob_prepass  = Float32[],
                trace_prob_infold   = Float32[],
            )
            writeArrow(fpath * PASS1_SIDECAR_SUFFIX, empty_df)
        end
        return return_predictions ? (scores = Float32[], targets = Bool[]) : nothing
    end

    fold_rows = _pass1_fold_row_indices(tbl.cv_fold)
    idx0 = fold_rows.fold0
    idx1 = fold_rows.fold1
    X0, X1 = _pass1_fold_feature_matrices(tbl, features, idx0, idx1)

    @inline function _predict_block(cls, rows::Vector{Int}, X_fold::Matrix{Float32})
        isempty(rows) && return Float32[]
        if cls isa NamedTuple && cls.kind == :constant
            return fill(clamp(Float32(cls.value), 1f-6, 1f0 - 1f-4), length(rows))
        end
        # Design A: Pass-3 runs SEQUENTIALLY on the main thread (see
        # _predict_pass1_files below), so the predict can safely use all cores.
        # libomp is only ever entered from the main thread (never a worker), so
        # no lock and no crash/corruption. Multi-threaded predict per file.
        raw = LightGBM.predict(cls, X_fold; num_threads=Threads.nthreads())
        s = ndims(raw) == 2 ? dropdims(raw; dims = 2) : raw
        return clamp.(Float32.(s), 1f-6, 1f0 - 1f-4)
    end

    if targets_out !== nothing
        target_c = tbl.target
        @inbounds for i in 1:n
            targets_out[i] = Bool(target_c[i])
        end
    end

    # OOF: opposite-fold booster on each subset.
    oof = scores_out === nothing ? Vector{Float32}(undef, n) : scores_out
    @inbounds oof[idx0] .= _predict_block(cls_trained_on.fold1, idx0, X0)
    @inbounds oof[idx1] .= _predict_block(cls_trained_on.fold0, idx1, X1)

    # In-fold: same-fold booster on each subset. During semi-supervised
    # iterations we only need OOF scores for q-value masks, so skip the NaN
    # placeholder unless a sidecar will actually be written.
    infold = if compute_infold
        v = Vector{Float32}(undef, n)
        @inbounds v[idx0] .= _predict_block(cls_trained_on.fold0, idx0, X0)
        @inbounds v[idx1] .= _predict_block(cls_trained_on.fold1, idx1, X1)
        v
    elseif write_sidecar
        fill(NaN32, n)
    else
        nothing
    end

    if write_sidecar
        side_df = DataFrame(
            precursor_idx       = collect(UInt32.(tbl.precursor_idx)),
            scan_idx            = collect(UInt32.(tbl.scan_idx)),
            trace_prob_prepass  = oof,
            trace_prob_infold   = infold,
        )
        writeArrow(fpath * PASS1_SIDECAR_SUFFIX, side_df)
    end
    if return_predictions
        targets = targets_out === nothing ? Vector{Bool}(tbl.target) : Vector{Bool}(targets_out)
        return (scores = scores_out === nothing ? oof : Vector{Float32}(scores_out), targets = targets)
    end
    return nothing
end

function _predict_pass1_files(
    file_paths::Vector{String},
    features::Vector{Symbol},
    cls_trained_on::NamedTuple,
    compute_infold::Bool;
    write_sidecars::Bool,
    scores_out = nothing,
    targets_out = nothing,
    file_offsets = nothing,
)
    direct_output = scores_out !== nothing || targets_out !== nothing
    collect_results = !direct_output && !write_sidecars
    result_type = NamedTuple{(:scores, :targets), Tuple{Vector{Float32}, Vector{Bool}}}
    results = collect_results ? Vector{result_type}(undef, length(file_paths)) : nothing
    # Design A: sequential on the main thread. The dominant cost (the LightGBM
    # predict) is recovered via multi-threaded predict per file (_predict_block),
    # which is only safe because it runs on the main thread here. Matrix build +
    # sidecar write lose across-file parallelism (acceptable if predict-bound).
    for f_idx in 1:length(file_paths)
        first_row = direct_output ? file_offsets[f_idx] + 1 : 1
        last_row = direct_output ? file_offsets[f_idx + 1] : 0
        score_slice = scores_out === nothing ? nothing : @view scores_out[first_row:last_row]
        target_slice = targets_out === nothing ? nothing : @view targets_out[first_row:last_row]
        result = _predict_pass1_to_sidecar(
            file_paths[f_idx],
            features,
            cls_trained_on,
            compute_infold,
            write_sidecars;
            scores_out = score_slice,
            targets_out = target_slice,
            return_predictions = collect_results,
        )
        collect_results && (results[f_idx] = result)
    end
    return results
end

"""
    train_and_predict_pass1_oom!(file_paths, precursors;
                                 features, compute_infold,
                                 lgbm_hp = SHARED_LGBM_HP,
                                 k_per_fold = PASS1_OOM_K_PER_FOLD,
                                 rng = MersenneTwister(1776))
        -> NamedTuple

Out-of-memory Pass-1 LightGBM training. Streams the per-file second-pass
Arrow tables (no in-memory `best_psms`), reservoir-samples up to `k_per_fold`
rows per cv_fold, trains one booster per fold, and writes per-file
`.pass1_sidecar.arrow` with `trace_prob_prepass` (OOF) and
`trace_prob_infold` (in-fold, if requested).

Returns a NamedTuple with diagnostics — including the last classifier
for downstream feature-importance reporting.
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
    t_total_start = time()

    # Pass 1: metadata.
    available = _resolve_available_features(file_paths, features)
    if isempty(available)
        error("train_and_predict_pass1_oom!: no requested features are present " *
              "in the per-file Arrow schema")
    end
    total_n0, total_n1 = _count_psm_rows_by_fold(file_paths)
    n_total = total_n0 + total_n1
    if n_total == 0
        @user_warn "train_and_predict_pass1_oom!: no PSM rows found across $(length(file_paths)) files"
        return (n_total = 0, last_classifier = nothing, available_features = available)
    end

    @debug_l1 "Pass-1 OOM training:"
    @debug_l1 "  rows: total=$n_total  fold0=$total_n0  fold1=$total_n1"
    @debug_l1 "  sampling: k_per_fold=$k_per_fold"
    @debug_l1 "  features in use: $(length(available)) / $(length(features))"

    row_counts, total_rows = _pass1_file_row_counts(file_paths)
    file_offsets = _pass1_file_offsets(row_counts)
    scores_buffer = Vector{Float32}(undef, total_rows)
    targets_buffer = Vector{Bool}(undef, total_rows)

    total_sample = 0.0
    total_fit = 0.0
    total_predict = 0.0

    function _run_iteration(iter_idx::Int, training_masks)
        n0, n1 = _count_psm_rows_by_fold(file_paths, training_masks)
        n_fit = n0 + n1
        n_fit == 0 && error("train_and_predict_pass1_oom!: no rows left for semi-supervised training")
        k0 = min(n0, k_per_fold)
        k1 = min(n1, k_per_fold)

        # Pass 2: single-pass sampler walks each file once, routing rows to the
        # appropriate per-fold reservoir.
        t_sample_start = time()
        X0, y0, X1, y1 = _sample_both_folds(
            file_paths,
            available,
            k_per_fold,
            n0,
            n1,
            rng,
            training_masks,
        )
        t_sample = time() - t_sample_start

        t_fit_start = time()
        cls_trained_on_fold0 = _fit_pass1_booster(X0, y0, lgbm_hp)
        cls_trained_on_fold1 = _fit_pass1_booster(X1, y1, lgbm_hp)
        t_fit = time() - t_fit_start

        n_train_targets = count(y0) + count(y1)
        n_train_decoys = length(y0) + length(y1) - n_train_targets

        # Free the sample buffers before per-file predict — peak memory drops
        # from two sample matrices down to one file's matrix at a time.
        X0 = Matrix{Float32}(undef, 0, 0)
        X1 = Matrix{Float32}(undef, 0, 0)
        y0 = Bool[]; y1 = Bool[]
        GC.gc()

        cls_trained_on = (fold0 = cls_trained_on_fold0, fold1 = cls_trained_on_fold1)
        t_predict_start = time()
        _predict_pass1_files(
            file_paths,
            available,
            cls_trained_on,
            false;
            write_sidecars = false,
            scores_out = scores_buffer,
            targets_out = iter_idx == 1 ? targets_buffer : nothing,
            file_offsets = file_offsets,
        )
        t_predict = time() - t_predict_start

        total_sample += t_sample
        total_fit += t_fit
        total_predict += t_predict

        metrics = _scoring_semisupervised_metrics_and_mask(
            scores_buffer,
            targets_buffer;
            train_q_threshold = semisupervised_train_q_threshold,
            stop_q_threshold = semisupervised_stop_q_threshold,
        )
        @debug_l1 "  ScoringSearch semi-supervised iter $iter_idx (Pass-1 OOM): " *
                   "train targets=$n_train_targets decoys=$n_train_decoys; " *
                   "fit_rows=$n_fit sample_fold0=$k0 sample_fold1=$k1; " *
                   "q≤.01 targets=$(metrics.target_q01) decoys=$(metrics.decoy_q01)"

        importance_classifier = _pass1_importance_classifier(cls_trained_on_fold0)
        if importance_classifier === nothing
            importance_classifier = _pass1_importance_classifier(cls_trained_on_fold1)
        end

        return (
            iter = iter_idx,
            cls_trained_on = cls_trained_on,
            last_classifier = importance_classifier,
            metrics = metrics,
            target_q01 = metrics.target_q01,
            decoy_q01 = metrics.decoy_q01,
            n_train_targets = n_train_targets,
            n_train_decoys = n_train_decoys,
        )
    end

    best_state = nothing
    training_masks = nothing
    previous_target_q01 = -1
    max_iterations = semisupervised ? semisupervised_max_iterations : 1
    for iter_idx in 1:max_iterations
        state = _run_iteration(iter_idx, training_masks)
        candidate_state = (
            iter = state.iter,
            cls_trained_on = state.cls_trained_on,
            last_classifier = state.last_classifier,
            target_q01 = state.target_q01,
            decoy_q01 = state.decoy_q01,
            n_train_targets = state.n_train_targets,
            n_train_decoys = state.n_train_decoys,
        )

        if semisupervised && iter_idx > 1 && !_scoring_target_gain_sufficient(
            previous_target_q01,
            state.target_q01;
            min_fraction = semisupervised_min_gain,
        )
            best_state = _scoring_better_iteration_state(best_state, candidate_state)
            @debug_l1 "  ScoringSearch semi-supervised stopping (Pass-1 OOM): " *
                       "iter $iter_idx q≤.01 targets=$(state.target_q01) did not improve by " *
                       "$(round(100 * semisupervised_min_gain, digits=2))% over $previous_target_q01; " *
                       "using iter $(best_state.iter) with q≤.01 targets=$(best_state.target_q01)"
            break
        end

        best_state = _scoring_better_iteration_state(best_state, candidate_state)
        if !semisupervised
            break
        elseif iter_idx == max_iterations
            @debug_l1 "  ScoringSearch semi-supervised stopping (Pass-1 OOM): " *
                       "hit max iterations $max_iterations; using iter $(best_state.iter) " *
                       "with q≤.01 targets=$(best_state.target_q01)"
            break
        end

        previous_target_q01 = state.target_q01
        training_masks = _split_scoring_train_masks(
            row_counts,
            state.metrics.training_mask,
        )
    end

    # Pass 3: per-file predict + sidecar write for the selected iteration.
    # Each file's matrix-build + LightGBM.predict + Arrow.write are
    # independent (read-only on the shared booster handles, disjoint output
    # paths). LightGBM's C API documents `LGBM_BoosterPredictForMat` as
    # thread-safe so long as no thread is updating the booster (we aren't).
    t_predict_start = time()
    _predict_pass1_files(
        file_paths,
        available,
        best_state.cls_trained_on,
        compute_infold;
        write_sidecars = true,
    )
    t_predict = time() - t_predict_start
    total_predict += t_predict

    t_total = time() - t_total_start
    @debug_l1 "Pass-1 OOM elapsed: total=$(round(t_total, digits = 2))s " *
               "(sample=$(round(total_sample, digits = 2))s, fit=$(round(total_fit, digits = 2))s, " *
               "predict=$(round(total_predict, digits = 2))s)"

    return (
        n_total            = n_total,
        n_fold0            = total_n0,
        n_fold1            = total_n1,
        k_per_fold         = k_per_fold,
        last_classifier    = best_state.last_classifier,   # for importance reporting
        cls_trained_on     = best_state.cls_trained_on,
        available_features = available,
        semisupervised_iter = best_state.iter,
        target_q01         = best_state.target_q01,
        decoy_q01          = best_state.decoy_q01,
        elapsed_s          = t_total,
    )
end
