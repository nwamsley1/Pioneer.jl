# Pass-1 LightGBM training, out-of-memory variant.
#
# Replaces the previous load-everything-into-best_psms approach with three
# streaming passes over the per-file Arrow tables:
#
#   Pass 1 (metadata): Count rows per cv_fold; pick the feature list. Reads
#                      only Arrow metadata, no column data.
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

# Suffix for the per-file Pass-1 OOF-score sidecar this trainer writes (read by score_psms.jl and
# two_round_scoring.jl). Previously lived in the now-deleted mbr_streaming.jl.
const PASS1_SIDECAR_SUFFIX = ".pass1_sidecar.arrow"

# Per-fold sampling cap. Resolved at call time (SCORING_LGBM_MAX_TRAIN isn't
# bound at parse time — MainSearch/scoring.jl loads AFTER this file in
# importScripts.jl).
default_pass1_oom_k_per_fold() = SCORING_LGBM_MAX_TRAIN

# Pass 1: count rows per cv_fold across all files. Returns (n_fold0, n_fold1).
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
# arrows produced by MainSearch share the same columns).
function _resolve_available_features(file_paths::Vector{String}, requested::Vector{Symbol})
    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        length(tbl.precursor_idx) == 0 && continue
        return filter(f -> hasproperty(tbl, f), requested)
    end
    return Symbol[]
end

# Pass 2: single-pass reservoir-sample both folds simultaneously.
# Walks each file ONCE (vs the previous two-pass version that walked every
# file twice — once per fold), routing each row to its fold-specific
# reservoir based on `cv_fold`. When `k_per_fold >= n_avail_in_fold` (the
# default regime — k_per_fold = 2.5M, typical n_avail = 0.5-2M), every
# row is accepted, no random draws happen, and the function is purely an
# I/O + feature-gather pass. When `k_per_fold < n_avail_in_fold`, standard
# reservoir replace-with-prob k/seen kicks in for that fold.
#
# Called ONCE per `train_and_predict_pass1_oom!` over the FULL population; the
# semi-supervised loop then re-uses this fixed sample by masking rows rather than
# re-gathering (see `_masked_sample_rows`). Rows land in [file, row] order within each
# fold, which is the order the semi-supervised mask is laid out in.
function _sample_both_folds(
    file_paths::Vector{String},
    features::Vector{Symbol},
    k_per_fold::Int,
    n0::Int,
    n1::Int,
    rng::AbstractRNG;
    inject_shadows::Bool = false,
    graft_cols = (),
)
    nfeat = length(features)
    n_files = length(file_paths)

    # Shadow injection: each sampled real row contributes one paired shadow, so the
    # per-fold budget splits evenly and the training set still tops out at k_per_fold
    # (matching the on-disk injection, whose reservoir drew k rows from an already
    # doubled population). A file with only one class cannot be paired at all; its rows
    # are excluded from the budget so no matrix slot is left unfilled.
    pairable = trues(n_files)
    n0_eff = n0
    n1_eff = n1
    if inject_shadows
        n0_eff = 0
        n1_eff = 0
        for (file_idx, fpath) in enumerate(file_paths)
            tbl = Arrow.Table(fpath)
            m = length(tbl.cv_fold)
            if m == 0
                pairable[file_idx] = false
                continue
            end
            tgt = tbl.target
            has_t = false
            has_d = false
            @inbounds for i in 1:m
                Bool(tgt[i]) ? (has_t = true) : (has_d = true)
                has_t && has_d && break
            end
            if !(has_t && has_d)
                pairable[file_idx] = false
                continue
            end
            fold_c = tbl.cv_fold
            @inbounds for i in 1:m
                UInt8(fold_c[i]) == UInt8(0) ? (n0_eff += 1) : (n1_eff += 1)
            end
        end
    end

    k_real = inject_shadows ? max(1, k_per_fold ÷ 2) : k_per_fold
    k0 = min(n0_eff, k_real)
    k1 = min(n1_eff, k_real)
    rows0 = inject_shadows ? 2 * k0 : k0
    rows1 = inject_shadows ? 2 * k1 : k1
    X0 = Matrix{Float32}(undef, rows0, nfeat)
    X1 = Matrix{Float32}(undef, rows1, nfeat)
    y0 = Vector{Bool}(undef, rows0)
    y1 = Vector{Bool}(undef, rows1)
    is_shadow0 = falses(rows0)
    is_shadow1 = falses(rows1)
    is_graft = Bool[f in graft_cols for f in features]

    # Phase A (serial): run the reservoir-sampling decision logic to build a
    # per-file plan of (src_row_in_file => dst_row_in_matrix) pairs. We also
    # fill y0/y1 here. Keeping this serial preserves the existing RNG
    # sequence so IDs are bit-stable vs the prior column-major loop.
    #
    # A shadow of the real row occupying slot `p` lands at slot `k + p`, so an eviction
    # replaces the pair together and the layout is [reals ; shadows] per fold. Its
    # non-graft columns come from the nearest-iRT opposite-class row (`shadow_src`),
    # its graft columns from the original row (`shadow_graft`) — that graft is what
    # gives every cross-run feature a 1:1 target/decoy marginal.
    plan0 = Vector{Vector{Tuple{Int32, Int32}}}(undef, n_files)
    plan1 = Vector{Vector{Tuple{Int32, Int32}}}(undef, n_files)
    shadow_src0 = Vector{Vector{Tuple{Int32, Int32}}}(undef, n_files)
    shadow_src1 = Vector{Vector{Tuple{Int32, Int32}}}(undef, n_files)
    shadow_graft0 = Vector{Vector{Tuple{Int32, Int32}}}(undef, n_files)
    shadow_graft1 = Vector{Vector{Tuple{Int32, Int32}}}(undef, n_files)
    all_feat_cols = Vector{Vector{AbstractVector}}(undef, n_files)
    seen0 = 0
    seen1 = 0
    for (file_idx, fpath) in enumerate(file_paths)
        tbl = Arrow.Table(fpath)
        n = length(tbl.cv_fold)
        plan0[file_idx] = Tuple{Int32, Int32}[]
        plan1[file_idx] = Tuple{Int32, Int32}[]
        shadow_src0[file_idx] = Tuple{Int32, Int32}[]
        shadow_src1[file_idx] = Tuple{Int32, Int32}[]
        shadow_graft0[file_idx] = Tuple{Int32, Int32}[]
        shadow_graft1[file_idx] = Tuple{Int32, Int32}[]
        all_feat_cols[file_idx] = AbstractVector[getproperty(tbl, f) for f in features]
        n == 0 && continue
        inject_shadows && !pairable[file_idx] && continue
        fold_c = tbl.cv_fold
        target_c = tbl.target
        partner = inject_shadows ?
            _shadow_partners(Vector{Bool}(target_c), Float32.(tbl.irt_obs)) : nothing
        @inbounds for i in 1:n
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
                if inject_shadows
                    s = partner[i]
                    sp = Int32(k0 + pos)
                    push!(shadow_src0[file_idx], (s, sp))
                    push!(shadow_graft0[file_idx], (Int32(i), sp))
                    y0[sp] = Bool(target_c[s])
                    is_shadow0[sp] = true
                end
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
                if inject_shadows
                    s = partner[i]
                    sp = Int32(k1 + pos)
                    push!(shadow_src1[file_idx], (s, sp))
                    push!(shadow_graft1[file_idx], (Int32(i), sp))
                    y1[sp] = Bool(target_c[s])
                    is_shadow1[sp] = true
                end
            end
        end
    end
    @assert seen0 == n0_eff "Sampler saw $seen0 rows for fold=0 but metadata counted $n0_eff"
    @assert seen1 == n1_eff "Sampler saw $seen1 rows for fold=1 but metadata counted $n1_eff"

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
            col = all_feat_cols[file_idx][j]
            _copy_sampled_column!(X0, X1, col, plan0[file_idx], plan1[file_idx], j)
            if inject_shadows
                if is_graft[j]
                    _copy_sampled_column!(X0, X1, col,
                                          shadow_graft0[file_idx], shadow_graft1[file_idx], j)
                else
                    _copy_sampled_column!(X0, X1, col,
                                          shadow_src0[file_idx], shadow_src1[file_idx], j)
                end
            end
        end
    end

    # Count the shadows actually occupying slots, not the number of pushes: under reservoir
    # rejection a slot is rewritten each time its real row is evicted and replaced.
    inject_shadows && @user_info "shadow-decoy injection (in-memory): " *
        "$(count(is_shadow0) + count(is_shadow1)) shadows paired into the training sample " *
        "($(count(!, pairable)) file(s) skipped for missing opposite class)"
    return X0, y0, is_shadow0, X1, y1, is_shadow1
end

# Nearest-iRT opposite-class partner for every row of one file. `partner[r]` is the row
# a shadow of `r` copies its non-graft columns (and its label) from. Returns all-zero
# when the file has only one class, in which case it cannot be paired at all.
function _shadow_partners(targets::Vector{Bool}, irt::Vector{Float32})
    n = length(targets)
    partner = zeros(Int32, n)
    target_idx = findall(targets)
    decoy_idx = findall(!, targets)
    (isempty(target_idx) || isempty(decoy_idx)) && return partner
    target_perm = sortperm(view(irt, target_idx))
    decoy_perm = sortperm(view(irt, decoy_idx))
    target_sorted = irt[target_idx[target_perm]]
    decoy_sorted = irt[decoy_idx[decoy_perm]]
    @inbounds for r in 1:n
        if targets[r]
            pos = _nearest_sorted_idx(decoy_sorted, irt[r])
            partner[r] = Int32(decoy_idx[decoy_perm[pos]])
        else
            pos = _nearest_sorted_idx(target_sorted, irt[r])
            partner[r] = Int32(target_idx[target_perm[pos]])
        end
    end
    return partner
end

# Rows of one fold's fixed sample that pass `mask`. The semi-supervised mask is
# sample-local and laid out [fold0 ; fold1], so fold 0 reads it at `offset = 0` and
# fold 1 at `offset = length(fold-0 sample)`. `nothing` (iteration 1, or the
# non-semi-supervised case) means "use the whole sample".
function _masked_sample_rows(
    mask::Union{Nothing, AbstractVector{Bool}},
    offset::Int,
    n::Int,
)
    mask === nothing && return nothing
    rows = Vector{Int32}(undef, n)
    m = 0
    @inbounds for d in 1:n
        if mask[offset + d]
            m += 1
            rows[m] = Int32(d)
        end
    end
    resize!(rows, m)
    return rows
end

# Score a feature matrix with one booster (or the constant sentinel a degenerate
# single-class fold produces). Clamped to match the per-file predict path.
function _predict_matrix(cls, X::Matrix{Float32})
    size(X, 1) == 0 && return Float32[]
    if cls isa NamedTuple && cls.kind == :constant
        return fill(clamp(Float32(cls.value), 1f-6, 1f0 - 1f-4), size(X, 1))
    end
    raw = LightGBM.predict(cls, X; num_threads = Threads.nthreads())
    s = ndims(raw) == 2 ? dropdims(raw; dims = 2) : raw
    return clamp.(Float32.(s), 1f-6, 1f0 - 1f-4)
end

# OOF-score the fixed sample in memory: fold-0 rows by the fold-1 booster and vice
# versa. Returned in [fold0 ; fold1] order so it aligns with the sample's target
# vector and with `_masked_sample_rows`' offsets.
function _predict_sample_oof(
    X0::Matrix{Float32}, X1::Matrix{Float32}, cls_trained_on::NamedTuple,
)
    return vcat(
        _predict_matrix(cls_trained_on.fold1, X0),
        _predict_matrix(cls_trained_on.fold0, X1),
    )
end

# Gather the masked rows of the fixed sample into a fresh training matrix.
# `rows === nothing` returns the full sample without copying.
function _subset_sample(
    X::Matrix{Float32}, y::Vector{Bool}, rows::Union{Nothing, Vector{Int32}},
)
    rows === nothing && return X, y
    m = length(rows)
    nfeat = size(X, 2)
    Xs = Matrix{Float32}(undef, m, nfeat)
    ys = Vector{Bool}(undef, m)
    # Column-major: outer loop over columns keeps both reads and writes sequential.
    @inbounds for j in 1:nfeat
        @simd for k in 1:m
            Xs[k, j] = X[rows[k], j]
        end
    end
    @inbounds for k in 1:m
        ys[k] = y[rows[k]]
    end
    return Xs, ys
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
# Sidecar layout matches what mbr_streaming.jl's downstream readers expect:
#   (:precursor_idx, :scan_idx, :trace_prob_prepass, :trace_prob_infold)
#
# Scores are clamped to [1e-6, 1 - 1e-4] (matches the in-memory pipeline's
# clamp at score_psms.jl).
function _predict_pass1_to_sidecar(
    fpath::String, features::Vector{Symbol},
    cls_trained_on::NamedTuple,
    compute_infold::Bool,
)
    tbl = Arrow.Table(fpath)
    n = length(tbl.precursor_idx)
    if n == 0
        # Empty file — write an empty sidecar to keep downstream invariants.
        empty_df = DataFrame(
            precursor_idx       = UInt32[],
            scan_idx            = UInt32[],
            trace_prob_prepass  = Float32[],
            trace_prob_infold   = Float32[],
        )
        writeArrow(fpath * PASS1_SIDECAR_SUFFIX, empty_df)
        return nothing
    end

    fold_rows = _pass1_fold_row_indices(tbl.cv_fold)
    idx0 = fold_rows.fold0
    idx1 = fold_rows.fold1
    X0, X1 = _pass1_fold_feature_matrices(tbl, features, idx0, idx1)

    # Design A: Pass-3 runs SEQUENTIALLY on the main thread (see _predict_pass1_files
    # below), so `_predict_matrix` can safely use all cores. libomp is only ever
    # entered from the main thread (never a worker), so no lock and no corruption.
    # OOF: opposite-fold booster on each subset.
    oof = Vector{Float32}(undef, n)
    @inbounds oof[idx0] .= _predict_matrix(cls_trained_on.fold1, X0)
    @inbounds oof[idx1] .= _predict_matrix(cls_trained_on.fold0, X1)

    # In-fold: same-fold booster on each subset.
    infold = if compute_infold
        v = Vector{Float32}(undef, n)
        @inbounds v[idx0] .= _predict_matrix(cls_trained_on.fold0, X0)
        @inbounds v[idx1] .= _predict_matrix(cls_trained_on.fold1, X1)
        v
    else
        fill(NaN32, n)
    end

    side_df = DataFrame(
        precursor_idx       = collect(UInt32.(tbl.precursor_idx)),
        scan_idx            = collect(UInt32.(tbl.scan_idx)),
        trace_prob_prepass  = oof,
        trace_prob_infold   = infold,
    )
    writeArrow(fpath * PASS1_SIDECAR_SUFFIX, side_df)
    return nothing
end

function _predict_pass1_files(
    file_paths::Vector{String},
    features::Vector{Symbol},
    cls_trained_on::NamedTuple,
    compute_infold::Bool,
)
    # Design A: sequential on the main thread. The dominant cost (the LightGBM
    # predict) is recovered via multi-threaded predict per file (_predict_matrix),
    # which is only safe because it runs on the main thread here. Matrix build +
    # sidecar write lose across-file parallelism (acceptable if predict-bound).
    for fpath in file_paths
        _predict_pass1_to_sidecar(fpath, features, cls_trained_on, compute_infold)
    end
    return nothing
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
    inject_shadows::Bool = false,
    graft_cols = (),
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

    total_sample = 0.0
    total_fit = 0.0
    total_predict = 0.0

    # Pass 2, ONCE: reservoir-sample the FULL population. Semi-supervised iterations
    # re-use this fixed sample (masking rows) instead of re-gathering it per iteration.
    #
    # When no reservoir rejection happens (k_per_fold >= rows-in-fold), masking the full
    # sample yields exactly the rows a mask-then-sample pass would have gathered, in the
    # same order — the training matrices are identical (verified bit-identical on KEAP1:
    # same fit_rows and q<=.01 counts every iteration, same 708,388 precursors / 42,213
    # protein groups).
    #
    # Once the reservoir DOES reject, the two diverge: sampling then masking trains on the
    # masked subset of k rows rather than a fresh k drawn from the masked population, so
    # later iterations see fewer rows. Measured on KEAP1 with k_per_fold forced to 250k
    # (a ~4-9x subsample, harsher than production): iterations 2+ trained on ~165k instead
    # of a refilled 250k, costing 0.36% of precursors and 0.48% of protein groups. Accepted
    # deliberately — one code path, and re-gathering is the cost that grows with experiment
    # size while the fit stays capped.
    t_sample_start = time()
    X0_all, y0_all, is_shadow0, X1_all, y1_all, is_shadow1 = _sample_both_folds(
        file_paths,
        available,
        k_per_fold,
        total_n0,
        total_n1,
        rng;
        inject_shadows = inject_shadows,
        graft_cols = graft_cols,
    )
    total_sample += time() - t_sample_start

    # The semi-supervised q-value sweep runs over the SAMPLE, not the full out-of-memory
    # population: the sample is what training consumes, and OOF-scoring it in memory costs
    # one matrix predict instead of re-reading and re-scoring every file each iteration.
    # When the reservoir accepts everything (k_per_fold >= rows-in-fold) the sample IS the
    # full population, so the q-values are unchanged. Only the final Pass-3 predict touches
    # files.
    n_sample0 = length(y0_all)
    n_sample1 = length(y1_all)
    sample_targets = vcat(y0_all, y1_all)
    sample_is_shadow = vcat(is_shadow0, is_shadow1)

    function _run_iteration(iter_idx::Int, training_mask)
        rows0 = _masked_sample_rows(training_mask, 0, n_sample0)
        rows1 = _masked_sample_rows(training_mask, n_sample0, n_sample1)
        k0 = rows0 === nothing ? n_sample0 : length(rows0)
        k1 = rows1 === nothing ? n_sample1 : length(rows1)
        n_fit = k0 + k1
        n_fit == 0 && error("train_and_predict_pass1_oom!: no rows left for semi-supervised training")

        # Fit one fold at a time so at most one subset copy is live alongside the
        # retained full sample.
        t_fit_start = time()
        X0, y0 = _subset_sample(X0_all, y0_all, rows0)
        cls_trained_on_fold0 = _fit_pass1_booster(X0, y0, lgbm_hp)
        n_train_targets = count(y0)
        n_train_decoys = length(y0) - n_train_targets
        rows0 === nothing || (X0 = Matrix{Float32}(undef, 0, 0); y0 = Bool[])
        X1, y1 = _subset_sample(X1_all, y1_all, rows1)
        cls_trained_on_fold1 = _fit_pass1_booster(X1, y1, lgbm_hp)
        n_targets1 = count(y1)
        n_train_targets += n_targets1
        n_train_decoys += length(y1) - n_targets1
        rows1 === nothing || (X1 = Matrix{Float32}(undef, 0, 0); y1 = Bool[])
        t_fit = time() - t_fit_start

        cls_trained_on = (fold0 = cls_trained_on_fold0, fold1 = cls_trained_on_fold1)
        t_predict_start = time()
        sample_scores = _predict_sample_oof(X0_all, X1_all, cls_trained_on)
        t_predict = time() - t_predict_start

        total_fit += t_fit
        total_predict += t_predict

        metrics = _scoring_semisupervised_metrics_and_mask(
            sample_scores,
            sample_targets;
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
            sample_scores = sample_scores,
            target_q01 = metrics.target_q01,
            decoy_q01 = metrics.decoy_q01,
            n_train_targets = n_train_targets,
            n_train_decoys = n_train_decoys,
        )
    end

    best_state = nothing
    training_mask = nothing
    previous_target_q01 = -1
    max_iterations = semisupervised ? semisupervised_max_iterations : 1
    for iter_idx in 1:max_iterations
        state = _run_iteration(iter_idx, training_mask)
        candidate_state = (
            iter = state.iter,
            cls_trained_on = state.cls_trained_on,
            last_classifier = state.last_classifier,
            sample_scores = state.sample_scores,
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
        training_mask = state.metrics.training_mask
    end

    # Release the retained full sample before the final per-file predict.
    X0_all = Matrix{Float32}(undef, 0, 0); y0_all = Bool[]
    X1_all = Matrix{Float32}(undef, 0, 0); y1_all = Bool[]
    GC.gc()

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
        compute_infold,
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
        sample_scores      = best_state.sample_scores,
        sample_targets     = sample_targets,
        sample_is_shadow   = sample_is_shadow,
        target_q01         = best_state.target_q01,
        decoy_q01          = best_state.decoy_q01,
        elapsed_s          = t_total,
    )
end
