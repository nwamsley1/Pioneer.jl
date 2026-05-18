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

# Per-fold sampling cap. Default matches MainSearch/scoring.jl's
# SHARED_LGBM_MAX_TRAIN — the in-memory trainer already subsamples each
# fold to this size before fitting, so the OOM reservoir is functionally
# equivalent. Resolved at call time (SHARED_LGBM_MAX_TRAIN isn't bound at
# parse time — MainSearch/scoring.jl loads AFTER this file in
# importScripts.jl).
default_pass1_oom_k_per_fold() = SHARED_LGBM_MAX_TRAIN

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

# Pass 2: reservoir-sample up to k rows for the given cv_fold into a
# preallocated feature matrix + target vector. Single linear scan of the
# files; feature gather is skipped for rows that the reservoir rejects.
function _reservoir_sample_fold(
    file_paths::Vector{String},
    features::Vector{Symbol},
    target_fold::UInt8,
    k::Int,
    n_avail_in_fold::Int,
    rng::AbstractRNG,
)
    nfeat = length(features)
    k_actual = min(k, n_avail_in_fold)
    X = Matrix{Float32}(undef, k_actual, nfeat)
    y = Vector{Bool}(undef, k_actual)
    seen = 0  # in-fold row counter (reservoir index space)

    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        n = length(tbl.cv_fold)
        n == 0 && continue
        fold_c = tbl.cv_fold
        target_c = tbl.target
        feat_cols = ntuple(j -> getproperty(tbl, features[j]), nfeat)

        @inbounds for i in 1:n
            UInt8(fold_c[i]) == target_fold || continue
            seen += 1
            if seen <= k_actual
                # Reservoir fill phase: always accept.
                pos = seen
            else
                # Replace phase: replace at random position with prob k/seen.
                draw = rand(rng, 1:seen)
                draw > k_actual && continue
                pos = draw
            end
            for j in 1:nfeat
                X[pos, j] = Float32(feat_cols[j][i])
            end
            y[pos] = Bool(target_c[i])
        end
    end
    @assert seen == n_avail_in_fold "Reservoir saw $seen rows but Pass 1 counted $n_avail_in_fold for fold=$target_fold"
    return X, y
end

# Fit a LightGBM classifier on the sampled (X, y). Returns the booster.
function _fit_pass1_booster(
    X::Matrix{Float32}, y::Vector{Bool}, lgbm_hp::NamedTuple,
)
    cls = build_lightgbm_classifier(; lgbm_hp...)
    LightGBM.fit!(cls, X, _prepare_labels(y); verbosity = -1)
    return cls
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

    # Build the full file's feature matrix once (one mmap-touch per col).
    nfeat = length(features)
    X = Matrix{Float32}(undef, n, nfeat)
    feat_cols = ntuple(j -> getproperty(tbl, features[j]), nfeat)
    @inbounds for j in 1:nfeat, i in 1:n
        X[i, j] = Float32(feat_cols[j][i])
    end

    fold_c = tbl.cv_fold
    idx0 = Int[]; idx1 = Int[]
    @inbounds for i in 1:n
        f = UInt8(fold_c[i])
        f == UInt8(0) ? push!(idx0, i) :
        f == UInt8(1) ? push!(idx1, i) : nothing
    end

    @inline function _predict_block(cls, rows::Vector{Int})
        isempty(rows) && return Float32[]
        raw = LightGBM.predict(cls, X[rows, :])
        s = ndims(raw) == 2 ? dropdims(raw; dims = 2) : raw
        return clamp.(Float32.(s), 1f-6, 1f0 - 1f-4)
    end

    # OOF: opposite-fold booster on each subset.
    oof = Vector{Float32}(undef, n)
    @inbounds oof[idx0] .= _predict_block(cls_trained_on.fold1, idx0)
    @inbounds oof[idx1] .= _predict_block(cls_trained_on.fold0, idx1)

    # In-fold: same-fold booster on each subset (NaN sentinel when not requested).
    infold = if compute_infold
        v = Vector{Float32}(undef, n)
        @inbounds v[idx0] .= _predict_block(cls_trained_on.fold0, idx0)
        @inbounds v[idx1] .= _predict_block(cls_trained_on.fold1, idx1)
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
)
    t_total_start = time()

    # Pass 1: metadata.
    available = _resolve_available_features(file_paths, features)
    if isempty(available)
        error("train_and_predict_pass1_oom!: no requested features are present " *
              "in the per-file Arrow schema")
    end
    n0, n1 = _count_psm_rows_by_fold(file_paths)
    n_total = n0 + n1
    if n_total == 0
        @user_warn "train_and_predict_pass1_oom!: no PSM rows found across $(length(file_paths)) files"
        return (n_total = 0, last_classifier = nothing, available_features = available)
    end

    k0 = min(n0, k_per_fold)
    k1 = min(n1, k_per_fold)
    @user_info "Pass-1 OOM training:"
    @user_info "  rows: total=$n_total  fold0=$n0  fold1=$n1"
    @user_info "  sampling: k_per_fold=$k_per_fold → sample_fold0=$k0  sample_fold1=$k1"
    @user_info "  features in use: $(length(available)) / $(length(features))"

    # Pass 2: reservoir sample + train.
    t_sample_start = time()
    X0, y0 = _reservoir_sample_fold(file_paths, available, UInt8(0), k_per_fold, n0, rng)
    X1, y1 = _reservoir_sample_fold(file_paths, available, UInt8(1), k_per_fold, n1, rng)
    t_sample = time() - t_sample_start

    t_fit_start = time()
    cls_trained_on_fold0 = _fit_pass1_booster(X0, y0, lgbm_hp)
    cls_trained_on_fold1 = _fit_pass1_booster(X1, y1, lgbm_hp)
    t_fit = time() - t_fit_start

    # Free the sample buffers before per-file predict — peak memory drops
    # from two ~50 MB matrices down to one file's matrix at a time.
    X0 = Matrix{Float32}(undef, 0, 0)
    X1 = Matrix{Float32}(undef, 0, 0)
    y0 = Bool[]; y1 = Bool[]
    GC.gc()

    # Pass 3: per-file predict + sidecar write.
    cls_trained_on = (fold0 = cls_trained_on_fold0, fold1 = cls_trained_on_fold1)
    t_predict_start = time()
    for fpath in file_paths
        _predict_pass1_to_sidecar(fpath, available, cls_trained_on, compute_infold)
    end
    t_predict = time() - t_predict_start

    t_total = time() - t_total_start
    @user_info "Pass-1 OOM elapsed: total=$(round(t_total, digits = 2))s " *
               "(sample=$(round(t_sample, digits = 2))s, fit=$(round(t_fit, digits = 2))s, " *
               "predict=$(round(t_predict, digits = 2))s)"

    return (
        n_total            = n_total,
        n_fold0            = n0,
        n_fold1            = n1,
        k_per_fold         = k_per_fold,
        last_classifier    = cls_trained_on_fold0,   # for importance reporting
        available_features = available,
        elapsed_s          = t_total,
    )
end
