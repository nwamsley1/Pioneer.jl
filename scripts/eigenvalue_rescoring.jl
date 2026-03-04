#!/usr/bin/env julia
#
# eigenvalue_rescoring.jl — Self-contained script to test whether fragment correlation
# eigenvalue features improve LightGBM rescoring at the global level.
#
# Usage: julia scripts/eigenvalue_rescoring.jl

using Arrow, DataFrames, Statistics, LinearAlgebra, LightGBM, Random, Printf

###############################################################################
# Configuration
###############################################################################

const DATA_DIR = "/Users/nathanwamsley/Data/For_Figures/OlsenEclipse/OlsenEclipse_FixSinglePointChrom/temp_data"
const PSM_DIR  = joinpath(DATA_DIR, "first_pass_psms")
const CORR_DIR = joinpath(DATA_DIR, "first_pass_corr")

const STANDARD_FEATURES = [
    :scribe, :spectral_contrast, :city_block, :entropy_score,
    :poisson, :y_count, :err_norm, :scan_count,
]

const EIGEN_FEATURES = [
    :pearson_lambda1, :pearson_lambda2, :pearson_lambda3,
    :pearson_lambda1_frac, :pearson_lambda12_frac, :pearson_lambda13_frac,
    :spearman_lambda1, :spearman_lambda2, :spearman_lambda3,
    :spearman_lambda1_frac, :spearman_lambda12_frac, :spearman_lambda13_frac,
]

const NARROW_EIGEN_FEATURES = [
    :narrow_pearson_lambda1, :narrow_pearson_lambda2, :narrow_pearson_lambda3,
    :narrow_pearson_lambda1_frac, :narrow_pearson_lambda12_frac, :narrow_pearson_lambda13_frac,
    :narrow_spearman_lambda1, :narrow_spearman_lambda2, :narrow_spearman_lambda3,
    :narrow_spearman_lambda1_frac, :narrow_spearman_lambda12_frac, :narrow_spearman_lambda13_frac,
]

const ALL_EIGEN_FEATURES = vcat(EIGEN_FEATURES, NARROW_EIGEN_FEATURES)

const CHROM_FEATURES = [
    :sum_frag_corr, :min_frag_corr, :mean_frag_corr,
    :shape_symmetry, :shape_apex_frac, :shape_kurtosis,
    :rank_stability, :top1_dominance,
]

const NARROW_CHROM_FEATURES = [Symbol("narrow_" * string(f)) for f in CHROM_FEATURES]

const ALL_CHROM_FEATURES = vcat(CHROM_FEATURES, NARROW_CHROM_FEATURES)

const CHROM_TOP3_FEATURES = [:mean_frag_corr, :sum_frag_corr, :rank_stability]

const EIGEN_TOP4_FEATURES = [
    :spearman_lambda1, :spearman_lambda1_frac,
    :spearman_lambda12_frac, :spearman_lambda13_frac,
]

const AUG_EIGEN_TOP4_FEATURES = vcat(STANDARD_FEATURES, EIGEN_TOP4_FEATURES)
const AUG_FULL_FEATURES   = vcat(STANDARD_FEATURES, EIGEN_FEATURES)
const AUG_NARROW_FEATURES = vcat(STANDARD_FEATURES, NARROW_EIGEN_FEATURES)
const AUG_BOTH_FEATURES   = vcat(STANDARD_FEATURES, EIGEN_FEATURES, NARROW_EIGEN_FEATURES)
const AUG_CHROM_TOP3_FEATURES = vcat(STANDARD_FEATURES, CHROM_TOP3_FEATURES)
const AUG_CHROM_FEATURES  = vcat(STANDARD_FEATURES, CHROM_FEATURES, NARROW_CHROM_FEATURES)
const AUG_ALL_FEATURES    = vcat(STANDARD_FEATURES, EIGEN_FEATURES, NARROW_EIGEN_FEATURES,
                                  CHROM_FEATURES, NARROW_CHROM_FEATURES)

const NUM_ITERATIONS_PASS1 = 100
const NUM_ITERATIONS_PASS2 = 200
const Q_VALUE_CUTOFF       = 0.01

###############################################################################
# 1. Load Data
###############################################################################

function load_data()
    psm_files  = filter(f -> endswith(f, ".arrow"), readdir(PSM_DIR))
    corr_files = filter(f -> endswith(f, ".arrow"), readdir(CORR_DIR))

    println("Loading $(length(psm_files)) PSM files and $(length(corr_files)) corr files...")

    psm_dfs  = DataFrame[]
    corr_dfs = DataFrame[]

    for (i, f) in enumerate(sort(psm_files))
        df = DataFrame(Arrow.Table(joinpath(PSM_DIR, f)))
        df[!, :ms_file_idx] .= UInt32(i)  # ensure consistent file index
        push!(psm_dfs, df)
    end

    for (i, f) in enumerate(sort(corr_files))
        df = DataFrame(Arrow.Table(joinpath(CORR_DIR, f)))
        df[!, :ms_file_idx] .= UInt32(i)
        push!(corr_dfs, df)
    end

    psm_all  = vcat(psm_dfs...)
    corr_all = vcat(corr_dfs...)

    # Convert narrow float types to Float64 for ML
    for col in names(psm_all)
        T = eltype(psm_all[!, col])
        if T <: Union{Float16, Float32}
            psm_all[!, col] = Float64.(psm_all[!, col])
        end
    end

    n_files = length(psm_files)
    println("  PSMs:  $(nrow(psm_all)) rows across $n_files files")
    println("  Corr:  $(nrow(corr_all)) rows")

    return psm_all, corr_all, n_files
end

###############################################################################
# 2. Compute Eigenvalue Features
###############################################################################

function pearson_corrmat(X::Matrix{Float64}; eps=1e-12)
    n, p = size(X)
    means = vec(mean(X; dims=1))
    stds  = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        ss = 0.0
        for i in 1:n
            ss += (X[i, j] - means[j])^2
        end
        stds[j] = sqrt(ss / max(n - 1, 1))
    end

    C = Matrix{Float64}(undef, p, p)
    @inbounds for j1 in 1:p
        C[j1, j1] = 1.0
        for j2 in (j1+1):p
            if stds[j1] > eps && stds[j2] > eps
                cov_val = 0.0
                for i in 1:n
                    cov_val += (X[i, j1] - means[j1]) * (X[i, j2] - means[j2])
                end
                c = cov_val / (max(n - 1, 1) * stds[j1] * stds[j2])
                C[j1, j2] = c
                C[j2, j1] = c
            else
                C[j1, j2] = 0.0
                C[j2, j1] = 0.0
            end
        end
    end
    return C
end

function rank_matrix(X::Matrix{Float64})
    n, p = size(X)
    R = Matrix{Float64}(undef, n, p)
    @inbounds for j in 1:p
        col = view(X, :, j)
        ord = sortperm(col)
        i = 1
        while i <= n
            val = col[ord[i]]
            k = i
            while k <= n && col[ord[k]] == val
                k += 1
            end
            avg_rank = (i + k - 1) / 2.0
            for t in i:(k-1)
                R[ord[t], j] = avg_rank
            end
            i = k
        end
    end
    return R
end

function eigenvalue_features(C::Matrix{Float64})
    vals = try
        sort(eigvals(Symmetric(C)); rev=true)
    catch
        return nothing
    end
    # Clamp negatives to zero (numerical noise)
    vals = max.(vals, 0.0)
    total = sum(vals)
    if total < 1e-12 || length(vals) < 3
        return nothing
    end
    λ1, λ2, λ3 = vals[1], vals[2], vals[3]
    return (
        lambda1      = λ1,
        lambda2      = λ2,
        lambda3      = λ3,
        lambda1_frac = λ1 / total,
        lambda12_frac = λ2 > 1e-12 ? λ1 / λ2 : λ1 * 1e6,
        lambda13_frac = λ3 > 1e-12 ? λ1 / λ3 : λ1 * 1e6,
    )
end

"""
    compute_eigenvalue_features(corr_all; window=nothing, prefix="")

Compute eigenvalue features from fragment intensity correlation matrices.

- `window=nothing`: use all scans (full chromatogram)
- `window=2`: use only ±2 scans from the apex (highest total intensity)
- `prefix`: prepended to column names (e.g. "narrow_")
"""
function compute_eigenvalue_features(corr_all::DataFrame;
                                     window::Union{Nothing,Int}=nothing,
                                     prefix::String="")
    win_desc = window === nothing ? "full window" : "±$window from apex"
    println("Computing eigenvalue features ($win_desc, prefix=\"$prefix\")...")
    intensity_cols = [:intensity_1, :intensity_2, :intensity_3,
                      :intensity_4, :intensity_5, :intensity_6]

    # Group by (precursor_idx, ms_file_idx)
    gdf = groupby(corr_all, [:precursor_idx, :ms_file_idx])

    n_groups = length(gdf)
    result_pidx   = Vector{UInt32}(undef, n_groups)
    result_fidx   = Vector{UInt32}(undef, n_groups)
    result_feats  = Vector{Union{Nothing, NamedTuple}}(undef, n_groups)

    Threads.@threads for gi in 1:n_groups
        g = gdf[gi]
        result_pidx[gi] = g.precursor_idx[1]
        result_fidx[gi] = g.ms_file_idx[1]

        n_scans = nrow(g)
        if n_scans < 3
            result_feats[gi] = nothing
            continue
        end

        # Build intensity matrix (all scans, sorted by RT)
        rt_order = sortperm(g.rt)
        X_full = Matrix{Float64}(undef, n_scans, 6)
        for (j, c) in enumerate(intensity_cols)
            col = g[!, c]
            for (ri, oi) in enumerate(rt_order)
                X_full[ri, j] = Float64(col[oi])
            end
        end

        # Apply window if requested: select ±window scans from apex
        if window !== nothing
            # Find apex: scan with highest total intensity
            total_int = vec(sum(X_full; dims=2))
            apex_pos = argmax(total_int)
            lo = max(1, apex_pos - window)
            hi = min(n_scans, apex_pos + window)
            X = X_full[lo:hi, :]
            if size(X, 1) < 3
                result_feats[gi] = nothing
                continue
            end
        else
            X = X_full
        end

        # Pearson
        Cp = pearson_corrmat(X)
        pf = eigenvalue_features(Cp)

        # Spearman
        Xr = rank_matrix(X)
        Cs = pearson_corrmat(Xr)
        sf = eigenvalue_features(Cs)

        if pf === nothing || sf === nothing
            result_feats[gi] = nothing
            continue
        end

        result_feats[gi] = (
            pearson_lambda1       = pf.lambda1,
            pearson_lambda2       = pf.lambda2,
            pearson_lambda3       = pf.lambda3,
            pearson_lambda1_frac  = pf.lambda1_frac,
            pearson_lambda12_frac = pf.lambda12_frac,
            pearson_lambda13_frac = pf.lambda13_frac,
            spearman_lambda1       = sf.lambda1,
            spearman_lambda2       = sf.lambda2,
            spearman_lambda3       = sf.lambda3,
            spearman_lambda1_frac  = sf.lambda1_frac,
            spearman_lambda12_frac = sf.lambda12_frac,
            spearman_lambda13_frac = sf.lambda13_frac,
        )
    end

    # Build DataFrame from non-nothing results
    valid = findall(!isnothing, result_feats)
    println("  $(length(valid)) / $n_groups groups with valid features (≥3 scans in window)")

    eigen_df = DataFrame(
        :precursor_idx                         => result_pidx[valid],
        :ms_file_idx                           => result_fidx[valid],
        Symbol(prefix * "pearson_lambda1")     => [result_feats[i].pearson_lambda1 for i in valid],
        Symbol(prefix * "pearson_lambda2")     => [result_feats[i].pearson_lambda2 for i in valid],
        Symbol(prefix * "pearson_lambda3")     => [result_feats[i].pearson_lambda3 for i in valid],
        Symbol(prefix * "pearson_lambda1_frac")  => [result_feats[i].pearson_lambda1_frac for i in valid],
        Symbol(prefix * "pearson_lambda12_frac") => [result_feats[i].pearson_lambda12_frac for i in valid],
        Symbol(prefix * "pearson_lambda13_frac") => [result_feats[i].pearson_lambda13_frac for i in valid],
        Symbol(prefix * "spearman_lambda1")    => [result_feats[i].spearman_lambda1 for i in valid],
        Symbol(prefix * "spearman_lambda2")    => [result_feats[i].spearman_lambda2 for i in valid],
        Symbol(prefix * "spearman_lambda3")    => [result_feats[i].spearman_lambda3 for i in valid],
        Symbol(prefix * "spearman_lambda1_frac")  => [result_feats[i].spearman_lambda1_frac for i in valid],
        Symbol(prefix * "spearman_lambda12_frac") => [result_feats[i].spearman_lambda12_frac for i in valid],
        Symbol(prefix * "spearman_lambda13_frac") => [result_feats[i].spearman_lambda13_frac for i in valid],
    )

    return eigen_df
end

###############################################################################
# 3. Compute Chromatogram Features (DIA-NN-inspired)
###############################################################################

"""
    compute_chromatogram_features(corr_all; window=nothing, prefix="")

Compute DIA-NN-inspired chromatogram scoring features from fragment intensity traces.

Features computed:
- Per-fragment elution correlation (sum/min/mean of Pearson r with total elution profile)
- Peak shape (symmetry, apex fraction, kurtosis)
- Fragment rank stability and top-1 dominance
"""
function compute_chromatogram_features(corr_all::DataFrame;
                                       window::Union{Nothing,Int}=nothing,
                                       prefix::String="")
    win_desc = window === nothing ? "full window" : "±$window from apex"
    println("Computing chromatogram features ($win_desc, prefix=\"$prefix\")...")

    intensity_cols = [:intensity_1, :intensity_2, :intensity_3,
                      :intensity_4, :intensity_5, :intensity_6]

    gdf = groupby(corr_all, [:precursor_idx, :ms_file_idx])
    n_groups = length(gdf)

    result_pidx  = Vector{UInt32}(undef, n_groups)
    result_fidx  = Vector{UInt32}(undef, n_groups)
    result_feats = Vector{Union{Nothing, NamedTuple}}(undef, n_groups)

    Threads.@threads for gi in 1:n_groups
        g = gdf[gi]
        result_pidx[gi] = g.precursor_idx[1]
        result_fidx[gi] = g.ms_file_idx[1]

        n_scans = nrow(g)
        if n_scans < 3
            result_feats[gi] = nothing
            continue
        end

        # Build intensity matrix sorted by RT
        rt_order = sortperm(g.rt)
        X = Matrix{Float64}(undef, n_scans, 6)
        for (j, c) in enumerate(intensity_cols)
            col = g[!, c]
            for (ri, oi) in enumerate(rt_order)
                X[ri, j] = Float64(col[oi])
            end
        end

        # Find apex on full chromatogram, then apply window
        elution_full = vec(sum(X; dims=2))
        apex_pos = argmax(elution_full)

        if window !== nothing
            lo = max(1, apex_pos - window)
            hi = min(n_scans, apex_pos + window)
            X = X[lo:hi, :]
            if size(X, 1) < 3
                result_feats[gi] = nothing
                continue
            end
            # Recompute apex position within windowed data
            elution = vec(sum(X; dims=2))
            apex_pos = argmax(elution)
        else
            elution = elution_full
        end

        n = size(X, 1)

        # ── Group 1: Per-fragment elution correlation ──
        # For each fragment, Pearson correlation with the total elution profile
        mean_elution = mean(elution)
        std_elution = std(elution; corrected=true)

        frag_corrs = Vector{Float64}(undef, 6)
        if std_elution > 1e-12
            for j in 1:6
                frag_j = @view X[:, j]
                mean_j = mean(frag_j)
                std_j = std(frag_j; corrected=true)
                if std_j > 1e-12
                    cov_val = 0.0
                    @inbounds for i in 1:n
                        cov_val += (frag_j[i] - mean_j) * (elution[i] - mean_elution)
                    end
                    frag_corrs[j] = cov_val / ((n - 1) * std_j * std_elution)
                else
                    frag_corrs[j] = 0.0
                end
            end
        else
            fill!(frag_corrs, 0.0)
        end

        sum_frag_corr  = sum(frag_corrs)
        min_frag_corr  = minimum(frag_corrs)
        mean_frag_corr = mean(frag_corrs)

        # ── Group 2: Peak shape features ──
        total_intensity = sum(elution)
        if total_intensity < 1e-12
            result_feats[gi] = nothing
            continue
        end

        # Shape symmetry: left-half vs right-half intensity
        left_sum  = sum(@view elution[1:apex_pos-1]; init=0.0)
        right_sum = sum(@view elution[apex_pos+1:n]; init=0.0)
        denom = left_sum + right_sum
        shape_symmetry = denom > 1e-12 ? min(left_sum, right_sum) / max(left_sum, right_sum) : 1.0

        # Apex fraction: what fraction of total intensity is in the apex scan
        shape_apex_frac = elution[apex_pos] / total_intensity

        # Kurtosis of the elution profile (excess kurtosis)
        mean_e = mean(elution)
        var_e  = 0.0
        m4_e   = 0.0
        @inbounds for i in 1:n
            d = elution[i] - mean_e
            d2 = d * d
            var_e += d2
            m4_e  += d2 * d2
        end
        var_e /= n
        m4_e  /= n
        shape_kurtosis = var_e > 1e-12 ? (m4_e / (var_e * var_e)) - 3.0 : 0.0

        # ── Group 3: Fragment rank stability ──
        # Fraction of scans where intensity_1 > intensity_2 > ... > intensity_6
        n_ordered = 0
        @inbounds for i in 1:n
            ordered = true
            for j in 1:5
                if X[i, j] < X[i, j+1]
                    ordered = false
                    break
                end
            end
            n_ordered += ordered
        end
        rank_stability = n_ordered / n

        # Top-1 dominance: mean fraction of total intensity from fragment 1
        top1_sum = 0.0
        @inbounds for i in 1:n
            if elution[i] > 1e-12
                top1_sum += X[i, 1] / elution[i]
            end
        end
        top1_dominance = top1_sum / n

        result_feats[gi] = (
            sum_frag_corr  = sum_frag_corr,
            min_frag_corr  = min_frag_corr,
            mean_frag_corr = mean_frag_corr,
            shape_symmetry  = shape_symmetry,
            shape_apex_frac = shape_apex_frac,
            shape_kurtosis  = shape_kurtosis,
            rank_stability  = rank_stability,
            top1_dominance  = top1_dominance,
        )
    end

    # Build DataFrame from non-nothing results
    valid = findall(!isnothing, result_feats)
    println("  $(length(valid)) / $n_groups groups with valid features (≥3 scans in window)")

    chrom_df = DataFrame(
        :precursor_idx                      => result_pidx[valid],
        :ms_file_idx                        => result_fidx[valid],
        Symbol(prefix * "sum_frag_corr")    => [result_feats[i].sum_frag_corr for i in valid],
        Symbol(prefix * "min_frag_corr")    => [result_feats[i].min_frag_corr for i in valid],
        Symbol(prefix * "mean_frag_corr")   => [result_feats[i].mean_frag_corr for i in valid],
        Symbol(prefix * "shape_symmetry")   => [result_feats[i].shape_symmetry for i in valid],
        Symbol(prefix * "shape_apex_frac")  => [result_feats[i].shape_apex_frac for i in valid],
        Symbol(prefix * "shape_kurtosis")   => [result_feats[i].shape_kurtosis for i in valid],
        Symbol(prefix * "rank_stability")   => [result_feats[i].rank_stability for i in valid],
        Symbol(prefix * "top1_dominance")   => [result_feats[i].top1_dominance for i in valid],
    )

    return chrom_df
end

###############################################################################
# 4. Build Feature DataFrame
###############################################################################

function build_feature_df(psm_all::DataFrame, eigen_full::DataFrame, eigen_narrow::DataFrame,
                          chrom_full::DataFrame, chrom_narrow::DataFrame)
    println("Joining PSM data with eigenvalue and chromatogram features...")

    # Join eigenvalue features
    df = leftjoin(psm_all, eigen_full; on=[:precursor_idx, :ms_file_idx])
    df = leftjoin(df, eigen_narrow; on=[:precursor_idx, :ms_file_idx])

    # Join chromatogram features
    df = leftjoin(df, chrom_full; on=[:precursor_idx, :ms_file_idx])
    df = leftjoin(df, chrom_narrow; on=[:precursor_idx, :ms_file_idx])

    # Fill missing features with 0
    for col in vcat(ALL_EIGEN_FEATURES, ALL_CHROM_FEATURES)
        col_data = df[!, col]
        if eltype(col_data) >: Missing
            df[!, col] = coalesce.(col_data, 0.0)
        end
    end

    println("  $(nrow(df)) PSMs total")
    println("  $(nrow(eigen_full)) groups with full-window eigen features")
    println("  $(nrow(eigen_narrow)) groups with narrow-window eigen features")
    println("  $(nrow(chrom_full)) groups with full-window chrom features")
    println("  $(nrow(chrom_narrow)) groups with narrow-window chrom features")
    return df
end

###############################################################################
# 5. LightGBM Training with 2-Fold CV
###############################################################################

function build_classifier(; num_iterations::Int=100)
    return LightGBM.LGBMClassification(
        objective        = "binary",
        metric           = ["binary_logloss"],
        learning_rate    = 0.05,
        num_iterations   = num_iterations,
        num_leaves       = 63,
        max_depth        = -1,
        feature_fraction = 0.5,
        bagging_fraction = 0.25,
        bagging_freq     = 1,
        min_data_in_leaf = 200,
        min_gain_to_split = 0.0,
        lambda_l2        = 0.0,
        max_bin          = 255,
        num_threads      = Threads.nthreads(),
        num_class        = 1,
        verbosity        = -1,
        is_unbalance     = false,
        seed             = 1776,
        deterministic    = true,
        force_row_wise   = true,
    )
end

function make_feature_matrix(df::DataFrame, features::Vector{Symbol})
    n = nrow(df)
    m = length(features)
    X = Matrix{Float32}(undef, n, m)
    for (j, feat) in enumerate(features)
        X[:, j] = Float32.(df[!, feat])
    end
    return X
end

function get_qvalues(probs::AbstractVector, targets::AbstractVector{Bool})
    order = sortperm(probs; rev=true)
    qvals = Vector{Float64}(undef, length(probs))
    cum_t = 0
    cum_d = 0
    @inbounds for i in order
        cum_t += targets[i]
        cum_d += !targets[i]
        qvals[i] = cum_t > 0 ? cum_d / cum_t : 1.0
    end
    # Monotonize: walk in reverse order, cap at running minimum
    fdr = Inf
    @inbounds for i in reverse(order)
        if qvals[i] > fdr
            qvals[i] = fdr
        else
            fdr = qvals[i]
        end
    end
    return qvals
end

function train_lgbm_cv(df::DataFrame, features::Vector{Symbol}; label::String="model")
    println("  Training LightGBM ($label) with $(length(features)) features, two-pass approach...")
    n = nrow(df)
    m = length(features)

    # CV fold assignment: even/odd precursor_idx
    folds = UInt8.((df.precursor_idx .% 2) .+ 1)

    # Output probability vector (final out-of-fold predictions from pass 2)
    probs_out = zeros(Float64, n)

    # Intermediate in-sample predictions from pass 1 (used to select pass-2 training set)
    pass1_probs = zeros(Float64, n)

    # Accumulate feature importances across folds (from pass 2)
    gain_accum = zeros(Float64, m)
    n_folds_with_gains = 0

    targets = df.target

    # ── Pass 1: Train on all targets + all decoys, predict in-sample ──
    println("    Pass 1: initial training (all targets + all decoys)...")
    for fold in UInt8[1, 2]
        train_idx = findall(folds .!= fold)

        X_train = make_feature_matrix(df[train_idx, :], features)
        y_train = Int.(targets[train_idx])

        model = build_classifier(; num_iterations=NUM_ITERATIONS_PASS1)
        LightGBM.fit!(model, X_train, y_train; verbosity=-1)

        # Predict in-sample on training data for pass-1 q-value computation
        raw = LightGBM.predict(model, X_train)
        preds = ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
        pass1_probs[train_idx] = Float64.(preds)
    end

    # Compute q-values from pass-1 in-sample predictions across all data
    pass1_qvals = get_qvalues(pass1_probs, targets)
    n_pass1_targets = count((pass1_qvals .<= Q_VALUE_CUTOFF) .& targets)
    println("    Pass 1 targets at $(Q_VALUE_CUTOFF*100)% FDR: $n_pass1_targets")

    # ── Pass 2: Retrain on 1% FDR targets + all decoys, predict out-of-fold ──
    println("    Pass 2: refined training ($(Q_VALUE_CUTOFF*100)% FDR targets + all decoys)...")
    for fold in UInt8[1, 2]
        test_idx  = findall(folds .== fold)
        train_idx = findall(folds .!= fold)

        # Filter training set: keep targets with pass-1 q-value ≤ cutoff + all decoys
        train_targets = targets[train_idx]
        train_qvals   = pass1_qvals[train_idx]
        keep = .!train_targets .| (train_qvals .<= Q_VALUE_CUTOFF)
        sel = train_idx[keep]

        X_train = make_feature_matrix(df[sel, :], features)
        y_train = Int.(targets[sel])

        model = build_classifier(; num_iterations=NUM_ITERATIONS_PASS2)
        LightGBM.fit!(model, X_train, y_train; verbosity=-1)

        # Predict on held-out test fold — final out-of-fold scores
        X_test = make_feature_matrix(df[test_idx, :], features)
        raw = LightGBM.predict(model, X_test)
        preds = ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
        probs_out[test_idx] = Float64.(preds)

        # Extract feature importances from pass-2 models
        try
            gains = LightGBM.gain_importance(model)
            gain_accum .+= gains
            n_folds_with_gains += 1
        catch e
            println("  (Could not extract feature importances fold $fold: $e)")
        end
    end

    # Print averaged feature importances across both folds
    if n_folds_with_gains > 0
        avg_gains = gain_accum ./ n_folds_with_gains
        println("\n  Feature importances ($label, averaged over $n_folds_with_gains folds):")
        feat_imp = collect(zip(features, avg_gains))
        sort!(feat_imp; by=x -> -x[2])
        for (f, g) in feat_imp
            @printf("    %-36s  %.1f\n", f, g)
        end
    end

    return probs_out
end

###############################################################################
# 6. Global PEP Reconstruction
###############################################################################

function logodds_combine(probs::Vector{Float64}, top_n::Int)::Float64
    isempty(probs) && return 0.0
    n = min(length(probs), top_n)
    sorted = sort(probs; rev=true)
    selected = @view sorted[1:n]
    eps = 1e-6
    lo = log.(clamp.(selected, 0.1, 1 - eps) ./ (1 .- clamp.(selected, 0.1, 1 - eps)))
    avg = sum(lo) / n
    return 1.0 / (1 + exp(-avg))
end

function weighted_pava(y::Vector{Float64}, w::Vector{Float64})
    n = length(y)
    v    = Vector{Float64}(undef, n)
    wt   = Vector{Float64}(undef, n)
    lenv = Vector{Int}(undef, n)
    m = 0
    for i in 1:n
        m += 1
        @inbounds begin
            v[m]    = y[i]
            wt[m]   = w[i]
            lenv[m] = 1
        end
        while m > 1 && v[m-1] > v[m]
            @inbounds begin
                new_w   = wt[m-1] + wt[m]
                new_v   = (v[m-1]*wt[m-1] + v[m]*wt[m]) / new_w
                new_len = lenv[m-1] + lenv[m]
                m -= 1
                v[m]    = new_v
                wt[m]   = new_w
                lenv[m] = new_len
            end
        end
    end
    result = Vector{Float64}(undef, n)
    idx = 1
    for j in 1:m
        @inbounds for _ in 1:lenv[j]
            result[idx] = v[j]
            idx += 1
        end
    end
    return result
end

function compute_global_pep(scores::Vector{Float64}, targets::Vector{Bool})
    N = length(scores)
    order = sortperm(scores; rev=true)

    # Prepare labels + weights with pseudo-observation
    labels  = Vector{Float64}(undef, N + 1)
    weights = Vector{Float64}(undef, N + 1)
    labels[1]  = 0.0  # pseudo observation
    weights[1] = 1.0

    @inbounds for j in 1:N
        idx = order[j]
        labels[j+1]  = targets[idx] ? 0.0 : 1.0  # decoy=1, target=0
        weights[j+1] = 1.0
    end

    fitted = weighted_pava(labels, weights)
    fitted = fitted[2:end]  # remove pseudo row

    pep = fitted ./ (1 .- fitted)
    pep = clamp.(pep, 0.0, 1.0)

    result = Vector{Float64}(undef, N)
    @inbounds for (j, idx) in enumerate(order)
        result[idx] = pep[j]
    end
    return result
end

function compute_global_probs(df::DataFrame, score_col::Symbol, n_files::Int)
    top_n = max(1, floor(Int, sqrt(n_files)))

    # Best score per (precursor_idx, ms_file_idx)
    prec_file_best = Dict{UInt32, Dict{UInt32, Float64}}()
    prec_target    = Dict{UInt32, Bool}()

    for i in 1:nrow(df)
        pidx = UInt32(df.precursor_idx[i])
        fidx = UInt32(df.ms_file_idx[i])
        s    = df[i, score_col]

        if !haskey(prec_file_best, pidx)
            prec_file_best[pidx] = Dict{UInt32, Float64}(fidx => s)
            prec_target[pidx]    = df.target[i]
        else
            fb = prec_file_best[pidx]
            fb[fidx] = max(get(fb, fidx, -Inf), s)
        end
    end

    # Log-odds combine across files
    n_prec = length(prec_file_best)
    out_pidx   = Vector{UInt32}(undef, n_prec)
    out_target = Vector{Bool}(undef, n_prec)
    out_gprob  = Vector{Float64}(undef, n_prec)

    for (i, (pidx, file_dict)) in enumerate(prec_file_best)
        probs_vec = collect(values(file_dict))
        out_pidx[i]   = pidx
        out_target[i] = prec_target[pidx]
        out_gprob[i]  = logodds_combine(probs_vec, top_n)
    end

    return DataFrame(
        precursor_idx = out_pidx,
        target        = out_target,
        global_prob   = out_gprob,
    )
end

###############################################################################
# 7. Hybrid Filter + Comparison
###############################################################################

function hybrid_filter(gp_df::DataFrame;
                       qvalue_threshold::Float64=0.15,
                       pep_threshold::Float64=0.5,
                       min_floor::Int=50_000)

    # Q-values on global probs
    qvals = get_qvalues(gp_df.global_prob, gp_df.target)

    # Global PEP via isotonic regression
    global_pep = compute_global_pep(gp_df.global_prob, gp_df.target)

    # Counts at FDR thresholds
    n_01 = count((qvals .<= 0.01) .& gp_df.target)
    n_05 = count((qvals .<= 0.05) .& gp_df.target)

    # Hybrid filter: sort by global_prob descending
    order = sortperm(gp_df.global_prob; rev=true)
    cum_t = 0
    cum_d = 0
    n_qvalue_based = 0
    for (i, idx) in enumerate(order)
        if gp_df.target[idx]
            cum_t += 1
        else
            cum_d += 1
        end
        if cum_t > 0 && cum_d / cum_t <= qvalue_threshold
            n_qvalue_based = i
        end
    end

    n_passing_threshold = count(global_pep .<= pep_threshold)
    n_min_floor = min(min_floor, nrow(gp_df))
    n_to_keep = max(n_passing_threshold, n_qvalue_based, n_min_floor)

    kept_indices = order[1:min(n_to_keep, length(order))]
    hybrid_targets = count(i -> gp_df.target[i], kept_indices)

    return (
        n_01_fdr       = n_01,
        n_05_fdr       = n_05,
        hybrid_targets = hybrid_targets,
        n_to_keep      = n_to_keep,
        n_qvalue_based = n_qvalue_based,
        n_passing_threshold = n_passing_threshold,
        n_min_floor    = n_min_floor,
    )
end

###############################################################################
# 8. Main
###############################################################################

function print_results_table(results::Vector{Tuple{String, NamedTuple}}, baseline_key::String)
    @printf("%-34s %8s %8s %8s\n", "", "1% FDR", "5% FDR", "Hybrid")
    println("-" ^ 64)

    for (name, res) in results
        @printf("%-34s %8d %8d %8d targets\n", name * ":", res.n_01_fdr, res.n_05_fdr, res.hybrid_targets)
    end
    println("-" ^ 64)

    # Find baseline for gain computation
    base_res = nothing
    for (name, res) in results
        if name == baseline_key
            base_res = res
            break
        end
    end
    base_res === nothing && return

    for (name, res) in results
        name == baseline_key && continue
        d01  = res.n_01_fdr - base_res.n_01_fdr
        d05  = res.n_05_fdr - base_res.n_05_fdr
        dhyb = res.hybrid_targets - base_res.hybrid_targets
        p01  = base_res.n_01_fdr > 0 ? 100.0 * d01 / base_res.n_01_fdr : 0.0
        p05  = base_res.n_05_fdr > 0 ? 100.0 * d05 / base_res.n_05_fdr : 0.0
        phyb = base_res.hybrid_targets > 0 ? 100.0 * dhyb / base_res.hybrid_targets : 0.0
        @printf("  Δ %-28s %+7d  %+7d  %+7d\n", name * ":", d01, d05, dhyb)
        @printf("  %30s %+6.1f%%  %+6.1f%%  %+6.1f%%\n", "", p01, p05, phyb)
    end
end

function main()
    println("=" ^ 70)
    println("Eigenvalue + Chromatogram Feature Rescoring Analysis")
    println("=" ^ 70)

    # Step 1: Load data
    psm_all, corr_all, n_files = load_data()

    # Step 2: Compute eigenvalue features — full window and narrow (±2 from apex)
    eigen_full   = compute_eigenvalue_features(corr_all; window=nothing, prefix="")
    eigen_narrow = compute_eigenvalue_features(corr_all; window=2,       prefix="narrow_")

    # Step 3: Compute chromatogram features — full window and narrow (±2 from apex)
    chrom_full   = compute_chromatogram_features(corr_all; window=nothing, prefix="")
    chrom_narrow = compute_chromatogram_features(corr_all; window=2,       prefix="narrow_")

    # Free corr data
    corr_all = nothing
    GC.gc()

    # Step 4: Build feature DataFrame with all feature sets
    df = build_feature_df(psm_all, eigen_full, eigen_narrow, chrom_full, chrom_narrow)
    eigen_full = nothing
    eigen_narrow = nothing
    chrom_full = nothing
    chrom_narrow = nothing
    GC.gc()

    # Pipeline baseline
    println("\n--- Pipeline Baseline (probit prob) ---")
    gp_pipeline = compute_global_probs(df, :prob, n_files)
    res_pipeline = hybrid_filter(gp_pipeline)

    # Step 5: Train LightGBM models — 6 variants
    models = [
        ("LightGBM baseline",       STANDARD_FEATURES,       "baseline"),
        ("LightGBM + eigen top4",   AUG_EIGEN_TOP4_FEATURES, "aug_eigen_top4"),
        ("LightGBM + both eigen",   AUG_BOTH_FEATURES,       "aug_both"),
        ("LightGBM + chrom top3",   AUG_CHROM_TOP3_FEATURES, "aug_chrom_top3"),
        ("LightGBM + chrom",        AUG_CHROM_FEATURES,      "aug_chrom"),
        ("LightGBM + all",          AUG_ALL_FEATURES,        "aug_all"),
    ]

    model_results = Tuple{String, NamedTuple}[]
    push!(model_results, ("Pipeline (probit)", res_pipeline))

    for (name, features, col_suffix) in models
        println("\n--- $name ---")
        probs = train_lgbm_cv(df, features; label=col_suffix)
        col_name = Symbol("lgbm_prob_" * col_suffix)
        df[!, col_name] = probs

        gp = compute_global_probs(df, col_name, n_files)
        res = hybrid_filter(gp)
        push!(model_results, (name, res))
    end

    # Step 6: Report
    println("\n" * "=" ^ 70)
    println("RESULTS")
    println("=" ^ 70)
    print_results_table(model_results, "LightGBM baseline")

    println("\n--- Hybrid Filter Details ---")
    for (name, res) in model_results
        @printf("  %-28s  n_qval=%6d  n_pep=%6d  n_floor=%6d  → keep=%6d  targets=%6d\n",
                name, res.n_qvalue_based, res.n_passing_threshold, res.n_min_floor,
                res.n_to_keep, res.hybrid_targets)
    end

    println("\nDone.")
end

main()
