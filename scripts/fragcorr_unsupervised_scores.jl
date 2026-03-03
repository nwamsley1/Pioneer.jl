#=
Fragment Correlation Unsupervised Score Analysis
================================================
Standalone script for analyzing FRAGCORR precursor-level data.
Computes unsupervised coherence/coelution scores and compares
target vs decoy distributions.

Usage (from Pioneer.jl root or Julia REPL):
    include("scripts/fragcorr_unsupervised_scores.jl")
=#

using Arrow, DataFrames, Tables
using Statistics, LinearAlgebra
using StatsBase: mad
using Plots

# ============================================================
# Data Loading
# ============================================================
const ARROW_PATH = "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/OlsenEclipse_fragcorr_test/temp_data/first_pass_corr/E45H50Y5_2.arrow"
const INTENSITY_COLS = [:intensity_1, :intensity_2, :intensity_3, :intensity_4, :intensity_5]

function load_fragcorr_data(path::String=ARROW_PATH)
    println("Loading $path ...")
    df = DataFrame(Tables.columntable(Arrow.Table(path)))
    sort!(df, [:precursor_idx, :scan_idx])
    gdf = groupby(df, :precursor_idx)
    println("  $(nrow(df)) rows, $(length(gdf)) precursor groups")
    return gdf
end

# ============================================================
# Matrix extraction & preprocessing
# ============================================================

"""Extract n_scans × 5 Float64 intensity matrix from a subtable."""
function get_intensity_matrix(sub::SubDataFrame)
    n = nrow(sub)
    X = Matrix{Float64}(undef, n, 5)
    @inbounds for (j, col) in enumerate(INTENSITY_COLS)
        copyto!(view(X, :, j), sub[!, col])
    end
    return X
end

"""Clamp negatives to zero."""
clamp_nn!(X) = (X .= max.(X, 0.0); X)

"""Per-column robust scaling: (x - median) / (MAD + ε), then clamp to [-clip, clip]."""
function robust_scale!(Xout, X; eps=1.0f0, clip=5.0)
    @inbounds for j in axes(X, 2)
        col = @view X[:, j]
        med = median(col)
        m   = mad(col; normalize=false)
        for i in axes(X, 1)
            Xout[i, j] = clamp((col[i] - med) / (m + eps), -clip, clip)
        end
    end
    return Xout
end

"""Per-column unit-area normalization (L1)."""
function unit_area_normalize!(Xout, X; eps=1e-12)
    @inbounds for j in axes(X, 2)
        s = 0.0
        for i in axes(X, 1); s += abs(X[i, j]); end
        inv_s = 1.0 / (s + eps)
        for i in axes(X, 1); Xout[i, j] = X[i, j] * inv_s; end
    end
    return Xout
end

"""Per-column unit-norm normalization (L2)."""
function unit_norm_normalize!(Xout, X; eps=1e-12)
    @inbounds for j in axes(X, 2)
        s = 0.0
        for i in axes(X, 1); s += X[i, j]^2; end
        inv_s = 1.0 / (sqrt(s) + eps)
        for i in axes(X, 1); Xout[i, j] = X[i, j] * inv_s; end
    end
    return Xout
end

# ============================================================
# Apex utilities
# ============================================================

"""Consensus apex index: argmax of mean across max-normalized columns."""
function consensus_apex_idx(X::Matrix{Float64})
    n, p = size(X)
    best_i, best_v = 1, -Inf
    @inbounds for i in 1:n
        v = 0.0
        for j in 1:p
            mx = 0.0
            for ii in 1:n
                mx = max(mx, X[ii, j])
            end
            v += mx > 0.0 ? X[i, j] / mx : 0.0
        end
        if v > best_v
            best_v = v
            best_i = i
        end
    end
    return best_i
end

"""Per-fragment apex scan indices."""
function fragment_apex_indices(X::Matrix{Float64})
    return [argmax(@view X[:, j]) for j in axes(X, 2)]
end

"""
    apex_center(X; half_window=5)

Extract a sub-matrix of ±half_window scans around the consensus apex.
Returns the windowed matrix (may be smaller than 2*half_window+1 at edges).
"""
function apex_center(X::Matrix{Float64}; half_window::Int=5)
    n = size(X, 1)
    apex = consensus_apex_idx(X)
    lo = max(1, apex - half_window)
    hi = min(n, apex + half_window)
    return X[lo:hi, :]
end

# ============================================================
# Score functions
# ============================================================

# 1. Rank-1 explained variance (SVD)
function score_ev1(X::Matrix{Float64}; eps=1e-12)
    all_zero = true
    @inbounds for v in X; v != 0.0 && (all_zero = false; break); end
    all_zero && return 0.0
    S = svdvals(X)
    total = sum(abs2, S)
    return total > eps ? S[1]^2 / total : 0.0
end

# 2. Non-negative rank-1 fit via power iteration
function score_nn_rank1(X::Matrix{Float64}; eps=1e-12, n_iter=5)
    n, p = size(X)
    nrm_sq = sum(abs2, X)
    nrm_sq < eps && return 0.0

    # Power iteration: X ≈ u * v'  (both non-negative)
    v = ones(p); v ./= sqrt(p)
    u = Vector{Float64}(undef, n)

    for _ in 1:n_iter
        mul!(u, X, v)
        u .= max.(u, 0.0)
        un = norm(u)
        un < eps && return 0.0
        u ./= un
        mul!(v, X', u)
        v .= max.(v, 0.0)
        vn = norm(v)
        vn < eps && return 0.0
        v ./= vn
    end

    # Final scale: σ = u'Xv
    sigma = dot(u, X * v)
    sigma = max(sigma, 0.0)
    # Residual fraction: 1 - ||X - σuv'||² / ||X||²
    # ||X - σuv'||² = ||X||² - σ²  (since u,v unit norm, orthogonal decomp)
    return sigma^2 / nrm_sq
end

# 3. Pairwise correlation eigenvalue scores
function score_corr_eigen(X::Matrix{Float64}; eps=1e-12)
    n, p = size(X)
    # Build correlation matrix
    C = Matrix{Float64}(undef, p, p)
    stds = Vector{Float64}(undef, p)
    means = Vector{Float64}(undef, p)

    @inbounds for j in 1:p
        s = 0.0; ss = 0.0
        for i in 1:n; s += X[i, j]; end
        means[j] = s / n
        for i in 1:n; ss += (X[i, j] - means[j])^2; end
        stds[j] = sqrt(ss / (n - 1))
    end

    @inbounds for j1 in 1:p
        C[j1, j1] = 1.0
        for j2 in (j1+1):p
            if stds[j1] > eps && stds[j2] > eps
                cov = 0.0
                for i in 1:n
                    cov += (X[i, j1] - means[j1]) * (X[i, j2] - means[j2])
                end
                c = cov / ((n - 1) * stds[j1] * stds[j2])
                C[j1, j2] = c
                C[j2, j1] = c
            else
                C[j1, j2] = 0.0
                C[j2, j1] = 0.0
            end
        end
    end

    evals = eigvals(Symmetric(C))
    sort!(evals, rev=true)
    tr = sum(evals)

    lambda1_frac = tr > eps ? evals[1] / tr : 0.0
    eigengap = p >= 2 ? evals[1] - evals[2] : 0.0

    # Median off-diagonal correlation
    corrs = Vector{Float64}(undef, p * (p - 1) ÷ 2)
    k = 0
    @inbounds for j1 in 1:p, j2 in (j1+1):p
        k += 1; corrs[k] = C[j1, j2]
    end

    return (lambda1_frac=lambda1_frac, eigengap=eigengap, median_corr=median(corrs))
end

# 4. Entropy of fragment contribution around apex
function score_entropy(X::Matrix{Float64}; half_window=3, eps=1e-12)
    n, p = size(X)
    apex = consensus_apex_idx(X)
    lo = max(1, apex - half_window)
    hi = min(n, apex + half_window)

    sum_ent = 0.0
    n_rows = 0
    max_ent = log(p)

    @inbounds for i in lo:hi
        s = 0.0
        for j in 1:p; s += X[i, j]; end
        if s > eps
            ent = 0.0
            for j in 1:p
                pj = X[i, j] / s
                if pj > eps
                    ent -= pj * log(pj)
                end
            end
            sum_ent += ent
        else
            sum_ent += max_ent
        end
        n_rows += 1
    end

    mean_ent = sum_ent / max(n_rows, 1)
    return (mean_entropy=mean_ent, normalized_entropy=mean_ent / max_ent)
end

# 5. Consensus-template projection score
function score_consensus(X::Matrix{Float64}; eps=1e-12)
    _, p = size(X)
    # Build unit-norm version for template
    Xn = similar(X)
    unit_norm_normalize!(Xn, X)

    # Consensus trace = mean across normalized fragments
    consensus = vec(mean(Xn, dims=2))
    c_dot = dot(consensus, consensus)
    c_dot < eps && return (median_template_corr=0.0, template_residual=1.0)

    # Per-fragment correlation with consensus
    corrs = Vector{Float64}(undef, p)
    c_std  = std(consensus)
    valid = 0
    @inbounds for j in 1:p
        col = @view Xn[:, j]
        s = std(col)
        if s > eps && c_std > eps
            valid += 1
            corrs[valid] = cor(col, consensus)
        end
    end
    med_corr = valid > 0 ? median(@view corrs[1:valid]) : 0.0

    # Residual: ||Xn - consensus * w'|| / ||Xn||
    w = Xn' * consensus ./ (c_dot + eps)
    nrm = norm(Xn)
    residual = nrm > eps ? norm(Xn .- consensus * w') / nrm : 1.0

    return (median_template_corr=med_corr, template_residual=residual)
end

# ============================================================
# DIA-NN-inspired score functions
# ============================================================

"""
    smooth3!(dst, src)

3-point weighted smoothing (DIA-NN style).
dst and src must be length-n vectors. Can be same array.
"""
function smooth3!(dst::AbstractVector, src::AbstractVector)
    n = length(src)
    n <= 1 && (copyto!(dst, src); return dst)
    # Store boundary values before overwriting
    first_val = (2.0/3.0) * src[1] + (1.0/3.0) * src[2]
    last_val  = (2.0/3.0) * src[n] + (1.0/3.0) * src[n-1]
    @inbounds for i in 2:n-1
        dst[i] = 0.5 * src[i] + 0.25 * (src[i-1] + src[i+1])
    end
    dst[1] = first_val
    dst[n] = last_val
    return dst
end

"""
    best_fragment_idx(X)

Find the fragment (column) with highest sum of pairwise Pearson correlations
to all other fragments. Returns (best_col_idx, corr_matrix).
"""
function best_fragment_idx(X::Matrix{Float64}; eps=1e-12)
    n, p = size(X)
    means = Vector{Float64}(undef, p)
    stds  = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        s = 0.0
        for i in 1:n; s += X[i, j]; end
        means[j] = s / n
        ss = 0.0
        for i in 1:n; ss += (X[i, j] - means[j])^2; end
        stds[j] = sqrt(ss / max(n - 1, 1))
    end

    # Build correlation matrix
    C = zeros(p, p)
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
            end
        end
    end

    # Fragment with highest row-sum (excluding diagonal)
    best_j = 1
    best_sum = -Inf
    @inbounds for j in 1:p
        s = 0.0
        for j2 in 1:p
            j2 == j && continue
            s += C[j, j2]
        end
        if s > best_sum
            best_sum = s
            best_j = j
        end
    end
    return best_j, C
end

"""
    score_best_frag_corr(X; eps=1e-12)

DIA-NN pTimeCorr analog: find best fragment, smooth it, correlate each other
fragment against the smoothed reference. Returns:
  - corr_sum: sum of correlations (like pTimeCorr)
  - corr_min: minimum per-fragment correlation
  - corr_mean: mean per-fragment correlation
"""
function score_best_frag_corr(X::Matrix{Float64}; eps=1e-12)
    n, p = size(X)
    n < 3 && return (corr_sum=0.0, corr_min=0.0, corr_mean=0.0)

    best_j, _ = best_fragment_idx(X)

    # Smooth the best fragment to create reference elution
    elution = Vector{Float64}(undef, n)
    smooth3!(elution, @view X[:, best_j])

    # Correlate each fragment with the smoothed reference
    el_std = std(elution)
    corrs = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        col = @view X[:, j]
        s = std(col)
        if s > eps && el_std > eps
            corrs[j] = cor(col, elution)
        else
            corrs[j] = 0.0
        end
    end

    return (corr_sum=sum(corrs), corr_min=minimum(corrs), corr_mean=mean(corrs))
end

"""
    score_local_corr(X; eps=1e-12)

DIA-NN pLocCorr analog: best-fragment correlation in a narrow window
around the apex (half the full width).
"""
function score_local_corr(X::Matrix{Float64}; eps=1e-12)
    n, p = size(X)
    n < 3 && return 0.0

    best_j, _ = best_fragment_idx(X)
    apex = argmax(@view X[:, best_j])

    # Half-width window
    hw = max(n ÷ 4, 1)
    lo = max(1, apex - hw)
    hi = min(n, apex + hw)
    w  = hi - lo + 1
    w < 3 && return 0.0

    Xw = X[lo:hi, :]
    elution = Vector{Float64}(undef, w)
    smooth3!(elution, @view Xw[:, best_j])

    el_std = std(elution)
    s = 0.0
    @inbounds for j in 1:p
        col = @view Xw[:, j]
        cs = std(col)
        if cs > eps && el_std > eps
            s += cor(col, elution)
        end
    end
    return s
end

"""
    min_filter!(dst, src)

Replace each point with min(self, left_neighbor, right_neighbor).
Boundary: min(self, single_neighbor).
"""
function min_filter!(dst::AbstractVector, src::AbstractVector)
    n = length(src)
    n == 0 && return dst
    n == 1 && (dst[1] = src[1]; return dst)
    dst[1] = min(src[1], src[2])
    dst[n] = min(src[n], src[n-1])
    @inbounds for i in 2:n-1
        dst[i] = min(src[i], src[i-1], src[i+1])
    end
    return dst
end

"""
    score_min_corr(X; eps=1e-12)

DIA-NN pMinCorr analog: correlate smoothed best-fragment with
min-filtered chromatograms. Suppresses single-scan spikes.
"""
function score_min_corr(X::Matrix{Float64}; eps=1e-12)
    n, p = size(X)
    n < 3 && return 0.0

    best_j, _ = best_fragment_idx(X)
    elution = Vector{Float64}(undef, n)
    smooth3!(elution, @view X[:, best_j])

    el_std = std(elution)
    el_std < eps && return 0.0

    filtered = Vector{Float64}(undef, n)
    s = 0.0
    @inbounds for j in 1:p
        min_filter!(filtered, @view X[:, j])
        fs = std(filtered)
        if fs > eps
            s += cor(filtered, elution)
        end
    end
    return s
end

"""
    score_peak_shape(X)

DIA-NN signal_level analog: given the consensus apex, measure whether
signal rises on the left and falls on the right for each fragment.
Returns mean across fragments (0 = noise, 1 = perfect peak shape).
"""
function score_peak_shape(X::Matrix{Float64})
    n, p = size(X)
    n < 3 && return 0.0

    apex = consensus_apex_idx(X)
    total_score = 0.0

    @inbounds for j in 1:p
        signal = 0.0
        noise  = 0.0
        # Left of apex: rising = signal, falling = noise
        for i in 1:apex-1
            delta = X[i+1, j] - X[i, j]
            change = abs(delta)
            if delta >= 0.0
                signal += change
            else
                noise += change
            end
        end
        # Right of apex: falling = signal, rising = noise
        for i in apex+1:n
            delta = X[i, j] - X[i-1, j]
            change = abs(delta)
            if delta <= 0.0
                signal += change
            else
                noise += change
            end
        end
        denom = signal + noise
        total_score += denom > 0.0 ? signal / denom : 0.0
    end

    return total_score / p
end

"""
    score_frag_signal_fractions(X)

DIA-NN pSig analog: for each fragment, compute its fractional contribution
to total signal. Returns (max_frac, std_frac) — max single-fragment
dominance and spread of fractions.
"""
function score_frag_signal_fractions(X::Matrix{Float64}; eps=1e-12)
    p = size(X, 2)
    totals = vec(sum(X, dims=1))  # sum per fragment
    grand  = sum(totals)
    grand < eps && return (max_frac=0.0, std_frac=0.0)
    fracs = totals ./ grand
    return (max_frac=maximum(fracs), std_frac=std(fracs))
end

"""
    score_peak_shape_bins(X; nbins=5)

DIA-NN pShape analog: divide the consensus elution curve into nbins
equal segments, compute fraction of total signal in each bin.
Returns the entropy of the bin distribution (low = peaked, high = flat).
"""
function score_peak_shape_bins(X::Matrix{Float64}; nbins::Int=5, eps=1e-12)
    n = size(X, 1)
    n < nbins && return 0.0

    # Use sum across fragments as the elution curve
    elution = vec(sum(X, dims=2))
    total = sum(elution)
    total < eps && return log(nbins)  # max entropy = no signal

    bin_sums = zeros(nbins)
    @inbounds for i in 1:n
        # Map scan index to bin (1-based)
        b = min(nbins, 1 + ((i - 1) * nbins) ÷ n)
        bin_sums[b] += elution[i]
    end

    # Entropy of bin distribution
    ent = 0.0
    @inbounds for b in 1:nbins
        p_b = bin_sums[b] / total
        if p_b > eps
            ent -= p_b * log(p_b)
        end
    end
    return ent
end

"""
    score_elution_cosine(X; eps=1e-12)

DIA-NN pCos analog (without library reference): at each scan, compute
fragment intensity vector and its cosine similarity to the mean fragment
profile across scans. Weight by total signal squared. Returns weighted mean.
"""
function score_elution_cosine(X::Matrix{Float64}; eps=1e-12)
    n, p = size(X)
    # Mean fragment profile across all scans
    ref = vec(mean(X, dims=1))
    ref_norm_sq = sum(abs2, ref)
    ref_norm_sq < eps && return 0.0
    ref_norm = sqrt(ref_norm_sq)

    weighted_cos = 0.0
    total_weight = 0.0
    @inbounds for i in 1:n
        row_norm_sq = 0.0
        dot_val = 0.0
        for j in 1:p
            row_norm_sq += X[i, j]^2
            dot_val += X[i, j] * ref[j]
        end
        row_norm_sq < eps && continue
        cos_val = dot_val / (sqrt(row_norm_sq) * ref_norm)
        weight = row_norm_sq  # weight by signal squared (like DIA-NN uses elution^2)
        weighted_cos += cos_val * weight
        total_weight += weight
    end
    return total_weight > eps ? weighted_cos / total_weight : 0.0
end

# Additional: apex dispersion (std of per-fragment apex positions)
function score_apex_dispersion(X::Matrix{Float64})
    positions = fragment_apex_indices(X)
    return std(Float64.(positions))
end

# Additional: amplitude (sum of top-k intensities)
function score_amplitude(X::Matrix{Float64}; k=10)
    vals = vec(X)
    partialsort!(vals, 1:min(k, length(vals)), rev=true)
    return sum(@view vals[1:min(k, length(vals))])
end

# ============================================================
# Combined scoring for one precursor group
# ============================================================

function compute_scores(sub::SubDataFrame; min_scans=3, apex_half_window=5)
    X = get_intensity_matrix(sub)
    n = size(X, 1)
    n < min_scans && return nothing

    Xnn = clamp_nn!(copy(X))

    # Pre-allocate workspace for normalized variants
    Xwork = similar(Xnn)

    # ==========================================================
    # Original scores (raw)
    # ==========================================================
    ev1_raw       = score_ev1(Xnn)
    nn_rank1      = score_nn_rank1(Xnn)
    corr_raw      = score_corr_eigen(Xnn)
    ent           = score_entropy(Xnn)
    cons_raw      = score_consensus(Xnn)
    apex_disp     = score_apex_dispersion(Xnn)
    amp           = score_amplitude(Xnn)

    # --- Unit-norm (L2) preprocessing ---
    unit_norm_normalize!(Xwork, Xnn)
    ev1_unitnorm      = score_ev1(Xwork)
    corr_unitnorm     = score_corr_eigen(Xwork)
    cons_unitnorm     = score_consensus(Xwork)

    # --- Robust scaling ---
    robust_scale!(Xwork, Xnn)
    ev1_robust        = score_ev1(Xwork)

    # --- Unit-area (L1) preprocessing ---
    unit_area_normalize!(Xwork, Xnn)
    ev1_unitarea          = score_ev1(Xwork)
    corr_unitarea         = score_corr_eigen(Xwork)
    cons_unitarea         = score_consensus(Xwork)

    # --- Apex-centered (raw, ±half_window scans around consensus apex) ---
    Xapex = apex_center(Xnn; half_window=apex_half_window)
    n_apex = size(Xapex, 1)
    if n_apex >= min_scans
        ev1_apex          = score_ev1(Xapex)
        nn_rank1_apex     = score_nn_rank1(Xapex)
        corr_apex         = score_corr_eigen(Xapex)
        cons_apex         = score_consensus(Xapex)
        ent_apex          = score_entropy(Xapex; half_window=min(3, n_apex ÷ 2))
    else
        ev1_apex          = 0.0
        nn_rank1_apex     = 0.0
        corr_apex         = (lambda1_frac=0.0, eigengap=0.0, median_corr=0.0)
        cons_apex         = (median_template_corr=0.0, template_residual=1.0)
        ent_apex          = (mean_entropy=log(5), normalized_entropy=1.0)
    end

    # ==========================================================
    # DIA-NN-inspired scores (on raw full-window matrix)
    # ==========================================================
    bfc_raw       = score_best_frag_corr(Xnn)
    loc_corr_raw  = score_local_corr(Xnn)
    min_corr_raw  = score_min_corr(Xnn)
    pk_shape_raw  = score_peak_shape(Xnn)
    frac_raw      = score_frag_signal_fractions(Xnn)
    shape_bins_raw = score_peak_shape_bins(Xnn)
    elut_cos_raw  = score_elution_cosine(Xnn)

    # DIA-NN-inspired scores (on apex-centered matrix)
    if n_apex >= min_scans
        bfc_apex       = score_best_frag_corr(Xapex)
        loc_corr_apex  = score_local_corr(Xapex)
        min_corr_apex  = score_min_corr(Xapex)
        pk_shape_apex  = score_peak_shape(Xapex)
        shape_bins_apex = score_peak_shape_bins(Xapex)
        elut_cos_apex  = score_elution_cosine(Xapex)
    else
        bfc_apex       = (corr_sum=0.0, corr_min=0.0, corr_mean=0.0)
        loc_corr_apex  = 0.0
        min_corr_apex  = 0.0
        pk_shape_apex  = 0.0
        shape_bins_apex = log(5)
        elut_cos_apex  = 0.0
    end

    return (
        precursor_idx               = first(sub.precursor_idx),
        target                      = first(sub.target),
        n_scans                     = n,
        # ======== Original raw ========
        ev1_raw                     = ev1_raw,
        nn_rank1_fit                = nn_rank1,
        lambda1_frac_raw            = corr_raw.lambda1_frac,
        eigengap_raw                = corr_raw.eigengap,
        median_corr_raw             = corr_raw.median_corr,
        mean_entropy                = ent.mean_entropy,
        normalized_entropy          = ent.normalized_entropy,
        median_template_corr_raw    = cons_raw.median_template_corr,
        template_residual_raw       = cons_raw.template_residual,
        apex_dispersion             = apex_disp,
        amplitude                   = amp,
        # ======== Unit-norm (L2) ========
        ev1_unitnorm                = ev1_unitnorm,
        lambda1_frac_unitnorm       = corr_unitnorm.lambda1_frac,
        eigengap_unitnorm           = corr_unitnorm.eigengap,
        median_corr_unitnorm        = corr_unitnorm.median_corr,
        median_template_corr_unitnorm = cons_unitnorm.median_template_corr,
        template_residual_unitnorm  = cons_unitnorm.template_residual,
        # ======== Robust ========
        ev1_robust                  = ev1_robust,
        # ======== Unit-area (L1) ========
        ev1_unitarea                = ev1_unitarea,
        lambda1_frac_unitarea       = corr_unitarea.lambda1_frac,
        eigengap_unitarea           = corr_unitarea.eigengap,
        median_corr_unitarea        = corr_unitarea.median_corr,
        median_template_corr_unitarea = cons_unitarea.median_template_corr,
        template_residual_unitarea  = cons_unitarea.template_residual,
        # ======== Apex-centered (original scores) ========
        ev1_apex                    = ev1_apex,
        nn_rank1_apex               = nn_rank1_apex,
        lambda1_frac_apex           = corr_apex.lambda1_frac,
        eigengap_apex               = corr_apex.eigengap,
        median_corr_apex            = corr_apex.median_corr,
        median_template_corr_apex   = cons_apex.median_template_corr,
        template_residual_apex      = cons_apex.template_residual,
        mean_entropy_apex           = ent_apex.mean_entropy,
        normalized_entropy_apex     = ent_apex.normalized_entropy,
        # ======== DIA-NN-inspired (raw) ========
        best_frag_corr_sum          = bfc_raw.corr_sum,
        best_frag_corr_min          = bfc_raw.corr_min,
        best_frag_corr_mean         = bfc_raw.corr_mean,
        local_corr                  = loc_corr_raw,
        min_corr                    = min_corr_raw,
        peak_shape                  = pk_shape_raw,
        frag_max_frac               = frac_raw.max_frac,
        frag_std_frac               = frac_raw.std_frac,
        peak_shape_bins_entropy     = shape_bins_raw,
        elution_cosine              = elut_cos_raw,
        # ======== DIA-NN-inspired (apex-centered) ========
        best_frag_corr_sum_apex     = bfc_apex.corr_sum,
        best_frag_corr_min_apex     = bfc_apex.corr_min,
        best_frag_corr_mean_apex    = bfc_apex.corr_mean,
        local_corr_apex             = loc_corr_apex,
        min_corr_apex               = min_corr_apex,
        peak_shape_apex             = pk_shape_apex,
        peak_shape_bins_entropy_apex = shape_bins_apex,
        elution_cosine_apex         = elut_cos_apex,
    )
end

# ============================================================
# Batch scoring
# ============================================================

function score_all_precursors(gdf; min_scans=3, report_every=50_000)
    n_groups = length(gdf)
    results = Vector{NamedTuple}()
    sizehint!(results, n_groups)

    println("Scoring $n_groups precursor groups...")
    t0 = time()

    for (i, sub) in enumerate(gdf)
        r = compute_scores(sub; min_scans)
        r !== nothing && push!(results, r)
        if i % report_every == 0
            elapsed = round(time() - t0, digits=1)
            println("  $i / $n_groups  ($(elapsed)s)")
        end
    end

    elapsed = round(time() - t0, digits=1)
    println("Done in $(elapsed)s — scored $(length(results)) precursors")

    return DataFrame(results)
end

# ============================================================
# Summary / comparison
# ============================================================

const SCORE_COLS = [
    # --- Original raw ---
    :ev1_raw, :nn_rank1_fit,
    :lambda1_frac_raw, :eigengap_raw, :median_corr_raw,
    :mean_entropy, :normalized_entropy,
    :median_template_corr_raw, :template_residual_raw,
    :apex_dispersion, :amplitude,
    # --- Unit-norm (L2) ---
    :ev1_unitnorm,
    :lambda1_frac_unitnorm, :eigengap_unitnorm, :median_corr_unitnorm,
    :median_template_corr_unitnorm, :template_residual_unitnorm,
    # --- Robust ---
    :ev1_robust,
    # --- Unit-area (L1) ---
    :ev1_unitarea,
    :lambda1_frac_unitarea, :eigengap_unitarea, :median_corr_unitarea,
    :median_template_corr_unitarea, :template_residual_unitarea,
    # --- Apex-centered (original scores) ---
    :ev1_apex, :nn_rank1_apex,
    :lambda1_frac_apex, :eigengap_apex, :median_corr_apex,
    :median_template_corr_apex, :template_residual_apex,
    :mean_entropy_apex, :normalized_entropy_apex,
    # --- DIA-NN-inspired (raw) ---
    :best_frag_corr_sum, :best_frag_corr_min, :best_frag_corr_mean,
    :local_corr, :min_corr,
    :peak_shape,
    :frag_max_frac, :frag_std_frac,
    :peak_shape_bins_entropy,
    :elution_cosine,
    # --- DIA-NN-inspired (apex-centered) ---
    :best_frag_corr_sum_apex, :best_frag_corr_min_apex, :best_frag_corr_mean_apex,
    :local_corr_apex, :min_corr_apex,
    :peak_shape_apex,
    :peak_shape_bins_entropy_apex,
    :elution_cosine_apex,
]

"""Print target vs decoy summary for each score."""
function print_summary(scores_df::DataFrame)
    targets = scores_df[scores_df.target, :]
    decoys  = scores_df[.!scores_df.target, :]
    println("\nTargets: $(nrow(targets))  |  Decoys: $(nrow(decoys))")
    println("─"^95)
    println(rpad("score", 36), rpad("T_median", 12), rpad("D_median", 12),
            rpad("T_mean", 12), rpad("D_mean", 12), "separation")
    println("─"^95)

    for col in SCORE_COLS
        t_vals = targets[!, col]
        d_vals = decoys[!, col]
        t_med  = median(t_vals)
        d_med  = median(d_vals)
        t_mean = mean(t_vals)
        d_mean = mean(d_vals)
        pooled = std(vcat(t_vals, d_vals))
        sep    = pooled > 0 ? abs(t_mean - d_mean) / pooled : 0.0

        println(
            rpad(string(col), 36),
            rpad(round(t_med,  digits=4), 12),
            rpad(round(d_med,  digits=4), 12),
            rpad(round(t_mean, digits=4), 12),
            rpad(round(d_mean, digits=4), 12),
            round(sep, digits=4),
        )
    end
    println("─"^95)
end

# ============================================================
# Q-value / FDR analysis
# ============================================================

# Scores where higher values indicate better targets
const HIGHER_IS_BETTER = Set([
    # Original
    :ev1_raw, :ev1_unitnorm, :ev1_robust, :ev1_unitarea, :ev1_apex,
    :nn_rank1_fit, :nn_rank1_apex,
    :lambda1_frac_raw, :eigengap_raw, :median_corr_raw,
    :lambda1_frac_unitnorm, :eigengap_unitnorm, :median_corr_unitnorm,
    :lambda1_frac_unitarea, :eigengap_unitarea, :median_corr_unitarea,
    :lambda1_frac_apex, :eigengap_apex, :median_corr_apex,
    :median_template_corr_raw, :median_template_corr_unitnorm,
    :median_template_corr_unitarea, :median_template_corr_apex,
    :amplitude,
    # DIA-NN-inspired
    :best_frag_corr_sum, :best_frag_corr_min, :best_frag_corr_mean,
    :best_frag_corr_sum_apex, :best_frag_corr_min_apex, :best_frag_corr_mean_apex,
    :local_corr, :local_corr_apex,
    :min_corr, :min_corr_apex,
    :peak_shape, :peak_shape_apex,
    :elution_cosine, :elution_cosine_apex,
])
# Scores where lower values indicate better targets
const LOWER_IS_BETTER = Set([
    # Original
    :mean_entropy, :normalized_entropy,
    :mean_entropy_apex, :normalized_entropy_apex,
    :template_residual_raw, :template_residual_unitnorm,
    :template_residual_unitarea, :template_residual_apex,
    :apex_dispersion,
    # DIA-NN-inspired
    :frag_max_frac, :frag_std_frac,  # high dominance by one frag = chimeric
    :peak_shape_bins_entropy, :peak_shape_bins_entropy_apex,
])

"""
    compute_qvalues(scores, is_target; higher_is_better=true)

Standard target-decoy competition q-value calculation.
Returns a vector of q-values (same order as input).
"""
function compute_qvalues(scores::AbstractVector{<:Real},
                         is_target::AbstractVector{Bool};
                         higher_is_better::Bool=true)
    n = length(scores)
    # Sort indices by score
    order = sortperm(scores, rev=higher_is_better)

    # Forward pass: cumulative FDR
    cum_targets = 0
    cum_decoys  = 0
    fdr = Vector{Float64}(undef, n)

    for idx in order
        if is_target[idx]
            cum_targets += 1
        else
            cum_decoys += 1
        end
        fdr[idx] = cum_targets > 0 ? cum_decoys / cum_targets : 1.0
    end

    # Backward pass: q-value = running minimum FDR from bottom up
    qvals = Vector{Float64}(undef, n)
    min_fdr = 1.0
    for idx in reverse(order)
        min_fdr = min(min_fdr, fdr[idx])
        qvals[idx] = min_fdr
    end

    return qvals
end

"""
    targets_at_qvalue(qvals, is_target, threshold)

Count targets with q-value ≤ threshold.
"""
function targets_at_qvalue(qvals, is_target, threshold)
    count = 0
    @inbounds for i in eachindex(qvals)
        count += (is_target[i] && qvals[i] <= threshold)
    end
    return count
end

"""
    qvalue_summary(scores_df; cols, thresholds)

Print table: how many targets pass at each q-value threshold for each score.
"""
function qvalue_summary(scores_df::DataFrame;
                        cols::Vector{Symbol}=SCORE_COLS,
                        thresholds=[0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20])
    is_target = scores_df.target
    n_targets = sum(is_target)
    n_decoys  = sum(.!is_target)
    println("\nTargets: $n_targets  |  Decoys: $n_decoys")

    # Header
    thresh_strs = [string(round(t * 100, digits=1)) * "%" for t in thresholds]
    header = rpad("score", 36) * join(rpad.(thresh_strs, 10))
    println("─"^(36 + 10 * length(thresholds)))
    println(header)
    println("─"^(36 + 10 * length(thresholds)))

    for col in cols
        hib = col in HIGHER_IS_BETTER
        qvals = compute_qvalues(scores_df[!, col], is_target; higher_is_better=hib)
        counts = [targets_at_qvalue(qvals, is_target, t) for t in thresholds]
        row = rpad(string(col), 36) * join(rpad.(string.(counts), 10))
        println(row)
    end
    println("─"^(36 + 10 * length(thresholds)))
end

"""
    plot_qvalue_curves(scores_df; cols, max_qval=0.10)

Plot number of targets passing vs q-value threshold for selected scores.
"""
function plot_qvalue_curves(scores_df::DataFrame;
                            cols::Vector{Symbol}=[
                                # Best original scores
                                :lambda1_frac_raw, :eigengap_raw, :template_residual_raw,
                                # DIA-NN-inspired
                                :best_frag_corr_sum, :best_frag_corr_mean,
                                :local_corr, :min_corr,
                                :peak_shape, :elution_cosine,
                            ],
                            max_qval::Float64=0.10,
                            n_points::Int=200)
    is_target = scores_df.target
    thresholds = range(0.0, max_qval, length=n_points)

    p = plot(xlabel="q-value", ylabel="targets passing",
             title="Targets passing at q-value threshold",
             legend=:topleft, size=(800, 500))

    for (i, col) in enumerate(cols)
        hib = col in HIGHER_IS_BETTER
        qvals = compute_qvalues(scores_df[!, col], is_target; higher_is_better=hib)
        counts = [targets_at_qvalue(qvals, is_target, t) for t in thresholds]
        plot!(p, collect(thresholds), counts, label=string(col), lw=2, color=i)
    end

    display(p)
    return p
end

# ============================================================
# Plotting
# ============================================================

"""
    plot_score_histograms(scores_df; cols=SCORE_COLS, nbins=80)

Side-by-side target/decoy histograms for each score.
Returns a vector of plot objects.
"""
function plot_score_histograms(scores_df::DataFrame;
                               cols=SCORE_COLS, nbins=80)
    targets = scores_df[scores_df.target, :]
    decoys  = scores_df[.!scores_df.target, :]
    plots = []

    for col in cols
        t_vals = filter(isfinite, targets[!, col])
        d_vals = filter(isfinite, decoys[!, col])

        lo = min(quantile(t_vals, 0.01), quantile(d_vals, 0.01))
        hi = max(quantile(t_vals, 0.99), quantile(d_vals, 0.99))

        p = histogram(d_vals, bins=range(lo, hi, length=nbins),
                      alpha=0.5, label="decoy", normalize=:pdf, color=:red)
        histogram!(p, t_vals, bins=range(lo, hi, length=nbins),
                   alpha=0.5, label="target", normalize=:pdf, color=:blue)
        title!(p, string(col))
        push!(plots, p)
    end
    return plots
end

"""
    plot_score_grid(scores_df; cols=SCORE_COLS, nbins=80)

All score histograms in a single multi-panel figure.
"""
function plot_score_grid(scores_df::DataFrame;
                         cols=SCORE_COLS, nbins=80)
    ps = plot_score_histograms(scores_df; cols, nbins)
    ncols = 3
    nrows = ceil(Int, length(ps) / ncols)
    p = plot(ps..., layout=(nrows, ncols), size=(1400, 350 * nrows),
             margin=5Plots.mm, legend=:topright)
    display(p)
    return p
end

"""Plot chromatogram traces for a single precursor subtable."""
function plot_chromatogram(sub; title_prefix="")
    rt = sub[!, :rt]
    p = plot(title="$(title_prefix)prec=$(first(sub.precursor_idx)) target=$(first(sub.target))",
             xlabel="RT (min)", ylabel="intensity")
    for (j, col) in enumerate(INTENSITY_COLS)
        vals = sub[!, col]
        plot!(p, rt, vals, label="frag $j", color=j, lw=1.5)
        scatter!(p, rt, vals, color=j, ms=2, label=nothing)
    end
    display(p)
    return p
end

# ============================================================
# Run if included directly
# ============================================================

if abspath(PROGRAM_FILE) == @__FILE__
    gdf = load_fragcorr_data()
    scores_df = score_all_precursors(gdf)
    print_summary(scores_df)
end

println("""

Loaded fragcorr_unsupervised_scores.jl
────────────────────────────────────────
Quick start:
    gdf = load_fragcorr_data()
    scores_df = score_all_precursors(gdf)
    print_summary(scores_df)

Q-value analysis:
    qvalue_summary(scores_df)
    qvalue_summary(scores_df; cols=[:lambda1_frac_raw, :eigengap_raw, :template_residual_raw])
    plot_qvalue_curves(scores_df)
    plot_qvalue_curves(scores_df; cols=[:lambda1_frac_raw, :eigengap_raw, :template_residual_raw])

Plotting:
    plot_score_grid(scores_df)
    plot_chromatogram(gdf[1000])
""")
