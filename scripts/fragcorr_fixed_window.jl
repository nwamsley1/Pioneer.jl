#=
Fragment Correlation: Fixed-Window Experiments
==============================================
Exploring whether centering on a data-driven apex with a fixed scan count
improves lambda1_frac and eigengap discrimination.

Problem: precursors can hit the fragment index at two separate RT regions,
creating overly wide windows that dilute the correlation signal. We want
to find the "main event" and crop tightly around it.

Usage:
    include("scripts/fragcorr_fixed_window.jl")
    # Uses load/qvalue functions from fragcorr_unsupervised_scores.jl
=#

# Load the base script for data loading, q-value, and plotting utilities
include(joinpath(@__DIR__, "fragcorr_unsupervised_scores.jl"))

# ============================================================
# Improved apex finding
# ============================================================

"""
    apex_sum_smooth(X; n_smooth=2)

Find apex as argmax of smoothed sum-of-all-fragments.
More robust than max-normalized consensus — the total signal across
all fragments is the strongest indicator of where the real peak is.
Applies `n_smooth` rounds of 3-point smoothing before argmax.
"""
function apex_sum_smooth(X::Matrix{Float64}; n_smooth::Int=2)
    trace = vec(sum(X, dims=2))
    buf = similar(trace)
    for _ in 1:n_smooth
        smooth3!(buf, trace)
        trace, buf = buf, trace
    end
    return argmax(trace)
end

"""
    apex_best_fragment(X; n_smooth=1)

Find apex using DIA-NN strategy: pick the fragment with highest
pairwise correlation sum, smooth it, take argmax. Falls back to
sum-smooth if all correlations are zero.
"""
function apex_best_fragment(X::Matrix{Float64}; n_smooth::Int=1)
    n = size(X, 1)
    n < 3 && return argmax(vec(sum(X, dims=2)))

    best_j, _ = best_fragment_idx(X)
    trace = X[:, best_j]
    buf = similar(trace)
    for _ in 1:n_smooth
        smooth3!(buf, trace)
        trace, buf = buf, trace
    end
    return argmax(trace)
end

"""
    apex_topk_smooth(X; k=3, n_smooth=2)

Find apex using sum of top-k fragments (by total intensity),
smoothed. Balances between single-fragment noise and dilution
from zero-signal fragments.
"""
function apex_topk_smooth(X::Matrix{Float64}; k::Int=3, n_smooth::Int=2)
    n, p = size(X)
    # Rank fragments by total intensity
    totals = vec(sum(X, dims=1))
    top_idx = partialsortperm(totals, 1:min(k, p), rev=true)

    trace = zeros(n)
    for j in top_idx
        @inbounds for i in 1:n
            trace[i] += X[i, j]
        end
    end

    buf = similar(trace)
    for _ in 1:n_smooth
        smooth3!(buf, trace)
        trace, buf = buf, trace
    end
    return argmax(trace)
end

# ============================================================
# Fixed-window extraction
# ============================================================

"""
    extract_fixed_window(X, apex_idx, half_window)

Extract exactly 2*half_window+1 scans centered on apex_idx.
Pads with zeros if near edges (so output is always the same size).
"""
function extract_fixed_window(X::Matrix{Float64}, apex_idx::Int, half_window::Int)
    n, p = size(X)
    w = 2 * half_window + 1
    Xw = zeros(w, p)
    for offset in -half_window:half_window
        src_i = apex_idx + offset
        dst_i = offset + half_window + 1
        if 1 <= src_i <= n
            @inbounds for j in 1:p
                Xw[dst_i, j] = X[src_i, j]
            end
        end
    end
    return Xw
end

# ============================================================
# Score computation for fixed windows
# ============================================================

"""Compute lambda1_frac and eigengap on a fixed window."""
function score_corr_eigen_only(X::Matrix{Float64}; eps=1e-12)
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

    evals = eigvals(Symmetric(C))
    sort!(evals, rev=true)
    tr = sum(evals)

    lambda1_frac = tr > eps ? evals[1] / tr : 0.0
    eigengap = p >= 2 ? evals[1] - evals[2] : 0.0
    return (lambda1_frac=lambda1_frac, eigengap=eigengap)
end

# ============================================================
# Sweep: apex method × window size
# ============================================================

const APEX_METHODS = Dict(
    :sum_smooth     => apex_sum_smooth,
    :best_fragment  => apex_best_fragment,
    :topk_smooth    => apex_topk_smooth,
)

const WINDOW_SIZES = [5, 8, 10, 12, 15, 20]  # half-windows → total 11,17,21,25,31,41 scans

"""
    sweep_fixed_windows(gdf; min_scans=3, report_every=50_000)

For each precursor group, compute lambda1_frac and eigengap using
every combination of apex method × window size, plus the raw baseline.
Returns a DataFrame with all results.
"""
function sweep_fixed_windows(gdf; min_scans=3, report_every=50_000)
    n_groups = length(gdf)

    # Build column names
    col_names = Symbol[]
    push!(col_names, :precursor_idx, :target, :n_scans)
    push!(col_names, :lambda1_frac_raw, :eigengap_raw)
    for method_name in keys(APEX_METHODS)
        for hw in WINDOW_SIZES
            push!(col_names, Symbol("lambda1_$(method_name)_w$(2hw+1)"))
            push!(col_names, Symbol("eigengap_$(method_name)_w$(2hw+1)"))
        end
    end

    results = Vector{NamedTuple}()
    sizehint!(results, n_groups)

    println("Sweeping $n_groups groups × $(length(APEX_METHODS)) apex methods × $(length(WINDOW_SIZES)) window sizes...")
    t0 = time()

    for (i, sub) in enumerate(gdf)
        X = get_intensity_matrix(sub)
        n = size(X, 1)
        n < min_scans && continue

        Xnn = clamp_nn!(copy(X))

        # Raw baseline
        raw = score_corr_eigen_only(Xnn)

        # Fixed-window variants
        vals = Dict{Symbol, Float64}()
        vals[:precursor_idx] = Float64(first(sub.precursor_idx))
        vals[:target] = first(sub.target) ? 1.0 : 0.0
        vals[:n_scans] = Float64(n)
        vals[:lambda1_frac_raw] = raw.lambda1_frac
        vals[:eigengap_raw] = raw.eigengap

        for (method_name, method_fn) in APEX_METHODS
            apex = method_fn(Xnn)
            for hw in WINDOW_SIZES
                Xw = extract_fixed_window(Xnn, apex, hw)
                sc = score_corr_eigen_only(Xw)
                vals[Symbol("lambda1_$(method_name)_w$(2hw+1)")] = sc.lambda1_frac
                vals[Symbol("eigengap_$(method_name)_w$(2hw+1)")] = sc.eigengap
            end
        end

        # Convert to NamedTuple (stable column order)
        nt = NamedTuple{Tuple(col_names)}(Tuple(vals[k] for k in col_names))
        push!(results, nt)

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
# Analysis helpers
# ============================================================

"""
    sweep_qvalue_summary(sweep_df; thresholds)

Print q-value table for all sweep columns. All scores are higher-is-better.
"""
function sweep_qvalue_summary(sweep_df::DataFrame;
                              thresholds=[0.001, 0.005, 0.01, 0.02, 0.05, 0.10])
    is_target = sweep_df.target .== 1.0
    n_targets = sum(is_target)
    n_decoys  = sum(.!is_target)
    println("\nTargets: $n_targets  |  Decoys: $n_decoys")

    # Collect all lambda1 and eigengap columns
    score_cols = Symbol[]
    push!(score_cols, :lambda1_frac_raw, :eigengap_raw)
    for method_name in sort(collect(keys(APEX_METHODS)))
        for hw in WINDOW_SIZES
            push!(score_cols, Symbol("lambda1_$(method_name)_w$(2hw+1)"))
        end
    end
    for method_name in sort(collect(keys(APEX_METHODS)))
        for hw in WINDOW_SIZES
            push!(score_cols, Symbol("eigengap_$(method_name)_w$(2hw+1)"))
        end
    end

    thresh_strs = [string(round(t * 100, digits=1)) * "%" for t in thresholds]
    w = 40
    header = rpad("score", w) * join(rpad.(thresh_strs, 10))
    println("─"^(w + 10 * length(thresholds)))
    println(header)
    println("─"^(w + 10 * length(thresholds)))

    for col in score_cols
        qvals = compute_qvalues(sweep_df[!, col], is_target; higher_is_better=true)
        counts = [targets_at_qvalue(qvals, is_target, t) for t in thresholds]
        row = rpad(string(col), w) * join(rpad.(string.(counts), 10))
        println(row)
    end
    println("─"^(w + 10 * length(thresholds)))
end

"""
    plot_sweep_curves(sweep_df; metric=:lambda1, max_qval=0.05)

Plot q-value curves for one metric across all apex methods and window sizes.
"""
function plot_sweep_curves(sweep_df::DataFrame;
                           metric::Symbol=:lambda1,
                           max_qval::Float64=0.05,
                           n_points::Int=200)
    is_target = sweep_df.target .== 1.0
    thresholds = range(0.0, max_qval, length=n_points)

    p = plot(xlabel="q-value", ylabel="targets passing",
             title="$(metric): apex method × window size",
             legend=:outertopright, size=(1000, 600))

    # Raw baseline
    raw_col = Symbol("$(metric)_frac_raw")
    if raw_col in propertynames(sweep_df)
        qvals = compute_qvalues(sweep_df[!, raw_col], is_target; higher_is_better=true)
        counts = [targets_at_qvalue(qvals, is_target, t) for t in thresholds]
        plot!(p, collect(thresholds), counts, label="raw", lw=3, color=:black, ls=:dash)
    end

    colors = Dict(:sum_smooth => :blue, :best_fragment => :red, :topk_smooth => :green)
    styles = Dict(5 => :dot, 8 => :dashdot, 10 => :dash, 12 => :solid, 15 => :solid, 20 => :solid)
    widths = Dict(5 => 1.0, 8 => 1.0, 10 => 1.5, 12 => 1.5, 15 => 2.0, 20 => 2.0)

    for method_name in sort(collect(keys(APEX_METHODS)))
        for hw in WINDOW_SIZES
            col = Symbol("$(metric)_$(method_name)_w$(2hw+1)")
            col in propertynames(sweep_df) || continue
            qvals = compute_qvalues(sweep_df[!, col], is_target; higher_is_better=true)
            counts = [targets_at_qvalue(qvals, is_target, t) for t in thresholds]
            plot!(p, collect(thresholds), counts,
                  label="$(method_name) w=$(2hw+1)",
                  color=get(colors, method_name, :gray),
                  ls=get(styles, hw, :solid),
                  lw=get(widths, hw, 1.5))
        end
    end

    display(p)
    return p
end

# ============================================================
# Quick start
# ============================================================

println("""

Loaded fragcorr_fixed_window.jl
────────────────────────────────────────
Quick start:
    gdf = load_fragcorr_data()
    sweep_df = sweep_fixed_windows(gdf)
    sweep_qvalue_summary(sweep_df)

Plotting:
    plot_sweep_curves(sweep_df; metric=:lambda1)
    plot_sweep_curves(sweep_df; metric=:eigengap)
""")
