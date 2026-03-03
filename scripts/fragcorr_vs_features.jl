#=
FRAGCORR vs Traditional Feature Q-value Comparison
===================================================
Compares unsupervised FRAGCORR scores (multi-scan fragment correlation)
against traditional per-scan scoring features (scribe, spectral_contrast, etc.)
using target-decoy q-value curves.

Usage:
    include("scripts/fragcorr_vs_features.jl")
=#

# Reuse all FRAGCORR scoring + q-value utilities from the base script
include(joinpath(@__DIR__, "fragcorr_unsupervised_scores.jl"))

const DATA_DIR = "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/OlsenEclipse_fragcorr_test/temp_data"
const PSM_DIR  = joinpath(DATA_DIR, "first_pass_psms")
const CORR_DIR = joinpath(DATA_DIR, "first_pass_corr")

# Traditional feature metadata: (column_name, higher_is_better)
const TRADITIONAL_FEATURES = [
    (:scribe,            true),
    (:spectral_contrast, true),
    (:city_block,        false),
    (:entropy_score,     true),
    (:poisson,           true),
    (:y_count,           true),
    (:err_norm,          false),
    (:score,             true),
    (:prob,              true),
]

# Best FRAGCORR scores to highlight in comparison plots
const BEST_FRAGCORR = [
    :lambda1_frac_raw, :eigengap_raw, :median_corr_raw,
    :best_frag_corr_sum, :best_frag_corr_mean,
    :local_corr, :min_corr,
    :peak_shape, :elution_cosine,
]

"""
    load_and_join(filename) → DataFrame

Load a PSM arrow + CORR arrow pair, compute FRAGCORR scores, and inner join
on precursor_idx. Returns one row per precursor with both feature sets.
"""
function load_and_join(filename::String)
    psm_path  = joinpath(PSM_DIR, filename)
    corr_path = joinpath(CORR_DIR, filename)

    # Load PSMs (one row per precursor)
    psm_df = DataFrame(Tables.columntable(Arrow.Table(psm_path)))

    # Load CORR data and compute FRAGCORR scores
    gdf = load_fragcorr_data(corr_path)
    fragcorr_df = score_all_precursors(gdf)

    # Inner join — keeps only precursors present in both
    combined = innerjoin(psm_df, fragcorr_df;
                         on=[:precursor_idx, :target],
                         makeunique=true)
    println("  Joined: $(nrow(psm_df)) PSMs × $(nrow(fragcorr_df)) FRAGCORR → $(nrow(combined)) combined")
    return combined
end

"""
    load_all_files() → DataFrame

Load and join all 6 raw file pairs, concatenate into a single DataFrame.
"""
function load_all_files()
    files = filter(f -> endswith(f, ".arrow"), readdir(PSM_DIR))
    # Only process files that exist in both directories
    corr_files = Set(readdir(CORR_DIR))
    files = filter(f -> f in corr_files, files)

    println("\nLoading $(length(files)) file pairs...")
    dfs = DataFrame[]
    for f in files
        println("\n── $f ──")
        push!(dfs, load_and_join(f))
    end

    combined = vcat(dfs...; cols=:union)
    n_target = sum(combined.target)
    n_decoy  = sum(.!combined.target)
    println("\n════════════════════════════════════════")
    println("Combined: $(nrow(combined)) rows ($n_target targets, $n_decoy decoys)")
    println("════════════════════════════════════════")
    return combined
end

"""
    qvalue_comparison(df) → nothing

Compute q-values for all traditional features and best FRAGCORR scores,
print a summary table, and plot comparative curves.
"""
function qvalue_comparison(df::DataFrame;
                           thresholds=[0.001, 0.005, 0.01, 0.02, 0.05, 0.10])
    is_target = df.target
    n_targets = sum(is_target)
    n_decoys  = sum(.!is_target)

    # Collect all scores to compare: (name, qvals)
    all_scores = Tuple{Symbol, Vector{Float64}}[]

    # Traditional features
    println("\n── Traditional Features ──")
    for (col, hib) in TRADITIONAL_FEATURES
        qvals = compute_qvalues(df[!, col], is_target; higher_is_better=hib)
        push!(all_scores, (col, qvals))
    end

    # FRAGCORR scores
    println("── FRAGCORR Scores ──")
    for col in BEST_FRAGCORR
        hib = col in HIGHER_IS_BETTER
        qvals = compute_qvalues(df[!, col], is_target; higher_is_better=hib)
        push!(all_scores, (col, qvals))
    end

    # Print summary table
    thresh_strs = [string(round(t * 100, digits=1)) * "%" for t in thresholds]
    w = 36 + 10 * length(thresholds)
    println("\nTargets: $n_targets  |  Decoys: $n_decoys")
    println("─"^w)
    println(rpad("feature", 36), join(rpad.(thresh_strs, 10)))
    println("─"^w)

    # Traditional section
    println("  Traditional:")
    for (col, _) in TRADITIONAL_FEATURES
        qvals = [q for (c, q) in all_scores if c == col][1]
        counts = [targets_at_qvalue(qvals, is_target, t) for t in thresholds]
        println("    ", rpad(string(col), 32), join(rpad.(string.(counts), 10)))
    end

    # FRAGCORR section
    println("  FRAGCORR:")
    for col in BEST_FRAGCORR
        qvals = [q for (c, q) in all_scores if c == col][1]
        counts = [targets_at_qvalue(qvals, is_target, t) for t in thresholds]
        println("    ", rpad(string(col), 32), join(rpad.(string.(counts), 10)))
    end
    println("─"^w)

    # ── Plot: combined q-value curves ──
    n_points = 200
    qv_range = range(0.0, 0.10, length=n_points)

    p_all = plot(xlabel="q-value", ylabel="targets passing",
                 title="Traditional vs FRAGCORR: targets at q-value threshold",
                 legend=:topleft, size=(900, 550))

    # Traditional features: solid lines
    for (i, (col, hib)) in enumerate(TRADITIONAL_FEATURES)
        qvals = [q for (c, q) in all_scores if c == col][1]
        counts = [targets_at_qvalue(qvals, is_target, t) for t in qv_range]
        plot!(p_all, collect(qv_range), counts,
              label="T: $(col)", lw=2, ls=:solid, color=i)
    end

    # FRAGCORR scores: dashed lines
    for (i, col) in enumerate(BEST_FRAGCORR)
        qvals = [q for (c, q) in all_scores if c == col][1]
        counts = [targets_at_qvalue(qvals, is_target, t) for t in qv_range]
        plot!(p_all, collect(qv_range), counts,
              label="F: $(col)", lw=2, ls=:dash,
              color=length(TRADITIONAL_FEATURES) + i)
    end

    display(p_all)

    # ── Separate panel plots for clarity ──
    p_trad = plot(xlabel="q-value", ylabel="targets passing",
                  title="Traditional Features Only",
                  legend=:topleft, size=(800, 500))
    for (i, (col, hib)) in enumerate(TRADITIONAL_FEATURES)
        qvals = [q for (c, q) in all_scores if c == col][1]
        counts = [targets_at_qvalue(qvals, is_target, t) for t in qv_range]
        plot!(p_trad, collect(qv_range), counts,
              label=string(col), lw=2, color=i)
    end

    p_frag = plot(xlabel="q-value", ylabel="targets passing",
                  title="FRAGCORR Scores Only",
                  legend=:topleft, size=(800, 500))
    for (i, col) in enumerate(BEST_FRAGCORR)
        qvals = [q for (c, q) in all_scores if c == col][1]
        counts = [targets_at_qvalue(qvals, is_target, t) for t in qv_range]
        plot!(p_frag, collect(qv_range), counts,
              label=string(col), lw=2, color=i)
    end

    p_panels = plot(p_trad, p_frag, layout=(1, 2), size=(1600, 500))
    display(p_panels)

    return (combined=p_all, traditional=p_trad, fragcorr=p_frag, panels=p_panels)
end

# ============================================================
# Main
# ============================================================

println("\n╔══════════════════════════════════════════════╗")
println("║  FRAGCORR vs Traditional Feature Comparison  ║")
println("╚══════════════════════════════════════════════╝")

combined_df = load_all_files()

# Check for NaN/missing issues
_has_nan(v) = eltype(v) <: AbstractFloat ? count(isnan, v) : 0
println("\nNaN check:")
for (col, _) in TRADITIONAL_FEATURES
    n_nan = _has_nan(combined_df[!, col])
    n_nan > 0 && println("  $col: $n_nan NaN values")
end
for col in BEST_FRAGCORR
    n_nan = _has_nan(combined_df[!, col])
    n_nan > 0 && println("  $col: $n_nan NaN values")
end
println("  (no output above = all clean)")

plots = qvalue_comparison(combined_df)

println("\nDone! Plots displayed.")
println("Access the combined DataFrame as `combined_df`.")
println("Re-run comparison: qvalue_comparison(combined_df)")
