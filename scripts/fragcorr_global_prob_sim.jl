#=
FRAGCORR Global Prob Simulation
================================
Simulates Pioneer's global filter (log-odds combine → PEP → hybrid threshold)
using predictions from different scoring models. Tests whether FRAGCORR-augmented
models recover more unique precursors — and whether the extras are DIA-NN validated.

Requires: run fragcorr_model_comparison.jl then fragcorr_diann_recovery.jl first
(provides mc_df with :is_diann, and results with .conditions Dict)

Usage:
    include("scripts/fragcorr_model_comparison.jl")
    include("scripts/fragcorr_diann_recovery.jl")
    include("scripts/fragcorr_global_prob_sim.jl")
=#

# ─────────────────────────────────────────────────────────────────────
# Step 1: Log-odds combine (copied from getBestPrecursorsAccrossRuns.jl)
# ─────────────────────────────────────────────────────────────────────

"""
    gp_logodds_combine(probs, top_n) → Float64

Log-odds average of top-N probabilities, converted back to probability space.
Mirrors Pioneer's `_logodds_combine` from FirstPassSearch.
"""
function gp_logodds_combine(probs::Vector{Float64}, top_n::Int)::Float64
    isempty(probs) && return 0.0
    n = min(length(probs), top_n)
    sorted = sort(probs; rev=true)
    selected = @view sorted[1:n]
    eps = 1e-6
    lo = log.(clamp.(selected, 0.1, 1 - eps) ./ (1 .- clamp.(selected, 0.1, 1 - eps)))
    avg = sum(lo) / n
    return 1.0 / (1 + exp(-avg))
end

# ─────────────────────────────────────────────────────────────────────
# Step 2: Compute global probs for each model condition
# ─────────────────────────────────────────────────────────────────────

"""
    compute_global_probs(mc_df, scores, n_files) → DataFrame

For a vector of per-PSM scores (same length as mc_df), compute one global
probability per unique precursor via:
  1. max score per (precursor_idx, ms_file_idx)  — best PSM per precursor per file
  2. log-odds combine across files              — global prob per precursor

Returns a DataFrame with columns: precursor_idx, target, is_diann, global_prob
"""
function compute_global_probs(mc_df::DataFrame, scores::Vector{Float64}, n_files::Int)
    top_n = max(1, floor(Int, sqrt(n_files)))

    # Step 2a: Best score per (precursor_idx, ms_file_idx)
    # Use a Dict of Dicts: precursor_idx → (ms_file_idx → max_score)
    prec_file_best = Dict{UInt32, Dict{UInt32, Float64}}()
    prec_target = Dict{UInt32, Bool}()
    prec_diann  = Dict{UInt32, Bool}()

    for i in 1:nrow(mc_df)
        pidx = UInt32(mc_df.precursor_idx[i])
        fidx = UInt32(mc_df.ms_file_idx[i])
        s    = scores[i]

        if !haskey(prec_file_best, pidx)
            prec_file_best[pidx] = Dict{UInt32, Float64}(fidx => s)
            prec_target[pidx] = mc_df.target[i]
            prec_diann[pidx]  = mc_df.is_diann[i]
        else
            fb = prec_file_best[pidx]
            fb[fidx] = max(get(fb, fidx, -Inf), s)
            # is_diann is true if ANY PSM for this precursor is DIA-NN validated
            if mc_df.is_diann[i]
                prec_diann[pidx] = true
            end
        end
    end

    # Step 2b: Log-odds combine across files for each precursor
    n_prec = length(prec_file_best)
    out_pidx   = Vector{UInt32}(undef, n_prec)
    out_target = Vector{Bool}(undef, n_prec)
    out_diann  = Vector{Bool}(undef, n_prec)
    out_gprob  = Vector{Float64}(undef, n_prec)

    for (i, (pidx, file_dict)) in enumerate(prec_file_best)
        probs_vec = collect(values(file_dict))
        out_pidx[i]   = pidx
        out_target[i] = prec_target[pidx]
        out_diann[i]  = prec_diann[pidx]
        out_gprob[i]  = gp_logodds_combine(probs_vec, top_n)
    end

    return DataFrame(
        precursor_idx = out_pidx,
        target        = out_target,
        is_diann      = out_diann,
        global_prob   = out_gprob,
    )
end

# ─────────────────────────────────────────────────────────────────────
# Step 3: Simulate Pioneer's hybrid filter
# ─────────────────────────────────────────────────────────────────────

"""
    simulate_global_filter(gp_df; global_pep_threshold=0.5, qvalue_threshold=0.15, min_floor=50_000)

Given a DataFrame with (precursor_idx, target, is_diann, global_prob),
simulate Pioneer's hybrid filter. Returns a NamedTuple with counts at
various q-value thresholds plus the hybrid filter outcome.
"""
function simulate_global_filter(gp_df::DataFrame;
                                 global_pep_threshold::Float64=0.5,
                                 qvalue_threshold::Float64=0.15,
                                 min_floor::Int=50_000)
    # Compute q-values on global probs (higher global_prob = better)
    qvals = compute_qvalues(gp_df.global_prob, gp_df.target; higher_is_better=true)

    # Counts at various q-value thresholds
    thresholds = [0.01, 0.05, 0.10, 0.15, 0.50]
    qval_counts = Dict{Float64, Int}()
    qval_diann  = Dict{Float64, Int}()
    for t in thresholds
        passing = qvals .<= t
        qval_counts[t] = count(passing .& gp_df.target)
        qval_diann[t]  = count(passing .& gp_df.target .& gp_df.is_diann)
    end

    # Raw count at global_prob ≥ threshold (analogous to PEP ≤ threshold)
    prob_passing = gp_df.global_prob .>= global_pep_threshold
    n_prob_pass  = count(prob_passing .& gp_df.target)

    # Hybrid filter simulation: sort by global_prob descending
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

    n_passing_threshold = count(gp_df.global_prob .>= global_pep_threshold)
    n_min_floor = min(min_floor, nrow(gp_df))
    n_to_keep = max(n_passing_threshold, n_qvalue_based, n_min_floor)

    # Count targets and DIA-NN among kept precursors
    kept_indices = order[1:min(n_to_keep, length(order))]
    hybrid_targets = count(i -> gp_df.target[i], kept_indices)
    hybrid_diann   = count(i -> gp_df.target[i] && gp_df.is_diann[i], kept_indices)

    return (
        qval_counts       = qval_counts,
        qval_diann        = qval_diann,
        n_prob_pass       = n_prob_pass,
        n_qvalue_based    = n_qvalue_based,
        n_passing_threshold = n_passing_threshold,
        n_min_floor       = n_min_floor,
        n_to_keep         = n_to_keep,
        hybrid_targets    = hybrid_targets,
        hybrid_diann      = hybrid_diann,
        qvals             = qvals,
        gp_df             = gp_df,
    )
end

# ─────────────────────────────────────────────────────────────────────
# Main execution
# ─────────────────────────────────────────────────────────────────────

println("\n╔══════════════════════════════════════════════════════════╗")
println("║  FRAGCORR Global Prob Simulation                        ║")
println("╚══════════════════════════════════════════════════════════╝")

# Determine number of unique files
n_files = length(unique(mc_df.ms_file_idx))
top_n = max(1, floor(Int, sqrt(n_files)))
println("\nFiles: $n_files  |  top_n for log-odds: $top_n")
println("Unique precursors in mc_df: $(length(unique(mc_df.precursor_idx)))")

# ── Compute global probs for each condition ──
ordered_names = [
    "LightGBM baseline", "LightGBM augmented",
    "Probit baseline", "Probit augmented",
    "Pipeline prob (ref)", "eigengap_raw (unsup)",
]

gp_results = Dict{String, NamedTuple}()

println("\nComputing global probs and simulating filter...")
for name in ordered_names
    scores = results.conditions[name]
    gp_df = compute_global_probs(mc_df, scores, n_files)
    sim = simulate_global_filter(gp_df)
    gp_results[name] = sim
    n_prec = nrow(gp_df)
    n_tgt  = count(gp_df.target)
    n_dec  = n_prec - n_tgt
    println("  $name: $n_prec unique precursors ($n_tgt T, $n_dec D) → hybrid keeps $(sim.n_to_keep) ($(sim.hybrid_targets) targets)")
end

# ─────────────────────────────────────────────────────────────────────
# Step 4: Output tables
# ─────────────────────────────────────────────────────────────────────

# Table 1: Unique target precursors passing at q-value thresholds
thresholds = [0.01, 0.05, 0.10, 0.15, 0.50]
thresh_strs = [string(round(t * 100, digits=0)) * "%" for t in thresholds]

println("\n\n" * "═"^80)
println("TABLE 1: Unique TARGET precursors passing global filter simulation")
println("═"^80)
println(rpad("Condition", 32), join(rpad.("q≤" .* thresh_strs, 10)), "hybrid")
println("─"^90)

for name in ordered_names
    sim = gp_results[name]
    counts = [sim.qval_counts[t] for t in thresholds]
    println(rpad(name, 32), join(rpad.(string.(counts), 10)), sim.hybrid_targets)
end
println("─"^90)

# Table 2: DIA-NN recovery at each threshold
println("\n\n" * "═"^80)
println("TABLE 2: DIA-NN validated targets at q-value thresholds")
println("═"^80)
println(rpad("Condition", 32), join(rpad.("q≤" .* thresh_strs, 10)))
println("─"^80)

for name in ordered_names
    sim = gp_results[name]
    counts = [sim.qval_diann[t] for t in thresholds]
    println(rpad(name, 32), join(rpad.(string.(counts), 10)))
end
println("─"^80)

# Table 3: DIA-NN overlap % at q ≤ 15% (matching Pioneer's D/T floor)
println("\n\n" * "═"^60)
println("TABLE 3: DIA-NN Overlap at q ≤ 15% (Pioneer's D/T floor)")
println("═"^60)
println(rpad("Condition", 32), rpad("Total", 10), rpad("DIA-NN", 10), "DIA-NN%")
println("─"^60)

for name in ordered_names
    sim = gp_results[name]
    total = sim.qval_counts[0.15]
    diann = sim.qval_diann[0.15]
    pct = total > 0 ? round(100.0 * diann / total, digits=1) : 0.0
    println(rpad(name, 32), rpad(string(total), 10), rpad(string(diann), 10), "$(pct)%")
end
println("─"^60)

# ─────────────────────────────────────────────────────────────────────
# Step 5: Venn decomposition (baseline vs augmented global probs)
# ─────────────────────────────────────────────────────────────────────

println("\n\n" * "═"^70)
println("TABLE 4: Venn Decomposition — Global probs at q ≤ 15%")
println("        (LightGBM baseline vs augmented)")
println("═"^70)

# Get the global prob DataFrames for baseline and augmented
gp_base = gp_results["LightGBM baseline"].gp_df
gp_aug  = gp_results["LightGBM augmented"].gp_df
qv_base = gp_results["LightGBM baseline"].qvals
qv_aug  = gp_results["LightGBM augmented"].qvals

# Build sets of passing target precursor_idx at q ≤ 0.15
pass_base_set = Set{UInt32}()
pass_aug_set  = Set{UInt32}()
diann_base_set = Set{UInt32}()
diann_aug_set  = Set{UInt32}()

for i in 1:nrow(gp_base)
    gp_base.target[i] || continue
    if qv_base[i] <= 0.15
        push!(pass_base_set, gp_base.precursor_idx[i])
        gp_base.is_diann[i] && push!(diann_base_set, gp_base.precursor_idx[i])
    end
end
for i in 1:nrow(gp_aug)
    gp_aug.target[i] || continue
    if qv_aug[i] <= 0.15
        push!(pass_aug_set, gp_aug.precursor_idx[i])
        gp_aug.is_diann[i] && push!(diann_aug_set, gp_aug.precursor_idx[i])
    end
end

both_set     = intersect(pass_base_set, pass_aug_set)
aug_only_set = setdiff(pass_aug_set, pass_base_set)
base_only_set = setdiff(pass_base_set, pass_aug_set)

n_both     = length(both_set)
n_aug_only = length(aug_only_set)
n_base_only = length(base_only_set)

diann_both     = length(intersect(both_set, diann_aug_set ∪ diann_base_set))
diann_aug_only = length(intersect(aug_only_set, diann_aug_set))
diann_base_only = length(intersect(base_only_set, diann_base_set))

pct_both     = n_both > 0 ? round(100.0 * diann_both / n_both, digits=1) : 0.0
pct_aug_only = n_aug_only > 0 ? round(100.0 * diann_aug_only / n_aug_only, digits=1) : 0.0
pct_base_only = n_base_only > 0 ? round(100.0 * diann_base_only / n_base_only, digits=1) : 0.0

println("Both models:          $(rpad(n_both, 8))  ($diann_both = $(pct_both)% DIA-NN validated)")
println("Aug only (FRAGCORR):  $(rpad(n_aug_only, 8))  ($diann_aug_only = $(pct_aug_only)% DIA-NN validated)")
println("Base only:            $(rpad(n_base_only, 8))  ($diann_base_only = $(pct_base_only)% DIA-NN validated)")
println("─"^70)

# Also do the Venn at q ≤ 1% for a stricter threshold
println("\n── Venn at q ≤ 1% ──")
pass_base_1pct = Set{UInt32}()
pass_aug_1pct  = Set{UInt32}()
diann_base_1pct = Set{UInt32}()
diann_aug_1pct  = Set{UInt32}()

for i in 1:nrow(gp_base)
    gp_base.target[i] || continue
    if qv_base[i] <= 0.01
        push!(pass_base_1pct, gp_base.precursor_idx[i])
        gp_base.is_diann[i] && push!(diann_base_1pct, gp_base.precursor_idx[i])
    end
end
for i in 1:nrow(gp_aug)
    gp_aug.target[i] || continue
    if qv_aug[i] <= 0.01
        push!(pass_aug_1pct, gp_aug.precursor_idx[i])
        gp_aug.is_diann[i] && push!(diann_aug_1pct, gp_aug.precursor_idx[i])
    end
end

both_1pct     = intersect(pass_base_1pct, pass_aug_1pct)
aug_only_1pct = setdiff(pass_aug_1pct, pass_base_1pct)
base_only_1pct = setdiff(pass_base_1pct, pass_aug_1pct)

n_both_1     = length(both_1pct)
n_aug_only_1 = length(aug_only_1pct)
n_base_only_1 = length(base_only_1pct)

diann_both_1     = length(intersect(both_1pct, diann_aug_1pct ∪ diann_base_1pct))
diann_aug_only_1 = length(intersect(aug_only_1pct, diann_aug_1pct))
diann_base_only_1 = length(intersect(base_only_1pct, diann_base_1pct))

pct_both_1     = n_both_1 > 0 ? round(100.0 * diann_both_1 / n_both_1, digits=1) : 0.0
pct_aug_only_1 = n_aug_only_1 > 0 ? round(100.0 * diann_aug_only_1 / n_aug_only_1, digits=1) : 0.0
pct_base_only_1 = n_base_only_1 > 0 ? round(100.0 * diann_base_only_1 / n_base_only_1, digits=1) : 0.0

println("Both models:          $(rpad(n_both_1, 8))  ($diann_both_1 = $(pct_both_1)% DIA-NN validated)")
println("Aug only (FRAGCORR):  $(rpad(n_aug_only_1, 8))  ($diann_aug_only_1 = $(pct_aug_only_1)% DIA-NN validated)")
println("Base only:            $(rpad(n_base_only_1, 8))  ($diann_base_only_1 = $(pct_base_only_1)% DIA-NN validated)")
println("─"^70)

# ─────────────────────────────────────────────────────────────────────
# Step 6: Plots
# ─────────────────────────────────────────────────────────────────────

println("\nGenerating global prob simulation plots...")

n_points = 200
qv_range = range(0.0, 0.20, length=n_points)

colors_map = Dict(
    "LightGBM baseline"    => :blue,
    "LightGBM augmented"   => :red,
    "Probit baseline"      => :cyan,
    "Probit augmented"     => :magenta,
    "Pipeline prob (ref)"  => :gray,
    "eigengap_raw (unsup)" => :orange,
)
styles_map = Dict(
    "LightGBM baseline"    => :solid,
    "LightGBM augmented"   => :solid,
    "Probit baseline"      => :dash,
    "Probit augmented"     => :dash,
    "Pipeline prob (ref)"  => :dot,
    "eigengap_raw (unsup)" => :dot,
)

# Plot 1: Unique target precursors passing global filter at q-value threshold
p_global_total = plot(xlabel="q-value (global)", ylabel="Unique target precursors",
                      title="Global Prob Simulation: unique targets at q-value threshold",
                      legend=:topleft, size=(900, 550))

for name in ordered_names
    sim = gp_results[name]
    gp_df = sim.gp_df
    qvals = sim.qvals
    counts = Int[]
    for t in qv_range
        passing = qvals .<= t
        push!(counts, count(passing .& gp_df.target))
    end
    plot!(p_global_total, collect(qv_range), counts;
          label=name, lw=2.5,
          color=colors_map[name], ls=styles_map[name])
end
display(p_global_total)

# Plot 2: DIA-NN validated targets at global q-value threshold
p_global_diann = plot(xlabel="q-value (global)", ylabel="DIA-NN validated unique targets",
                      title="Global Prob Simulation: DIA-NN validated at q-value threshold",
                      legend=:topleft, size=(900, 550))

for name in ordered_names
    sim = gp_results[name]
    gp_df = sim.gp_df
    qvals = sim.qvals
    counts = Int[]
    for t in qv_range
        passing = qvals .<= t
        push!(counts, count(passing .& gp_df.target .& gp_df.is_diann))
    end
    plot!(p_global_diann, collect(qv_range), counts;
          label=name, lw=2.5,
          color=colors_map[name], ls=styles_map[name])
end
display(p_global_diann)

# Plot 3: LightGBM baseline vs augmented — total + DIA-NN overlay
p_global_lgbm = plot(xlabel="q-value (global)", ylabel="Unique target precursors",
                     title="Global Prob: LightGBM baseline vs augmented (total + DIA-NN)",
                     legend=:topleft, size=(900, 550))

for name in ["LightGBM baseline", "LightGBM augmented"]
    sim = gp_results[name]
    gp_df = sim.gp_df
    qvals = sim.qvals
    col = colors_map[name]

    total_counts = [count((qvals .<= t) .& gp_df.target) for t in qv_range]
    diann_counts = [count((qvals .<= t) .& gp_df.target .& gp_df.is_diann) for t in qv_range]

    plot!(p_global_lgbm, collect(qv_range), total_counts;
          label="$name (total)", lw=2.5, color=col, ls=:solid)
    plot!(p_global_lgbm, collect(qv_range), diann_counts;
          label="$name (DIA-NN)", lw=2.0, color=col, ls=:dash)
end
display(p_global_lgbm)

# Plot 4: Comparison — PSM-level q-values vs global-level q-values (LightGBM augmented)
p_psm_vs_global = plot(xlabel="q-value", ylabel="Targets passing",
                       title="PSM-level vs Global-level: LightGBM augmented",
                       legend=:topleft, size=(900, 550))

# PSM-level (from model comparison results)
psm_qvals = results.qvals["LightGBM augmented"]
psm_is_target = mc_df.target
psm_counts = [targets_at_qvalue(psm_qvals, psm_is_target, t) for t in qv_range]
plot!(p_psm_vs_global, collect(qv_range), psm_counts;
      label="PSM-level (pooled PSMs)", lw=2.5, color=:red, ls=:dot)

# Global-level
sim = gp_results["LightGBM augmented"]
global_counts = [count((sim.qvals .<= t) .& sim.gp_df.target) for t in qv_range]
plot!(p_psm_vs_global, collect(qv_range), global_counts;
      label="Global-level (unique precursors)", lw=2.5, color=:red, ls=:solid)

# Also add baseline for reference
psm_qvals_base = results.qvals["LightGBM baseline"]
psm_counts_base = [targets_at_qvalue(psm_qvals_base, psm_is_target, t) for t in qv_range]
plot!(p_psm_vs_global, collect(qv_range), psm_counts_base;
      label="PSM-level baseline (pooled)", lw=2.0, color=:blue, ls=:dot)

sim_base = gp_results["LightGBM baseline"]
global_counts_base = [count((sim_base.qvals .<= t) .& sim_base.gp_df.target) for t in qv_range]
plot!(p_psm_vs_global, collect(qv_range), global_counts_base;
      label="Global-level baseline (unique)", lw=2.0, color=:blue, ls=:solid)

display(p_psm_vs_global)

# ─────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────

println("\n\nPlots available:")
println("  p_global_total     — Unique target precursors at global q-value threshold")
println("  p_global_diann     — DIA-NN validated targets at global q-value threshold")
println("  p_global_lgbm      — LightGBM baseline vs augmented (total + DIA-NN)")
println("  p_psm_vs_global    — PSM-level vs global-level q-value curves")
println("\nAccess results as `gp_results` Dict (keyed by condition name).")
println("  gp_results[\"LightGBM augmented\"].qval_counts  — targets at each q-val threshold")
println("  gp_results[\"LightGBM augmented\"].gp_df        — per-precursor global probs DataFrame")
println("  gp_results[\"LightGBM augmented\"].hybrid_targets — targets kept by hybrid filter")

println("\nDone! Global prob simulation complete.")
