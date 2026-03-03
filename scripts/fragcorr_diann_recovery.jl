#=
FRAGCORR DIA-NN Recovery Analysis
==================================
Validates whether the extra targets found by LightGBM augmented (with FRAGCORR
features) are real by checking overlap with DIA-NN's independently identified
precursors.

Requires: run fragcorr_model_comparison.jl first (provides mc_df, results)

Usage:
    include("scripts/fragcorr_model_comparison.jl")
    include("scripts/fragcorr_diann_recovery.jl")
=#

using CSV

# ─────────────────────────────────────────────────────────────────────
# Step 1: Canonical key matching infrastructure
# ─────────────────────────────────────────────────────────────────────
# (Copied from diann_vs_pioneer_all_experiments.jl)

const DR_ModTuple = Tuple{Int, Int}  # (position, unimod_id)
const DR_CanonicalKey = Tuple{String, Vector{DR_ModTuple}, Int}

function dr_parse_diann_modified_seq(mod_seq::AbstractString, charge::Int)::DR_CanonicalKey
    stripped = IOBuffer()
    mods = DR_ModTuple[]
    pos = 0
    i = 1
    while i <= ncodeunits(mod_seq)
        c = mod_seq[i]
        if c == '('
            close_idx = findnext(')', mod_seq, i)
            mod_str = mod_seq[i+1:close_idx-1]
            colon_idx = findlast(':', mod_str)
            unimod_id = parse(Int, mod_str[colon_idx+1:end])
            push!(mods, (pos, unimod_id))
            i = close_idx + 1
        else
            pos += 1
            write(stripped, c)
            i += 1
        end
    end
    sort!(mods)
    return (String(take!(stripped)), mods, charge)
end

function dr_parse_pioneer_structural_mods(mods_str::AbstractString)::Vector{DR_ModTuple}
    isempty(mods_str) && return DR_ModTuple[]
    mods = DR_ModTuple[]
    for m in eachmatch(r"\((\d+),\w,(\w+:\d+)\)", mods_str)
        pos = parse(Int, m.captures[1])
        mod_id_str = m.captures[2]
        colon_idx = findlast(':', mod_id_str)
        unimod_id = parse(Int, mod_id_str[colon_idx+1:end])
        push!(mods, (pos, unimod_id))
    end
    sort!(mods)
    return mods
end

function dr_pioneer_canonical_key(seq::AbstractString, structural_mods::AbstractString, charge::Integer)::DR_CanonicalKey
    mods = dr_parse_pioneer_structural_mods(structural_mods)
    return (String(seq), mods, Int(charge))
end

# ─────────────────────────────────────────────────────────────────────
# Step 2: Build key maps
# ─────────────────────────────────────────────────────────────────────

const DR_LIBRARY_PATH = "/Users/nathanwamsley/Data/SPEC_LIBS/3P_rank_based.poin"
const DR_DIANN_PR_MATRIX = joinpath("/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse",
    "DIA-NN_results", "OlsenExplorisThreeProteome500ng-11-24-2025-report.pr_matrix.tsv")
const DR_PIONEER_FILES = ["E45H50Y5_2", "E45H50Y5_3", "E45H50Y5_4",
                          "E5H50Y45_1", "E5H50Y45_2", "E5H50Y45_4"]

println("\n╔══════════════════════════════════════════════════════════╗")
println("║  FRAGCORR DIA-NN Recovery Analysis                      ║")
println("╚══════════════════════════════════════════════════════════╝")

# --- Pioneer library: precursor_idx → canonical key ---
println("\nLoading Pioneer library...")
lib_df = DataFrame(Arrow.Table(joinpath(DR_LIBRARY_PATH, "precursors_table.arrow")))

# Only build keys for precursor_idx values present in mc_df (targets only)
mc_pidx_set = Set{UInt32}(mc_df[mc_df.target, :precursor_idx])
println("  Unique target precursor_idx in mc_df: $(length(mc_pidx_set))")

pidx_to_key = Dict{UInt32, DR_CanonicalKey}()
for i in 1:nrow(lib_df)
    pidx = UInt32(i)
    pidx ∈ mc_pidx_set || continue
    lib_df.is_decoy[i] && continue
    seq = lib_df.sequence[i]
    smods_raw = lib_df.structural_mods[i]
    smods = ismissing(smods_raw) ? "" : smods_raw
    charge = Int(lib_df.prec_charge[i])
    pidx_to_key[pidx] = dr_pioneer_canonical_key(seq, smods, charge)
end
println("  Built $(length(pidx_to_key)) precursor_idx → key mappings")

# --- DIA-NN fair key set ---
println("\nLoading DIA-NN pr_matrix...")
diann_df = CSV.read(DR_DIANN_PR_MATRIX, DataFrame; delim='\t')
n_diann = nrow(diann_df)
println("  Loaded $n_diann precursors")

# Identify run columns (same logic as diann_vs_pioneer_all_experiments.jl)
metadata_cols = Set(["Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                     "First.Protein.Description", "Proteotypic",
                     "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge", "Precursor.Id"])
run_cols = [c for c in names(diann_df) if c ∉ metadata_cols]

# Filter to Pioneer's 6 files
function dr_diann_col_to_pioneer_file(col::AbstractString)::Union{String, Nothing}
    m = match(r"(E\d+H\d+Y\d+)_\d+SPD_DIA_(\d+)", col)
    m === nothing && return nothing
    return "$(m.captures[1])_$(m.captures[2])"
end

pioneer_run_cols = String[]
for col in run_cols
    pf = dr_diann_col_to_pioneer_file(col)
    pf !== nothing && pf ∈ Set(DR_PIONEER_FILES) && push!(pioneer_run_cols, col)
end
println("  Pioneer run columns matched: $(length(pioneer_run_cols))")

# Parse DIA-NN canonical keys and build fair set (detected in at least one Pioneer run)
diann_fair_keys = Set{DR_CanonicalKey}()
diann_all_keys = Set{DR_CanonicalKey}()
for i in 1:n_diann
    key = dr_parse_diann_modified_seq(
        diann_df[i, "Modified.Sequence"],
        Int(diann_df[i, "Precursor.Charge"])
    )
    push!(diann_all_keys, key)

    # Check if detected in any Pioneer run
    for col in pioneer_run_cols
        val = diann_df[i, col]
        if !ismissing(val) && val isa Number && val > 0
            push!(diann_fair_keys, key)
            break
        end
    end
end
println("  DIA-NN total unique keys: $(length(diann_all_keys))")
println("  DIA-NN fair keys (detected in Pioneer runs): $(length(diann_fair_keys))")

# ─────────────────────────────────────────────────────────────────────
# Step 3: Annotate mc_df with is_diann
# ─────────────────────────────────────────────────────────────────────

println("\nAnnotating mc_df...")
is_diann = falses(nrow(mc_df))
let n_mapped = 0, n_matched = 0
    for i in 1:nrow(mc_df)
        mc_df.target[i] || continue
        pidx = UInt32(mc_df.precursor_idx[i])
        key = get(pidx_to_key, pidx, nothing)
        key === nothing && continue
        n_mapped += 1
        if key ∈ diann_fair_keys
            is_diann[i] = true
            n_matched += 1
        end
    end
    println("  Targets mapped to keys: $n_mapped")
    println("  Targets matching DIA-NN fair set: $n_matched")
end
mc_df.is_diann = is_diann

# ─────────────────────────────────────────────────────────────────────
# Step 4: Recovery analysis
# ─────────────────────────────────────────────────────────────────────

println("\n" * "═"^70)
println("DIA-NN OVERLAP AT MULTIPLE Q-VALUE THRESHOLDS")
println("═"^70)

thresholds = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10]
ordered_names = [
    "LightGBM baseline", "LightGBM augmented",
    "Probit baseline", "Probit augmented",
    "Pipeline prob (ref)", "eigengap_raw (unsup)",
]

is_target = mc_df.target
is_diann_col = mc_df.is_diann

# Table 1: Per-condition DIA-NN overlap
thresh_strs = [string(round(t * 100, digits=1)) * "%" for t in thresholds]
w = 32 + 18 * length(thresholds)
println("\n" * rpad("Condition", 32) * join([rpad("$t (total/diann/%)", 18) for t in thresh_strs]))
println("─"^w)

# Store results for later use
recovery_data = Dict{String, Dict{Float64, NamedTuple}}()

for name in ordered_names
    qvals = results.qvals[name]
    recovery_data[name] = Dict{Float64, NamedTuple}()

    row_parts = String[rpad(name, 32)]
    for t in thresholds
        passing = qvals .<= t
        total_targets = count(passing .& is_target)
        diann_targets = count(passing .& is_target .& is_diann_col)
        pct = total_targets > 0 ? 100.0 * diann_targets / total_targets : 0.0
        push!(row_parts, rpad("$total_targets/$diann_targets/$(round(pct, digits=1))%", 18))
        recovery_data[name][t] = (total=total_targets, diann=diann_targets, pct=pct)
    end
    println(join(row_parts))
end
println("─"^w)

# Clean summary table at 1% FDR
println("\n\n" * "═"^60)
println("TABLE 1: DIA-NN Overlap at 1% FDR")
println("═"^60)
println(rpad("Condition", 32) * rpad("Total", 10) * rpad("DIA-NN", 10) * "DIA-NN%")
println("─"^60)
for name in ordered_names
    r = recovery_data[name][0.01]
    println(rpad(name, 32) * rpad(string(r.total), 10) * rpad(string(r.diann), 10) *
            "$(round(r.pct, digits=1))%")
end
println("─"^60)

# ─────────────────────────────────────────────────────────────────────
# Step 5: Venn decomposition (LightGBM baseline vs augmented at 1% FDR)
# ─────────────────────────────────────────────────────────────────────

println("\n\n" * "═"^60)
println("TABLE 2: Venn Decomposition (LightGBM baseline vs augmented, 1% FDR)")
println("═"^60)

qv_base = results.qvals["LightGBM baseline"]
qv_aug  = results.qvals["LightGBM augmented"]

pass_base = (qv_base .<= 0.01) .& is_target
pass_aug  = (qv_aug  .<= 0.01) .& is_target

both_models = pass_base .& pass_aug
aug_only    = pass_aug .& .!pass_base
base_only   = pass_base .& .!pass_aug

n_both      = count(both_models)
n_aug_only  = count(aug_only)
n_base_only = count(base_only)

diann_both      = count(both_models .& is_diann_col)
diann_aug_only  = count(aug_only .& is_diann_col)
diann_base_only = count(base_only .& is_diann_col)

pct_aug_only_diann  = n_aug_only > 0 ? 100.0 * diann_aug_only / n_aug_only : 0.0
pct_base_only_diann = n_base_only > 0 ? 100.0 * diann_base_only / n_base_only : 0.0
pct_both_diann      = n_both > 0 ? 100.0 * diann_both / n_both : 0.0

println("Both models:          $(rpad(n_both, 8))  ($diann_both = $(round(pct_both_diann, digits=1))% are DIA-NN validated)")
println("Aug only (FRAGCORR):  $(rpad(n_aug_only, 8))  ($diann_aug_only = $(round(pct_aug_only_diann, digits=1))% are DIA-NN validated)")
println("Base only:            $(rpad(n_base_only, 8))  ($diann_base_only = $(round(pct_base_only_diann, digits=1))% are DIA-NN validated)")
println("─"^60)

# ─────────────────────────────────────────────────────────────────────
# Step 6: Q-value curves with DIA-NN overlay
# ─────────────────────────────────────────────────────────────────────

println("\nGenerating q-value curves with DIA-NN overlay...")

n_points = 200
qv_range = range(0.0, 0.10, length=n_points)

# Colors and styles
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

# Plot: total targets and DIA-NN targets for each condition
p_total = plot(xlabel="q-value", ylabel="Targets passing",
               title="Total targets at q-value threshold",
               legend=:topleft, size=(900, 550))

p_diann = plot(xlabel="q-value", ylabel="DIA-NN validated targets",
               title="DIA-NN validated targets at q-value threshold",
               legend=:topleft, size=(900, 550))

p_pct = plot(xlabel="q-value", ylabel="DIA-NN overlap %",
             title="DIA-NN overlap percentage at q-value threshold",
             legend=:bottomleft, size=(900, 550), ylims=(0, 100))

for name in ordered_names
    qvals = results.qvals[name]

    total_counts = Int[]
    diann_counts = Int[]
    pct_vals = Float64[]

    for t in qv_range
        passing = qvals .<= t
        total = count(passing .& is_target)
        diann = count(passing .& is_target .& is_diann_col)
        pct = total > 0 ? 100.0 * diann / total : 0.0
        push!(total_counts, total)
        push!(diann_counts, diann)
        push!(pct_vals, pct)
    end

    col = colors_map[name]
    ls = styles_map[name]

    plot!(p_total, collect(qv_range), total_counts;
          label=name, lw=2.5, color=col, ls=ls)
    plot!(p_diann, collect(qv_range), diann_counts;
          label=name, lw=2.5, color=col, ls=ls)
    plot!(p_pct, collect(qv_range), pct_vals;
          label=name, lw=2.5, color=col, ls=ls)
end

# Combined 2-panel: baseline vs augmented (total + DIA-NN)
p_combined = plot(xlabel="q-value", ylabel="Targets passing",
                  title="LightGBM: Total vs DIA-NN validated targets",
                  legend=:topleft, size=(900, 550))

for name in ["LightGBM baseline", "LightGBM augmented"]
    qvals = results.qvals[name]
    col = colors_map[name]

    total_counts = [count((qvals .<= t) .& is_target) for t in qv_range]
    diann_counts = [count((qvals .<= t) .& is_target .& is_diann_col) for t in qv_range]

    plot!(p_combined, collect(qv_range), total_counts;
          label="$name (total)", lw=2.5, color=col, ls=:solid)
    plot!(p_combined, collect(qv_range), diann_counts;
          label="$name (DIA-NN)", lw=2.0, color=col, ls=:dash)
end

display(p_combined)

println("\nPlots available:")
println("  p_total    — Total targets at q-value threshold (all conditions)")
println("  p_diann    — DIA-NN validated targets (all conditions)")
println("  p_pct      — DIA-NN overlap percentage (all conditions)")
println("  p_combined — LightGBM baseline vs augmented with DIA-NN overlay")

println("\nDone! DIA-NN recovery analysis complete.")
println("Access annotated DataFrame as `mc_df` (now has :is_diann column).")
