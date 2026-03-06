#=
DIA-NN Recovery Analysis
========================
Compares Pioneer second-pass PSMs against DIA-NN identifications to measure
how many DIA-NN precursors Pioneer recovers at various FDR thresholds.

Usage:
    include("scripts/diann_recovery_analysis.jl")

    # Run from a Pioneer config JSON (reads paths.results, writes report there)
    results = diann_recovery_analysis("data/OlsenEclipse_nofilters_noMBR_20a3d62c.json")

    # Run with explicit results directory
    results = diann_recovery_analysis(
        pioneer_results_dir = "/path/to/results"
    )

    # Print a previously computed result
    print_recovery_table(results)
=#

using Arrow, DataFrames, Tables, Parquet2, Printf, JSON

# ── Default paths ──
const DIANN_LIBRARY_PATH = "/Users/nathanwamsley/Data/SPEC_LIBS/3P_rank_based.poin"
const DIANN_REPORT_PATH  = "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/DIA-NN_results/OlsenExplorisThreeProteome500ng-11-24-2025-report.parquet"
const DIANN_PIONEER_FILES = ["E45H50Y5_2", "E45H50Y5_3", "E45H50Y5_4",
                             "E5H50Y45_1", "E5H50Y45_2", "E5H50Y45_4"]
const DIANN_FDR_THRESHOLDS = [0.01, 0.05, 0.10, 0.15, 0.20, 0.50, 1.00]

# ── Canonical key: (stripped_seq, [(pos,unimod_id)...], charge) ──
const _DRA_CKey = Tuple{String, Vector{Tuple{Int,Int}}, Int}

function _dra_parse_diann_seq(mod_seq::AbstractString, charge::Int)::_DRA_CKey
    stripped = IOBuffer()
    mods = Tuple{Int,Int}[]
    pos = 0; i = 1
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

function _dra_pioneer_key(seq::AbstractString, structural_mods::AbstractString, charge::Integer)::_DRA_CKey
    mods = Tuple{Int,Int}[]
    for m in eachmatch(r"\((\d+),\w,(\w+:\d+)\)", structural_mods)
        pos = parse(Int, m.captures[1])
        mod_str = m.captures[2]
        colon_idx = findlast(':', mod_str)
        unimod_id = parse(Int, mod_str[colon_idx+1:end])
        push!(mods, (pos, unimod_id))
    end
    sort!(mods)
    return (String(seq), mods, Int(charge))
end

function _dra_diann_run_to_pioneer(run::AbstractString)::Union{String, Nothing}
    m = match(r"(E\d+H\d+Y\d+)_\d+SPD_DIA_(\d+)", run)
    m === nothing && return nothing
    return "$(m.captures[1])_$(m.captures[2])"
end

"""
    compute_tdc_qvalues(probs, targets) -> Vector{Float64}

Target-decoy competition q-values: sort by score descending, compute
running FDR = #decoys / #targets, then monotonize from worst to best.
"""
function compute_tdc_qvalues(probs::AbstractVector{<:Real}, targets::AbstractVector{Bool})
    n = length(probs)
    order = sortperm(probs, rev=true)
    qvals = Vector{Float64}(undef, n)
    n_target = 0; n_decoy = 0
    for i in order
        if targets[i]
            n_target += 1
        else
            n_decoy += 1
        end
        qvals[i] = n_target > 0 ? n_decoy / n_target : 1.0
    end
    # Monotonize from worst to best
    min_q = 1.0
    for idx in reverse(order)
        qvals[idx] = min(qvals[idx], min_q)
        min_q = qvals[idx]
    end
    return qvals
end

struct RecoveryResult
    file::String
    n_diann::Int
    # Per FDR threshold: (total_passing, diann_recovered, diann_missing, recovery_pct)
    thresholds::Vector{Float64}
    total_passing::Vector{Int}
    diann_recovered::Vector{Int}
    diann_missing::Vector{Int}
    recovery_pct::Vector{Float64}
end

"""
    diann_recovery_analysis(config_path::String; kwargs...) -> Vector{RecoveryResult}

Run recovery analysis using a Pioneer config JSON. Reads `paths.results` from the
config and writes `diann_recovery_report.txt` to that directory.
"""
function diann_recovery_analysis(config_path::String; kwargs...)
    config = JSON.parsefile(config_path)
    results_dir = config["paths"]["results"]
    return diann_recovery_analysis(; pioneer_results_dir=results_dir, kwargs...)
end

"""
    diann_recovery_analysis(; kwargs...) -> Vector{RecoveryResult}

Compare Pioneer second-pass PSMs against DIA-NN identifications.

Keyword arguments:
- `pioneer_results_dir`: Path to Pioneer results (must contain temp_data/second_pass_psms/)
- `library_path`: Path to .poin spectral library
- `diann_report`: Path to DIA-NN report parquet
- `pioneer_files`: Vector of file stems to analyze
- `fdr_thresholds`: Vector of FDR thresholds to evaluate
- `quiet`: Suppress per-file output (default: false)
"""
function diann_recovery_analysis(;
    pioneer_results_dir::String = "/Users/nathanwamsley/Data/For_Figures/OlsenEclipse/OlsenEclipse_nofilters_noMBR_20a3d62c",
    library_path::String = DIANN_LIBRARY_PATH,
    diann_report::String = DIANN_REPORT_PATH,
    pioneer_files::Vector{String} = DIANN_PIONEER_FILES,
    fdr_thresholds::Vector{Float64} = DIANN_FDR_THRESHOLDS,
    quiet::Bool = false,
)
    second_pass_dir = joinpath(pioneer_results_dir, "temp_data", "second_pass_psms")
    if !isdir(second_pass_dir)
        error("Second-pass PSM directory not found: $second_pass_dir")
    end

    # ── 1. Load library ──
    !quiet && println("Loading library...")
    lib_df = DataFrame(Arrow.Table(joinpath(library_path, "precursors_table.arrow")))
    pidx_to_key = Dict{UInt32, _DRA_CKey}()
    for i in 1:nrow(lib_df)
        lib_df.is_decoy[i] && continue
        smods = ismissing(lib_df.structural_mods[i]) ? "" : lib_df.structural_mods[i]
        pidx_to_key[UInt32(i)] = _dra_pioneer_key(lib_df.sequence[i], smods, lib_df.prec_charge[i])
    end
    !quiet && println("  $(length(pidx_to_key)) target precursors in library")

    # ── 2. Load DIA-NN report ──
    !quiet && println("Loading DIA-NN report...")
    diann_df = DataFrame(Parquet2.Dataset(diann_report); copycols=false)
    pioneer_file_set = Set(pioneer_files)

    diann_per_file = Dict{String, Set{_DRA_CKey}}()
    for f in pioneer_files
        diann_per_file[f] = Set{_DRA_CKey}()
    end

    for i in 1:nrow(diann_df)
        pf = _dra_diann_run_to_pioneer(diann_df[i, "Run"])
        pf === nothing && continue
        pf in pioneer_file_set || continue
        qty = diann_df[i, "Precursor.Quantity"]
        (ismissing(qty) || !(qty isa Number) || qty <= 0) && continue
        key = _dra_parse_diann_seq(diann_df[i, "Modified.Sequence"],
                                    Int(diann_df[i, "Precursor.Charge"]))
        push!(diann_per_file[pf], key)
    end

    if !quiet
        for f in pioneer_files
            println("  DIA-NN $f: $(length(diann_per_file[f])) unique precursors")
        end
        diann_union = union(values(diann_per_file)...)
        println("  DIA-NN union: $(length(diann_union)) unique precursors")
    end

    # ── 3. Analyze each Pioneer file ──
    results = RecoveryResult[]

    if !quiet
        println("\n" * "="^90)
        println("FDR RECOVERY ANALYSIS: $(basename(pioneer_results_dir))")
        println("="^90)
    end

    for f in pioneer_files
        diann_keys = diann_per_file[f]
        n_diann = length(diann_keys)

        psms = DataFrame(Arrow.Table(joinpath(second_pass_dir, "$f.arrow")))

        # Best PSM per precursor_idx (highest prec_prob)
        best = combine(groupby(psms, :precursor_idx)) do sub
            idx = argmax(sub.prec_prob)
            return (prec_prob=sub.prec_prob[idx], target=sub.target[idx])
        end

        qvals = compute_tdc_qvalues(best.prec_prob, best.target)

        # Map precursor_idx -> canonical key
        best_keys = Vector{Union{_DRA_CKey, Nothing}}(nothing, nrow(best))
        for i in 1:nrow(best)
            best_keys[i] = get(pidx_to_key, UInt32(best.precursor_idx[i]), nothing)
        end

        total_vec = Int[]
        recov_vec = Int[]
        miss_vec  = Int[]
        pct_vec   = Float64[]

        for threshold in fdr_thresholds
            passing_keys = Set{_DRA_CKey}()
            n_passing = 0
            for i in 1:nrow(best)
                best.target[i] || continue
                qvals[i] <= threshold || continue
                n_passing += 1
                k = best_keys[i]
                k !== nothing && push!(passing_keys, k)
            end

            recovered = length(intersect(passing_keys, diann_keys))
            missing_count = n_diann - recovered
            recov_pct = 100.0 * recovered / n_diann

            push!(total_vec, n_passing)
            push!(recov_vec, recovered)
            push!(miss_vec, missing_count)
            push!(pct_vec, recov_pct)
        end

        push!(results, RecoveryResult(f, n_diann, copy(fdr_thresholds),
                                       total_vec, recov_vec, miss_vec, pct_vec))
    end

    !quiet && print_recovery_table(results)

    # Write report to results directory
    report_path = joinpath(pioneer_results_dir, "diann_recovery_report.txt")
    open(report_path, "w") do io
        write_recovery_report(io, results, pioneer_results_dir)
    end
    !quiet && println("Report written to: $report_path")

    return results
end

"""
    print_recovery_table(results::Vector{RecoveryResult})

Pretty-print the recovery analysis results.
"""
function print_recovery_table(results::Vector{RecoveryResult})
    println("\n" * "="^90)
    println("FDR RECOVERY ANALYSIS")
    println("="^90)

    for r in results
        println("\n── $(r.file) ──")
        println("  DIA-NN unique precursors: $(r.n_diann)")
        @printf("  %-8s  %8s  %8s  %8s  %8s\n", "FDR", "Total", "DIA-NN\u2191", "Missing", "Recov%")
        println("  " * "-"^50)
        for i in eachindex(r.thresholds)
            @printf("  %-8s  %8d  %8d  %8d  %7.1f%%\n",
                    "$(Int(r.thresholds[i]*100))%",
                    r.total_passing[i], r.diann_recovered[i],
                    r.diann_missing[i], r.recovery_pct[i])
        end
    end

    # Summary row: averages across files
    n_thresh = length(results[1].thresholds)
    println("\n── AVERAGE across $(length(results)) files ──")
    avg_diann = mean(r.n_diann for r in results)
    @printf("  Avg DIA-NN unique precursors: %.0f\n", avg_diann)
    @printf("  %-8s  %8s  %8s  %8s  %8s\n", "FDR", "Total", "DIA-NN\u2191", "Missing", "Recov%")
    println("  " * "-"^50)
    for i in 1:n_thresh
        avg_total = mean(r.total_passing[i] for r in results)
        avg_recov = mean(r.diann_recovered[i] for r in results)
        avg_miss  = mean(r.diann_missing[i] for r in results)
        avg_pct   = mean(r.recovery_pct[i] for r in results)
        @printf("  %-8s  %8.0f  %8.0f  %8.0f  %7.1f%%\n",
                "$(Int(results[1].thresholds[i]*100))%",
                avg_total, avg_recov, avg_miss, avg_pct)
    end
    println("\n" * "="^90)
end

"""
    write_recovery_report(io::IO, results, results_dir)

Write the recovery table to an IO stream (file or stdout).
"""
function write_recovery_report(io::IO, results::Vector{RecoveryResult}, results_dir::String)
    println(io, "DIA-NN Recovery Analysis")
    println(io, "Run: $(basename(results_dir))")
    println(io, "Date: $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM"))")
    println(io, "="^90)

    for r in results
        println(io, "\n-- $(r.file) --")
        println(io, "  DIA-NN unique precursors: $(r.n_diann)")
        @printf(io, "  %-8s  %8s  %8s  %8s  %8s\n", "FDR", "Total", "DIA-NN+", "Missing", "Recov%")
        println(io, "  " * "-"^50)
        for i in eachindex(r.thresholds)
            @printf(io, "  %-8s  %8d  %8d  %8d  %7.1f%%\n",
                    "$(Int(r.thresholds[i]*100))%",
                    r.total_passing[i], r.diann_recovered[i],
                    r.diann_missing[i], r.recovery_pct[i])
        end
    end

    n_thresh = length(results[1].thresholds)
    println(io, "\n-- AVERAGE across $(length(results)) files --")
    avg_diann = mean(r.n_diann for r in results)
    @printf(io, "  Avg DIA-NN unique precursors: %.0f\n", avg_diann)
    @printf(io, "  %-8s  %8s  %8s  %8s  %8s\n", "FDR", "Total", "DIA-NN+", "Missing", "Recov%")
    println(io, "  " * "-"^50)
    for i in 1:n_thresh
        avg_total = mean(r.total_passing[i] for r in results)
        avg_recov = mean(r.diann_recovered[i] for r in results)
        avg_miss  = mean(r.diann_missing[i] for r in results)
        avg_pct   = mean(r.recovery_pct[i] for r in results)
        @printf(io, "  %-8s  %8.0f  %8.0f  %8.0f  %7.1f%%\n",
                "$(Int(results[1].thresholds[i]*100))%",
                avg_total, avg_recov, avg_miss, avg_pct)
    end
    println(io, "\n" * "="^90)
end

# Import mean and Dates
using Statistics, Dates
