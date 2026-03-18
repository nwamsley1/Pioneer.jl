#=
Precursor-in-MS2 Feature Analysis
==================================
Post-hoc computation and analysis of the `log2_prec_ms2_intensity` feature
for ALL second-pass PSMs (not just recovered ones). Compares distributions
across score tiers and cross-references DIA-NN identifications.

Usage:
    include("scripts/analyze_prec_ms2_feature.jl")
    analyze_prec_ms2_feature()   # uses defaults from OlsenEclipse
=#

using Arrow, DataFrames, Tables, Statistics, Printf, Parquet2

# ── Default paths ──
const _PM_LIBRARY_PATH = "/Users/nathanwamsley/Data/SPEC_LIBS/3P_rank_based.poin"
const _PM_MS_DATA_DIR  = "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse"
const _PM_DIANN_REPORT = "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/DIA-NN_results/OlsenExplorisThreeProteome500ng-11-24-2025-report.parquet"

const _PM_FILE_MAP = Dict(
    "E45H50Y5_2" => "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E45H50Y5_30SPD_DIA_2.arrow",
    "E45H50Y5_3" => "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E45H50Y5_30SPD_DIA_3.arrow",
    "E45H50Y5_4" => "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E45H50Y5_30SPD_DIA_4.arrow",
    "E5H50Y45_1" => "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.arrow",
    "E5H50Y45_2" => "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_2.arrow",
    "E5H50Y45_4" => "20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_4.arrow",
)

# ── Precursor key type (matches diann_recovery_analysis.jl) ──
const _PM_CKey = Tuple{String, Vector{Tuple{Int,Int}}, Int}

function _pm_parse_diann_seq(mod_seq::AbstractString, charge::Int)::_PM_CKey
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

function _pm_diann_run_to_pioneer(run::AbstractString)::Union{String, Nothing}
    m = match(r"(E\d+H\d+Y\d+)_\d+SPD_DIA_(\d+)", run)
    m === nothing && return nothing
    return "$(m.captures[1])_$(m.captures[2])"
end

"""
    estimate_mass_correction(mz_arrays, int_arrays, psms, lib_mzs; top_n=5000, tol_ppm=20.0)

Estimate the systematic mass offset (ppm) from the highest-scoring target PSMs.
Finds the closest peak to each precursor m/z and computes the median ppm error.
"""
function estimate_mass_correction(mz_arrays, int_arrays, psms::DataFrame, lib_mzs::Vector{Float32};
                                   top_n::Int=5000, tol_ppm::Float64=20.0)
    targets = filter(r -> r.target, psms)
    top = first(sort(targets, :prec_prob, rev=true), min(top_n, nrow(targets)))

    ppm_errors = Float64[]
    for row in eachrow(top)
        prec_mz = lib_mzs[row.precursor_idx]
        mz_arr = mz_arrays[row.scan_idx]
        tol = prec_mz * tol_ppm / 1f6
        lo = searchsortedfirst(mz_arr, prec_mz - tol)
        hi = searchsortedlast(mz_arr, prec_mz + tol)
        lo > hi && continue
        lo > length(mz_arr) && continue

        best_j = lo
        best_dist = abs(coalesce(mz_arr[lo], Inf32) - prec_mz)
        for j in (lo+1):hi
            d = abs(coalesce(mz_arr[j], Inf32) - prec_mz)
            if d < best_dist
                best_dist = d
                best_j = j
            end
        end
        obs_mz = Float64(coalesce(mz_arr[best_j], NaN32))
        isnan(obs_mz) && continue
        obs_int = Float64(coalesce(int_arrays[row.scan_idx][best_j], 0f0))
        obs_int <= 0 && continue
        push!(ppm_errors, (obs_mz - prec_mz) / prec_mz * 1e6)
    end

    if isempty(ppm_errors)
        @warn "Could not estimate mass correction, defaulting to 0 ppm"
        return 0.0
    end
    med = median(ppm_errors)
    @info @sprintf("  Mass correction estimate: %.2f ppm (from %d peaks, IQR: %.2f to %.2f)",
                   med, length(ppm_errors),
                   quantile(ppm_errors, 0.25), quantile(ppm_errors, 0.75))
    return med
end

"""
    compute_prec_ms2_intensity(mz_arrays, int_arrays, psms, lib_mzs, mass_offset_ppm, tol_ppm)

Compute log2_prec_ms2_intensity for every PSM, applying mass correction
to raw peaks before matching (matching the corrected pipeline logic).
"""
function compute_prec_ms2_intensity(mz_arrays, int_arrays, psms::DataFrame,
                                     lib_mzs::Vector{Float32},
                                     mass_offset_ppm::Float64, tol_ppm::Float64)
    N = nrow(psms)
    log2_int = zeros(Float32, N)

    for i in 1:N
        prec_mz = lib_mzs[psms.precursor_idx[i]]
        mz_arr = mz_arrays[psms.scan_idx[i]]
        int_arr = int_arrays[psms.scan_idx[i]]

        # Corrected bounds (what the corrected peak should fall within)
        ppm = prec_mz / 1f6
        low = prec_mz - tol_ppm * ppm
        high = prec_mz + tol_ppm * ppm

        # Widen binary search to cover raw (uncorrected) m/z values
        offset = abs(mass_offset_ppm) * ppm
        lo_idx = searchsortedfirst(mz_arr, low - offset)
        hi_idx = searchsortedlast(mz_arr, high + offset)

        lo_idx > hi_idx && continue
        lo_idx > length(mz_arr) && continue

        best_int = Float32(0)
        for j in lo_idx:hi_idx
            raw_mz = Float64(coalesce(mz_arr[j], NaN32))
            isnan(raw_mz) && continue
            # Apply mass correction: corrected = raw - offset*ppm
            corrected_mz = raw_mz - mass_offset_ppm * (raw_mz / 1e6)
            if corrected_mz >= low && corrected_mz <= high
                peak_int = Float32(coalesce(int_arr[j], 0f0))
                best_int = max(best_int, peak_int)
            end
        end
        if best_int > 0
            log2_int[i] = log2(best_int)
        end
    end
    return log2_int
end

"""
    compute_tdc_qvalues(probs, targets)

Target-decoy competition q-values.
"""
function compute_tdc_qvalues(probs::AbstractVector{<:Real}, targets::AbstractVector{Bool})
    n = length(probs)
    order = sortperm(probs, rev=true)
    qvals = Vector{Float64}(undef, n)
    n_target = 0; n_decoy = 0
    for i in order
        if targets[i]; n_target += 1 else n_decoy += 1 end
        qvals[i] = n_target > 0 ? n_decoy / n_target : 1.0
    end
    min_q = 1.0
    for idx in reverse(order)
        qvals[idx] = min(qvals[idx], min_q)
        min_q = qvals[idx]
    end
    return qvals
end

"""
    summarize_group(name, vals)

Print summary statistics for a named group of log2_prec_ms2_intensity values.
"""
function summarize_group(name::String, vals::Vector{Float32})
    n = length(vals)
    n_found = count(>(0), vals)
    pct = 100.0 * n_found / max(n, 1)
    found_vals = vals[vals .> 0]
    if isempty(found_vals)
        @printf("  %-35s  n=%6d  found=%5d (%5.1f%%)\n", name, n, n_found, pct)
    else
        @printf("  %-35s  n=%6d  found=%5d (%5.1f%%)  med=%.1f  mean=%.1f  [%.1f, %.1f]\n",
                name, n, n_found, pct,
                median(found_vals), mean(found_vals),
                minimum(found_vals), maximum(found_vals))
    end
end

"""
    analyze_prec_ms2_feature(; kwargs...)

Compute and analyze log2_prec_ms2_intensity for all second-pass PSMs.

Groups analyzed:
- Targets at q<1%  (high confidence)
- Targets at 1%<q<5%
- Targets at 5%<q<20%
- Targets at q>20% (low confidence)
- Decoys
- DIA-NN-only IDs (found by DIA-NN but missed by Pioneer at 1% FDR)
- DIA-NN-recovered IDs (found by both)
"""
function analyze_prec_ms2_feature(;
    results_dir::String = "/Users/nathanwamsley/Data/For_Figures/OlsenEclipse/OlsenEclipse_FixSinglePointChrom",
    library_path::String = _PM_LIBRARY_PATH,
    ms_data_dir::String  = _PM_MS_DATA_DIR,
    diann_report::String = _PM_DIANN_REPORT,
    tol_ppm::Float64 = 10.0,
    pioneer_files::Vector{String} = ["E45H50Y5_2", "E45H50Y5_3", "E45H50Y5_4",
                                      "E5H50Y45_1", "E5H50Y45_2", "E5H50Y45_4"],
)
    second_pass_dir = joinpath(results_dir, "temp_data", "second_pass_psms")
    @assert isdir(second_pass_dir) "Not found: $second_pass_dir"

    # ── 1. Load library precursor m/z ──
    println("Loading library precursors...")
    lib_df = DataFrame(Arrow.Table(joinpath(library_path, "precursors_table.arrow")))
    lib_mzs = Vector{Float32}(lib_df.mz)

    # Build precursor_idx -> canonical key mapping (targets only)
    pidx_to_key = Dict{UInt32, _PM_CKey}()
    for i in 1:nrow(lib_df)
        lib_df.is_decoy[i] && continue
        smods = ismissing(lib_df.structural_mods[i]) ? "" : lib_df.structural_mods[i]
        mods = Tuple{Int,Int}[]
        for m in eachmatch(r"\((\d+),\w,(\w+:\d+)\)", smods)
            pos = parse(Int, m.captures[1])
            mod_str = m.captures[2]
            colon_idx = findlast(':', mod_str)
            unimod_id = parse(Int, mod_str[colon_idx+1:end])
            push!(mods, (pos, unimod_id))
        end
        sort!(mods)
        pidx_to_key[UInt32(i)] = (String(lib_df.sequence[i]), mods, Int(lib_df.prec_charge[i]))
    end
    println("  $(length(pidx_to_key)) target precursors in library")

    # ── 2. Load DIA-NN report ──
    println("Loading DIA-NN report...")
    diann_df = DataFrame(Parquet2.Dataset(diann_report); copycols=false)
    diann_per_file = Dict{String, Set{_PM_CKey}}()
    for f in pioneer_files
        diann_per_file[f] = Set{_PM_CKey}()
    end
    pioneer_set = Set(pioneer_files)
    for i in 1:nrow(diann_df)
        pf = _pm_diann_run_to_pioneer(diann_df[i, "Run"])
        pf === nothing && continue
        pf in pioneer_set || continue
        qty = diann_df[i, "Precursor.Quantity"]
        (ismissing(qty) || !(qty isa Number) || qty <= 0) && continue
        key = _pm_parse_diann_seq(diann_df[i, "Modified.Sequence"],
                                   Int(diann_df[i, "Precursor.Charge"]))
        push!(diann_per_file[pf], key)
    end
    for f in pioneer_files
        println("  DIA-NN $f: $(length(diann_per_file[f])) precursors")
    end

    # ── 3. Process each file ──
    # Collect pooled data for aggregate analysis
    all_results = DataFrame[]

    for f in pioneer_files
        println("\n" * "="^80)
        println("Processing $f")
        println("="^80)

        # Load second-pass PSMs
        psms = DataFrame(Arrow.Table(joinpath(second_pass_dir, "$f.arrow")))
        println("  $(nrow(psms)) total PSMs ($(count(psms.target)) targets, $(count(.!psms.target)) decoys)")

        # Load raw MS data
        ms_file = joinpath(ms_data_dir, _PM_FILE_MAP[f])
        println("  Loading MS data: $(basename(ms_file))")
        ms_table = Arrow.Table(ms_file)
        mz_arrays = ms_table[:mz_array]
        int_arrays = ms_table[:intensity_array]

        # Estimate mass correction from top targets
        mass_offset = estimate_mass_correction(mz_arrays, int_arrays, psms, lib_mzs)

        # Compute feature for ALL PSMs
        println("  Computing log2_prec_ms2_intensity...")
        log2_int = compute_prec_ms2_intensity(mz_arrays, int_arrays, psms, lib_mzs, mass_offset, tol_ppm)
        psms[!, :log2_prec_ms2_intensity] = log2_int

        n_found = count(>(0), log2_int)
        println("  Overall: $(n_found)/$(nrow(psms)) ($(round(100*n_found/nrow(psms), digits=1))%) precursor peaks found")

        # Compute q-values for best PSM per precursor
        best = combine(groupby(psms, :precursor_idx)) do sub
            idx = argmax(sub.prec_prob)
            return (prec_prob = sub.prec_prob[idx],
                    target = sub.target[idx],
                    log2_prec_ms2_intensity = sub.log2_prec_ms2_intensity[idx],
                    scan_idx = sub.scan_idx[idx])
        end
        qvals = compute_tdc_qvalues(best.prec_prob, best.target)
        best[!, :qval] = qvals

        # Map to canonical keys for DIA-NN comparison
        best_keys = Vector{Union{_PM_CKey, Nothing}}(nothing, nrow(best))
        for i in 1:nrow(best)
            best_keys[i] = get(pidx_to_key, UInt32(best.precursor_idx[i]), nothing)
        end

        # Pioneer passing at 1% FDR (target keys)
        pioneer_1pct = Set{_PM_CKey}()
        for i in 1:nrow(best)
            best.target[i] && qvals[i] <= 0.01 && best_keys[i] !== nothing && push!(pioneer_1pct, best_keys[i])
        end
        diann_keys = diann_per_file[f]
        diann_only = setdiff(diann_keys, pioneer_1pct)
        diann_recovered = intersect(diann_keys, pioneer_1pct)

        println("  Pioneer 1% FDR: $(length(pioneer_1pct)) | DIA-NN: $(length(diann_keys))")
        println("  DIA-NN recovered: $(length(diann_recovered)) | DIA-NN only: $(length(diann_only))")

        # ── Group PSMs by category ──
        # Use best PSM per precursor for the tier analysis
        targets_q01  = best[best.target .& (qvals .<= 0.01), :log2_prec_ms2_intensity]
        targets_q05  = best[best.target .& (qvals .> 0.01) .& (qvals .<= 0.05), :log2_prec_ms2_intensity]
        targets_q20  = best[best.target .& (qvals .> 0.05) .& (qvals .<= 0.20), :log2_prec_ms2_intensity]
        targets_rest = best[best.target .& (qvals .> 0.20), :log2_prec_ms2_intensity]
        decoys_all   = best[.!best.target, :log2_prec_ms2_intensity]

        # DIA-NN subgroups: look up feature values for DIA-NN-only and DIA-NN-recovered precursors
        diann_only_vals = Float32[]
        diann_recov_vals = Float32[]
        for i in 1:nrow(best)
            best.target[i] || continue
            k = best_keys[i]
            k === nothing && continue
            if k in diann_only
                push!(diann_only_vals, best.log2_prec_ms2_intensity[i])
            elseif k in diann_recovered
                push!(diann_recov_vals, best.log2_prec_ms2_intensity[i])
            end
        end

        println("\n  log2_prec_ms2_intensity by group (best PSM per precursor):")
        summarize_group("Targets q<1%", Vector{Float32}(targets_q01))
        summarize_group("Targets 1%<q<5%", Vector{Float32}(targets_q05))
        summarize_group("Targets 5%<q<20%", Vector{Float32}(targets_q20))
        summarize_group("Targets q>20%", Vector{Float32}(targets_rest))
        summarize_group("Decoys", Vector{Float32}(decoys_all))
        summarize_group("DIA-NN recovered (by Pioneer)", diann_recov_vals)
        summarize_group("DIA-NN only (missed by Pioneer)", diann_only_vals)

        # Also analyze ALL PSMs (not deduplicated) for scan-level view
        println("\n  Scan-level (all PSMs, not deduplicated):")
        all_targets = psms[psms.target, :log2_prec_ms2_intensity]
        all_decoys  = psms[.!psms.target, :log2_prec_ms2_intensity]
        summarize_group("All target PSMs", Vector{Float32}(all_targets))
        summarize_group("All decoy PSMs", Vector{Float32}(all_decoys))

        # Store for pooled analysis
        best[!, :file] .= f
        best[!, :is_diann] .= false
        best[!, :is_diann_only] .= false
        for i in 1:nrow(best)
            k = best_keys[i]
            k === nothing && continue
            if k in diann_keys
                best[i, :is_diann] = true
            end
            if k in diann_only
                best[i, :is_diann_only] = true
            end
        end
        push!(all_results, best)
    end

    # ── 4. Pooled summary ──
    pooled = vcat(all_results...)
    pooled_qvals = compute_tdc_qvalues(pooled.prec_prob, pooled.target)

    println("\n" * "="^80)
    println("POOLED ANALYSIS ($(length(pioneer_files)) files, $(nrow(pooled)) unique precursor-file pairs)")
    println("="^80)

    p_q01  = pooled[pooled.target .& (pooled_qvals .<= 0.01), :log2_prec_ms2_intensity]
    p_q05  = pooled[pooled.target .& (pooled_qvals .> 0.01) .& (pooled_qvals .<= 0.05), :log2_prec_ms2_intensity]
    p_q20  = pooled[pooled.target .& (pooled_qvals .> 0.05) .& (pooled_qvals .<= 0.20), :log2_prec_ms2_intensity]
    p_rest = pooled[pooled.target .& (pooled_qvals .> 0.20), :log2_prec_ms2_intensity]
    p_dec  = pooled[.!pooled.target, :log2_prec_ms2_intensity]
    p_diann_recov = pooled[pooled.target .& pooled.is_diann .& .!pooled.is_diann_only, :log2_prec_ms2_intensity]
    p_diann_only  = pooled[pooled.target .& pooled.is_diann_only, :log2_prec_ms2_intensity]

    println("\n  log2_prec_ms2_intensity by group:")
    summarize_group("Targets q<1%", Vector{Float32}(p_q01))
    summarize_group("Targets 1%<q<5%", Vector{Float32}(p_q05))
    summarize_group("Targets 5%<q<20%", Vector{Float32}(p_q20))
    summarize_group("Targets q>20%", Vector{Float32}(p_rest))
    summarize_group("Decoys", Vector{Float32}(p_dec))
    summarize_group("DIA-NN recovered (by Pioneer)", Vector{Float32}(p_diann_recov))
    summarize_group("DIA-NN only (missed by Pioneer)", Vector{Float32}(p_diann_only))

    # ── 5. Discriminative power summary ──
    println("\n  Discriminative power (found rate):")
    for (name, vals) in [
        ("Targets q<1%", p_q01), ("Targets q>20%", p_rest), ("Decoys", p_dec),
        ("DIA-NN only", p_diann_only),
    ]
        n = length(vals)
        n_found = count(>(0), vals)
        @printf("    %-35s  %5.1f%% found\n", name, 100.0 * n_found / max(n, 1))
    end

    return pooled
end
