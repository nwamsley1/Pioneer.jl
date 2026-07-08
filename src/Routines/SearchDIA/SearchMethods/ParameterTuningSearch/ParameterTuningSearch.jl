# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
    ParameterTuningSearch

Search method for tuning mass error and retention time parameters.

This search:
1. Estimates mass error distributions for fragment matching
2. Fits retention time alignment between library and empirical data
3. Generates QC plots for parameter estimation quality
4. Stores models in SearchContext for use by other search methods

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :isotope_err_bounds => (3, 1),
    :presearch_params => Dict(
        "frag_tol_ppm" => 30.0,
        "frag_err_quantile" => 0.01,
        "min_samples" => 1000,
        "max_qval" => 0.01,
        "max_presearch_iters" => 10,
        "min_index_search_score" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "n_frag_isotopes" => 2,
        "max_frag_rank" => 10,
        "abreviate_precursor_calc" => false
    ),
    :irt_mapping_params => Dict(
        "n_sigma_tol" => 3
    )
)

# Execute search
results = execute_search(ParameterTuningSearch(), search_context, params)
```
"""
# Type definitions moved to types.jl


#==========================================================
State Management Functions
==========================================================#

#==========================================================
Results Access Methods
==========================================================#
getMassErrorModel(ptsr::ParameterTuningSearchResults) = ptsr.mass_err_model[]
getRtToIrtModel(ptsr::ParameterTuningSearchResults) = ptsr.rt_to_irt_model[]
getQcPlotsFolder(ptsr::ParameterTuningSearchResults) = ptsr.qc_plots_folder_path
getDiagnostics(ptsr::ParameterTuningSearchResults) = ptsr.diagnostics
getParameterHistory(ptsr::ParameterTuningSearchResults) = ptsr.parameter_history

function set_rt_to_irt_model!(
    ptsr::ParameterTuningSearchResults,
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64,
    model::Tuple{RtConversionModel, Vector{Float32}, Vector{Float32}, Float32}
) where {P<:ParameterTuningSearchParameters}
    
    ptsr.rt_to_irt_model[] = model[1]
    resize!(ptsr.irt, 0)
    resize!(ptsr.rt, 0)
    append!(ptsr.rt, model[2])
    append!(ptsr.irt, model[3])
    
    # CRITICAL: Store RT model in SearchContext for downstream methods
    setRtIrtMap!(search_context, model[1], ms_file_idx)
    
    #parsed_fname = getParsedFileName(search_context, ms_file_idx)
    getIrtErrors(search_context)[ms_file_idx] = model[4] * TUNING_IRT_TOL_SIGMA
end


#==========================================================
Interface Implementation
==========================================================#

get_parameters(::ParameterTuningSearch, params::Any) = ParameterTuningSearchParameters(params)

function init_search_results(::ParameterTuningSearchParameters, search_context::SearchContext)
    out_dir = getDataOutDir(search_context)
    qc_dir = joinpath(out_dir, "qc_plots")
    !isdir(qc_dir) && mkdir(qc_dir)
    rt_alingment_plots = joinpath(qc_dir, "rt_alignment_plots")
    !isdir(rt_alingment_plots) && mkdir(rt_alingment_plots)
    mass_error_plots = joinpath(qc_dir, "mass_error_plots")
    !isdir(mass_error_plots) && mkdir(mass_error_plots)
    ms1_mass_error_plots = joinpath(qc_dir, "ms1_mass_error_plots")
    !isdir(ms1_mass_error_plots ) && mkdir(ms1_mass_error_plots )

    # ParameterTuning is the first search method to run — verify the library
    # is m/z-sorted per precursor before any fused kernel calls. MainSearch's
    # verify_mz_sorted is then a safe no-op repeat (also on shared library).
    t = time()
    lookup = getFragmentLookupTable(getSpecLib(search_context))
    verify_mz_sorted(getFragments(lookup), lookup.prec_frag_ranges)
    @debug_l1 "  ParameterTuning: verify_mz_sorted OK " *
               "($(round(time() - t, digits=2))s)"
    return ParameterTuningSearchResults(
        Ref{AbstractMassErrorModel}(MassErrorModel(0.0f0, (30.0f0, 30.0f0))),
        Ref{RtConversionModel}(),
        Vector{Float32}(),
        Vector{Float32}(),
        Vector{Float32}(),
        Vector{Float32}(),  # frag_mzs
        Vector{UInt8}[],  # rt_plots (legacy)
        Vector{UInt8}[],  # mass_plots (legacy)
        Any[],            # rt_plot_objects
        Any[],            # mass_plot_objects
        Any[],            # nce_plot_objects
        qc_dir,
        ParameterTuningDiagnostics(),
        ParameterHistory(),
        Ref{Union{Nothing, IterationState}}(nothing)  # current_iteration_state
    )
end

#==========================================================
Helper Functions for Refactored process_file!
==========================================================#

# Collect raw PSMs (library search + columns) without scoring or FDR filtering.
# Returns unscored DataFrame ready for probit scoring.
function collect_raw_psms(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64;
    scan_indices::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    max_peaks::Int = 0
) where {P<:ParameterTuningSearchParameters}
    psms = library_search(spectra, search_context, params, ms_file_idx;
                          scan_indices = scan_indices, max_peaks = max_peaks)
    if !iszero(size(psms, 1))
        precursors = getPrecursors(getSpecLib(search_context))
        add_tuning_search_columns!(
            psms, spectra,
            getIsDecoy(precursors), getIrt(precursors),
            getCharge(precursors), getRetentionTimes(spectra),
            getTICs(spectra))
    end
    return psms
end

"""
    initialize_models!(search_context, ms_file_idx, params)

Initialize mass error and quad transmission models for file.
"""
function initialize_models!(search_context, ms_file_idx, params)
    # Seed the per-file mass-error and quad-transmission models. The mass
    # error tolerance here is a placeholder — the wide-scout phase
    # (~WIDE_SCOUT_TOL_PPM) immediately overwrites it before any search
    # is run on this file.
    setMassErrorModel!(search_context, ms_file_idx, MassErrorModel(
        0.0f0,
        (WIDE_SCOUT_TOL_PPM, WIDE_SCOUT_TOL_PPM)
    ))

    # Quad transmission: trust the stated isolation width with a hard
    # square cutoff during tuning, before any file-specific Razo fit.
    setQuadTransmissionModel!(search_context, ms_file_idx, SquareQuadModel(0.0f0))
end


"""
    store_final_results!(results, search_context, params, ms_file_idx, converged, n_attempts, final_psm_count, iteration_state)

Store tuning results and diagnostic information, including best attempt fallback status.
"""
function store_final_results!(results, search_context, params, ms_file_idx, 
                              converged, n_attempts, final_psm_count, iteration_state::IterationState)
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    
    # Determine convergence type
    convergence_type = if converged
        "CONVERGED"
    elseif iteration_state.best_psm_count > 0
        "BEST_ATTEMPT ($(iteration_state.best_psm_count) PSMs)"
    else
        "FAILED"
    end
    
    # Record tuning results
    tuning_result = TuningResults(
        getMassOffset(getMassErrorModel(search_context, ms_file_idx)),
        (getLeftTol(getMassErrorModel(search_context, ms_file_idx)), 
         getRightTol(getMassErrorModel(search_context, ms_file_idx))),
        converged,
        final_psm_count,
        n_attempts
    )
    store_tuning_results!(results.parameter_history, ms_file_idx, tuning_result)
    
    # Record diagnostic status
    status = ParameterTuningStatus(
        ms_file_idx,
        parsed_fname,
        converged,
        !converged,  # used_fallback
        !converged ? "Failed to converge after $n_attempts attempts" : "",
        n_attempts,
        final_psm_count,
        getMassOffset(getMassErrorModel(search_context, ms_file_idx)),
        (getLeftTol(getMassErrorModel(search_context, ms_file_idx)), 
         getRightTol(getMassErrorModel(search_context, ms_file_idx)))
    )
    record_tuning_status!(results.diagnostics, status)
    
    # Store final mass error model in results
    results.mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)
end


#==========================================================
Main process_file! Function
==========================================================#

"""
    collect_and_converge!(results, params, search_context, ms_file_idx,
        spectra, scan_priority, initial_scans, max_scans, iteration_state)

Collect PSMs with 3× scan growth. Raw PSMs are accumulated across batches
and jointly re-scored with probit + FDR each iteration (better model with
more data). Returns (converged::Bool, scored_psms::DataFrame, n_scans_used::Int, rate::Float64).
"""
function accumulate_psms!(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::ParameterTuningSearchParameters,
    ms_file_idx::Int64,
    scan_priority::AbstractVector{<:Integer};
    mass_model::AbstractMassErrorModel,
    target_psms::Int64,
    initial_scans::Int64,
    max_peaks::Int = 0,
    score_tiers = TUNING_SCORE_TIERS,
    n_required_top::Int = TUNING_N_REQUIRED_TOP,
    fdr_threshold::Float16 = Float16(0.01),
    label::String = "accumulate"
)
    all_scan_indices = scan_priority[1:min(length(scan_priority), length(scan_priority))]
    max_scans = length(all_scan_indices)
    fdr_scale = getLibraryFdrScaleFactor(search_context)

    scored_psms = DataFrame()
    n_passing = 0
    total_scans_used = 0

    for (tier_idx, score) in enumerate(score_tiers)
        setMassErrorModel!(search_context, ms_file_idx, mass_model)
        setCurrentMinScore!(params, score)
        lut = make_top_n_required_lut(n_required_top, Int(score))
        setBitVecFilter!(search_context, ms_file_idx, lut)

        t_tier = time()
        raw_psms = DataFrame()
        prev = 0
        n_passing = 0
        # First tier: start at initial_scans. Subsequent: start at max_scans.
        scan_target = tier_idx == 1 ? min(initial_scans, max_scans) : max_scans

        while prev < max_scans
            batch_indices = all_scan_indices[(prev+1):scan_target]
            batch = collect_raw_psms(
                spectra, search_context, params, ms_file_idx;
                scan_indices = batch_indices,
                max_peaks = max_peaks)
            if size(batch, 1) > 0
                if isempty(raw_psms)
                    raw_psms = batch
                else
                    append!(raw_psms, batch)
                end
            end
            prev = scan_target

            # Score all accumulated PSMs jointly
            n_raw = size(raw_psms, 1)
            if n_raw >= 50
                scored_tmp = copy(raw_psms)
                score_presearch!(scored_tmp)
                scored_tmp[!, :q_value] = zeros(Float16, nrow(scored_tmp))
                get_qvalues!(scored_tmp[!,:prob], scored_tmp[!,:target],
                             scored_tmp[!,:q_value]; fdr_scale_factor=fdr_scale)
                filter!(row -> row.q_value::Float16 <= fdr_threshold, scored_tmp)
                filter!(row -> row.target::Bool, scored_tmp)
                n_passing = nrow(scored_tmp)
                scored_psms = scored_tmp
            end

            @debug_l1 "  $(label) (score≥$(score)): $(prev) scans, $(n_raw) raw, " *
                       "$(n_passing) at $(round(Float64(fdr_threshold)*100, digits=1))% FDR"

            if n_passing >= target_psms
                break
            end

            # Estimate additional scans from rate
            scans_remaining = max_scans - prev
            scans_remaining <= 0 && break
            rate = n_passing / max(prev, 1)
            additional = if rate > 0
                remaining = target_psms - n_passing
                clamp(ceil(Int, remaining / rate * 1.5), 1, scans_remaining)
            else
                min(prev, scans_remaining)
            end
            scan_target = min(prev + additional, max_scans)
            scan_target <= prev && break
        end

        total_scans_used = prev
        if n_passing >= target_psms
            break
        end
        @debug_l1 "  $(label): score≥$(score) insufficient ($(n_passing) < $(target_psms)), " *
                   "backing off ($(round(time()-t_tier, digits=2))s)"
    end

    # Clean up LUT filter
    delete!(search_context.bitvec_filter, ms_file_idx)

    converged = n_passing >= target_psms
    rate = n_passing / max(total_scans_used, 1)
    return converged, scored_psms, total_scans_used, rate
end

"""
NCE sweep: reuse existing (scan, precursor) pairs to find optimal NCE per precursor.
Skips fragment index — builds PerScanPrecursorIndex directly from PSMs.
Returns NCE model or nothing if insufficient data.
"""
function fit_nce_from_psms!(
    search_context::SearchContext,
    params::ParameterTuningSearchParameters,
    ms_file_idx::Int64,
    spectra::MassSpecData,
    psms::DataFrame;
    nce_grid::AbstractVector{Float32} = LinRange{Float32}(21.0f0, 40.0f0, 20)
)
    nrow(psms) < 50 && return nothing

    spec_lib = getSpecLib(search_context)
    precursors = getPrecursors(spec_lib)
    ion_list = getFragmentLookupTable(spec_lib)
    search_data = getSearchData(search_context)
    qtm = getQuadTransmissionModel(search_context, ms_file_idx)
    mem = getMassErrorModel(search_context, ms_file_idx)
    rt_to_irt = getRtIrtModel(search_context, ms_file_idx)
    irt_tol = get_irt_tolerance(search_context, params, ms_file_idx)

    # Build PerScanPrecursorIndex from existing PSMs
    scan_idxs = psms[!, :scan_idx]
    prec_idxs = psms[!, :precursor_idx]
    sorted = sortperm(scan_idxs)
    scan_idxs_sorted = scan_idxs[sorted]
    prec_idxs_sorted = UInt32.(prec_idxs[sorted])

    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(missing, length(spectra))
    i = 1
    while i <= length(scan_idxs_sorted)
        si = scan_idxs_sorted[i]
        j = i
        while j <= length(scan_idxs_sorted) && scan_idxs_sorted[j] == si
            j += 1
        end
        scan_to_prec_idx[si] = i:(j-1)
        i = j
    end
    prec_index = PerScanPrecursorIndex(scan_to_prec_idx, prec_idxs_sorted)

    # Build thread tasks from unique scans
    unique_scans = unique(Int.(scan_idxs_sorted))
    n_threads = Threads.nthreads()
    thread_tasks = [(i, Int[]) for i in 1:n_threads]
    for (idx, si) in enumerate(unique_scans)
        push!(thread_tasks[mod1(idx, n_threads)][2], si)
    end
    filter!(tt -> !isempty(last(tt)), thread_tasks)

    # Sweep NCE grid (fused — same dispatch as library_search).
    t_nce = time()
    all_results = map(nce_grid) do nce_val
        nce_model = PiecewiseNceModel(nce_val)
        tasks = map(thread_tasks) do thread_task
            Threads.@spawn process_scans_fused!(
                last(thread_task), spectra, prec_index,
                ms_file_idx,
                search_data[first(thread_task)], params, precursors, ion_list,
                nce_model, qtm, mem, rt_to_irt, irt_tol)
        end
        result = vcat(fetch.(tasks)...)
        if !isempty(result)
            result[!, :nce] .= nce_val
        end
        return result
    end
    nce_psms = vcat(all_results...)
    dt_nce = round(time() - t_nce, digits=2)

    if nrow(nce_psms) < 50
        @debug_l1 "NCE sweep: too few PSMs ($(nrow(nce_psms)))"
        return nothing
    end

    # Add precursor metadata needed for NCE fit
    prec_mzs = getMz(precursors)
    prec_charges = getCharge(precursors)
    nce_psms[!, :prec_mz] = Float32[prec_mzs[pid] for pid in nce_psms[!, :precursor_idx]]
    nce_psms[!, :charge] = UInt8[prec_charges[pid] for pid in nce_psms[!, :precursor_idx]]

    # Keep best NCE per precursor (highest gof — deconvolution goodness of fit)
    sort!(nce_psms, :gof, rev=true)
    best_nce = combine(groupby(nce_psms, :precursor_idx), first)

    # Fit NCE model
    nce_model = fit_binned_median_nce(
        best_nce[!, :prec_mz],
        best_nce[!, :nce],
        best_nce[!, :charge],
        Float32(median(nce_grid)))

    setNceModel!(search_context, ms_file_idx, nce_model)

    n_precs = nrow(best_nce)
    charges = sort(unique(best_nce[!, :charge]))
    @debug_l1 "NCE: $(n_precs) precursors, $(length(charges)) charges, $(length(nce_grid)) grid pts ($(dt_nce)s)"

    # Generate per-charge diagnostic plots
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    nce_plots = Plots.Plot[]
    for charge in charges
        mask = best_nce[!, :charge] .== charge
        n_c = count(mask)
        n_c < 10 && continue
        charge_mz = best_nce[mask, :prec_mz]
        charge_nce = best_nce[mask, :nce]
        ci = Int(UInt8(charge))

        # Get the bin edges from the fitted model
        has_bins = ci >= 1 && ci <= 6 && nce_model.offsets[ci] != 0x00
        if has_bins
            nb = Int(nce_model.n_bins[ci])
            bw = Float64(nce_model.bin_width[ci])
            mz_lo = Float64(nce_model.mz_min[ci])
            bin_edges = [mz_lo + (b - 1) * bw for b in 1:nb+1]
            bin_medians = [Float64(nce_model.medians[Int(nce_model.offsets[ci]) + b - 1]) for b in 1:nb]
        else
            nb = 1
            mz_lo_f, mz_hi_f = extrema(charge_mz)
            bin_edges = [Float64(mz_lo_f), Float64(mz_hi_f) + 1.0]
            bin_medians = [Float64(nce_model(median(charge_mz), charge))]
        end

        p = Plots.plot(
            xlabel = "Precursor m/z", ylabel = "Best NCE",
            title = _split_title(parsed_fname, "NCE +$(charge)") *
                    "\nn=$n_c, $(nb) bins, $(length(nce_grid)) grid pts",
            size = (900, 900), topmargin = 15Plots.mm,
            ylim = (minimum(nce_grid) - 1, maximum(nce_grid) + 1))

        for b in 1:nb
            lo = bin_edges[b]
            hi = bin_edges[b + 1]
            in_bin = (b == nb) ? (charge_mz .>= Float32(lo)) : ((charge_mz .>= Float32(lo)) .& (charge_mz .< Float32(hi)))
            bin_mz = charge_mz[in_bin]
            bin_nce_vals = charge_nce[in_bin]
            isempty(bin_mz) && continue
            center = (lo + hi) / 2
            half_w = (hi - lo) / 2

            # Jittered raw points
            jittered_x = [center + half_w * 0.8 * (2 * rand() - 1) for _ in eachindex(bin_mz)]
            Plots.scatter!(p, jittered_x, Float64.(bin_nce_vals),
                alpha = 0.12, markersize = 1.5, color = :steelblue, label = (b == 1 ? "data" : nothing))

            # Boxplot whiskers and box
            nce_sorted = sort(Float64.(bin_nce_vals))
            q1 = quantile(nce_sorted, 0.25)
            q3 = quantile(nce_sorted, 0.75)
            iqr = q3 - q1
            wlo = max(minimum(nce_sorted), q1 - 1.5 * iqr)
            whi = min(maximum(nce_sorted), q3 + 1.5 * iqr)
            box_hw = half_w * 0.35
            # Box
            Plots.plot!(p, [center - box_hw, center + box_hw, center + box_hw, center - box_hw, center - box_hw],
                [q1, q1, q3, q3, q1], lw = 1.5, color = :gray30, fillalpha = 0.15, fill = true,
                label = nothing)
            # Whiskers
            Plots.plot!(p, [center, center], [wlo, q1], lw = 1, color = :gray30, label = nothing)
            Plots.plot!(p, [center, center], [q3, whi], lw = 1, color = :gray30, label = nothing)
            Plots.plot!(p, [center - box_hw * 0.5, center + box_hw * 0.5], [wlo, wlo], lw = 1, color = :gray30, label = nothing)
            Plots.plot!(p, [center - box_hw * 0.5, center + box_hw * 0.5], [whi, whi], lw = 1, color = :gray30, label = nothing)

            # Median line
            Plots.plot!(p, [lo, hi], [bin_medians[b], bin_medians[b]], lw = 3, color = :red,
                label = (b == 1 ? "bin median" : nothing))
        end

        for b in 1:nb+1
            Plots.vline!(p, [bin_edges[b]], lw = 0.5, ls = :dot, color = :gray60, label = nothing)
        end

        push!(nce_plots, p)
    end

    return nce_model, nce_plots
end

function generate_wide_scout_plot(wide_frags, scout_model, parsed_fname)
    (scout_model === nothing || length(wide_frags) < 20) && return nothing
    frag_mzs = Float64[s.theoretical_mz for s in wide_frags]
    da_errs = Float64[s.observed_mz - s.theoretical_mz for s in wide_frags]
    mz_range = range(minimum(frag_mzs), maximum(frag_mzs), length=200)
    tol_mda = Float64(scout_model.tolerance_da * 1e3)
    bias_mda = [Float64(_scout_mz_bias_da(scout_model, Float32(m))) * 1e3 for m in mz_range]

    p = Plots.scatter(frag_mzs, da_errs .* 1e3,
        alpha=0.1, markersize=1.5, color=:steelblue, label=nothing,
        xlabel="Fragment m/z", ylabel="Raw error (mDa)",
        title="$(parsed_fname)\nWide scout m/z bias (n=$(length(wide_frags)), tol=±$(round(tol_mda, digits=1)) mDa)",
        size=(600, 600), topmargin=10Plots.mm)
    Plots.plot!(p, mz_range, bias_mda, lw=2.5, color=:red, label="m/z bias (robust linear)")
    Plots.plot!(p, mz_range, bias_mda .+ tol_mda, lw=1.5, ls=:dash, color=:red,
        label="collection tol: ±$(round(tol_mda, digits=1)) mDa")
    Plots.plot!(p, mz_range, bias_mda .- tol_mda, lw=1.5, ls=:dash, color=:red, label=nothing)
    Plots.hline!(p, [0.0], color=:black, lw=1, ls=:dot, label=nothing)
    return p
end

"""
Process a single MS file to determine optimal mass error and RT parameters.
Uses a scout-then-collect strategy:
1. Scout 500 scans with strict bitmask filter (score≥8) to estimate signal and pick bias
2. Collect PSMs with informed scan count estimate and doubling
3. Fall back to score≥7 if score≥8 is insufficient
"""
function process_file!(
    results::ParameterTuningSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:ParameterTuningSearchParameters}

    file_name = try
        getFileIdToName(getMSData(search_context), ms_file_idx)
    catch
        "file_$ms_file_idx"
    end
    @debug_l1 "ParameterTuning file $ms_file_idx ($file_name)"

    converged = false
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    final_psm_count = 0
    iteration_state = IterationState()

    try
        initialize_models!(search_context, ms_file_idx, params)
        scan_priority = get_ms2_scan_priority_order(spectra)
        total_ms2 = length(scan_priority)
        # Sciex ZT scanning DIA: top-TIC scans are the center bins of the most
        # abundant precursors, and their adjacent Q1 bins are near-duplicate
        # (overlapping ion content). Selecting calibration scans by TIC therefore
        # clusters the sample on a few precursors / RT regions ("blips") and
        # mis-calibrates the rest of the run. Shuffle (seeded per file) so the
        # calibration draws a representative sample across the whole run — same
        # rationale as the bitvec random-sampling fix.
        if something(tryparse(Int, get(ENV, "PIONEER_ZT_METASCAN_K", "")), 0) > 0
            shuffle!(MersenneTwister(1_800_017 + total_ms2), scan_priority)
        end
        if total_ms2 == 0
            iteration_state.failed_with_exception = true
            throw(ErrorException("No usable MS2 scans"))
        end
        min_psms_needed = getMinPsms(params)

        # Two-phase loop: Phase 1 (wide scout) discovers bias → Phase 2 (collection) refines
        phases = (
            (label = "Wide scout",
             mass_model = MassErrorModel(0.0f0, (WIDE_SCOUT_TOL_PPM, WIDE_SCOUT_TOL_PPM)),
             target_psms = WIDE_SCOUT_TARGET_PSMS,
             initial_scans = WIDE_SCOUT_INITIAL_SCANS,
             max_peaks = Int(TUNING_TOPN_PEAKS)),
            (label = "Collection",
             mass_model = nothing,  # filled from scout model after phase 1
             target_psms = Int64(min_psms_needed),
             initial_scans = Int64(TUNING_MIN_COLLECT_SCANS),
             max_peaks = 0),
        )

        scored_psms = DataFrame()
        for (phase_idx, phase) in enumerate(phases)
            t_phase = time()
            model = phase.mass_model !== nothing ? phase.mass_model :
                getMassErrorModel(search_context, ms_file_idx)

            converged, scored_psms, n_scans, _ = accumulate_psms!(
                spectra, search_context, params, ms_file_idx, scan_priority;
                mass_model = model,
                target_psms = phase.target_psms,
                initial_scans = phase.initial_scans,
                max_peaks = phase.max_peaks,
                label = phase.label)
            n_passing = nrow(scored_psms)
            @debug_l1 "  $(phase.label): $(n_scans) scans → $(n_passing) PSMs ($(round(time()-t_phase, digits=2))s)"

            # Extract fragments and fit model
            frags = n_passing > 0 ?
                get_matched_fragments(spectra, scored_psms, search_context, params, ms_file_idx) :
                MassErrSample[]
            n_frags = length(frags)

            if phase_idx == 1
                # Phase 1: fit scout calibration model
                scout_model = n_frags >= SCOUT_MIN_FRAGS ?
                    fit_scout_calibrated_model(frags;
                        with_intensity = n_passing >= SCOUT_INTENSITY_THRESHOLD) : nothing
                if scout_model !== nothing
                    setMassErrorModel!(search_context, ms_file_idx, scout_model)
                else
                    setMassErrorModel!(search_context, ms_file_idx,
                        MassErrorModel(0.0f0, (WIDE_SCOUT_FALLBACK_TOL_PPM, WIDE_SCOUT_FALLBACK_TOL_PPM)))
                    @debug_l1 "  Scout: <$(SCOUT_MIN_FRAGS) frags, fallback ±$(WIDE_SCOUT_FALLBACK_TOL_PPM) ppm"
                end
                iteration_state.wide_scout_plot = generate_wide_scout_plot(frags, scout_model, parsed_fname)

            else
                # Phase 2: fit RT model + final mass error model
                if !converged && n_passing > 0
                    @debug_l1 "Did not converge ($(n_passing) PSMs, need $(min_psms_needed))"
                end

                if n_passing > 0
                    rt_psms = filter_top_psms_per_precursor(scored_psms, 3)
                    if n_passing >= MIN_PSMS_FOR_RT
                        rt_model_data = fit_irt_model(params, rt_psms)
                        set_rt_to_irt_model!(results, search_context, params, ms_file_idx, rt_model_data)
                        @debug_l1 "  RT model: $(n_passing) PSMs"
                    else
                        setRtIrtMap!(search_context, IdentityModel(), ms_file_idx)
                        results.rt_to_irt_model[] = IdentityModel()
                        getIrtErrors(search_context)[ms_file_idx] = typemax(Float32)
                        @debug_l1 "  RT: insufficient PSMs ($(n_passing) < $(MIN_PSMS_FOR_RT))"
                    end

                    iteration_state.best_fragments = frags
                    if n_frags >= MIN_FRAGS_FOR_INTENSITY_MODEL
                        k_val = Float32(quantile(Normal(), (1.0 + TUNING_GAUSSIAN_COVERAGE) / 2.0))
                        fit_and_install_intensity_model!(search_context, ms_file_idx, frags; k=k_val)
                        @debug_l1 "  IntensityMassErrorModel: $(n_frags) frags, k=$(round(k_val, digits=3))"
                    elseif n_frags > 0
                        mass_err_model, _, _ = fit_models_from_fragments(params, frags)
                        mass_err_model !== nothing && setMassErrorModel!(search_context, ms_file_idx, mass_err_model)
                        @debug_l1 "  SimpleMassErrorModel: $(n_frags) frags"
                    end

                    # MS1 mass-error fit (median bias + 5·MAD tolerance) from
                    # MS2-accepted PSMs. Installs into SearchContext for use by
                    # MainSearch MS1 features and IntegrateChromatogramsSearch.
                    try
                        ms1_residuals = collect_ms1_residuals(spectra, scored_psms, search_context, ms_file_idx)
                        parsed_fname_ms1 = getParsedFileName(search_context, ms_file_idx)
                        ms1_dir = joinpath(getDataOutDir(search_context), "qc_plots", "ms1_mass_error_plots")
                        isdir(ms1_dir) || mkpath(ms1_dir)
                        generate_ms1_residual_histogram(ms1_residuals, parsed_fname_ms1,
                            joinpath(ms1_dir, "$(parsed_fname_ms1).png"))
                        fit = fit_ms1_model_from_residuals(ms1_residuals)
                        if fit !== nothing
                            ms1_model, ms1_med, ms1_mad = fit
                            setMs1MassErrorModel!(search_context, ms_file_idx, ms1_model)
                            @debug_l1 "  MS1 model: bias=$(round(ms1_med, digits=2)) ppm, " *
                                      "tol=±$(round(MS1_DIAG_TOL_K_MAD*ms1_mad, digits=2)) ppm " *
                                      "($(length(ms1_residuals)) residuals)"
                        else
                            @debug_l1 "  MS1 model: insufficient residuals ($(length(ms1_residuals)))"
                        end
                    catch ms1_err
                        @debug_l1 "  MS1 diag failed: $(sprint(showerror, ms1_err))"
                    end

                    # NCE sweep: reuse collected PSMs (skip fragment index).
                    # Plots are accumulated for the combined NCE PDF written
                    # by summarize_results!; the redundant per-file PDF was
                    # dropped 2026-06-26 (same rationale as the mass-error PDF).
                    nce_result = fit_nce_from_psms!(search_context, params, ms_file_idx, spectra, scored_psms)
                    if nce_result !== nothing
                        append!(results.nce_plot_objects, nce_result[2])
                    end

                    iteration_state.best_psms = rt_psms
                    iteration_state.best_psm_count = n_passing
                    results.mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)
                    converged = true
                end
            end
        end
        
    catch e
        # Get the actual file name for clear error reporting
        file_name = try
            getFileIdToName(getMSData(search_context), ms_file_idx)
        catch
            "file_$ms_file_idx"
        end
        
        # Don't mark as failed - just continue with defaults
        @user_warn "Parameter tuning failed for MS data file: $file_name. Error type: $(typeof(e)). Using conservative default parameters to continue analysis."
        # Log full stacktrace for diagnostics
        try
            bt = catch_backtrace()
            @user_error sprint(showerror, e, bt)
        catch
        end

        # Set conservative defaults to allow pipeline to continue
        converged = false
        iteration_state.failed_with_exception = true
        
        # Set default IRT error to infinite tolerance (no IRT filtering)
        getIrtErrors(search_context)[ms_file_idx] = typemax(Float32)
        
        # Clear any partial state
    end
    
    # Store results and handle fallback if needed
    if !converged
        # Check if we have a best attempt to use
        if iteration_state.best_mass_error_model !== nothing
            left_tol = round(getLeftTol(iteration_state.best_mass_error_model), digits = 1)
            right_tol = round(getRightTol(iteration_state.best_mass_error_model), digits = 1)
            @debug_l1 "Failed to converge for file $ms_file_idx. " *
                  "Using best attempt: score≥$(iteration_state.best_score), " *
                  "$(iteration_state.best_psm_count) PSMs, " *
                  "$(iteration_state.best_scan_count) scans, " *
                  "offset=$(round(getMassOffset(iteration_state.best_mass_error_model), digits=1)) ppm, " *
                  "tolerance -$left_tol/+$right_tol ppm"
            
            # Apply best attempt models
            setMassErrorModel!(search_context, ms_file_idx, iteration_state.best_mass_error_model)

            # Fit RT model from best PSMs (deferred fitting)
            if iteration_state.best_psms !== nothing
                @debug_l1 "Fitting RT model from best attempt PSMs ($(iteration_state.best_psm_count) PSMs)"
                rt_model_data = fit_irt_model(params, iteration_state.best_psms)
                set_rt_to_irt_model!(results, search_context, params, ms_file_idx, rt_model_data)
            else
                # If no PSMs available, use identity
                @debug_l1 "No PSMs found for RT model, using IdentityModel"
                setRtIrtMap!(search_context, IdentityModel(), ms_file_idx)
                results.rt_to_irt_model[] = IdentityModel()
            end
            
            # Build detailed warning message
            left_tol = round(getLeftTol(iteration_state.best_mass_error_model), digits = 1)
            right_tol = round(getRightTol(iteration_state.best_mass_error_model), digits = 1)
            final_psm_count = iteration_state.best_psm_count
            
        else
            # Determine why we failed to provide appropriate messaging
            if iteration_state.failed_with_exception
                file_name = try
                    getFileIdToName(getMSData(search_context), ms_file_idx)
                catch
                    "file_$ms_file_idx"
                end
                @debug_l1 "Processing failed with exception for MS data file: $file_name. Using conservative default parameters."
            else
                # No valid attempts found any PSMs - use conservative defaults or borrow
                @debug_l1 "Failed to converge for file $ms_file_idx after $(iteration_state.scan_attempt) attempts. " *
                          "No attempts found sufficient PSMs (best was $(iteration_state.best_psm_count)), using fallback strategy"
            end
            
            fallback_mass_err, fallback_rt_model, borrowed_from = get_fallback_parameters(
                search_context, ms_file_idx
            )

            setMassErrorModel!(search_context, ms_file_idx, fallback_mass_err)
            setRtIrtMap!(search_context, fallback_rt_model, ms_file_idx)

            # Seed irt_errors with infinite tolerance so downstream library
            # search has a usable iRT tolerance for sparse files.
            if !haskey(getIrtErrors(search_context), ms_file_idx)
                getIrtErrors(search_context)[ms_file_idx] = typemax(Float32)
            end

            # Update results with fallback
            results.mass_err_model[] = fallback_mass_err
            results.rt_to_irt_model[] = fallback_rt_model
            
            if borrowed_from === nothing
                @debug_l1 "CONSERVATIVE_FALLBACK: No valid attempts, used defaults"
            end
        end
    end
    
    # Record tuning status
    store_final_results!(
        results, search_context, params, ms_file_idx,
        converged, iteration_state.scan_attempt, final_psm_count, iteration_state
    )
    
    # Add to diagnostics
    record_file_status!(
        results.diagnostics, ms_file_idx, parsed_fname,
        converged, !converged, iteration_state
    )

    results.current_iteration_state[] = iteration_state

    # Clear the ParameterTuning LUT filter so downstream tuning methods
    # (NceTuning, QuadTuning) use their own CountFilter, not our LUT.
    delete!(search_context.bitvec_filter, ms_file_idx)

    return results
end

"""
Process search results and generate QC plots for a single MS file.
"""
function process_search_results!(
    results::ParameterTuningSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:ParameterTuningSearchParameters}
    try
        rt_alignment_folder = getRtAlignPlotFolder(search_context)
        mass_error_folder = getMassErrPlotFolder(search_context)
        parsed_fname = getParsedFileName(search_context, ms_file_idx)
        
        # Get iteration_state from results
        iteration_state = results.current_iteration_state[]
        # Note: No mass-error buffer is applied. Plots reflect the fitted model.
        
        # (Plot objects collected into file_rt_plots and file_mass_plots below)

        # Generate RT alignment plot (as Plot object for PDF output)
        file_rt_plots = Plots.Plot[]
        if length(results.rt) > 0
            irt_tol = getIrtErrors(search_context)[ms_file_idx]
            rt_plot = generate_rt_plot(results, parsed_fname; irt_tol=irt_tol)
            push!(file_rt_plots, rt_plot)
        elseif iteration_state !== nothing && iteration_state.best_psms !== nothing &&
               hasproperty(iteration_state.best_psms, :rt) &&
               nrow(iteration_state.best_psms) > 0
            psms = iteration_state.best_psms
            precursors = getPrecursors(getSpecLib(search_context))
            irts = Float32[getIrt(precursors)[pid] for pid in psms[!, :precursor_idx]]
            rts = Float32[getRetentionTimes(spectra)[sid] for sid in psms[!, :scan_idx]]
            p = Plots.scatter(rts, irts,
                alpha=0.1, markersize=2, color=:steelblue, label=nothing,
                xlabel="Retention Time RT (min)",
                ylabel="Indexed Retention Time iRT (min)",
                title=parsed_fname * "\n⚠️ Insufficient PSMs for RT model " *
                      "($(iteration_state.best_psm_count) PSMs, need 750)",
                size=(600, 600))
            push!(file_rt_plots, p)
        else
            fallback_plot = generate_fallback_rt_plot_in_memory(results, parsed_fname, search_context, ms_file_idx)
            if fallback_plot !== nothing
                push!(file_rt_plots, fallback_plot)
            end
        end

        # Generate mass error plot (only store in memory, no individual files)
        current_model = getMassErrorModel(search_context, ms_file_idx)
        has_intensity_model = current_model isa IntensityMassErrorModel

        has_fragments = iteration_state !== nothing &&
                        iteration_state.best_fragments !== nothing &&
                        !isempty(iteration_state.best_fragments)

        # Collect Plot objects for per-file PDF (display-based, correct orientation)
        file_mass_plots = Plots.Plot[]

        # Include wide scout plot if stored in iteration_state
        if iteration_state !== nothing && iteration_state.wide_scout_plot !== nothing
            push!(file_mass_plots, iteration_state.wide_scout_plot)
            iteration_state.wide_scout_plot = nothing
        end

        if has_fragments
            frag_data = extract_fragment_plot_data(iteration_state.best_fragments)
            iteration_state.best_fragments = nothing
            n_frags = length(frag_data.da_errs)

            mda_plots = generate_mass_error_plot_mda(frag_data, current_model, parsed_fname;
)
            if mda_plots !== nothing
                append!(file_mass_plots, mda_plots)
            end

            if has_intensity_model
                @debug_l1 "IntensityMassErrorModel plots: $n_frags fragments, " *
                    "model α=$(current_model.mz_spread_α), β=$(current_model.mz_spread_β), γ=$(current_model.mz_spread_γ)"
                intensity_plots = generate_intensity_model_plots(frag_data, current_model, parsed_fname)
                @debug_l2 "Generated $(length(intensity_plots)) intensity model plots"
                append!(file_mass_plots, intensity_plots)
            end
        else
            fallback_plot = generate_fallback_mass_error_plot_in_memory(results, parsed_fname, search_context, ms_file_idx)
            if fallback_plot !== nothing
                push!(file_mass_plots, fallback_plot)
            end
        end

        # Accumulate Plot objects for the combined mass-error PDF written by
        # summarize_results!. The per-file mass-error PDF was dropped 2026-06-26:
        # writing it took ~2.4 s per file (Plots.jl PDF backend; ~15 s of the
        # ~44 s warm Param Tuning stage on 6-file Olsen). The combined PDF
        # already paginates per file, so no diagnostic info is lost.
        if !isempty(file_mass_plots)
            append!(results.mass_plot_objects, file_mass_plots)
        end

        # Accumulate RT plots for the combined RT-alignment PDF written by
        # summarize_results!. The per-file RT PDF was dropped 2026-06-26 (same
        # rationale as the mass-error PDF: the combined PDF already paginates
        # per file, so no diagnostic info is lost).
        if !isempty(file_rt_plots)
            append!(results.rt_plot_objects, file_rt_plots)
        end

        # Update models in search context
        setMassErrorModel!(search_context, ms_file_idx, getMassErrorModel(results))
        setRtIrtMap!(search_context, getRtToIrtModel(results), ms_file_idx)

        # Clear plotting data to save memory
        resize!(results.rt, 0)
        resize!(results.irt, 0)
        resize!(results.ppm_errs, 0)
        resize!(results.frag_mzs, 0)
        results.current_iteration_state[] = nothing  # Clear iteration state after use
    catch e
        # Plot-generation failures are non-fatal and should not mark the file as failed.
        bt = catch_backtrace()
        @user_error "PLOT GENERATION FAILED for file $ms_file_idx: $(typeof(e))"
        @user_error sprint(showerror, e, bt)
    end
end

"""
Reset results state between files.
"""
function reset_results!(results::ParameterTuningSearchResults)
    # Clear data vectors
    resize!(results.irt, 0)
    resize!(results.rt, 0)
    resize!(results.ppm_errs, 0)
    
    # Models are per-file, so they get reset at the start of each file
    # No need to reset them here
end

# Deprecated: apply_final_mass_error_buffer! removed – buffer is applied per-file before plotting

"""
Summarize results across all MS files.
"""
function summarize_results!(results::ParameterTuningSearchResults, params::P, search_context::SearchContext) where {P<:ParameterTuningSearchParameters}
    # Combine individual plots into merged PDFs using save_multipage_pdf
    
    try
        rt_plots_folder = getRtAlignPlotFolder(search_context)
        mass_error_plots_folder = getMassErrPlotFolder(search_context)
        
        # Combined RT alignment PDF from all files
        if !isempty(results.rt_plot_objects)
            rt_combined_path = joinpath(rt_plots_folder, "rt_alignment_plots.pdf")
            save_multipage_pdf(Plots.Plot[p for p in results.rt_plot_objects], rt_combined_path)
            empty!(results.rt_plot_objects)
        end

        # Combined mass error PDF from all files
        if !isempty(results.mass_plot_objects)
            mass_combined_path = joinpath(mass_error_plots_folder, "mass_error_plots.pdf")
            save_multipage_pdf(Plots.Plot[p for p in results.mass_plot_objects], mass_combined_path)
            empty!(results.mass_plot_objects)
        end

        # Combined NCE alignment PDF from all files
        if !isempty(results.nce_plot_objects)
            nce_dir = joinpath(getDataOutDir(search_context), "qc_plots", "collision_energy_alignment")
            mkpath(nce_dir)
            nce_combined_path = joinpath(nce_dir, "nce_alignment_plots.pdf")
            save_multipage_pdf(Plots.Plot[p for p in results.nce_plot_objects], nce_combined_path)
            empty!(results.nce_plot_objects)
        end
        empty!(results.rt_plots)
        empty!(results.mass_plots)
        
        # Generate summary report
        # TODO: Implement generate_summary_report if detailed report needed
        # For now, the diagnostic summary below provides the key information
        
    catch e
        @user_warn "Failed to merge QC plots" exception=(e, catch_backtrace())
    end
    
    # Apply buffer to all mass error models AFTER plots are generated
    # This ensures plots show the actual fitted values, not buffered ones
    # Do not apply an additional buffer here; models were buffered once per file
    
    # Log diagnostic summary
    diagnostics = getDiagnostics(results)
    # Fixed: use values() to iterate over dictionary values
    file_statuses = values(diagnostics.file_statuses)
    converged_count = sum(s.converged for s in file_statuses)
    fallback_count = sum(s.used_fallback for s in file_statuses)
    total_files = length(file_statuses)
    
    @debug_l1 "Parameter Tuning Summary: $total_files files, $converged_count converged, $fallback_count fallback"
    
end
