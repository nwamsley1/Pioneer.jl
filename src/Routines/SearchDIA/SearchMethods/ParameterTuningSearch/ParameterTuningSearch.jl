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
    if nrow(psms) < 50
        record_calibration_qc!(search_context.calibration_qc, :nce, ms_file_idx,
            assess_calibration_qc(:nce,nrow(psms),(NaN,NaN,NaN,NaN); min_support=50))
        select_calibration_plot!(search_context.calibration_qc, :nce, ms_file_idx) &&
            calibration_notice!(search_context, :nce, ms_file_idx)
        return nothing
    end

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
        intensity_model = prepare_fragment_intensity_model(ion_list, nce_model)
        tasks = map(thread_tasks) do thread_task
            Threads.@spawn process_scans_fused!(
                last(thread_task), spectra, prec_index,
                ms_file_idx,
                search_data[first(thread_task)], params, precursors, ion_list,
                intensity_model, qtm, mem, rt_to_irt, irt_tol)
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
        record_calibration_qc!(search_context.calibration_qc, :nce, ms_file_idx,
            assess_calibration_qc(:nce,nrow(nce_psms),(NaN,NaN,NaN,NaN); min_support=50))
        select_calibration_plot!(search_context.calibration_qc, :nce, ms_file_idx) &&
            calibration_notice!(search_context, :nce, ms_file_idx)
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

    weak = 0
    for group in groupby(nce_psms, :precursor_idx)
        best_score = group.gof[1]
        alternative = findfirst(!=(group.nce[1]), group.nce)
        weak += alternative === nothing || !isfinite(best_score) ||
            best_score - group.gof[alternative] <= 0.01 * max(abs(best_score), eps(Float32))
    end
    endpoints = count(x -> x == first(nce_grid) || x == last(nce_grid), best_nce.nce)
    supported = count(c -> 1 <= c <= 6 && nce_model.offsets[Int(c)] != 0, best_nce.charge)
    record_calibration_qc!(search_context.calibration_qc, :nce, ms_file_idx,
        assess_calibration_qc(:nce,n_precs,(supported/n_precs,NaN,weak/n_precs,endpoints/n_precs);
            min_support=50, fallback=supported == 0))
    if select_calibration_plot!(search_context.calibration_qc, :nce, ms_file_idx)
        render_calibration_safely(search_context, :nce, ms_file_idx) do
            if any(charge -> count(==(charge), best_nce.charge) >= 10, charges)
                plot_nce_calibration!(search_context, ms_file_idx, best_nce, nce_grid, nce_model, charges)
            else
                calibration_notice!(search_context, :nce, ms_file_idx)
            end
        end
    end
    return nce_model
end

function plot_nce_calibration!(search_context, ms_file_idx, best_nce, nce_grid, nce_model, charges)
    # Generate per-charge diagnostic plots
    parsed_fname = calibration_qc_title(search_context, :nce, ms_file_idx, getParsedFileName(search_context, ms_file_idx))
    plot_rng = MersenneTwister(1844 + ms_file_idx)
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
            jittered_x = [center + half_w * 0.8 * (2 * rand(plot_rng) - 1) for _ in eachindex(bin_mz)]
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

        write_calibration_page!(search_context, :nce, p)
    end

    return nothing
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
                        ms1_coordinates = (Float32[], Float32[])
                        ms1_residuals = collect_ms1_residuals(spectra, scored_psms, search_context, ms_file_idx;
                            qc_coordinates=ms1_coordinates)
                        parsed_fname_ms1 = getParsedFileName(search_context, ms_file_idx)
                        ms1_dir = joinpath(getDataOutDir(search_context), "qc_plots", "ms1_mass_error_plots")
                        isdir(ms1_dir) || mkpath(ms1_dir)
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
                        if fit !== nothing || any(i -> getMsOrder(spectra, i) == 1, 1:length(spectra))
                            trend = if fit === nothing || !all(isfinite, (fit[2],fit[3]))
                                NaN
                            else
                                corrected = ms1_residuals .- fit[2]
                                maximum(_qc_binned_bias(xs, corrected) for xs in ms1_coordinates)
                            end
                            record_calibration_qc!(search_context.calibration_qc, :ms1_mass, ms_file_idx,
                                assess_calibration_qc(:ms1_mass,length(ms1_residuals),
                                    (NaN,fit === nothing ? NaN : fit[3],trend,NaN);
                                    min_support=10, failed=fit !== nothing && !all(isfinite, (fit[2],fit[3]))))
                            if select_calibration_plot!(search_context.calibration_qc, :ms1_mass, ms_file_idx)
                                if fit === nothing
                                    calibration_notice!(search_context, :ms1_mass, ms_file_idx)
                                else
                                    render_calibration_safely(search_context, :ms1_mass, ms_file_idx) do
                                        generate_ms1_residual_histogram(ms1_residuals,
                                            calibration_qc_title(search_context, :ms1_mass, ms_file_idx, parsed_fname_ms1),
                                            joinpath(ms1_dir, "$(parsed_fname_ms1).png"))
                                    end
                                end
                            end
                        end
                    catch ms1_err
                        record_calibration_qc!(search_context.calibration_qc, :ms1_mass, ms_file_idx,
                            assess_calibration_qc(:ms1_mass,0,(NaN,NaN,NaN,NaN); failed=true))
                        select_calibration_plot!(search_context.calibration_qc, :ms1_mass, ms_file_idx) &&
                            calibration_notice!(search_context, :ms1_mass, ms_file_idx)
                        @debug_l1 "  MS1 diag failed: $(sprint(showerror, ms1_err))"
                    end

                    try
                        fit_nce_from_psms!(search_context, params, ms_file_idx, spectra, scored_psms)
                    catch err
                        record_calibration_qc!(search_context.calibration_qc, :nce, ms_file_idx,
                            assess_calibration_qc(:nce,0,(NaN,NaN,NaN,NaN); failed=true))
                        select_calibration_plot!(search_context.calibration_qc, :nce, ms_file_idx) &&
                            calibration_notice!(search_context, :nce, ms_file_idx)
                        rethrow()
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
    state = results.current_iteration_state[]
    name = getParsedFileName(search_context, ms_file_idx)
    try
        rt_record = rt_calibration_qc(results.rt, results.irt, getRtToIrtModel(results),
            getRetentionTimes(spectra); min_support=MIN_PSMS_FOR_RT)
        record_calibration_qc!(search_context.calibration_qc, :rt, ms_file_idx, rt_record)
        if select_calibration_plot!(search_context.calibration_qc, :rt, ms_file_idx)
            if isempty(results.rt)
                calibration_notice!(search_context, :rt, ms_file_idx)
            else
                render_calibration_safely(search_context, :rt, ms_file_idx) do
                    p = generate_rt_plot(results, calibration_qc_title(search_context, :rt, ms_file_idx, name);
                        irt_tol=get(getIrtErrors(search_context), ms_file_idx, Inf32))
                    write_calibration_page!(search_context, :rt, p)
                end
            end
        end
        fragments = state === nothing ? nothing : state.best_fragments
        model = getMassErrorModel(search_context, ms_file_idx)
        fallback = results.diagnostics.file_statuses[ms_file_idx].used_fallback
        record_calibration_qc!(search_context.calibration_qc, :ms2_mass, ms_file_idx,
            ms2_calibration_qc(fragments, model; fallback))
        if select_calibration_plot!(search_context.calibration_qc, :ms2_mass, ms_file_idx)
            if fragments === nothing || isempty(fragments)
                calibration_notice!(search_context, :ms2_mass, ms_file_idx)
            else
                render_calibration_safely(search_context, :ms2_mass, ms_file_idx) do
                    data = extract_fragment_plot_data(fragments)
                    title = calibration_qc_title(search_context, :ms2_mass, ms_file_idx, name)
                    plots = generate_mass_error_plot_mda(data, model, title)
                    plots !== nothing && foreach(p -> write_calibration_page!(search_context, :ms2_mass, p), plots)
                    if model isa IntensityMassErrorModel
                        foreach(p -> write_calibration_page!(search_context, :ms2_mass, p),
                            generate_intensity_model_plots(data, model, title))
                    end
                end
            end
        end
    finally
        empty!(results.rt); empty!(results.irt); empty!(results.ppm_errs); empty!(results.frag_mzs)
        results.current_iteration_state[] = nothing
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
    root = joinpath(getDataOutDir(search_context), "qc_plots")
    for (stage, path) in (
        (:rt, joinpath(getRtAlignPlotFolder(search_context), "rt_alignment_plots.pdf")),
        (:ms2_mass, joinpath(getMassErrPlotFolder(search_context), "mass_error_plots.pdf")),
        (:ms1_mass, joinpath(root, "ms1_mass_error_plots", "ms1_calibration_notices.pdf")),
        (:nce, joinpath(root, "collision_energy_alignment", "nce_alignment_plots.pdf")),
    )
        finish_calibration_report!(search_context, stage, path)
    end

    # Log diagnostic summary
    diagnostics = getDiagnostics(results)
    # Fixed: use values() to iterate over dictionary values
    file_statuses = values(diagnostics.file_statuses)
    converged_count = sum(s.converged for s in file_statuses)
    fallback_count = sum(s.used_fallback for s in file_statuses)
    total_files = length(file_statuses)
    
    @debug_l1 "Parameter Tuning Summary: $total_files files, $converged_count converged, $fallback_count fallback"
    
end
