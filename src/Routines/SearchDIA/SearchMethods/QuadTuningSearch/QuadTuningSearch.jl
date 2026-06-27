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
    QuadTuningSearch

Search method for optimizing quadrupole transmission models.

This search:
1. Uses fitted MassErrorModel from ParameterTuningSearch for accurate fragment matching
2. Collects PSMs with extended precursor isotope patterns
3. Performs deconvolution to estimate relative isotope abundances
4. Fits transmission model based on isotope ratio deviations
5. Stores optimized models in SearchContext for other methods

Note: Mass tolerances are determined by the fitted MassErrorModel from ParameterTuningSearch,
not by parameter settings. This ensures data-driven, instrument-specific tolerances.

# Configuration
Quadrupole tuning parameters are configured in the parameter_tuning.quad_tuning section:

```json
{
    "parameter_tuning": {
        "quad_tuning": {
            "min_psms_per_thompson": 250,  // PSMs required per Thomson isolation width
            "min_fragments": 3,             // Minimum fragments for PSM acceptance
            "initial_percent": 2.5,         // Initial sampling percentage
            "min_initial_scans": 5000       // Minimum scans for initial sample
        }
    }
}
```

The search automatically:
- Calculates dynamic PSM requirements based on isolation window width
- Uses progressive sampling starting at max(initial_percent, min_initial_scans/total_scans)
- Prioritizes scans by m/z distribution for optimal spectral coverage
"""
struct QuadTuningSearch <: TuningMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for quadrupole tuning search.
"""
struct QuadTuningSearchResults <: SearchResults
    tuning_results::Vector{Vector{@NamedTuple{
        precursor_idx::UInt32,
        scan_idx::UInt32,
        weight::Float32,
        iso_idx::UInt8,
        center_mz::Float32,
        n_matches::UInt8
    }}}
    quad_model::Base.Ref{QuadTransmissionModel}
    quad_plot_dir::String
    quad_plot_objects::Vector{Any}
    # Per-file (file_name, fitted_model, window_width) for the Razo LM overlay.
    per_file_models::Vector{Tuple{String, QuadTransmissionModel, Float64}}
end

"""
Parameters for quadrupole tuning search.
Configures deconvolution and quad model fitting.
"""
struct QuadTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # Search parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_index_search_score::UInt8
    min_log2_matched_ratio::Float32
    min_spectral_contrast::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::UInt8
    max_frag_rank::UInt8
    n_frag_isotopes::Int64
    irt_tol::Float32
    spec_order::Set{Int64}
    relative_improvement_threshold::Float32
    
    # Deconvolution solver
    deconvolution_solver::DeconvolutionSolver
    max_iter_outer::Int64
    max_diff::Float32

    # Quad tuning specific parameters
    min_quad_tuning_fragments::Int64
    min_quad_tuning_psms_per_thompson::Int64
    initial_percent::Float32
    prec_estimation::P

    function QuadTuningSearchParameters(params::PioneerParameters)
        prec_estimation = PartialPrecCapture()

        new{typeof(prec_estimation)}(
            (UInt8(0), UInt8(0)),                                          # isotope_err_bounds
            UInt8(maximum(TUNING_MIN_SCORE)),                              # min_index_search_score (formerly derived from fragment_index_search.min_score=15 → max(TUNING_MIN_SCORE))
            typemin(Float32),                                              # min_log2_matched_ratio
            TUNING_MIN_SPECTRAL_CONTRAST,                                  # min_spectral_contrast
            (Int64(first(TUNING_MIN_TOPN_OF_M)), Int64(last(TUNING_MIN_TOPN_OF_M))), # min_topn_of_m
            UInt8(1),                                                      # max_best_rank
            TUNING_MAX_FRAG_RANK,                                          # max_frag_rank
            Int64(1),                                                      # n_frag_isotopes
            typemax(Float32),                                              # irt_tol
            Set{Int64}([2]),                                               # spec_order
            TUNING_RELATIVE_IMPROVEMENT_THRESHOLD,                         # relative_improvement_threshold
            PoissonMMSolver(),                                             # deconvolution_solver — hardcoded
            DECONV_MAX_ITER,                                               # max_iter_outer
            DECONV_CONVERGENCE_TOL,                                        # max_diff
            QUAD_TUNING_MIN_FRAGMENTS,                                     # min_quad_tuning_fragments
            QUAD_TUNING_MIN_PSMS_PER_THOMPSON,                             # min_quad_tuning_psms_per_thompson
            QUAD_TUNING_INITIAL_PERCENT,                                   # initial_percent
            prec_estimation                                                # prec_estimation
        )
    end
end

#==========================================================
Interface Implementation
==========================================================#

# Getters
getQuadModel(q::QuadTuningSearchResults) = q.quad_model[]

# Setters
function setQuadModel(q::QuadTuningSearchResults, model::Q) where {Q<:QuadTransmissionModel}
    q.quad_model[] = model
end

get_parameters(::QuadTuningSearch, params::Any) = QuadTuningSearchParameters(params)

function init_search_results(::QuadTuningSearchParameters, search_context::SearchContext)
    # Initialize empty tuning results vector in each search data structure
    out_dir = getDataOutDir(search_context)
    qc_dir = joinpath(out_dir, "qc_plots")
    !isdir(qc_dir) && mkdir(qc_dir)

    qpp = joinpath(qc_dir, "quad_transmission_model")
    mkpath(qpp)
    temp_data = Vector{Vector{@NamedTuple{
            precursor_idx::UInt32,
            scan_idx::UInt32,
            weight::Float32,
            iso_idx::UInt8,
            center_mz::Float32,
            n_matches::UInt8
        }}}()#(undef, length(getSearchData(search_context)))
    for i in range(1, length(getSearchData(search_context)))
        push!(temp_data, Vector{@NamedTuple{
            precursor_idx::UInt32,
            scan_idx::UInt32,
            weight::Float32,
            iso_idx::UInt8,
            center_mz::Float32,
            n_matches::UInt8
        }}())
    end
    return QuadTuningSearchResults(
        temp_data,
        Ref{QuadTransmissionModel}(),
        qpp,
        Any[],
        Tuple{String, QuadTransmissionModel, Float64}[]
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Main file processing method for quad tuning search.
"""
function process_file!(
    results::QuadTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData) where {P<:QuadTuningSearchParameters}

    setQuadTransmissionModel!(search_context, ms_file_idx, SquareQuadModel(0.5f0))

    # Get file name for debugging
    file_name = try
        getFileIdToName(getMSData(search_context), ms_file_idx)
    catch
        "file_$ms_file_idx"
    end

    # Log fitted mass error model from ParameterTuningSearch
    fitted_model = getMassErrorModel(search_context, ms_file_idx)

    # (Removed try/catch — let errors propagate so they can be diagnosed.)
    begin
        t_file_start = time()
        # Check if file has any scans
        if length(spectra) == 0
            @user_warn "Skipping quad tuning for $file_name — file contains no scans"
            setQuadModel(results, SquareQuadModel(0.0f0))
            return results
        end

        # Check for inconsistent array types (data quality issue)
        try
            # Test array access - this will fail if there are type mismatches
            test_scan_idx = findfirst(i -> getMsOrder(spectra, i) == 2, 1:length(spectra))
            if test_scan_idx !== nothing
                _ = getMzArray(spectra, test_scan_idx)
                _ = getIntensityArray(spectra, test_scan_idx)
            end
        catch type_error
            if isa(type_error, MethodError) || contains(string(type_error), "SubArray")
                @user_warn "Data type inconsistency detected in $file_name. Array types don't match expected schema. This may be due to dummy/test data with inconsistent typing. Skipping quad tuning."
                setQuadModel(results, SquareQuadModel(0.0f0))
                return results
            else
                rethrow(type_error)
            end
        end

        # Adjust arrays for isotope variants
        adjust_precursor_arrays!(search_context)

        # Check window widths
        window_widths = check_window_widths(spectra)
        if length(window_widths) != 1
            @user_warn "Multiple window sizes detected: $(join(collect(window_widths), ';'))"
            setQuadModel(results, SquareQuadModel(0.0f0))
            return results
        end
        window_width = first(window_widths)
        # Build scan priority index (metadata only, no peak data)
        scan_index = build_quad_scan_priority_index(spectra)

        # Two-tier acceptance: try for 2× the per-Thompson rate; fall back to
        # using whatever we collected as long as it clears the per-Thompson
        # minimum, else SquareQuadModel.
        target_psms_int   = round(Int, 2 * params.min_quad_tuning_psms_per_thompson * window_width)
        min_psms_int      = round(Int,     params.min_quad_tuning_psms_per_thompson * window_width)

        initial_psms, converged, _ = progressive_quad_psm_collection!(
            scan_index,
            spectra,
            search_context,
            params,
            ms_file_idx;
            target_psms=target_psms_int,
            fallback_min_psms=min_psms_int
        )
        n_collected = nrow(initial_psms)

        # Models used downstream — Razo LM (active) or SquareQuad fallback.
        active_model::QuadTransmissionModel = SquareQuadModel(0.0f0)
        razo_initial_model::Union{Nothing, RazoQuadModel} = nothing
        total_psms = DataFrame()
        use_fallback = true

        if !converged
            @user_warn "QuadTuning [$file_name]: only $n_collected PSMs collected (target=$target_psms_int, fallback_min=$min_psms_int); using SquareQuadModel fallback"
        else
            total_psms = process_quad_pipeline(initial_psms, spectra, search_context, results, params, ms_file_idx, window_width)
            if nrow(total_psms) >= 50
                fitted_params, initial_params = fit_quad_model(total_psms, window_width)
                fitted_model = RazoQuadModel(fitted_params)
                razo_initial_model = RazoQuadModel(initial_params)
                active_model = fitted_model
                use_fallback = false
            end
        end

        setQuadModel(results, active_model)

        # Per-file QC plots. Fallback files still get the SquareQuad transmission
        # plot; the scatter/median plots are skipped when there's no data.
        file_plots = Plots.Plot[]
        fname = getFileIdToName(getMSData(search_context), ms_file_idx)
        razo_lm_for_plot = use_fallback ? nothing :
            (@isdefined(fitted_model) ? fitted_model : nothing)
        if !isempty(total_psms)
            push!(file_plots, plot_charge_distributions(total_psms, results, fname;
                                                         quad_model=razo_lm_for_plot,
                                                         window_width=window_width))
        end
        push!(file_plots, plot_quad_model(active_model, window_width, results, fname;
                                           initial_model=razo_initial_model))
        if !isempty(total_psms)
            push!(file_plots, plot_sliding_median_smoother(total_psms, fname;
                                                           window_width=window_width,
                                                           quad_model=razo_lm_for_plot,
                                                           initial_model=razo_initial_model))
        end

        # Accumulate plot objects for the combined quad-transmission PDF
        # written by summarize_results!. The per-file PDF was dropped
        # 2026-06-26 (same rationale as the ParameterTuningSearch mass-error
        # PDF: the combined PDF paginates per file, so no diagnostic info is
        # lost, and writing it doubled the QuadTuning plot I/O cost).
        parsed_fname = getParsedFileName(search_context, ms_file_idx)
        append!(results.quad_plot_objects, file_plots)
        push!(results.per_file_models, (parsed_fname, active_model, Float64(window_width)))

        t_wall = time() - t_file_start
        if !use_fallback
            p = fitted_model.params
            @debug_l1 "  QuadTuning [$file_name]: $n_collected PSMs, $(nrow(total_psms)) deconv pts, wall=$(round(t_wall,digits=2))s; Razo LM: al=$(round(p.al,digits=2)) ar=$(round(p.ar,digits=2)) bl=$(round(p.bl,digits=1)) br=$(round(p.br,digits=1))"
        end
    end

    return results
end

function process_search_results!(
    results::QuadTuningSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    ::MassSpecData
) where {P<:QuadTuningSearchParameters}

    setQuadTransmissionModel!(search_context, ms_file_idx, getQuadModel(results))
end

function summarize_results!(
    results::QuadTuningSearchResults,
    ::P,
    search_context::SearchContext
) where {P<:QuadTuningSearchParameters}
    
    # Cross-file overlay plots: one per model variant so we can see how
    # consistent each model's fit is across files.
    if !isempty(results.per_file_models)
        max_window = maximum(w for (_, _, w) in results.per_file_models)
        half_width = 2.0 + max_window / 2
        plot_bins = LinRange(-half_width, half_width, 200)

        # Razo LM overlay — split into pages of 12 files so the legend stays readable.
        files_per_page = 12
        n_files = length(results.per_file_models)
        n_pages = cld(n_files, files_per_page)
        for page in 1:n_pages
            lo = (page - 1) * files_per_page + 1
            hi = min(page * files_per_page, n_files)
            title_suffix = n_pages > 1 ? " ($(lo)-$(hi) of $n_files)" : ""
            overlay = plot(title="Per-file quad transmission — Razo LM$title_suffix",
                           xlabel="m/z offset", ylabel="transmission",
                           legend=:outertopright, size=(800, 500))
            for (name, model, _) in results.per_file_models[lo:hi]
                f = getQuadTransmissionFunction(model, 0.0f0, 2.0f0)
                plot!(overlay, plot_bins, f.(plot_bins), lw=1.5, alpha=0.6, label=name)
            end
            push!(results.quad_plot_objects, overlay)
        end
    end

    if !isempty(results.quad_plot_objects)
        combined_path = joinpath(results.quad_plot_dir, "quad_transmission_plots.pdf")
        save_multipage_pdf(Plots.Plot[p for p in results.quad_plot_objects], combined_path)
        empty!(results.quad_plot_objects)
    end
    empty!(results.per_file_models)

    reset_precursor_arrays!(search_context)
    return nothing
end

function reset_results!(results::QuadTuningSearchResults)
    for r in results.tuning_results
        empty!(r)
    end
    return nothing
end
