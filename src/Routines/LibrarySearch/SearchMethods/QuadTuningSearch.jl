"""
    QuadTuningSearch

Search method for optimizing quadrupole transmission models.

This search:
1. Collects PSMs with extended precursor isotope patterns
2. Performs deconvolution to estimate relative isotope abundances
3. Fits transmission model based on isotope ratio deviations
4. Stores optimized models in SearchContext for other methods

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :isotope_err_bounds => (0, 0),  # Fixed for quad tuning
    :presearch_params => Dict(
        "frag_tol_ppm" => 30.0,
        "min_index_search_score" => 3,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "n_frag_isotopes" => 2,
        "max_frag_rank" => 10,
        "quad_tuning_sample_rate" => 0.1,
        "min_quad_tuning_fragments" => 3,
        "min_quad_tuning_psms" => 1000,
        "abreviate_precursor_calc" => false
    ),
    :deconvolution_params => Dict(
        "max_iter_newton" => 100,
        "max_iter_bisection" => 100,
        "max_iter_outer" => 100,
        "accuracy_newton" => 1e-5,
        "accuracy_bisection" => 1e-4,
        "max_diff" => 1e-5
    )
)

# Execute search
results = execute_search(QuadTuningSearch(), search_context, params)
```
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
end

"""
Parameters for quadrupole tuning search.
Configures deconvolution and quad model fitting.
"""
struct QuadTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # Search parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    frag_tol_ppm::Float32
    min_index_search_score::UInt8
    min_log2_matched_ratio::Float32
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::UInt8
    max_frag_rank::UInt8
    n_frag_isotopes::Int64
    sample_rate::Float32
    irt_tol::Float32
    spec_order::Set{Int64}
    
    # Deconvolution parameters
    max_iter_newton::Int64
    max_iter_bisection::Int64
    max_iter_outer::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    max_diff::Float32
    
    # Quad tuning specific parameters
    min_quad_tuning_fragments::Int64
    min_quad_tuning_psms::Int64
    prec_estimation::P

    function QuadTuningSearchParameters(params::Any)
        pp = params[:presearch_params]
        dp = params[:deconvolution_params]
        prec_estimation = PartialPrecCapture() #pp["abreviate_precursor_calc"] ? FullPrecCapture() : PartialPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(0), UInt8(0)),  # Fixed for quad tuning
            Float32(pp["frag_tol_ppm"]),
            UInt8(pp["min_index_search_score"]),
            typemin(Float32),
            Int64(pp["min_frag_count"]),
            Float32(pp["min_spectral_contrast"]),
            (Int64(first(pp["min_topn_of_m"])), Int64(last(pp["min_topn_of_m"]))),
            UInt8(pp["max_best_rank"]),
            UInt8(pp["max_frag_rank"]),
            one(Int64),
            Float32(pp["quad_tuning_sample_rate"]),
            typemax(Float32),
            Set(2),

            Int64(dp["max_iter_newton"]),
            Int64(dp["max_iter_bisection"]),
            Int64(dp["max_iter_outer"]),
            Float32(dp["accuracy_newton"]),
            Float32(dp["accuracy_bisection"]),
            Float32(dp["max_diff"]),

            Int64(pp["min_quad_tuning_fragments"]),
            Int64(pp["min_quad_tuning_psms"]),
            prec_estimation
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
    !isdir(qpp) && mkdir(qpp)
    quad_data = joinpath(qpp, "quad_data")
    !isdir(quad_data) && mkdir(quad_data)
    quad_models = joinpath(qpp, "quad_models")
    !isdir(quad_models) && mkdir(quad_models)
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
        qpp
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
    spectra::Arrow.Table) where {P<:QuadTuningSearchParameters}

    setQuadTransmissionModel!(search_context, ms_file_idx, SquareQuadModel(1.0f0))
    
    try
        setNceModel!(
            getFragmentLookupTable(getSpecLib(search_context)), 
            getNceModelModel(search_context, ms_file_idx)
        )
        # Adjust arrays for isotope variants
        adjustPrecursorArrays!(search_context)
        
        # Check window widths
        window_widths = check_window_widths(spectra)
        if length(window_widths) != 1
            @warn "Multiple window sizes detected: $(join(collect(window_widths), ';'))"
            setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
            return results
        end
        
        # Collect and process PSMs
        total_psms = collect_psms(spectra, search_context, results, params, ms_file_idx)
        
        if nrow(total_psms) == 0
            setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
            return results
        end

        # Plot charge states
        plot_charge_distributions(total_psms, results, getFileIdToName(getMSData(search_context), ms_file_idx))
        
        # Fit quad model
        window_width = parse(Float64, first(window_widths))

        fitted_model = fit_quad_model(total_psms, window_width)
        setQuadModel(results, RazoQuadModel(fitted_model))
        
    catch e
        throw(e)
        @warn "Quad transmission function fit failed" exception=e
        setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
    end
    
    return results
end

function process_search_results!(
    results::QuadTuningSearchResults,
    ::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    ::Arrow.Table
) where {P<:QuadTuningSearchParameters}
    
    setQuadTransmissionModel!(search_context, ms_file_idx, getQuadModel(results))
end

function summarize_results!(
    results::QuadTuningSearchResults,
    ::P,
    search_context::SearchContext
) where {P<:QuadTuningSearchParameters}
    
    plot_bins = LinRange(0-3, 0+3, 100)
    for (ms_file_idx, quad_model) in pairs(search_context.quad_transmission_model)
        fname = getFileIdToName(getMSData(search_context), ms_file_idx)
        quad_func = getQuadTransmissionFunction(quad_model, 0.0f0, 2.0f0)
        p = plot(plot_bins, quad_func.(plot_bins), lw = 2, alpha = 0.5, title = "$fname")
        savefig(p, joinpath(results.quad_plot_dir, "quad_models", fname*".pdf"))
    end

    qmp = [x for x in readdir(joinpath(results.quad_plot_dir, "quad_models"), join=true) if endswith(x, ".pdf")]
    if !isempty(qmp)
        merge_pdfs(qmp, 
                  joinpath(results.quad_plot_dir, "quad_models", "quad_model_plots.pdf"), 
                  cleanup=true)
    end

    qmp = [x for x in readdir(joinpath(results.quad_plot_dir, "quad_data"), join=true) if endswith(x, ".pdf")]
    if !isempty(qmp)
        merge_pdfs(qmp, 
                  joinpath(results.quad_plot_dir, "quad_data", "quad_data_plots.pdf"), 
                  cleanup=true)
    end

    resetPrecursorArrays!(search_context)
    return nothing
end

function reset_results!(results::QuadTuningSearchResults)
    for r in results.tuning_results
        empty!(r)
    end
    return nothing
end

#==========================================================
Helper Methods
==========================================================#

"""
Check MS2 window widths in spectra.
"""
function check_window_widths(spectra::Arrow.Table)
    window_widths = Set{String}()
    for i in 1:length(spectra[:isolationWidthMz])
        if spectra[:msOrder][i] == 2 && !ismissing(spectra[:isolationWidthMz][i])
            push!(window_widths, string(spectra[:isolationWidthMz][i]))
        end
    end
    return window_widths
end

"""
Collect and process PSMs for quad tuning.
"""
function collect_psms(
    spectra::Arrow.Table,
    search_context::SearchContext,
    results::QuadTuningSearchResults,
    params::QuadTuningSearchParameters,
    ms_file_idx::Int64
)
    total_psms = DataFrame()
    n = 0
    unique_precursors = Set{UInt32}()
    function getCharges(prec_charges::AbstractVector{UInt8}, precursor_idx::AbstractVector{UInt32})
        charges = zeros(UInt8, length(precursor_idx))
        for i in range(1, length(precursor_idx))
            charges[i] = prec_charges[precursor_idx[i]]
        end
        return charges
    end
    while n < 5
        # Get initial PSMs
        psms = library_search(spectra, search_context, params, ms_file_idx)
        isempty(psms) && return total_psms
        psms = psms[[pid∉unique_precursors for pid in psms[!,:precursor_idx]],:]
        unique_precursors = union(unique_precursors, Set(psms[!,:precursor_idx]))

        # Process PSMs
        processed_psms = process_initial_psms(psms, spectra, search_context)
        
        # Get scan mapping and perform quad search
        scan_idx_to_prec_idx = getScanToPrecIdx(
            processed_psms[!, :scan_idx],
            processed_psms[!, :precursor_idx],
            spectra[:centerMz],
            spectra[:isolationWidthMz]
        )
        
        quad_psms = performQuadTransmissionSearch(
            spectra,
            results,
            scan_idx_to_prec_idx,
            search_context,
            params,
            ms_file_idx
        )

        
        # Filter and process results
        quad_psms[!,:charge] = getCharges(getCharge(getPrecursors(getSpecLib(search_context))), quad_psms[!,:precursor_idx])
        quad_psms = quad_psms[
            filter_quad_psms(
                quad_psms[!,:iso_idx],
                quad_psms[!,:n_matches],
                quad_psms[!,:weight],
                quad_psms[!,:charge],
                params
            ),
            :
        ]
        
        processed_psms = process_quad_results(
            quad_psms,
            getPrecursors(getSpecLib(search_context)),
            getIsoSplines(first(getSearchData(search_context)))
        )
        processed_psms[!,:half_width_mz] = zeros(Float32, size(processed_psms, 1))
        for (i, scan_idx) in enumerate(processed_psms[!,:scan_idx])
            processed_psms[i,:half_width_mz] = spectra[:isolationWidthMz][scan_idx]/2
        end
        keep_data = zeros(Bool, size(processed_psms, 1))
        for i in range(1, size(processed_psms, 1))
            x0 = processed_psms[i,:x0]::Float32
            hw = processed_psms[i,:half_width_mz]::Float32
            if x0 > zero(Float32)
                if (x0 - hw) < (NEUTRON/4 + 0.1)
                    keep_data[i] = true
                end
            else
                if (abs(x0) - hw) < (NEUTRON/2 + 0.1)
                    keep_data[i] = true
                end
            end
        end
        processed_psms = processed_psms[keep_data,:]

        append!(total_psms, processed_psms)
        if size(total_psms, 1) > params.min_quad_tuning_psms
            break
        else
            n += 1
        end
    end
    return total_psms
end

"""
Process initial PSMs from library search.
"""
function process_initial_psms(
    psms::DataFrame,
    spectra::Arrow.Table,
    search_context::SearchContext
)
    addPreSearchColumns!(
        psms,
        spectra,
        getIsDecoy(getPrecursors(getSpecLib(search_context))),#[:is_decoy],
        getIrt(getPrecursors(getSpecLib(search_context))),#[:irt],
        getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge],
        spectra[:retentionTime],
        spectra[:TIC]
    )
    
    scorePresearch!(psms)
    getQvalues!(psms[!, :prob], psms[!, :target], psms[!, :q_value])
    
    filter!(:q_value => x -> x <= 0.01, psms)
    filter!(:target => identity, psms)
    
    psms[!, :best_psms] .= false
    for group in groupby(psms, :precursor_idx)
        best_idx = argmax(group.prob)
        group[best_idx, :best_psms] = true
    end
    
    filter!(row -> row.best_psms, psms)
    return psms
end

"""
Process quad search results.
"""
function process_quad_results(
    psms::DataFrame,
    precursors::BasicLibraryPrecursors,
    iso_splines::IsotopeSplineModel
)
    sort!(psms, [:scan_idx, :precursor_idx, :iso_idx])
    
    processed = hcat(psms, addColumns(
        psms[!, :precursor_idx],
        getMz(precursors),
        getCharge(precursors),#[:prec_charge],
        getSulfurCount(precursors),#[:sulfur_count],
        psms[!, :iso_idx],
        psms[!, :center_mz],
        iso_splines
    ))

    sort!(processed, [:scan_idx, :precursor_idx, :iso_idx])
    combined = combine(groupby(processed, [:scan_idx, :precursor_idx])) do group
        summarizePrecursor(
            group[!, :iso_idx],
            group[!, :center_mz],
            group[!, :iso_mz],
            group[!, :prec_charge],
            group[!, :weight],
            group[!, :δ]
        )
    end
    
    postprocess_combined_results!(combined)
    return combined
end

"""
Additional processing of combined results.
"""
function postprocess_combined_results!(combined::DataFrame)
    filter!(row -> !ismissing(row.yt), combined)
    combined[!, :prec_charge] = UInt8.(combined[!, :prec_charge])
    combined[!, :x0] = Float32.(combined[!, :x0])
    combined[!, :yt] = Float32.(combined[!, :yt])
    combined[!, :x1] = Float32.(combined[!, :x1])
    #filter!(row -> row.prec_charge < 3, combined)
    return combined
end

"""
Plot charge state distributions.
"""
function plot_charge_distributions(psms::DataFrame, results::QuadTuningSearchResults, fname::String)
    p = plot(title = "Quad Model Data for $fname")
    for charge in 2:3
        mask = psms[!, :prec_charge] .== charge
        plot!(p, 
            psms[mask, :x0],
            psms[mask, :yt],
            seriestype=:scatter,
            alpha=0.1,
            label="Charge $charge",
            xlabel = "m/z offset of M0",
            ylabel = L"log(\delta_i\frac{x_0}{x_1})",
            ylim = (-5, 5)
        )
    end
    savefig(p, joinpath(results.quad_plot_dir, "quad_data", fname*".pdf"))
end

"""
Fit quad model to binned PSM data.
"""
function fit_quad_model(psms::DataFrame, window_width::Float64)
    binned_psms = MergeBins(
        psms,
        (-(window_width + 1.0), window_width + 1.0),
        min_bin_size=20,
        min_bin_width=0.1
    )

    return fitRazoQuadModel(
        window_width,
        binned_psms[!, :median_x0],
        binned_psms[!, :median_x1],
        binned_psms[!, :median_yt],
        λ0=1e-2,
        ϵ1=1e-5,
        ϵ2=1e-4,
        ϵ3=1e-5
    )
end

#==========================================================
Array Management Methods
==========================================================#

"""
Adjust arrays to accommodate isotope variants.
"""
function adjustPrecursorArrays!(search_context::SearchContext)
    target_size = length(getPrecursors(getSpecLib(search_context))) * 3 + 1

    # Only adjust if not already at target size
    if length(first(getSearchData(search_context)).precursor_weights) == length(getPrecursors(getSpecLib(search_context)))
        for search_data in getSearchData(search_context)
            search_data.id_to_col = ArrayDict(UInt32, UInt16, target_size)
            
            if length(search_data.precursor_weights) == length(getPrecursors(getSpecLib(search_context)))
                resize!(search_data.precursor_weights, target_size)
                search_data.precursor_weights[length(search_data.precursor_weights):end] .= zero(Float32)
            end
        end
    end
    return search_context
end

"""
Reset arrays to original size.
"""
function resetPrecursorArrays!(search_context::SearchContext)
    original_size = length(getPrecursors(getSpecLib(search_context)))

    if length(first(getSearchData(search_context)).precursor_weights) != original_size
        for search_data in getSearchData(search_context)
            search_data.id_to_col = ArrayDict(UInt32, UInt16, original_size)
            resize!(search_data.precursor_weights, original_size)
        end
    end
    return search_context
end

#==========================================================
Filtering Methods
==========================================================#

"""
Filter quad PSMs based on criteria:
- M0/M1 isotopes only
- Minimum number of matches
- Non-zero abundance
"""
function filter_quad_psms(
    iso_idx::AbstractVector{UInt8},
    n_matches::AbstractVector{UInt8},
    weight::AbstractVector{Float32}, 
    charge::AbstractVector{UInt8},
    params::QuadTuningSearchParameters
)
    n = length(iso_idx)
    mask = Vector{Bool}(undef, n)
    @inbounds for i in 1:n
        mask[i] = (
            (iso_idx[i] < 3) &&  # M0/M1 only
            (n_matches[i] >= params.min_quad_tuning_fragments) &&  # Min matches
            (weight[i] > 0)  &&
            (charge[i] == 2)# Non-zero abundance
        )
    end
    return mask
end

#==========================================================
Quad Transmission Search
==========================================================#

"""
Perform quadrupole transmission search on MS data.
"""
function performQuadTransmissionSearch(
    spectra::Arrow.Table,
    results::QuadTuningSearchResults,
    scan_idx_to_prec_idx::Dictionary{UInt32, Vector{UInt32}},
    search_context::SearchContext,
    params::QuadTuningSearchParameters,
    ms_file_idx::Int64
)
    thread_tasks = partition_scans(spectra, Threads.nthreads())
    scan_idxs = Set(keys(scan_idx_to_prec_idx))

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            # Get thread-specific data structures
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]
            tuning_results = results.tuning_results[thread_id]
            # Get working arrays
            Hs = getHs(search_data)
            weights = getTempWeights(search_data)
            precursor_weights = getPrecursorWeights(search_data)
            residuals = getResiduals(search_data)

            # Process each scan
            for scan_idx in last(thread_task)
                process_scan!(
                    scan_idx,
                    scan_idxs,
                    spectra,
                    tuning_results,
                    search_context,
                    search_data,
                    params,
                    ms_file_idx,
                    Hs,
                    weights,
                    precursor_weights,
                    residuals,
                    scan_idx_to_prec_idx
                )
            end
            
            return DataFrame(tuning_results)
        end
    end
    
    return vcat(fetch.(tasks)...)
end

"""
Process a single scan for quad transmission search.
"""
function process_scan!(
    scan_idx::Int,
    scan_idxs::Set{UInt32},
    spectra::Arrow.Table,
    tuning_results::Vector{@NamedTuple{
        precursor_idx::UInt32,
        scan_idx::UInt32,
        weight::Float32,
        iso_idx::UInt8,
        center_mz::Float32,
        n_matches::UInt8
    }},
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::QuadTuningSearchParameters,
    ms_file_idx::Int64,
    Hs::SparseArray,
    weights::Vector{Float32},
    precursor_weights::Vector{Float32},
    residuals::Vector{Float32},
    scan_idx_to_prec_idx::Dictionary{UInt32, Vector{UInt32}}  # Added this parameter
)
    scan_idx ∉ scan_idxs && return
    
    msn = spectra[:msOrder][scan_idx]
    msn ∉ params.spec_order && return
    
    # Select transitions
    ion_idx, _ = selectTransitions!(
        getIonTemplates(search_data),
        QuadEstimationTransitionSelection(),
        PartialPrecCapture(),
        getFragmentLookupTable(getSpecLib(search_context)),
        scan_idx_to_prec_idx[scan_idx],
        getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
        getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge],
        getSulfurCount(getPrecursors(getSpecLib(search_context))),#[:sulfur_count],
        getIsoSplines(search_data),
        getPrecursorTransmission(search_data),
        getIsotopes(search_data),
        (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
        block_size = 10000
    )
    
    # Match peaks
    nmatches, nmisses = matchPeaks!(
        getIonMatches(search_data),
        getIonMisses(search_data),
        getIonTemplates(search_data),
        ion_idx,
        spectra[:mz_array][scan_idx],
        spectra[:intensity_array][scan_idx],
        getMassErrorModel(search_context, ms_file_idx),
        spectra[:highMz][scan_idx],
        UInt32(scan_idx),
        UInt32(ms_file_idx)
    )
    
    nmatches ≤ 2 && return
    
    # Sort matches
    sort!(@view(getIonMatches(search_data)[1:nmatches]), 
          by = x->(x.peak_ind, x.prec_id),
          alg=QuickSort)
    
    # Build design matrix and deconvolve
    perform_deconvolution!(
        Hs,
        weights,
        precursor_weights,
        residuals,
        search_data,
        nmatches,
        nmisses,
        params
    )
    
    # Record results
    record_scan_results!(
        search_data,
        tuning_results,
        weights,
        Hs,
        scan_idx,
        spectra[:centerMz][scan_idx]
    )
    
    # Reset for next scan
    reset!(getIdToCol(search_data))
    reset!(Hs)
end

"""
Perform deconvolution for a single scan.
"""
function perform_deconvolution!(
    Hs::SparseArray,
    weights::Vector{Float32},
    precursor_weights::Vector{Float32},
    residuals::Vector{Float32},
    search_data::SearchDataStructures,
    nmatches::Int,
    nmisses::Int,
    params::QuadTuningSearchParameters
)
    buildDesignMatrix!(
        Hs,
        getIonMatches(search_data),
        getIonMisses(search_data),
        nmatches,
        nmisses,
        getIdToCol(search_data)
    )
    
    # Resize arrays if needed
    if getIdToCol(search_data).size > length(weights)
        new_entries = getIdToCol(search_data).size - length(weights) + 1000
        resize!(weights, length(weights) + new_entries)
        resize!(getSpectralScores(search_data), length(getSpectralScores(search_data)) + new_entries)
        append!(getUnscoredPsms(search_data), [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])
    end
    
    # Initialize weights
    for i in 1:getIdToCol(search_data).size
        weights[getIdToCol(search_data)[getIdToCol(search_data).keys[i]]] = 
            precursor_weights[getIdToCol(search_data).keys[i]]
    end
    
    # Solve deconvolution problem
    initResiduals!(residuals, Hs, weights)
    solveHuber!(
        Hs,
        residuals,
        weights,
        Float32(100000.0f0),
        0.0f0,
        params.max_iter_newton,
        params.max_iter_bisection,
        params.max_iter_outer,
        params.accuracy_newton,
        params.accuracy_bisection,
        10.0,
        params.max_diff
    )
end

"""
Record results for a single scan.
"""
function record_scan_results!(
    search_data::SearchDataStructures,
    tuning_results::Vector{@NamedTuple{
        precursor_idx::UInt32,
        scan_idx::UInt32,
        weight::Float32,
        iso_idx::UInt8,
        center_mz::Float32,
        n_matches::UInt8
    }},
    weights::Vector{Float32},
    Hs::SparseArray,
    scan_idx::Int,
    center_mz::Float32
)
    #tuning_results = getTuningResults(search_data)
    
    for i in 1:getIdToCol(search_data).size
        id = getIdToCol(search_data).keys[i]
        colid = getIdToCol(search_data)[id]
        
        # Update precursor weights
        search_data.precursor_weights[id] = weights[colid]
        
        # Calculate indices
        isotope_idx = UInt8(((id - 1) % 3) + 1)
        pid = UInt32(((id - 1) ÷ 3) + 1)
        
        # Count matches
        n_matches = sum(Hs.matched[j] 
            for j in Hs.colptr[colid]:(Hs.colptr[colid+1] - 1))
        
        # Record result
        push!(tuning_results, (
            precursor_idx = pid,
            scan_idx = scan_idx,
            weight = weights[colid],
            iso_idx = isotope_idx,
            center_mz = center_mz,
            n_matches = n_matches
        ))
    end
    #Arrow.write("/Users/n.t.wamsley/Desktop/test.arrow", DataFrame(tuning_results))
    #ttable = DataFrame(Tables.columntable(Arrow.Table("/Users/n.t.wamsley/Desktop/test.arrow")))

end