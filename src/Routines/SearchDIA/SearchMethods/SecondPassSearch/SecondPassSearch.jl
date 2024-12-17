"""
    SecondPassSearch

Second pass search method using optimized parameters from initial searches.

This search:
1. Uses optimized parameters from first pass search
2. Performs PSM identification with previously calculated Huber delta
3. Tracks retention time windows for efficient searching
4. Records chromatogram information for later analysis

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :second_pass_params => Dict(
        "min_y_count" => 1,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "n_frag_isotopes" => 2,
        "max_frag_rank" => 10
    )
)

# Execute search
results = execute_search(SecondPassSearch(), search_context, params)
```
"""
struct SecondPassSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for second pass search.
"""
struct SecondPassSearchResults <: SearchResults
    psms::Base.Ref{DataFrame}          # PSMs for each file
end

"""
Parameters for second pass search.
"""
struct SecondPassSearchParameters{P<:PrecEstimation,I<:IsotopeTraceType} <: FragmentIndexSearchParameters
    # Core parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    sample_rate::Float32
    spec_order::Set{Int64}

    # Deconvolution parameters
    lambda::Float32
    max_iter_newton::Int64
    max_iter_bisection::Int64
    max_iter_outer::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    max_diff::Float32

    # PSM filtering
    min_y_count::Int64
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::Int64
    
    # Precursor estimation strategy
    isotope_tracetype::I
    prec_estimation::P

    function SecondPassSearchParameters(params::Any) 
        sp = params[:quant_search_params]
        dp = params[:deconvolution_params]
        _ISOTOPE_TRACE_TYPE_ = nothing
        if params[:quant_search_params]["combine_isotope_traces"]
            _ISOTOPE_TRACE_TYPE_ = CombineTraces(Float32(params_[:quant_search_params]["min_fraction_transmitted"]))
            @warn "Combine Traces"
        else
            _ISOTOPE_TRACE_TYPE_ = SeperateTraces()
            @warn "Seperate Traces"
        end

        new{typeof(PartialPrecCapture()),typeof(_ISOTOPE_TRACE_TYPE_)}(
            (UInt8(3), UInt8(0)),  # Fixed isotope bounds
            Int64(sp["n_frag_isotopes"]),
            UInt8(sp["max_frag_rank"]),
            1.0f0,
            Set(2),
            Float32(dp["lambda"]),
            Int64(dp["max_iter_newton"]),
            Int64(dp["max_iter_bisection"]),
            Int64(dp["max_iter_outer"]),
            Float32(dp["accuracy_newton"]),
            Float32(dp["accuracy_bisection"]),
            Float32(dp["max_diff"]),
            Int64(sp["min_y_count"]),
            Int64(sp["min_frag_count"]),
            Float32(sp["min_spectral_contrast"]),
            Float32(sp["min_log2_matched_ratio"]),
            (Int64(first(sp["min_topn_of_m"])), Int64(last(sp["min_topn_of_m"]))),
            Int64(sp["max_best_rank"]),
            _ISOTOPE_TRACE_TYPE_,
            PartialPrecCapture()
        )
    end
end

getIsotopeTraceType(p::SecondPassSearchParameters) = p.isotope_tracetype
#==========================================================
Interface Implementation
==========================================================#

get_parameters(::SecondPassSearch, params::Any) = SecondPassSearchParameters(params)

function init_search_results(::P, search_context::SearchContext) where {P<:SecondPassSearchParameters}
    second_pass_psms = joinpath(getDataOutDir(search_context), "second_pass_psms")
    !isdir(second_pass_psms) && mkdir(second_pass_psms)
    return SecondPassSearchResults(
        DataFrame()
    )
end


#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single file for second pass search.
"""
function process_file!(
    results::SecondPassSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::Arrow.Table
) where {P<:SecondPassSearchParameters}

    try
        setNceModel!(
            getFragmentLookupTable(getSpecLib(search_context)), 
            getNceModelModel(search_context, ms_file_idx)
        )
        # Get RT index
        rt_index = buildRtIndex(
            DataFrame(Arrow.Table(getRtIndex(getMSData(search_context), ms_file_idx))),
            bin_rt_size = 0.1)

        # Perform second pass search
        psms = perform_second_pass_search(
            spectra,
            rt_index,
            search_context,
            params,
            ms_file_idx
        )

        results.psms[] = psms

    catch e
        @warn "Second pass search failed" ms_file_idx exception=e
        rethrow(e)
    end

    return results
end

function process_search_results!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::Arrow.Table
) where {P<:SecondPassSearchParameters}

    try
        # Get PSMs from results container
        psms = results.psms[]
        
        # Add basic search columns (RT, charge, target/decoy status)
        add_second_search_columns!(psms, 
            spectra[:retentionTime],
            getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge], 
            getIsDecoy(getPrecursors(getSpecLib(search_context))),#[:is_decoy],
            getPrecursors(getSpecLib(search_context))
            );

        # Determine which precursor isotopes are captured in each scan's isolation window
        get_isotopes_captured!(
            psms,
            getIsotopeTraceType(params),#.isotope_tracetype,
            getQuadTransmissionModel(search_context, ms_file_idx),
            psms[!, :scan_idx],
            getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge],
            getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
            spectra[:centerMz],
            spectra[:isolationWidthMz]
        )

        # Remove PSMs where only M2+ isotopes are captured (expect poor quantification)
        filter!(row -> first(row.isotopes_captured) < 2, psms)

        # Initialize columns for best scan selection and summary statistics
        psms[!,:best_scan] = zeros(Bool, size(psms, 1));
        init_summary_columns!(psms);

        # Calculate summary scores for each PSM group
        for (key, gpsms) in pairs(groupby(psms, getPsmGroupbyCols(getIsotopeTraceType(params))))
            get_summary_scores!(
                gpsms, 
                gpsms[!,:weight],
                gpsms[!,:gof],
                gpsms[!,:matched_ratio],
                gpsms[!,:fitted_manhattan_distance],
                gpsms[!,:fitted_spectral_contrast],
                gpsms[!,:y_count]
            );
        end

        # Keep only apex scans for each PSM group
        filter!(x->x.best_scan, psms);

        # Add additional features for final analysis
        add_features!(
            psms,
            search_context,
            spectra[:TIC],
            spectra[:mz_array],
            ms_file_idx,
            getRtIrtModel(search_context, ms_file_idx),
            getPrecursorDict(search_context)
        )

        # Initialize probability scores (will be calculated later)
        psms[!,:prob], psms[!,:max_prob], psms[!,:mean_prob], psms[!,:min_prob] = zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1))
        
        # Save processed results
        temp_path = joinpath(
            getDataOutDir(search_context),
            "second_pass_psms",
            getParsedFileName(search_context, ms_file_idx) * ".arrow"
        )
        Arrow.write(temp_path, psms)
        setSecondPassPsms!(getMSData(search_context), ms_file_idx, temp_path)
    catch e
        @warn "Failed to process search results" ms_file_idx exception=e
        rethrow(e)
    end

    return nothing
end


"""
Reset results containers.
"""
function reset_results!(results::SecondPassSearchResults)
    results.psms[] = DataFrame()
end

function summarize_results!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:SecondPassSearchParameters}

    return nothing
end
