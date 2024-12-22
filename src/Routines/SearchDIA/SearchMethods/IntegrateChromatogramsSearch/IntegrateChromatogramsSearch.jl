"""
    IntegrateChromatogramSearch

Search method for analyzing chromatograms to get quantitative information.

This search:
1. Uses precursor and trace information from previous searches
2. Builds chromatograms for each precursor
3. Integrates areas for quantification
4. Incorporates isotope pattern information
"""
struct IntegrateChromatogramSearch <: SearchMethod end

#==========================================================
Type Definitions  
==========================================================#

"""
Results container for chromatogram integration search.
"""
struct IntegrateChromatogramSearchResults <: SearchResults
    psms::Base.Ref{DataFrame}  # Chromatogram data per file
end

"""
Parameters for chromatogram integration search.
"""
struct IntegrateChromatogramSearchParameters{P<:PrecEstimation, I<:IsotopeTraceType} <: FragmentIndexSearchParameters
    # Core parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    sample_rate::Float32
    spec_order::Set{Int64}

    # Chromatogram parameters
    wh_smoothing_strength::Float32
    n_pad::Int64
    max_apex_offset::Int64
    
    # Deconvolution parameters
    lambda::Float32
    max_iter_newton::Int64
    max_iter_bisection::Int64
    max_iter_outer::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    max_diff::Float32

    # Analysis strategies
    isotope_tracetype::I
    prec_estimation::P

    function IntegrateChromatogramSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        global_params = params.global_settings
        quant_params = params.quant_search
        frag_params = quant_params.fragment_settings
        chrom_params = quant_params.chromatogram
        deconv_params = params.optimization.deconvolution
        
        # Determine isotope trace type
        isotope_trace_type = if haskey(global_params.isotope_settings, :combine_traces) && 
                               global_params.isotope_settings.combine_traces
            CombineTraces(0.0f0)  # Default min_fraction_transmitted
        else
            SeperateTraces()
        end

        # Always use partial precursor capture for integrate chromatogram
        prec_estimation = global_params.isotope_settings.partial_capture ? PartialPrecCapture() : FullPrecCapture()

        new{typeof(prec_estimation), typeof(isotope_trace_type)}(
            (UInt8(3), UInt8(0)),  # Fixed isotope bounds
            Int64(frag_params.n_isotopes),
            UInt8(frag_params.max_rank),
            1.0f0,  # Full sampling
            Set{Int64}([2]),
            
            Float32(chrom_params.smoothing_strength),
            Int64(chrom_params.padding),
            Int64(chrom_params.max_apex_offset),
            
            Float32(deconv_params.lambda),
            Int64(deconv_params.newton_iters),
            Int64(deconv_params.newton_iters),
            Int64(deconv_params.newton_iters),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.max_diff),
            
            isotope_trace_type,
            prec_estimation
        )
    end
end


#==========================================================
Interface Implementation
==========================================================#

get_parameters(::IntegrateChromatogramSearch, params::Any) = IntegrateChromatogramSearchParameters(params)

function init_search_results(::IntegrateChromatogramSearchParameters, search_context::SearchContext)
    return IntegrateChromatogramSearchResults(
        Ref(DataFrame())
    )
end

"""
Process a single file for chromatogram integration.
"""
function process_file!(
    results::IntegrateChromatogramSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData) where {P<:IntegrateChromatogramSearchParameters}

    try
        # Set the NCE model from the search context for fragment matching
        setNceModel!(
            getFragmentLookupTable(getSpecLib(search_context)), 
            getNceModelModel(search_context, ms_file_idx)
        )

        # Build retention time index for efficient precursor lookup
        # Groups RTs into bins of 0.1 minute width
        rt_index = buildRtIndex(
            DataFrame(Arrow.Table(getRtIndex(getMSData(search_context), ms_file_idx))),
            bin_rt_size = 0.1)

        # Load PSMs that passed previous filtering steps
        # Convert to DataFrame for processing
        passing_psms = DataFrame(Tables.columntable(Arrow.Table(getPassingPsms(getMSData(search_context), ms_file_idx))))#load_passing_psms(search_context, parsed_fname)
        
        # Keep only target (non-decoy) PSMs
        filter!(row -> row.target, passing_psms)

        # Initialize columns to store integration results
        # peak_area: Integrated area of chromatographic peak
        # new_best_scan: Updated apex scan after refinement
        passing_psms[!, :peak_area] = zeros(Float32, nrow(passing_psms))
        passing_psms[!, :new_best_scan] = zeros(UInt32, nrow(passing_psms))

        # Extract chromatograms for all passing PSMs
        # Builds chromatograms using parallel processing across scan ranges
        chromatograms = extract_chromatograms(
            spectra,
            passing_psms,
            rt_index,
            search_context,
            params,
            ms_file_idx
        )

        # Determine which isotopes are captured in each isolation window
        # Uses quadrupole transmission model to check isotope coverage
        get_isotopes_captured!(
            chromatograms,
            params.isotope_tracetype,
            getQuadTransmissionModel(search_context, ms_file_idx),
            chromatograms[!, :scan_idx],
            getCharge(getPrecursors(getSpecLib(search_context))),
            getMz(getPrecursors(getSpecLib(search_context))),
            getCenterMzs(spectra),
            getIsolationWidthMzs(spectra)
        )

        # Remove chromatograms where M2+ isotopes are captured
        # These can interfere with accurate quantification
        filter!(row -> first(row.isotopes_captured) < 2, chromatograms)

        # Integrate chromatographic peaks for each precursor
        # Updates peak_area and new_best_scan in passing_psms        
        integrate_precursors(
            chromatograms,
            params.isotope_tracetype,
            passing_psms[!, :precursor_idx],
            passing_psms[!, :isotopes_captured],
            passing_psms[!, :scan_idx],
            passing_psms[!, :peak_area],
            passing_psms[!, :new_best_scan],
            ms_file_idx,
            λ = params.wh_smoothing_strength,
            n_pad = params.n_pad,
            max_apex_offset = params.max_apex_offset
        )

        # Clear chromatograms to free memory
        chromatograms = nothing

        # Store processed PSMs in results
        results.psms[] = passing_psms
    catch e
        # Log error and re-throw for debugging
        @warn "Chromatogram integration failed" ms_file_idx exception=e
        rethrow(e)
    end

    return results
end

function process_search_results!(
    results::IntegrateChromatogramSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:IntegrateChromatogramSearchParameters}
    try
        passing_psms = results.psms[]
        parsed_fname = getFileIdToName(getMSData(search_context), ms_file_idx)
        # Process final PSMs
        process_final_psms!(
            passing_psms,
            search_context,
            parsed_fname,
            ms_file_idx
        )
        # Save results
        
        Arrow.write(
            getPassingPsms(getMSData(search_context))[ms_file_idx],
            passing_psms)
        
    catch e
        @warn "Chromatogram processing failed" ms_file_idx exception=e
        rethrow(e)
    end
    return nothing
end

function reset_results!(results::IntegrateChromatogramSearchResults)
    results.psms[] = DataFrame()
    GC.gc()
    return nothing
end

function summarize_results!(
    ::IntegrateChromatogramSearchResults,
    ::P,
    ::SearchContext
) where {P<:IntegrateChromatogramSearchParameters}
    return nothing
end
