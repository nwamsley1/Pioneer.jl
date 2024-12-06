struct NceTuningSearch <: TuningMethod end

struct NceTuningSearchResults <: SearchResults
    nce_models::Dict{Int64, NceModel}  # Store NCE models for each file
    nce_psms::DataFrame  # Store PSMs from NCE grid search
end

struct NceTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    frag_tol_ppm::Float32
    min_index_search_score::UInt8
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::UInt8
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    sample_rate::Float32
    spec_order::Set{Int64}
    nce_grid::LinRange{Float32, Int64}
    nce_breakpoint::Float32
    max_q_val::Float32
    prec_estimation::P

    function NceTuningSearchParameters(params::Any)
        pp = params[:presearch_params]
        prec_estimation = pp["abreviate_precursor_calc"] ? FullPrecCapture() : PartialPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(first(params[:isotope_err_bounds])), UInt8(last(params[:isotope_err_bounds]))),
            Float32(pp["frag_tol_ppm"]),
            UInt8(pp["min_index_search_score"]),
            Int64(pp["min_frag_count"]),
            Float32(pp["min_spectral_contrast"]),
            Float32(pp["min_log2_matched_ratio"]),
            (Int64(first(pp["min_topn_of_m"])), Int64(last(pp["min_topn_of_m"]))),
            UInt8(pp["max_best_rank"]),
            Int64(pp["n_frag_isotopes"]),
            UInt8(pp["max_frag_rank"]),
            Float32(pp["sample_rate"]),
            Set(2),
            LinRange(21.0f0, 40.0f0, 15),  # Default NCE grid
            NCE_MODEL_BREAKPOINT,  # Global constant for NCE model
            0.01f0,  # Default q-value threshold
            prec_estimation
        )
    end
end

get_parameters(search_type::NceTuningSearch, params::Any) = NceTuningSearchParameters(params)

function init_search_results(search_parameters::NceTuningSearchParameters, search_context::SearchContext, ms_file_idx::Int64)
    return NceTuningSearchResults(
        Dict{Int64, NceModel}(),
        DataFrame()
    )
end

function process_file!(
    results::NceTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::Arrow.Table) where {P<:NceTuningSearchParameters}

    try
        # Get models from context
        #rt_to_irt_model = getRtIrtModel(search_context, ms_file_idx)
        ##mass_error_model = getMassErrorModel(search_context, ms_file_idx)
        #irt_err = getIrtErrs(search_context)[getParsedFileName(search_context, ms_file_idx)]
        #quad_model = getQuadTransmissionModel(search_context, ms_file_idx)

            
        psms = library_search(spectra, search_context, params, ms_file_idx)

        # Add necessary columns
        addPreSearchColumns!(
            psms,
            spectra,
            getPrecursors(getSpecLib(search_context))[:is_decoy],
            getPrecursors(getSpecLib(search_context))[:irt],
            getPrecursors(getSpecLib(search_context))[:prec_charge],
            spectra[:retentionTime],
            spectra[:TIC]
        )

        # Score and filter PSMs
        scorePresearch!(psms)
        getQvalues!(psms[!, :prob], psms[!, :target], psms[!, :q_value])

        # Get best PSMs per precursor and scan
        spsms = combine(groupby(psms, [:precursor_idx, :scan_idx])) do group
            max_idx = argmax(group[!, :scribe])
            return group[max_idx:max_idx, :]
        end

        # Filter by q-value and target
        filter!(row -> row.target && row.q_value <= params.max_q_val, spsms)
        
        # Get passing precursors and filter original PSMs
        passing_precs = Set(spsms[!, :precursor_idx])
        filter!(row -> row.precursor_idx ∈ passing_precs, psms)

        # Add precursor m/z and find best PSMs
        psms[!, :prec_mz] = [getPrecursors(getSpecLib(search_context))[:mz][pid] for pid in psms[!, :precursor_idx]]
        psms[!, :best_psms] .= false
        
        for group in groupby(psms, :precursor_idx)
            best_idx = argmax(group[!, :scribe])
            group[best_idx, :best_psms] = true
        end
        
        filter!(row -> row.best_psms, psms)

        # Fit NCE model and store results
        nce_model = fit_nce_model(
            PiecewiseNceModel(0.0f0),
            psms[!, :prec_mz],
            psms[!, :nce],
            psms[!, :charge],
            params.nce_breakpoint
        )
        
        results.nce_models[ms_file_idx] = nce_model
        append!(results.nce_psms, psms)

    catch e
        @warn "NCE tuning search failed for file index $ms_file_idx" exception=(e, catch_backtrace())
        rethrow(e)
    end

    return results
end

function process_search_results!(results::NceTuningSearchResults, params::P, search_context::SearchContext, ms_file_idx::Int64) where {P<:NceTuningSearchParameters}
    # Store NCE model in search context
    setNceModel!(search_context, ms_file_idx, results.nce_models[ms_file_idx])
end

function summarize_results!(results::NceTuningSearchResults, params::P, search_context::SearchContext) where {P<:NceTuningSearchParameters}
    # Could add summary statistics or plots about NCE optimization if desired
    return nothing
end

function reset_results!(results::NceTuningSearchResults)
    empty!(results.nce_models)
    empty!(results.nce_psms)
end