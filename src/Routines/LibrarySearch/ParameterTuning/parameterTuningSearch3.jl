struct ParameterTuningSearch <: TuningMethod end

struct ParameterTuningSearchResults <: SearchResults
    rt_to_irt_map::Dict
    frag_err_dist::Dict
    irt_errs::Dict
end

struct PresearchParameters{P<:PrecEstimation} <: SearchParameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_index_search_score::UInt8
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::UInt8
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    sample_rate::Float32
    irt_tol::Float32
    spec_order::Set{Int64}
    prec_estimation::P

    function PresearchParameters(params::Any)
        pp = params[:presearch_params]
        prec_estimation = pp["abreviate_precursor_calc"] ? FullPrecCapture() : PartialPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(first(params[:isotope_err_bounds])), UInt8(last(params[:isotope_err_bounds]))),
            UInt8(pp["min_index_search_score"]),
            Int64(pp["min_frag_count"]),
            Float32(pp["min_spectral_contrast"]),
            Float32(pp["min_log2_matched_ratio"]),
            (Int64(first(pp["min_topn_of_m"])), Int64(last(pp["min_topn_of_m"]))),
            UInt8(pp["max_best_rank"]),
            Int64(pp["n_frag_isotopes"]),
            UInt8(pp["max_frag_rank"]),
            Float32(pp["sample_rate"]),
            typemax(Float32),
            Set(2),
            prec_estimation
        )
    end
end

get_parameters(search_type::ParameterTuningSearch, params::Any) = PresearchParameters(params)


function process_file!(
    results::R,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::Arrow.Table
) where {R<:SearchResults, P<:PresearchParameters}

    try 

        setMassErrorModel!(
            search_context, ms_file_idx,
            MassErrorModel(zero(Float32), (getFragTolPpm(params), getFragTolPpm(params)))
        )
        psms = DataFrame()
        for i in range(1, getMaxPresearchIters(params))
            psms = vcat(LibrarySearch(
                    spectra,
                    UInt32(ms_file_idx),
                    getPresearchFragmentIndex(getSpecLib(search_context)),
                    getSpecLib(search_context),
                    getSearchData(search_context),
                    getQuadTransmissionModel(search_context, ms_file_idx),
                    getMassErrorModel(search_context, ms_file_idx),
                    getRtIrtModel(search_context, ms_file_idx),
                    search_parameters,
                ), psms)
            filterAndScorePsms!(params,psms) >= getMinPsms(params) && break
        end

        rt_model, irt_mad = fitIrtModel(
            search_parameters,
            psms
        )

        fragments = getMatchedFragments(search_parameters, psms)



    catch e 
        throw(e)
    end
    println("size ", size(psms))
end

# Core algorithm implementations using multiple dispatch
function fitIrtModel(params::P, psms::DataFrame) where {P<:PresearchParameters}
    # Initial fit
    rt_to_irt_map = UniformSpline(
        psms[!,:irt_predicted],
        psms[!,:rt],
        getSplineDegree(params),
        getSplineNKnots(spline_knots),
    )
    
    # Calculate residuals
    psms[!,:irt_observed] = rt_to_irt_map.(psms.rt::Vector{Float32})
    residuals::Vector{Float32} = psms[!,:irt_observed]::Vector{Float32} .- psms[!,:irt_predicted]::Vector{Float32}
    irt_mad = mad(residuals)::Float32
    
    # Remove outliers and refit
    valid_psms = @view(psms[abs.(residuals) .< (irt_mad * getOutlierThreshold(params)), :])
    
    final_map = UniformSpline(
        valid_psms[!,:irt_predicted],
        valid_psms[!,:rt],
        getSplineDegree(params),
        getSplineNKnots(spline_knots),
    )
    
    return final_map, irt_mad
end
