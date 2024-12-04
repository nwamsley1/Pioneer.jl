struct ParameterTuningSearch <: TuningMethod end

struct ParameterTuningSearchResults <: SearchResults 
    mass_err_model::Base.Ref{MassErrorModel}
    irt_to_rt_model::Base.Ref{SplineRtConversionModel}
    qc_plots_folder_path::String
end

# Getters
getMassErrorModel(ptsr::ParameterTuningSearchResults) = ptsr.mass_err_model[]
getIrtToRtModel(ptsr::ParameterTuningSearchResults) = ptsr.irt_to_rt_model[]
getQcPlotsFolder(ptsr::ParameterTuningSearchResults) = ptsr.qc_plots_folder_path
# Setters
function set_mass_err_model!(ptsr::ParameterTuningSearchResults, model::MassErrorModel)
    ptsr.mass_err_model[] = model
end

function set_irt_to_rt_model!(ptsr::ParameterTuningSearchResults, model::SplineRtConversionModel)
    ptsr.irt_to_rt_model[] = model
end

struct PresearchParameters{P<:PrecEstimation} <: SearchParameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_psms::Int64
    max_q_val::Float32
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
    spline_degree::Int64
    spline_n_knots::Int64,
    spline_fit_outlier_sd::Int64
    prec_estimation::P

    function PresearchParameters(params::Any)
        pp = params[:presearch_params]
        prec_estimation = pp["abreviate_precursor_calc"] ? FullPrecCapture() : PartialPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(first(params[:isotope_err_bounds])), UInt8(last(params[:isotope_err_bounds]))),
            Int64(pp["min_samples"]),
            Float3(pp["max_qval"]),
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
            3,
            5,
            10,
            prec_estimation
        )
    end
end

get_parameters(search_type::ParameterTuningSearch, params::Any) = PresearchParameters(params)

init_search_results(search_parameters::P, search_context::SearchContext, ms_file_idx::Int64) = ParameterTuningSearchResults()

function process_file!(
    results::ParameterTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::Arrow.Table) where {P<:PresearchParameters}

    try 

        setMassErrorModel!(
            search_context, ms_file_idx,
            MassErrorModel(zero(Float32), (getFragTolPpm(params), getFragTolPpm(params)))
        )
        psms = DataFrame()
        for i in range(1, getMaxPresearchIters(params))
            new_psms = library_search(spectra, search_context, search_parameters, ms_file_idx)
            iszero(size(new_psms, 1)) && continue
            #Add new psms that pass the scoring thresholds 
            add_columns_and_concat!(pams, new_psms, spectra, getPrecursors(spec_lib), params)
            try 
                filter_and_score_psms!(params,psms) >= getMinPsms(params) && break
            catch
                continue
            end
        end

        #Model conversion between library (irt) and empirical (rt) retention times
        set_irt_to_rt_model!(results, fit_irt_model(search_parameters, psms))

        fragments = get_matched_fragments!(search_parameters, psms)



    catch e 
        throw(e)
    end
    return results
end


#=
function plot_rt_alignment!(results::ParameterTuningSearchResults, psms::DataFrame)
    qc_plot_folder = getQcPlotsFolder(results)
    rt_alignment_folder = joinpath(qc_plot_folder, "rt_alignment")
    if !isdir(rt_alignment_folder)
        mkdir(rt_alignment_folder)
    end
    plotRTAlign(

    getQcPlotsFolder(results)
    )
end
=#
function add_columns_and_concat!(psms::DataFrame, new_psms::DataFrame, spectra::Arrow.Table, precursors::Arrow.Table, params::P) where {P<:PresearchParameters}
    addPreSearchColumns!(psms, spectra, precursors[:is_decoy], precursors[:irt], precursors[:prec_charge], spectra[:retentionTime], spectra[:TIC])
    if new_psms === nothing
        return nothing
    else
        append!(psms, new_psms)
        return nothing
    end
end

function merge_new_psms!(psms::DataFrame, new_psms::DataFrame)
end

function filter_and_score_psms!(psms::DataFrame, params::P) where {P<:PresearchParameters}
        scorePresearch(psms)
        getQvalues!(psms[!,:prob], psms[!,:target], psms[!,:q_value])
        max_q_val, min_psms = getMaxQVal(params), getMinPsms(params)
        n_pssing_psms = sum(rtpsms[!,:q_value].<=max_q_val)
        if n_passing_psms >= min_psms
            filter!(row -> row.q_value::Float32 <= min_psms, psms)
            psms[!,:best_psms] = zeros(Bool, size(psms, 1))
            grouped_psms = groupby(rtpsms,:precursor_idx)
            for sub_psms in grouped_psms
                best_idx = argmax(sub_psms.prob::Float32)
                psms[best_idx,:best_psms] = true
            end
            filter!(x->x.best_psms::Bool, psms)
            return n_passing_psms 
        end
    return -1
end

# Core algorithm implementations using multiple dispatch
function fit_irt_model!(params::P, psms::DataFrame) where {P<:PresearchParameters}
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
    #By defaualt the `mad` function returns an adjusted median absolute deviation
    #that estimates the standard deviation.
    irt_mad = mad(residuals)::Float32
    
    # Remove extreme outliers and refit
    valid_psms = @view(psms[abs.(residuals) .< (irt_mad * getOutlierThreshold(params)), :])
    
    final_map = UniformSpline(
        valid_psms[!,:irt_predicted],
        valid_psms[!,:rt],
        getSplineDegree(params),
        getSplineNKnots(spline_knots),
    )
    
    return final_map
end
