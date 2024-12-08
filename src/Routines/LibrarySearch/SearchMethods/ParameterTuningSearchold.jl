struct ParameterTuningSearch <: TuningMethod end

struct ParameterTuningSearchResults <: SearchResults 
    mass_err_model::Base.Ref{MassErrorModel}
    irt_to_rt_model::Base.Ref{SplineRtConversionModel}
    irt::Vector{Float32}
    rt::Vector{Float32}
    ppm_errs::Vector{Float32}
    qc_plots_folder_path::String
end

# Getters
getMassErrorModel(ptsr::ParameterTuningSearchResults) = ptsr.mass_err_model[]
getIrtToRtModel(ptsr::ParameterTuningSearchResults) = ptsr.irt_to_rt_model[]
getQcPlotsFolder(ptsr::ParameterTuningSearchResults) = ptsr.qc_plots_folder_path
# Setters
function set_mass_err_model!(ptsr::ParameterTuningSearchResults, model::Tuple{MassErrorModel, Vector{Float32}})
    ptsr.mass_err_model[] = model[1]
    append!(ptsr.ppm_errs, model[2])

end

function reset_results!(ptsr::ParameterTuningSearchResults)
    resize!(ptsr.irt, 0)
    resize!(ptsr.rt, 0)
end

struct ParameterTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    frag_tol_ppm::Float32
    frag_err_quantile::Float32
    min_psms::Int64
    max_q_val::Float32
    max_presearch_iters::Int64
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
    spline_n_knots::Int64
    spline_fit_outlier_sd::Int64
    irt_tol_sd::Int64
    prec_estimation::P

    function ParameterTuningSearchParameters(params::Any)
        pp = params[:presearch_params]
        prec_estimation = pp["abreviate_precursor_calc"] ? FullPrecCapture() : PartialPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(first(params[:isotope_err_bounds])), UInt8(last(params[:isotope_err_bounds]))),
            Float32(pp["frag_tol_ppm"]),
            Float32(pp["frag_err_quantile"]),
            Int64(pp["min_samples"]),
            Float32(pp["max_qval"]),
            Int64(pp["max_presearch_iters"]),
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
            Int64(params[:irt_mapping_params]["n_sigma_tol"]),
            prec_estimation
        )
    end
end
function set_irt_to_rt_model!(ptsr::ParameterTuningSearchResults, search_context::SearchContext, params::P, ms_file_idx::Int64, model::Tuple{SplineRtConversionModel, Vector{Float32}, Vector{Float32}, Float32}) where {P<:ParameterTuningSearchParameters}
    ptsr.irt_to_rt_model[] = model[1]
    if length(ptsr.irt) > 0
        resize!(ptsr.irt, 0)
    end 
    if length(ptsr.rt) > 0
        resize!(ptsr.rt, 0)
    end 
    append!(ptsr.rt, model[2])
    append!(ptsr.irt, model[3])
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    getIrtErrors(search_context)[parsed_fname] =  model[4]*params.irt_tol_sd
end

get_parameters(search_type::ParameterTuningSearch, params::Any) = ParameterTuningSearchParameters(params)
function init_search_results(search_parameters::ParameterTuningSearchParameters, search_context::SearchContext, ms_file_idx::Int64)
    out_dir = getDataOutDir(search_context)
    qc_dir = joinpath(out_dir, "qc_plots")
    if !isdir(qc_dir)
        mkdir(qc_dir)
    end
    return ParameterTuningSearchResults(
        getMassErrorModel(search_context, ms_file_idx),
        Base.Ref{SplineRtConversionModel}(),
        Vector{Float32}(),
        Vector{Float32}(),
        Vector{Float32}(),
        qc_dir,
    )
end
function process_file!(
    results::ParameterTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::Arrow.Table) where {P<:ParameterTuningSearchParameters}

    try 
        setMassErrorModel!(
            search_context, ms_file_idx,
            MassErrorModel(zero(Float32), (getFragTolPpm(params), getFragTolPpm(params)))
        )
        psms = DataFrame()
        for i in range(1, getMaxPresearchIters(params))
            new_psms = library_search(spectra, search_context, params, ms_file_idx)
            iszero(size(new_psms, 1)) && continue
            #Add new psms that pass the scoring thresholds 
            add_columns_and_concat!(psms, new_psms, spectra, getPrecursors(getSpecLib(search_context)), params)
            try 
                filter_and_score_psms!(psms, params) >= getMinPsms(params) && break
            catch e
                throw(e)
                continue
            end
        end
        #Model conversion between library (irt) and empirical (rt) retention times
        set_irt_to_rt_model!(results, search_context, params, ms_file_idx, fit_irt_model(params, psms))

        fragments = get_matched_fragments(spectra, psms, search_context, params, ms_file_idx)

        set_mass_err_model!(results, fit_mass_err_model(params, fragments))

    catch e 
        throw(e)
    end
    return results
end

function process_search_results!(results::ParameterTuningSearchResults, params::P, search_context::SearchContext, ms_file_idx::Int64,
    ::Arrow.Table) where {P<:ParameterTuningSearchParameters}
    rt_alignment_folder = getRtAlignPlotFolder(search_context)
    mass_error_folder = getMassErrPlotFolder(search_context)
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    rt_plot_path = joinpath(rt_alignment_folder, parsed_fname*".pdf")
    mass_plot_path = joinpath(mass_error_folder, parsed_fname*".pdf")
    n = length(results.rt)
    p = Plots.plot(
        results.rt,
        results.irt,
        seriestype=:scatter, title = parsed_fname*"\n n = $n",                        
        xlabel = "Retention Time RT (min)",
        ylabel ="Indexed Retention Time iRT (min)",
        label = nothing,
        alpha = 0.1,
        size = 100*[13.3, 7.5]
    )
    pbins = LinRange(minimum(results.rt), maximum(results.rt), 100)
    Plots.plot!(pbins, getIrtToRtModel(results).(pbins), lw = 3, label = nothing)
    savefig(p, rt_plot_path)
    p = histogram(results.ppm_errs)
    savefig(p, mass_plot_path)
    setMassErrorModel!(search_context, ms_file_idx, getMassErrorModel(results))
    setRtIrtModel!(search_context, ms_file_idx, getIrtToRtModel(results))
end
#function set_irt_to_rt_model!(results::ParameterTuningSearchResults, model::RtConversionModel)#
#
#end
function fit_mass_err_model(search_parameters::P, fragments::Vector{FragmentMatch{Float32}}) where {P<:FragmentIndexSearchParameters}
    
    function _getPPM(a::T, b::T) where {T<:AbstractFloat}
        Float32((b - a)/(a/1e6))
    end
    frag_err_quantile = getFragErrQuantile(search_parameters)
    ppm_errs = [_getPPM(match.theoretical_mz, match.match_mz) for match in fragments];

    bins = LinRange(minimum(ppm_errs), maximum(ppm_errs), 100)
    mass_err = median(ppm_errs)
    ppm_errs = ppm_errs .- mass_err
    l_bound = quantile(ppm_errs, frag_err_quantile)
    r_bound = quantile(ppm_errs, 1 - frag_err_quantile)
    return MassErrorModel(
        Float32(mass_err),
        (Float32(abs(l_bound)), Float32(abs(r_bound)))
        ), ppm_errs
end

function get_matched_fragments(spectra::Arrow.Table, psms::DataFrame, search_context::SearchContext, search_parameters::P, ms_file_idx::Int64) where {P<:FragmentIndexSearchParameters}
    return vcat(massErrorSearch(
        spectra,
        psms[!,:scan_idx],
        psms[!,:precursor_idx],
        UInt32(ms_file_idx),
        getSpecLib(search_context),
        getSearchData(search_context),
        getMassErrorModel(search_context, ms_file_idx),
        search_parameters,
    )...)
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
function add_columns_and_concat!(psms::DataFrame, new_psms::DataFrame, spectra::Arrow.Table, precursors::Arrow.Table, params::P) where {P<:ParameterTuningSearchParameters}
    addPreSearchColumns!(new_psms, spectra, precursors[:is_decoy], precursors[:irt], precursors[:prec_charge], spectra[:retentionTime], spectra[:TIC])
    if new_psms === nothing
        return nothing
    else
        append!(psms, new_psms)
        return nothing
    end
end


function filter_and_score_psms!(psms::DataFrame, params::P) where {P<:ParameterTuningSearchParameters}
        scorePresearch!(psms)
        getQvalues!(psms[!,:prob], psms[!,:target], psms[!,:q_value])
        max_q_val, min_psms = getMaxQVal(params), getMinPsms(params)
        n_passing_psms = sum(psms[!,:q_value].<=max_q_val)
        if n_passing_psms >= min_psms
            filter!(row -> row.q_value::Float16 <= max_q_val, psms)
            psms[!,:best_psms] = zeros(Bool, size(psms, 1))
            grouped_psms = groupby(psms,:precursor_idx)
            for sub_psms in grouped_psms
                best_idx = argmax(sub_psms.prob::AbstractVector{Float32})
                sub_psms[best_idx,:best_psms] = true
            end
            filter!(x->x.best_psms::Bool, psms)
            return n_passing_psms 
        end
    return -1
end

# Core algorithm implementations using multiple dispatch
function fit_irt_model(params::P, psms::DataFrame) where {P<:ParameterTuningSearchParameters}
    # Initial fit
    rt_to_irt_map = UniformSpline(
        psms[!,:irt_predicted],
        psms[!,:rt],
        getSplineDegree(params),
        getSplineNKnots(params),
    )
    
    # Calculate residuals
    psms[!,:irt_observed] = rt_to_irt_map.(psms.rt::Vector{Float32})
    residuals::Vector{Float32} = psms[!,:irt_observed]::Vector{Float32} .- psms[!,:irt_predicted]::Vector{Float32}
    #By defaualt the `mad` function returns an adjusted median absolute deviation
    #that estimates the standard deviation.
    irt_mad = mad(residuals)::Float32
    
    # Remove extreme outliers and refit
    valid_psms = psms[abs.(residuals) .< (irt_mad * getOutlierThreshold(params)), :]
    return (SplineRtConversionModel(UniformSpline(
        valid_psms[!,:irt_predicted],
        valid_psms[!,:rt],
        getSplineDegree(params),
        getSplineNKnots(params),
    )), valid_psms[!,:rt], valid_psms[!,:irt_predicted], irt_mad)
end

function summarize_results!(results::ParameterTuningSearchResults, params::P, search_context::SearchContext) where {P<:ParameterTuningSearchParameters}
    @info "Merging QC plots..."
    
    # Get plot folders from search context
    rt_alignment_folder = getRtAlignPlotFolder(search_context)
    mass_error_folder = getMassErrPlotFolder(search_context)
    
    # Merge retention time alignment plots
    rt_plots = [joinpath(rt_alignment_folder, x) for x in readdir(rt_alignment_folder) if endswith(x, ".pdf")]
    if !isempty(rt_plots)
        merge_pdfs(rt_plots, 
                  joinpath(rt_alignment_folder, "rt_alignment_plots.pdf"), 
                  cleanup=true)
    end
    
    # Merge mass error estimation plots
    mass_plots = [joinpath(mass_error_folder, x) for x in readdir(mass_error_folder) if endswith(x, ".pdf")]
    if !isempty(mass_plots)
        merge_pdfs(mass_plots, 
                  joinpath(mass_error_folder, "mass_error_plots.pdf"), 
                  cleanup=true)
    end

    @info "QC plot merging complete"
    return nothing
end
