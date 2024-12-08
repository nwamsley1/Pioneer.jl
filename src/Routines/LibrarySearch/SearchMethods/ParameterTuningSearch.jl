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
    :isotope_err_bounds => (0, 2),
    :presearch_params => Dict(
        "frag_tol_ppm" => 30.0,
        "frag_err_quantile" => 0.01,
        "min_samples" => 1000,
        "max_qval" => 0.01,
        "max_presearch_iters" => 10,
        "min_index_search_score" => 3,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "n_frag_isotopes" => 2,
        "max_frag_rank" => 10,
        "sample_rate" => 0.1,
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
struct ParameterTuningSearch <: TuningMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for parameter tuning search.
Holds mass error models, RT alignment models, and associated data.
"""
struct ParameterTuningSearchResults <: SearchResults 
    mass_err_model::Base.Ref{<:MassErrorModel}
    irt_to_rt_model::Base.Ref{<:RtConversionModel}
    irt::Vector{Float32}
    rt::Vector{Float32}
    ppm_errs::Vector{Float32}
    qc_plots_folder_path::String
end

"""
Parameters for parameter tuning search.
Configures fragment matching, RT alignment, and general search behavior.
"""
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
            3,  # Spline degree
            5,  # Number of knots
            10, # Outlier threshold
            Int64(params[:irt_mapping_params]["n_sigma_tol"]),
            prec_estimation
        )
    end
end

#==========================================================
Interface Implementation
==========================================================#

get_parameters(::ParameterTuningSearch, params::Any) = ParameterTuningSearchParameters(params)

function init_search_results(::ParameterTuningSearchParameters, search_context::SearchContext)
    out_dir = getDataOutDir(search_context)
    qc_dir = joinpath(out_dir, "qc_plots")
    !isdir(qc_dir) && mkdir(qc_dir)
    
    return ParameterTuningSearchResults(
        Base.Ref{MassErrorModel}(),
        Ref{SplineRtConversionModel}(),
        Vector{Float32}(),
        Vector{Float32}(),
        Vector{Float32}(),
        qc_dir
    )
end

#==========================================================
Results Access Methods
==========================================================#

getMassErrorModel(ptsr::ParameterTuningSearchResults) = ptsr.mass_err_model[]
getIrtToRtModel(ptsr::ParameterTuningSearchResults) = ptsr.irt_to_rt_model[]
getQcPlotsFolder(ptsr::ParameterTuningSearchResults) = ptsr.qc_plots_folder_path

function set_mass_err_model!(ptsr::ParameterTuningSearchResults, model::Tuple{MassErrorModel, Vector{Float32}})
    ptsr.mass_err_model[] = model[1]
    append!(ptsr.ppm_errs, model[2])
end

function set_irt_to_rt_model!(
    ptsr::ParameterTuningSearchResults, 
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64,
    model::Tuple{SplineRtConversionModel, Vector{Float32}, Vector{Float32}, Float32}
) where {P<:ParameterTuningSearchParameters}
    
    ptsr.irt_to_rt_model[] = model[1]
    resize!(ptsr.irt, 0)
    resize!(ptsr.rt, 0)
    append!(ptsr.rt, model[2])
    append!(ptsr.irt, model[3])
    
    #parsed_fname = getParsedFileName(search_context, ms_file_idx)
    getIrtErrors(search_context)[ms_file_idx] = model[4] * params.irt_tol_sd
end

function reset_results!(ptsr::ParameterTuningSearchResults)
    resize!(ptsr.irt, 0)
    resize!(ptsr.rt, 0)
    resize!(ptsr.ppm_errs, 0)
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Main file processing method for parameter tuning search.
"""
function process_file!(
    results::ParameterTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::Arrow.Table
) where {P<:ParameterTuningSearchParameters}

    try
        results.mass_err_model[] =  MassErrorModel(zero(Float32), (getFragTolPpm(params), getFragTolPpm(params)))
        setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
        setQuadTransmissionModel!(search_context, ms_file_idx, GeneralGaussModel(5.0f0, 0.0f0))

        # Collect PSMs through iterations
        psms = collect_psms(spectra, search_context, params, ms_file_idx)
        
        # Fit RT alignment model
        set_irt_to_rt_model!(results, search_context, params, ms_file_idx, 
                            fit_irt_model(params, psms))

        # Get fragments and fit mass error model
        fragments = get_matched_fragments(spectra, psms, results,search_context, params, ms_file_idx)
        mass_err_model, ppm_errs = fit_mass_err_model(params, fragments)
        results.mass_err_model[] = mass_err_model
        append!(results.ppm_errs, ppm_errs)
    catch e 
        throw(e)
    end
    
    return results
end

"""
Process search results and generate QC plots for a single MS file.
"""
function process_search_results!(
    results::ParameterTuningSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64
) where {P<:ParameterTuningSearchParameters}
    
    rt_alignment_folder = getRtAlignPlotFolder(search_context)
    mass_error_folder = getMassErrPlotFolder(search_context)
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    
    # Generate and save RT alignment plot
    rt_plot_path = joinpath(rt_alignment_folder, parsed_fname*".pdf")
    generate_rt_plot(results, rt_plot_path, parsed_fname)
    
    # Generate and save mass error plot
    mass_plot_path = joinpath(mass_error_folder, parsed_fname*".pdf")
    generate_mass_error_plot(results, mass_plot_path)
    
    # Update models in search context
    setMassErrorModel!(search_context, ms_file_idx, getMassErrorModel(results))
    setRtIrtModel!(search_context, ms_file_idx, getIrtToRtModel(results))
end

"""
Summarize results across all files and merge QC plots.
"""
function summarize_results!(
    results::ParameterTuningSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:ParameterTuningSearchParameters}
    
    @info "Merging QC plots..."
    
    # Merge RT alignment plots
    rt_alignment_folder = getRtAlignPlotFolder(search_context)
    rt_plots = [joinpath(rt_alignment_folder, x) for x in readdir(rt_alignment_folder) 
                if endswith(x, ".pdf")]
    if !isempty(rt_plots)
        merge_pdfs(rt_plots, 
                  joinpath(rt_alignment_folder, "rt_alignment_plots.pdf"), 
                  cleanup=true)
    end
    
    # Merge mass error plots
    mass_error_folder = getMassErrPlotFolder(search_context)
    mass_plots = [joinpath(mass_error_folder, x) for x in readdir(mass_error_folder) 
                 if endswith(x, ".pdf")]
    if !isempty(mass_plots)
        merge_pdfs(mass_plots, 
                  joinpath(mass_error_folder, "mass_error_plots.pdf"), 
                  cleanup=true)
    end
    
    @info "QC plot merging complete"
end

#==========================================================
PSM Collection and Processing
==========================================================#

"""
Collect PSMs through multiple iterations until sufficient high-quality PSMs found.
"""
function collect_psms(
    spectra::Arrow.Table,
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64
) where {P<:ParameterTuningSearchParameters}
    
    psms = DataFrame()
    for i in 1:getMaxPresearchIters(params)
        new_psms = library_search(spectra, search_context, params, ms_file_idx)
        iszero(size(new_psms, 1)) && continue
        
        add_columns_and_concat!(psms, new_psms, spectra, 
                              getPrecursors(getSpecLib(search_context)), params)
        
        try 
            filter_and_score_psms!(psms, params) >= getMinPsms(params) && break
        catch e
            throw(e)
        end
    end
    return psms
end

"""
Add columns and concatenate new PSMs to existing DataFrame.
"""
function add_columns_and_concat!(
    psms::DataFrame,
    new_psms::DataFrame,
    spectra::Arrow.Table,
    precursors::Arrow.Table,
    params::P
) where {P<:ParameterTuningSearchParameters}
    
    addPreSearchColumns!(
        new_psms,
        spectra,
        precursors[:is_decoy],
        precursors[:irt],
        precursors[:prec_charge],
        spectra[:retentionTime],
        spectra[:TIC]
    )
    
    if new_psms !== nothing
        append!(psms, new_psms)
    end
end

"""
Filter PSMs based on score and select best matches.
"""
function filter_and_score_psms!(
    psms::DataFrame,
    params::P
) where {P<:ParameterTuningSearchParameters}
    
    scorePresearch!(psms)
    getQvalues!(psms[!,:prob], psms[!,:target], psms[!,:q_value])
    
    max_q_val, min_psms = getMaxQVal(params), getMinPsms(params)
    n_passing_psms = sum(psms[!,:q_value] .<= max_q_val)
    
    if n_passing_psms >= min_psms
        filter!(row -> row.q_value::Float16 <= max_q_val, psms)
        psms[!,:best_psms] = zeros(Bool, size(psms, 1))
        
        # Select best PSM per precursor
        for sub_psms in groupby(psms, :precursor_idx)
            best_idx = argmax(sub_psms.prob::AbstractVector{Float32})
            sub_psms[best_idx,:best_psms] = true
        end
        
        filter!(x -> x.best_psms::Bool, psms)
        return n_passing_psms
    end
    
    return -1
end

#==========================================================
Model Fitting Methods
==========================================================#

"""
Fit retention time alignment model between library and empirical RTs.
"""
function fit_irt_model(
    params::P,
    psms::DataFrame
) where {P<:ParameterTuningSearchParameters}
    
    # Initial spline fit
    rt_to_irt_map = UniformSpline(
        psms[!,:irt_predicted],
        psms[!,:rt],
        getSplineDegree(params),
        getSplineNKnots(params)
    )
    
    # Calculate residuals
    psms[!,:irt_observed] = rt_to_irt_map.(psms.rt::Vector{Float32})
    residuals = psms[!,:irt_observed] .- psms[!,:irt_predicted]
    irt_mad = mad(residuals)::Float32
    
    # Remove outliers and refit
    valid_psms = psms[abs.(residuals) .< (irt_mad * getOutlierThreshold(params)), :]
    
    final_model = SplineRtConversionModel(UniformSpline(
        valid_psms[!,:irt_predicted],
        valid_psms[!,:rt],
        getSplineDegree(params),
        getSplineNKnots(params)
    ))
    
    return (final_model, valid_psms[!,:rt], valid_psms[!,:irt_predicted], irt_mad)
end

"""
Fit mass error model from fragment matches.
"""
function fit_mass_err_model(
    params::P,
    fragments::Vector{FragmentMatch{Float32}}
) where {P<:FragmentIndexSearchParameters}
    
    # Calculate PPM errors
    ppm_errs = [calc_ppm_error(match.theoretical_mz, match.match_mz) for match in fragments]
    mass_err = median(ppm_errs)
    ppm_errs .-= mass_err
    
    # Calculate error bounds
    frag_err_quantile = getFragErrQuantile(params)
    l_bound = quantile(ppm_errs, frag_err_quantile)
    r_bound = quantile(ppm_errs, 1 - frag_err_quantile)
    
    return MassErrorModel(
        Float32(mass_err),
        (Float32(abs(l_bound)), Float32(abs(r_bound)))
    ), ppm_errs
end

"""
Get matched fragments for mass error estimation.
"""
function get_matched_fragments(
    spectra::Arrow.Table,
    psms::DataFrame,
    results::ParameterTuningSearchResults,
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64
) where {P<:FragmentIndexSearchParameters}
    
    return vcat(massErrorSearch(
        spectra,
        psms[!,:scan_idx],
        psms[!,:precursor_idx],
        UInt32(ms_file_idx),
        getSpecLib(search_context),
        getSearchData(search_context),
        getMassErrorModel(results),
        params
    )...)
end

#==========================================================
Plotting Helpers
==========================================================#

"""
Generate retention time alignment plot.
"""
function generate_rt_plot(
    results::ParameterTuningSearchResults,
    plot_path::String,
    title::String
)
    n = length(results.rt)
    p = Plots.plot(
        results.rt,
        results.irt,
        seriestype=:scatter,
        title = title*"\n n = $n",
        xlabel = "Retention Time RT (min)",
        ylabel = "Indexed Retention Time iRT (min)",
        label = nothing,
        alpha = 0.1,
        size = 100*[13.3, 7.5]
    )
    
    pbins = LinRange(minimum(results.rt), maximum(results.rt), 100)
    Plots.plot!(pbins, getIrtToRtModel(results).(pbins), lw=3, label=nothing)
    savefig(p, plot_path)
end

"""
Generate mass error distribution plot.
"""
function generate_mass_error_plot(
    results::ParameterTuningSearchResults,
    plot_path::String
)
    p = histogram(results.ppm_errs)
    savefig(p, plot_path)
end

#==========================================================
Utility Functions
==========================================================#

"""
Calculate PPM error between theoretical and observed mass.
"""
function calc_ppm_error(theoretical::T, observed::T) where T<:AbstractFloat
    return Float32((observed - theoretical)/(theoretical/1e6))
end