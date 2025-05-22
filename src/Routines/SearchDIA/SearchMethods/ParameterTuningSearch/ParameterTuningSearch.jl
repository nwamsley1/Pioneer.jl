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
    ms1_mass_err_model::Base.Ref{<:MassErrorModel}
    rt_to_irt_model::Base.Ref{<:RtConversionModel}
    irt::Vector{Float32}
    rt::Vector{Float32}
    ppm_errs::Vector{Float32}
    ms1_ppm_errs::Vector{Float32}
    qc_plots_folder_path::String
end

"""
Parameters for parameter tuning search.
Configures fragment matching, RT alignment, and general search behavior.
"""
struct ParameterTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # Core parameters from the original struct
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_fraction_transmitted::Float32
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
    relative_improvement_threshold::Float32
    spline_degree::Int64
    spline_n_knots::Int64
    spline_fit_outlier_sd::Int64
    irt_tol_sd::Int64
    prec_estimation::P

    function ParameterTuningSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        global_params = params.global_settings
        tuning_params = params.parameter_tuning
        frag_params = tuning_params.fragment_settings
        search_params = tuning_params.search_settings
        rt_params = params.rt_alignment
        
        # Convert isotope error bounds
        isotope_bounds = global_params.isotope_settings.err_bounds_first_pass
        
        # Create precursor estimation type
        prec_estimation = global_params.isotope_settings.partial_capture ? PartialPrecCapture() : FullPrecCapture()
        
        # Construct with appropriate type conversions
        new{typeof(prec_estimation)}(
            # Core parameters
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            Float32(global_params.isotope_settings.min_fraction_transmitted),
            Float32(frag_params.tol_ppm),
            Float32(search_params.frag_err_quantile),
            Int64(search_params.min_samples),
            Float32(global_params.scoring.q_value_threshold),
            Int64(search_params.max_presearch_iters),
            UInt8(frag_params.min_score),
            Int64(frag_params.min_count),
            Float32(frag_params.min_spectral_contrast),
            Float32(frag_params.min_log2_ratio),
            (Int64(first(frag_params.min_top_n)), Int64(last(frag_params.min_top_n))),
            UInt8(1), # max_best_rank default
            Int64(frag_params.n_isotopes),
            UInt8(frag_params.max_rank),
            Float32(search_params.sample_rate),
            typemax(Float32), # irt_tol default
            Set{Int64}([2]), # spec_order default
            Float32(frag_params.relative_improvement_threshold),
            3,  # spline_degree default
            5,  # spline_n_knots default
            5,  # spline_fit_outlier_sd default
            Int64(rt_params.sigma_tolerance),
            prec_estimation
        )
    end
end

#==========================================================
Results Access Methods
==========================================================#
getMassErrorModel(ptsr::ParameterTuningSearchResults) = ptsr.mass_err_model[]
getMs1MassErrorModel(ptsr::ParameterTuningSearchResults) = ptsr.ms1_mass_err_model[]
getRtToIrtModel(ptsr::ParameterTuningSearchResults) = ptsr.rt_to_irt_model[]
getQcPlotsFolder(ptsr::ParameterTuningSearchResults) = ptsr.qc_plots_folder_path

function set_rt_to_irt_model!(
    ptsr::ParameterTuningSearchResults, 
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64,
    model::Tuple{SplineRtConversionModel, Vector{Float32}, Vector{Float32}, Float32}
) where {P<:ParameterTuningSearchParameters}
    
    ptsr.rt_to_irt_model[] = model[1]
    resize!(ptsr.irt, 0)
    resize!(ptsr.rt, 0)
    append!(ptsr.rt, model[2])
    append!(ptsr.irt, model[3])
    
    #parsed_fname = getParsedFileName(search_context, ms_file_idx)
    getIrtErrors(search_context)[ms_file_idx] = model[4] * params.irt_tol_sd
end

#==========================================================
Interface Implementation
==========================================================#

get_parameters(::ParameterTuningSearch, params::PioneerParameters) = ParameterTuningSearchParameters(params)

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
    return ParameterTuningSearchResults(
        Base.Ref{MassErrorModel}(),
        Base.Ref{MassErrorModel}(),
        Ref{SplineRtConversionModel}(),
        Vector{Float32}(),
        Vector{Float32}(),
        Vector{Float32}(),
        Vector{Float32}(),
        qc_dir
    )
end

"""
Main file processing method for parameter tuning search.
"""
function process_file!(
    results::ParameterTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:ParameterTuningSearchParameters}

    """
    Collect PSMs through multiple iterations until sufficient high-quality PSMs found.
    """
    function collect_psms(
        spectra::MassSpecData,
        search_context::SearchContext,
        params::P,
        ms_file_idx::Int64
    ) where {P<:ParameterTuningSearchParameters}
        
        """
        Add columns and concatenate new PSMs to existing DataFrame.
        """
        function add_columns_and_concat!(
                psms::DataFrame,
                new_psms::DataFrame,
                spectra::MassSpecData,
                precursors::LibraryPrecursors,
                params::P
            ) where {P<:ParameterTuningSearchParameters}
                
                add_tuning_search_columns!(
                    new_psms,
                    spectra,
                    getIsDecoy(precursors),#[:is_decoy],
                    getIrt(precursors),#[:irt],
                    getCharge(precursors),#[:prec_charge],
                    getRetentionTimes(spectra),
                    getTICs(spectra)
                )
                
                if new_psms !== nothing
                    append!(psms, new_psms)
                end
        end

        psms = DataFrame()
        for i in 1:getMaxPresearchIters(params)
            new_psms = library_search(spectra, search_context, params, ms_file_idx)
            iszero(size(new_psms, 1)) && continue
            
            add_columns_and_concat!(psms, new_psms, spectra, 
                                getPrecursors(getSpecLib(search_context)), params)
            #Arrow.write("C:\\Users\\n.t.wamsley\\Desktop\\testpsms.arrow", psms)
            #[println(x) for x in names(psms)]
            try 
                filter_and_score_psms!(psms, params) >= getMinPsms(params) && break
            catch e
                throw(e)
            end
        end
        return psms
    end

    try
        mass_err_passing = false
        results.mass_err_model[] =  MassErrorModel(
            0.0f0, 
            (getFragTolPpm(params), getFragTolPpm(params)))
        ppm_errs = nothing
        ms1_ppm_errs = nothing
        prev_mass_err = 0.0f0
        for i in range(0, 1)
            setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
            setQuadTransmissionModel!(search_context, ms_file_idx, GeneralGaussModel(5.0f0, 0.0f0))

            # Collect PSMs through iterations
            psms = collect_psms(spectra, search_context, params, ms_file_idx)
            
            # Fit RT alignment model
            set_rt_to_irt_model!(results, search_context, params, ms_file_idx, 
                                fit_irt_model(params, psms))

            # Get fragments and fit mass error model
            fragments = get_matched_fragments(spectra, psms, results,search_context, params, ms_file_idx)
            mass_err_model, ppm_errs = fit_mass_err_model(params, fragments)
            mass_errs = get_matched_precursors(spectra, psms, results,search_context, params, ms_file_idx)
            ms1_mass_err_model, ms1_ppm_errs = mass_err_ms1(mass_errs, params)
            #If the mass offset is much larger than expected, need to research with an updated mass offset estimate 
            if abs(getMassOffset(mass_err_model))>(getFragTolPpm(params)/4)
                prev_mass_err = getMassOffset(mass_err_model)
                results.mass_err_model[] = MassErrorModel(getMassOffset(mass_err_model), (getFragTolPpm(params), getFragTolPpm(params)))
                results.ms1_mass_err_model[] = ms1_mass_err_model
            else
                results.mass_err_model[] = MassErrorModel(getMassOffset(mass_err_model) + prev_mass_err, (getLeftTol(mass_err_model), getRightTol(mass_err_model)))
                results.ms1_mass_err_model[] = ms1_mass_err_model
                break
            end
        end

        append!(results.ppm_errs, ppm_errs)
        append!(results.ms1_ppm_errs, ms1_ppm_errs)
    catch e
        setFailedIndicator!(getMSData(search_context), ms_file_idx, true)
        @warn "Could not tune parameters for $ms_file_idx"
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
    ms_file_idx::Int64,
    ::MassSpecData
) where {P<:ParameterTuningSearchParameters}
    try
    rt_alignment_folder = getRtAlignPlotFolder(search_context)
    mass_error_folder = getMassErrPlotFolder(search_context)
    ms1_mass_error_folder = getMs1MassErrPlotFolder(search_context)
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    
    # Generate and save RT alignment plot
    rt_plot_path = joinpath(rt_alignment_folder, parsed_fname*".pdf")
    generate_rt_plot(results, rt_plot_path, parsed_fname)
    
    # Generate and save mass error plot
    mass_plot_path = joinpath(mass_error_folder, parsed_fname*".pdf")
    generate_mass_error_plot(results, parsed_fname, mass_plot_path)
    
    # Generate and save mass error plot
    mass_plot_path = joinpath(ms1_mass_error_folder, parsed_fname*".pdf")
    generate_ms1_mass_error_plot(results, parsed_fname, mass_plot_path)
    
    # Update models in search context
    setMassErrorModel!(search_context, ms_file_idx, getMassErrorModel(results))

    # Update models in search context
    setMs1MassErrorModel!(search_context, ms_file_idx, getMs1MassErrorModel(results))
    
    setRtIrtMap!(search_context, getRtToIrtModel(results), ms_file_idx)
    catch
        setFailedIndicator!(getMSData(search_context), ms_file_idx, true)
        nothing
    end
end

function reset_results!(ptsr::ParameterTuningSearchResults)
    resize!(ptsr.irt, 0)
    resize!(ptsr.rt, 0)
    resize!(ptsr.ppm_errs, 0)
    resize!(ptsr.ms1_ppm_errs, 0)
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
    output_path = joinpath(rt_alignment_folder, "rt_alignment_plots.pdf")
    try
        if isfile(output_path)
            rm(output_path)
        end
    catch e
        @warn "Could not clear existing file: $e"
    end
    rt_plots = [joinpath(rt_alignment_folder, x) for x in readdir(rt_alignment_folder) 
    if endswith(x, ".pdf")]
    
    if !isempty(rt_plots)
        merge_pdfs(rt_plots, 
                    output_path, 
                  cleanup=true)
    end
    
    # Merge mass error plots
    mass_error_folder = getMassErrPlotFolder(search_context)
    output_path = joinpath(mass_error_folder, "mass_error_plots.pdf")
    try
        if isfile(output_path)
            rm(output_path)
        end
    catch e
        @warn "Could not clear existing file: $e"
    end
    mass_plots = [joinpath(mass_error_folder, x) for x in readdir(mass_error_folder) 
                    if endswith(x, ".pdf")]

    if !isempty(mass_plots)
        merge_pdfs(mass_plots, 
                  output_path, 
                  cleanup=true)
    end
    
    # Merge mass error plots
    ms1_mass_error_folder = getMs1MassErrPlotFolder(search_context)
    output_path = joinpath(ms1_mass_error_folder, "ms1_mass_error_plots.pdf")
    try
        if isfile(output_path)
            rm(output_path)
        end
    catch e
        @warn "Could not clear existing file: $e"
    end
    mass_plots = [joinpath(ms1_mass_error_folder, x) for x in readdir(ms1_mass_error_folder) 
                    if endswith(x, ".pdf")]

    if !isempty(mass_plots)
        merge_pdfs(mass_plots, 
                    output_path, 
                    cleanup=true)
    end

    @info "QC plot merging complete"
end


