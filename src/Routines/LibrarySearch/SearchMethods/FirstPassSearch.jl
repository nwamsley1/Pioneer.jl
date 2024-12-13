"""
    FirstPassSearch

Initial search method to identify PSMs and establish retention time calibration.

This search:
1. Performs initial PSM identification with learned scoring
2. Calculates retention time indices and FWHM statistics
3. Maps between library and empirical retention times
4. Generates RT calibration curves for subsequent searches

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :isotope_err_bounds => (0, 2),
    :first_search_params => Dict(
        "n_train_rounds_probit" => 10,
        "max_iter_probit" => 100,
        "max_q_value_probit_rescore" => 0.01,
        "min_index_search_score" => 3,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "max_precursors_passing" => 5000,
        "abreviate_precursor_calc" => false
    ),
    :summarize_first_search_params => Dict(
        "min_inference_points" => 1000,
        "max_q_val_for_irt" => 0.01,
        "max_precursors" => 10000,
        "max_irt_bin_size" => 0.1,
        "max_prob_to_impute" => 0.99
    ),
    :irt_mapping_params => Dict(
        "min_prob" => 0.9
    )
)

# Execute search
results = execute_search(FirstPassSearch(), search_context, params)
```
"""
struct FirstPassSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for first pass search.
Holds FWHM statistics, PSM file paths, and model updates.
"""
struct FirstPassSearchResults <: SearchResults
    fwhms::Dictionary{Int64, @NamedTuple{median_fwhm::Float32,mad_fwhm::Float32}}
    psms::Base.Ref{DataFrame}
end

"""
Parameters for first pass search.
Configures PSM identification, scoring, and RT calibration.
"""
struct FirstPassSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # Search parameters
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
    
    # Scoring parameters
    n_train_rounds_probit::Int64
    max_iter_probit::Int64
    max_q_value_probit_rescore::Float32
    
    # RT parameters
    max_precursors_passing::Int64
    min_inference_points::Int64
    max_q_val_for_irt::Float32
    min_prob_for_irt_mapping::Float32
    max_precursors::Int64
    max_irt_bin_size::Float32
    max_prob_to_impute::Float32
    fwhm_nstd::Float32
    irt_nstd::Float32
    prec_estimation::P

    function FirstPassSearchParameters(params::Any)
        fp = params[:first_search_params]
        sp = params[:summarize_first_search_params]
        prec_estimation = fp["abreviate_precursor_calc"] ? FullPrecCapture() : PartialPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(first(params[:isotope_err_bounds])), UInt8(last(params[:isotope_err_bounds]))),
            0.0f0,  # No fragment tolerance for first pass
            UInt8(fp["min_index_search_score"]),
            Int64(fp["min_frag_count"]),
            Float32(fp["min_spectral_contrast"]),
            Float32(fp["min_log2_matched_ratio"]),
            (Int64(first(fp["min_topn_of_m"])), Int64(last(fp["min_topn_of_m"]))),
            UInt8(fp["max_best_rank"]),
            Int64(fp["n_frag_isotopes"]),
            UInt8(fp["max_frag_rank"]),
            1.0f0,  # Full sampling for first pass
            Set(2),
            Int64(fp["n_train_rounds_probit"]),
            Int64(fp["max_iter_probit"]),
            Float32(fp["max_q_value_probit_rescore"]),
            Int64(fp["max_precursors_passing"]),
            Int64(sp["min_inference_points"]),
            Float32(sp["max_q_val_for_irt"]),
            Float32(params[:irt_mapping_params]["min_prob"]),
            Int64(sp["max_precursors"]),
            Float32(sp["max_irt_bin_size"]),
            Float32(sp["max_prob_to_impute"]),
            Float32(sp["fwhm_nstd"]),
            Float32(sp["irt_nstd"]),
            prec_estimation
        )
    end
end

#==========================================================
Interface Implementation
==========================================================#

get_parameters(::FirstPassSearch, params::Any) = FirstPassSearchParameters(params)

function init_search_results(
    ::FirstPassSearchParameters,
    search_context::SearchContext
)
    temp_folder = joinpath(getDataOutDir(search_context), "first_pass_psms")
    !isdir(temp_folder) && mkdir(temp_folder)
    
    return FirstPassSearchResults(
        Dictionary{Int64, NamedTuple{(:median_fwhm, :mad_fwhm), Tuple{Float32, Float32}}}(),
        Base.Ref{DataFrame}()
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single MS file in the first pass search.
"""
function process_file!(
    results::FirstPassSearchResults,
    params::P, 
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::Arrow.Table
) where {P<:FirstPassSearchParameters}

    try
        # Get models and update fragment lookup table
        psms = perform_library_search(spectra, search_context, params, ms_file_idx)
        results.psms[] = process_psms!(psms, spectra, search_context, params, ms_file_idx)
    catch e
        @warn "First pass search failed" ms_file_idx exception=e
        rethrow(e)
    end

    return results
end

"""
Initial file processing complete, no additional processing needed.
"""
function process_search_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    ::Arrow.Table
) where {P<:FirstPassSearchParameters}
    psms = results.psms[]
    fwhms = skipmissing(psms[!, :fwhm])
    fwhm_points = count(!ismissing, fwhms)
    if fwhm_points >= 1#params.min_inference_points
        insert!(results.fwhms, ms_file_idx, (
            median_fwhm = median(fwhms),
            mad_fwhm = mad(fwhms, normalize=true)))
    else
        @warn "Insuficient fwhm_points to estimate for $ms_file_idx"
        insert!(results.fwhms, ms_file_idx, (
            median_fwhm = 0.2f0,
            mad_fwhm = 0.2f0))
    end
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    temp_path = joinpath(getDataOutDir(search_context), "first_pass_psms", parsed_fname * ".arrow")
    psms[!, :ms_file_idx] .= UInt32(ms_file_idx)
    Arrow.write(
        temp_path,
        select!(psms, [:ms_file_idx, :scan_idx, :precursor_idx, :rt,
            :irt_predicted, :q_value, :score, :prob, :scan_count])
    )
    setFirstPassPsms!(getMSData(search_context), ms_file_idx, temp_path)
end

"""
Summarize results across all files.
"""
function summarize_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:FirstPassSearchParameters}
    
    #isempty(results.psms_paths) && return nothing
    @info "Summarizing first pass search results..."
    # Map retention times
    map_retention_times!(search_context, results, params)
    # Process precursors
    precursor_dict = get_best_precursors_accross_runs!(search_context, results, params)
    setPrecursorDict!(search_context, precursor_dict)
    # Calculate RT indices
    create_rt_indices!(search_context, results, precursor_dict, params)
    
    @info "Search results summarization complete"
end


"""
No cleanup needed between files.
"""
function reset_results!(results::FirstPassSearchResults)
    empty!(results.psms[])
    return nothing
end

#==========================================================
Helper Methods
==========================================================#

"""
Perform library search with current parameters.
"""
function perform_library_search(
    spectra::Arrow.Table,
    search_context::SearchContext,
    params::FirstPassSearchParameters,
    ms_file_idx::Int64
)
    setNceModel!(
        getFragmentLookupTable(getSpecLib(search_context)), 
        getNceModelModel(search_context, ms_file_idx)
    )
    
    return library_search(spectra, search_context, params, ms_file_idx)
end

"""
Process PSMs from library search.
"""
function process_psms!(
    psms::DataFrame,
    spectra::Arrow.Table,
    search_context::SearchContext,
    params::FirstPassSearchParameters,
    ms_file_idx::Int64
)
    rt_model = getRtIrtModel(search_context, ms_file_idx)
    
    # Add columns
    add_psm_columns!(psms, spectra, search_context, rt_model, ms_file_idx)
    
    # Score PSMs
    score_psms!(psms, params)
    
    # Get best PSMs
    select_best_psms!(
        psms,
        getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
        params
    )
    
    return psms
end

"""
Add necessary columns to PSM DataFrame.
"""
function add_psm_columns!(
    psms::DataFrame,
    spectra::Arrow.Table,
    search_context::SearchContext,
    rt_model::RtConversionModel,
    ms_file_idx::Int64
)
    addMainSearchColumns!(
        psms,
        getModel(rt_model),
        getStructuralMods(getPrecursors(getSpecLib(search_context))),#,
        getMissedCleavages(getPrecursors(getSpecLib(search_context))),#[:missed_cleavages],
        getIsDecoy(getPrecursors(getSpecLib(search_context))),#[:is_decoy],
        getIrt(getPrecursors(getSpecLib(search_context))),#[:irt],
        getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge],
        spectra[:retentionTime],
        spectra[:TIC],
        spectra[:mz_array]
    )
    
    # Calculate RT values
    psms[!, :irt_observed] = rt_model.(psms[!, :rt])
    psms[!, :irt_error] = Float16.(abs.(psms[!, :irt_observed] .- psms[!, :irt_predicted]))
    psms[!, :charge2] = UInt8.(psms[!, :charge] .== 2)
    psms[!, :ms_file_idx] .= UInt32(ms_file_idx)
end

"""
Score PSMs using probit model.
"""
function score_psms!(
    psms::DataFrame,
    params::FirstPassSearchParameters
)
    column_names = [
        :spectral_contrast, :city_block, :entropy_score, :scribe,
        :charge2, :poisson, :irt_error, :missed_cleavage, :Mox,
        :charge, :TIC, :y_count, :err_norm, :spectrum_peak_count, :intercept
    ]
    
    # Select scoring columns
    select!(psms, vcat(column_names, [:ms_file_idx, :score, :precursor_idx, :scan_idx,
        :q_value, :log2_summed_intensity, :irt, :rt, :irt_predicted, :target]))

    # Score PSMs
    scoreMainSearchpsms!(
        psms,
        column_names,
        n_train_rounds=params.n_train_rounds_probit,
        max_iter_per_round=params.max_iter_probit,
        max_q_value=Float64(params.max_q_value_probit_rescore)
    )

    # Process scores
    select!(psms, [:ms_file_idx, :score, :precursor_idx, :scan_idx,
        :q_value, :log2_summed_intensity, :irt, :rt, :irt_predicted])
    getProbs!(psms)
    sort!(psms, :irt)
end

"""
Select best PSMs based on criteria.
"""
function select_best_psms!(
    psms::DataFrame,
    precursor_mzs::AbstractVector,
    params::FirstPassSearchParameters
)
    getBestPSMs!(
        psms,
        precursor_mzs,
        max_q_val=params.max_q_val_for_irt,
        max_psms=params.max_precursors_passing
    )
end

"""
Map retention times between library and empirical scales.
"""
function map_retention_times!(
    search_context::SearchContext,
    results::FirstPassSearchResults,
    params::FirstPassSearchParameters
)
    @info "Mapping library to empirical retention times..."
    
    for (ms_file_idx, psms_path) in enumerate(getFirstPassPsms(getMSData(search_context)))
        psms = Arrow.Table(psms_path)
        best_hits = psms[:prob].>params.min_prob_for_irt_mapping#Map rts using only the best psms
        if sum(best_hits) > 100
            best_rts = psms[:rt][best_hits]
            best_irts = psms[:irt_predicted][best_hits]
            irt_to_rt_spline = UniformSpline(
                                        best_rts,
                                        best_irts,
                                        3, 
                                        5
            )
            rt_to_irt_spline = UniformSpline(
                best_irts,
                best_rts,
                3, 
                5
            )
            #Build rt=>irt and irt=> rt mappings for the file and add to the dictionaries 
            setRtIrtMap!(search_context, SplineRtConversionModel(rt_to_irt_spline), ms_file_idx)
            setIrtRtMap!(search_context, SplineRtConversionModel(rt_to_irt_spline), ms_file_idx)
        else
            throw("add a default option here...")
            #sensible default here?
            continue
        end
    end
    return nothing
end

"""
Process precursors and calculate iRT errors.
"""
function get_best_precursors_accross_runs!(
    search_context::SearchContext,
    results::FirstPassSearchResults,
    params::FirstPassSearchParameters
)
    @info "Finding best precursors across runs..."
    
    # Get best precursors
    return get_best_precursors_accross_runs(
        getFirstPassPsms(getMSData(search_context)),
        getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
        getRtIrtMap(search_context),
        max_q_val=params.max_q_val_for_irt,
        max_precursors=params.max_precursors
    )
end

PrecToIrtType = Dictionary{UInt32, 
    NamedTuple{
        (:best_prob, :best_ms_file_idx, :best_scan_idx, :best_irt, :mean_irt, :var_irt, :n, :mz), 
        Tuple{Float32, UInt32, UInt32, Float32, Union{Missing, Float32}, Union{Missing, Float32}, Union{Missing, UInt16}, Float32}
    }
}

"""
Create retention time indices.
"""
function create_rt_indices!(
    search_context::SearchContext,
    results::FirstPassSearchResults,
    precursor_dict::PrecToIrtType,
    params::FirstPassSearchParameters
)


    # Calculate iRT errors
    @info "Calculating iRT errors..."
    irt_errs = get_irt_errs(results.fwhms, precursor_dict, params)


    setIrtErrors!(search_context, irt_errs)

    @info "Creating RT indices..."
    # Create precursor to iRT mapping
    prec_to_irt = map(x -> (irt=x[:best_irt], mz=x[:mz]), 
                      precursor_dict)

    # Set up indices folder
    rt_indices_folder = joinpath(getDataOutDir(search_context), "rt_indices")
    !isdir(rt_indices_folder) && mkdir(rt_indices_folder)

    # Make RT indices
    rt_index_paths = makeRTIndices(
        rt_indices_folder,
        getFirstPassPsms(getMSData(search_context)),
        prec_to_irt,
        getRtIrtMap(search_context),
        min_prob=params.max_prob_to_impute
    )
    for (ms_file_idx, rt_index_path) in pairs(rt_index_paths)
        setRtIndex!(getMSData(search_context), ms_file_idx, rt_index_path)
    end
end

function get_irt_errs(
    fwhms::Dictionary{Int64, 
                        @NamedTuple{
                            median_fwhm::Float32,
                            mad_fwhm::Float32
                        }},
    prec_to_irt::Dictionary{UInt32, 
    @NamedTuple{best_prob::Float32, 
                best_ms_file_idx::UInt32, 
                best_scan_idx::UInt32, 
                best_irt::Float32, 
                mean_irt::Union{Missing, Float32}, 
                var_irt::Union{Missing, Float32}, 
                n::Union{Missing, UInt16}, 
                mz::Float32}}
    ,
    params::FirstPassSearchParameters
)
    #Get upper bound on peak fwhm. Use median + n*standard_deviation
    #estimate standard deviation by the median absolute deviation. 
    #n is a user-defined paramter. 
    fwhms = map(x->x[:median_fwhm] + params.fwhm_nstd*x[:mad_fwhm],
    fwhms)

    #Get variance in irt of apex accross runs. Only consider precursor identified below q-value threshold
    #in more than two runs .
    irt_std = median(
                skipmissing(map(x-> (x[:n] > 2) ? sqrt(x[:var_irt]/(x[:n] - 1)) : missing, prec_to_irt))
                )
    #Number of standard deviations to cover 
    irt_std *= params.irt_nstd
    #dictionary maping file name to irt tolerance. 
    return map(x->Float32(x+irt_std)::Float32, fwhms)::Dictionary{Int64, Float32}
end