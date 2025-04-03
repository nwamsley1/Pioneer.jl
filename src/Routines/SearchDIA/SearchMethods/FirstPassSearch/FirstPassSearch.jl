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
    # Core parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_fraction_transmitted::Float32
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
    match_between_runs::Bool
    
    # Scoring parameters
    n_train_rounds_probit::Int64
    max_iter_probit::Int64
    max_q_value_probit_rescore::Float32
    
    # RT parameters
    min_inference_points::Int64
    max_q_val_for_irt::Float32
    min_prob_for_irt_mapping::Float32
    max_irt_bin_size::Float32
    max_prob_to_impute::Float32
    fwhm_nstd::Float32
    irt_nstd::Float32
    prec_estimation::P

    function FirstPassSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        global_params = params.global_settings
        first_params = params.first_search
        frag_params = first_params.fragment_settings
        score_params = first_params.scoring_settings
        rt_params = params.rt_alignment
        irt_mapping_params = first_params.irt_mapping
        # Convert isotope error bounds
        isotope_bounds = global_params.isotope_settings.err_bounds_first_pass
        # Determine precursor estimation strategy
        prec_estimation = global_params.isotope_settings.partial_capture ? PartialPrecCapture() : FullPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            0.0f0,  # No transmission threshold for first pass
            0.0f0,  # No fragment tolerance for first pass
            UInt8(frag_params.min_score),
            Int64(frag_params.min_count),
            Float32(frag_params.min_spectral_contrast),
            Float32(frag_params.min_log2_ratio),
            (Int64(first(frag_params.min_top_n)), Int64(last(frag_params.min_top_n))),
            UInt8(1), # max_best_rank
            Int64(frag_params.n_isotopes),
            UInt8(frag_params.max_rank),
            1.0f0,  # Full sampling for first pass
            Set{Int64}([2]),
            global_params.match_between_runs,
            
            Int64(score_params.n_train_rounds),
            Int64(score_params.max_iterations),
            Float32(score_params.max_q_value),
            
            Int64(1000), # Default min_inference_points
            Float32(rt_params.min_probability),
            Float32(rt_params.min_probability),
            Float32(0.1), # Default max_irt_bin_size
            Float32(irt_mapping_params.max_prob_to_impute_irt),  # Default max_prob_to_impute
            Float32(irt_mapping_params.fwhm_nstd),   # Default fwhm_nstd
            Float32(irt_mapping_params.irt_nstd),   # Default irt_nstd
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
    temp_folder = joinpath(getDataOutDir(search_context), "temp_data", "first_pass_psms")
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
    spectra::MassSpecData
) where {P<:FirstPassSearchParameters}

    """
    Perform library search with current parameters.
    """
    function perform_library_search(
        spectra::MassSpecData,
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
        spectra::MassSpecData,
        search_context::SearchContext,
        params::FirstPassSearchParameters,
        ms_file_idx::Int64
    )

        """
        Select best PSMs based on criteria.
        """
        function select_best_psms!(
            psms::DataFrame,
            precursor_mzs::AbstractVector,
            params::FirstPassSearchParameters
        )
            get_best_psms!(
                psms,
                precursor_mzs,
                max_q_val=params.max_q_val_for_irt
            )
        end

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
        spectra::MassSpecData,
        search_context::SearchContext,
        rt_model::RtConversionModel,
        ms_file_idx::Int64
    )
        add_main_search_columns!(
            psms,
            getModel(rt_model),
            getStructuralMods(getPrecursors(getSpecLib(search_context))),
            getMissedCleavages(getPrecursors(getSpecLib(search_context))),
            getIsDecoy(getPrecursors(getSpecLib(search_context))),
            getIrt(getPrecursors(getSpecLib(search_context))),
            getCharge(getPrecursors(getSpecLib(search_context))),
            getRetentionTimes(spectra),
            getTICs(spectra),
            getMzArrays(spectra)
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
            :charge2, :poisson, :irt_error, 
            :missed_cleavage, 
            :Mox,
            #:charge, Only works with charge 2 if at least 3 charge states presence. otherwise singular error
            :TIC, :y_count, :err_norm, :spectrum_peak_count, :intercept
        ]

        # Select scoring columns
        select!(psms, vcat(column_names, [:ms_file_idx, :score, :precursor_idx, :scan_idx,
            :q_value, :log2_summed_intensity, :irt, :rt, :irt_predicted, :target]))
        # Score PSMs
        try
            score_main_search_psms!(
                psms,
                column_names,
                n_train_rounds=params.n_train_rounds_probit,
                max_iter_per_round=params.max_iter_probit,
                max_q_value=Float64(params.max_q_value_probit_rescore)
            )
        catch
            column_names = [
            :spectral_contrast, :city_block, :entropy_score, :scribe,
            :charge2, :poisson, :irt_error, :TIC, :y_count, :err_norm, :spectrum_peak_count, :intercept
            ]
            score_main_search_psms!(
                psms,
                column_names,
                n_train_rounds=params.n_train_rounds_probit,
                max_iter_per_round=params.max_iter_probit,
                max_q_value=Float64(params.max_q_value_probit_rescore)
            )
        end
        # Process scores
       
        select!(psms, [:ms_file_idx, :score, :precursor_idx, :scan_idx,
            :q_value, :log2_summed_intensity, :irt, :rt, :irt_predicted])
        get_probs!(psms, psms[!,:score])
        sort!(psms, :rt)
    end

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
    ::MassSpecData
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
    temp_path = joinpath(getDataOutDir(search_context), "temp_data", "first_pass_psms", parsed_fname * ".arrow")
    psms[!, :ms_file_idx] .= UInt32(ms_file_idx)
    Arrow.write(
        temp_path,
        select!(psms, [:ms_file_idx, :scan_idx, :precursor_idx, :rt,
            :irt_predicted, :q_value, :score, :prob, :scan_count])
    )
    setFirstPassPsms!(getMSData(search_context), ms_file_idx, temp_path)
end

"""
No cleanup needed between files.
"""
function reset_results!(results::FirstPassSearchResults)
    empty!(results.psms[])
    return nothing
end

"""
Summarize results across all files.
"""
function summarize_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:FirstPassSearchParameters}
    
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
            getRtIrtModel(search_context),
            max_q_val=params.max_q_val_for_irt
        )
    end
    #isempty(results.psms_paths) && return nothing
    @info "Summarizing first pass search results..."
    # Map retention times
    map_retention_times!(search_context, results, params)
    # Process precursors
    precursor_dict = get_best_precursors_accross_runs!(search_context, results, params)

    if params.match_between_runs==true
        #######
        #Each target has a corresponding decoy and vice versa
        #Add the complement targets/decoys to the precursor dict 
        #if the `sibling_peptide_scores` parameter is set to true
        #In the target/decoy scoring (see SearchMethods/ScoringSearch)
        #the maximum score for each target/decoy pair is shared accross runs
        #in an iterative training scheme. 
        precursors = getPrecursors(getSpecLib(search_context))
        for (pid, val) in pairs(precursor_dict)
            partner_pid = getPartnerPrecursorIdx(precursors)[pid]
            if ismissing(partner_pid)
                continue
            end
            if !haskey(precursor_dict, partner_pid)
                insert!(precursor_dict, partner_pid, val)
            end
        end
    end
    setPrecursorDict!(search_context, precursor_dict)
    # Calculate RT indices
    create_rt_indices!(search_context, results, precursor_dict, params)
    
    @info "Search results summarization complete"
end

