struct FirstPassSearch <: SearchMethod end

struct FirstPassSearchResults <: SearchResults
    peak_fwhms::Dict{String, NamedTuple{(:median_fwhm, :mad_fwhm), Tuple{Float32, Float32}}}
    psms_paths::Dict{String, String}
    model_updates::Dict{Int64, Any}  # Store any model updates from the search
end

struct FirstPassSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    frag_tol_ppm::Float32
    n_train_rounds_probit::Int64
    max_iter_probit::Int64
    max_q_value_probit_rescore::Float32
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
    max_precursors_passing::Int64
    min_inference_points::Int64
    max_q_val_for_irt::Float32
    min_prob_for_irt_mapping::Float32
    max_precursors::Int64
    max_irt_bin_size::Float32
    max_prob_to_impute::Float32
    prec_estimation::P

    function FirstPassSearchParameters(params::Any)
        fp = params[:first_search_params]
        sp = params[:summarize_first_search_params]
        prec_estimation = fp["abreviate_precursor_calc"] ? FullPrecCapture() : PartialPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(first(params[:isotope_err_bounds])), UInt8(last(params[:isotope_err_bounds]))),
            0.0f0,#Float32(fp["frag_tol_ppm"]),
            Int64(fp["n_train_rounds_probit"]),
            Int64(fp["max_iter_probit"]),
            Float32(fp["max_q_value_probit_rescore"]),
            UInt8(fp["min_index_search_score"]),
            Int64(fp["min_frag_count"]),
            Float32(fp["min_spectral_contrast"]),
            Float32(fp["min_log2_matched_ratio"]),
            (Int64(first(fp["min_topn_of_m"])), Int64(last(fp["min_topn_of_m"]))),
            UInt8(fp["max_best_rank"]),
            Int64(fp["n_frag_isotopes"]),
            UInt8(fp["max_frag_rank"]),
            1.0f0,#Float32(fp["sample_rate"]),
            Set(2),
            Int64(fp["max_precursors_passing"]),
            Int64(sp["min_inference_points"]),
            Float32(sp["max_q_val_for_irt"]),
            Float32(params[:irt_mapping_params]["min_prob"]),
            Int64(sp["max_precursors"]),
            Float32(sp["max_irt_bin_size"]),
            Float32(sp["max_prob_to_impute"]),
            prec_estimation
        )
    end
end

# Implement required interface methods
get_parameters(search_type::FirstPassSearch, params::Any) = FirstPassSearchParameters(params)

function init_search_results(search_parameters::FirstPassSearchParameters, search_context::SearchContext, ms_file_idx::Int64)
    temp_folder = joinpath(getDataOutDir(search_context), "temp_psms")
    if !isdir(temp_folder)
        mkdir(temp_folder)
    end
    
    return FirstPassSearchResults(
        Dict{String, NamedTuple{(:median_fwhm, :mad_fwhm), Tuple{Float32, Float32}}}(),
        Dict{String, String}(),
        Dict{Int64, Any}()
    )
end

function process_file!(
    results::FirstPassSearchResults,
    params::P, 
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::Arrow.Table) where {P<:FirstPassSearchParameters}

    try
        # Get models from context
        rt_to_irt_model = getRtIrtModel(search_context, ms_file_idx)

        # Update fragment lookup table with NCE model
        setNceModel!(
            getFragmentLookupTable(getSpecLib(search_context)), 
            getNceModelModel(search_context, ms_file_idx)
        )
        # Get file name for looking up irt_err
        parsed_fname = getParsedFileName(search_context, ms_file_idx)
        
        # Perform library search
        psms = library_search(spectra, search_context, params, ms_file_idx)
        # Add columns and process results
        addMainSearchColumns!(
            psms,
            rt_to_irt_model.model,
            getPrecursors(getSpecLib(search_context))[:structural_mods],
            getPrecursors(getSpecLib(search_context))[:missed_cleavages],
            getPrecursors(getSpecLib(search_context))[:is_decoy],
            getPrecursors(getSpecLib(search_context))[:irt],
            getPrecursors(getSpecLib(search_context))[:prec_charge],
            spectra[:retentionTime],
            spectra[:TIC],
            spectra[:mz_array]
        )
        # Calculate iRT values
        psms[!, :irt_observed] = rt_to_irt_model.(psms[!, :rt])
        psms[!, :irt_error] = Float16.(abs.(psms[!, :irt_observed] .- psms[!, :irt_predicted]))
        psms[!, :charge2] = UInt8.(psms[!, :charge] .== 2)

        # Add MS file index
        psms[!, :ms_file_idx] .= UInt32(ms_file_idx)

        # Select and score columns
        column_names = [
            :spectral_contrast, :city_block, :entropy_score, :scribe,
            :charge2, :poisson, :irt_error, :missed_cleavage, :Mox,
            :charge, :TIC, :y_count, :err_norm, :spectrum_peak_count, :intercept
        ]
        
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

        # Further processing
        select!(psms, [:ms_file_idx, :score, :precursor_idx, :scan_idx,
            :q_value, :log2_summed_intensity, :irt, :rt, :irt_predicted])
        
        getProbs!(psms)
        sort!(psms, :irt)
        
        # Get best PSMs
        getBestPSMs!(
            psms,
            getPrecursors(getSpecLib(search_context))[:mz],
            max_q_val=params.max_q_val_for_irt,
            max_psms=params.max_precursors_passing
        )

        # Calculate FWHM statistics
        fwhms = skipmissing(psms[!, :fwhm])
        fwhm_points = count(!ismissing, fwhms)

        parsed_fname = getParsedFileName(search_context, ms_file_idx)
        
        if fwhm_points >= params.min_inference_points
            results.peak_fwhms[parsed_fname] = (
                median_fwhm = median(fwhms),
                mad_fwhm = mad(fwhms, normalize=true)
            )
        end

        # Save results
        temp_path = joinpath(getDataOutDir(search_context), "temp_psms", parsed_fname * ".arrow")
        psms[!, :ms_file_idx] .= UInt32(ms_file_idx)
        
        Arrow.write(
            temp_path,
            select!(psms, [:ms_file_idx, :scan_idx, :precursor_idx, :rt,
                :irt_predicted, :q_value, :score, :prob, :scan_count])
        )

        results.psms_paths[parsed_fname] = temp_path

    catch e
        @warn "First pass search failed for file index $ms_file_idx" exception=(e, catch_backtrace())
        rethrow(e)
    end

    return results
end

function process_search_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64) where {P<:FirstPassSearchParameters}
    
    # No additional processing needed after file processing
    # All results are already saved to files and dictionaries
    return nothing
end

function reset_results!(results::FirstPassSearchResults)
    #empty!(results.peak_fwhms)
    #empty!(results.psms_paths)
    #empty!(results.model_updates)
    return nothing
end




#need to implement
function summarize_results!(results::FirstPassSearchResults, params::P, search_context::SearchContext) where {P<:FirstPassSearchParameters}
    # First validate that we have results to process
    if isempty(results.psms_paths)
        @warn "No PSM results found to summarize"
        return nothing
    end

    @info "Summarizing search results..."

    # Write PSM counts
    #=
    first_psms = DataFrame(Arrow.Table(collect(values(results.psms_paths))))
    filter!(:q_value => x -> x <= 0.01, first_psms)
    psms_counts = combine(groupby(first_psms, :ms_file_idx), nrow)
    println("psms_counts ", psms_counts)
    println("typeof(psms_counts) ", typeof(psms_counts))
    csv_path = joinpath(getDataOutDir(search_context), "first_search_psms_counts.csv")
    CSV.write(csv_path, psms_counts)
    println("passed?")
    =#
    # Map library retention times to empirical retention times
    @info "Mapping library to empirical retention times..."
    irt_rt, rt_irt = mapLibraryToEmpiricalRT(
        collect(values(results.psms_paths)),
        search_context.rt_to_irt_model,
        min_prob=params.min_prob_for_irt_mapping
    )
    setIrtRtMap!(search_context, irt_rt)
    setRtIrtMap!(search_context, rt_irt)

    # Get best precursors across runs
    @info "Finding best precursors across runs..."
    precursor_dict = getBestPrecursorsAccrossRuns(
        collect(values(results.psms_paths)),
        getPrecursors(getSpecLib(search_context))[:mz],
        rt_irt,
        max_q_val=params.max_q_val_for_irt,
        max_precursors=params.max_precursors
    )
    setPrecursorDict!(search_context, precursor_dict)

    # Calculate iRT errors
    @info "Calculating iRT errors..."
    if length(keys(results.peak_fwhms)) > 1
        irt_errs = getIrtErrs(
            results.peak_fwhms,
            precursor_dict,
            params
        )
    else
        irt_errs = Dict(zip(keys(results.peak_fwhms), zeros(Float32, length(keys(results.peak_fwhms)))))
    end
    setIrtErrors!(search_context, irt_errs)

    # Create RT indices
    @info "Creating RT indices..."
    prec_to_irt = map(x -> (irt=x[:best_irt], mz=x[:mz]), precursor_dict)

    rt_indices_folder = joinpath(getDataOutDir(search_context), "rt_indices")
    !isdir(rt_indices_folder) && mkdir(rt_indices_folder)

    rt_index_paths = makeRTIndices(
        rt_indices_folder,
        collect(values(results.psms_paths)),
        prec_to_irt,
        rt_irt,
        min_prob=params.max_prob_to_impute
    )
    setRtIndexPaths!(search_context, rt_index_paths)

    @info "Search results summarization complete"
    return nothing
end