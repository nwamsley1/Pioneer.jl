
#=
function convertToQuadTuningContext!(search_context::SearchContext) 
    iso_splines = getIsoSplines(search_context.temp_structures[1])
    
    # Create new quad structures with correct type
    quad_structures = initQuadTuningSearchContexts(
        iso_splines,
        getNPrecursors(search_context),
        getNThreads(search_context),
        getBufferSize(search_context)
    )
    
    # Empty existing vector
    empty!(search_context.temp_structures)
    
    # Append each new structure individually
    for quad_struct in quad_structures
        push!(search_context.temp_structures, quad_struct)
    end
    
    return search_context
end
=#
# Main search types
struct QuadTuningSearch <: TuningMethod end

struct QuadTuningSearchResults <: SearchResults
    quad_models::Base.Ref{Dict{Int64, QuadTransmissionModel}}
    tuning_psms::Base.Ref{DataFrame}
    window_widths::Base.Ref{Set{String}}
end

struct QuadTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
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
    irt_tol::Float32
    spec_order::Set{Int64}
    
    # Deconvolution parameters
    max_iter_newton::Int64
    max_iter_bisection::Int64
    max_iter_outer::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    max_diff::Float32
    
    # Quad tuning specific parameters
    min_quad_tuning_fragments::Int64
    min_quad_tuning_psms::Int64
    quad_tuning_sample_rate::Float32
    prec_estimation::P

    function QuadTuningSearchParameters(params::Any)
        pp = params[:presearch_params]
        dp = params[:deconvolution_params]
        prec_estimation = pp["abreviate_precursor_calc"] ? FullPrecCapture() : PartialPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(0), UInt8(0)),  # isotope_err_bounds fixed for quad tuning
            Float32(pp["frag_tol_ppm"]),
            UInt8(pp["min_index_search_score"]),
            Int64(pp["min_frag_count"]),
            Float32(pp["min_spectral_contrast"]),
            Float32(pp["min_log2_matched_ratio"]),
            (Int64(first(pp["min_topn_of_m"])), Int64(last(pp["min_topn_of_m"]))),
            UInt8(pp["max_best_rank"]),
            Int64(pp["n_frag_isotopes"]),
            UInt8(pp["max_frag_rank"]),
            Float32(pp["quad_tuning_sample_rate"]),
            typemax(Float32),
            Set(2),
            Int64(dp["max_iter_newton"]),
            Int64(dp["max_iter_bisection"]),
            Int64(dp["max_iter_outer"]),
            Float32(dp["accuracy_newton"]),
            Float32(dp["accuracy_bisection"]),
            Float32(dp["max_diff"]),
            Int64(pp["min_quad_tuning_fragments"]),
            Int64(pp["min_quad_tuning_psms"]),
            Float32(pp["quad_tuning_sample_rate"]),
            prec_estimation
        )
    end
end
# And modify the constructor to:
function init_search_results(search_parameters::QuadTuningSearchParameters, search_context::SearchContext, ms_file_idx::Int64)
    return QuadTuningSearchResults(
        Ref(Dict{Int64, QuadTransmissionModel}()),
        Ref(DataFrame()),
        Ref(Set{String}())
    )
end

# Interface implementations
get_parameters(search_type::QuadTuningSearch, params::Any) = QuadTuningSearchParameters(params)

function adjustPrecursorArrays!(search_context::SearchContext)
    target_size = getNPrecursors(search_context) * 3 + 1

    # Only adjust if not already the right size
    if length(first(search_context.temp_structures).precursor_weights) == getNPrecursors(search_context)
        for search_data in search_context.temp_structures
            # Create new IDtoCol of expanded size
            search_data.id_to_col = ArrayDict(UInt32, UInt16, target_size)
            
            # Expand weights
            if length(search_data.precursor_weights) == getNPrecursors(search_context)
                resize!(search_data.precursor_weights, target_size)
                search_data.precursor_weights[length(search_data.precursor_weights):end] .= zero(Float32)
            end
        end
    end
    return search_context
end

function resetPrecursorArrays!(search_context::SearchContext)
    original_size = getNPrecursors(search_context)

    # Only reset if currently expanded
    if length(first(search_context.temp_structure).precursor_weights) != original_size
        for search_data in search_context.temp_structures
            # Reset IDtoCol to original size
            search_data.id_to_col = ArrayDict(UInt32, UInt16, original_size)
            
            # Reset weights to original size
            resize!(search_data.precursor_weights, original_size)
        end
    end
    return search_context
end

function process_file!(
    results::QuadTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::Arrow.Table) where {P<:QuadTuningSearchParameters}

    #convertToQuadTuningContext!(search_context)
    println("TEST")
    setQuadTransmissionModel!(search_context, ms_file_idx, SquareQuadModel(1.0f0))
    try
        adjustPrecursorArrays!(search_context)
        #convertToQuadTuningContext!(search_context)
        # Check for multiple window widths
        window_widths = Set{String}()
        for i in 1:length(spectra[:isolationWidthMz])
            if spectra[:msOrder][i] == 2
                if !ismissing(spectra[:isolationWidthMz][i])
                    push!(window_widths, string(spectra[:isolationWidthMz][i]))
                end
            end
        end
        
        if length(window_widths) != 1
            @warn "Multiple window sizes detected: $(join(collect(window_widths), ';'))"
            results.quad_models[ms_file_idx] = GeneralGaussModel(5.0f0, 0.0f0)
            return results
        end
        
        #$results.window_widths = window_widths
        
        # Collect PSMs through multiple iterations
        total_psms = DataFrame()
        n_iterations = 0
        
        while n_iterations < 1
            # Get initial PSMs using library search
            psms = library_search(spectra, search_context, params, ms_file_idx)
            
            if isempty(psms)
                n_iterations += 1
                continue
            end
            
            # Add columns and filter PSMs
            addPreSearchColumns!(
                psms,
                spectra,
                getPrecursors(getSpecLib(search_context))[:is_decoy],
                getPrecursors(getSpecLib(search_context))[:irt],
                getPrecursors(getSpecLib(search_context))[:prec_charge],
                spectra[:retentionTime],
                spectra[:TIC]
            )
            Arrow.write("/Users/n.t.wamsley/Desktop/og_psms.arrow", psms)
            # Score and filter PSMs
            scorePresearch!(psms)
            getQvalues!(psms[!, :prob], psms[!, :target], psms[!, :q_value])
            println("a size(psms) ", size(psms))
            filter!(:q_value => x -> x <= 0.01, psms)
            println("b size(psms) ", size(psms))
            filter!(:target => identity, psms)
            
            println("c size(psms) ", size(psms))
            # Get best PSMs per precursor
            psms[!, :best_psms] .= false
            for group in groupby(psms, :precursor_idx)
                best_idx = argmax(group.prob)
                group[best_idx, :best_psms] = true
            end
            filter!(row -> row.best_psms, psms)
            println("d size(psms) ", size(psms))
            # Create scan to precursor mapping
            scan_idx_to_prec_idx = getScanToPrecIdx(
                psms[!, :scan_idx],
                psms[!, :precursor_idx],
                spectra[:centerMz],
                spectra[:isolationWidthMz]
            )
            println("total ", sum([abs(length(x)) for x in values(scan_idx_to_prec_idx)]))
            # Perform quad transmission search
            quad_psms = performQuadTransmissionSearch(
                spectra,
                scan_idx_to_prec_idx,
                search_context,
                params,
                ms_file_idx
            )
            
            # Filter and process quad search results
            println("size(quad_psms ) ", size(quad_psms))
            quad_psms = quad_psms[filter_quad_psms(quad_psms[!,:iso_idx], quad_psms[!,:n_matches], quad_psms[!,:weight], params),:]
            println("size(quad_psms) ", size(quad_psms))
            Arrow.write("/Users/n.t.wamsley/Desktop/tpsms.arrow", quad_psms)
            processed_psms = process_quad_results(
                quad_psms,
                getPrecursors(getSpecLib(search_context)),
                getIsoSplines(first(getSearchData(search_context)))
            )
            
            append!(total_psms, processed_psms)
            
            if nrow(total_psms) > 0#params.min_quad_tuning_psms
                break
            end
            
            n_iterations += 1
        end
        
        if nrow(total_psms) < 0#params.min_quad_tuning_psms
            @warn "Insufficient PSMs for quad model estimation"
            results.quad_models[ms_file_idx] = GeneralGaussModel(5.0f0, 0.0f0)
            return results
        end
        Arrow.write("/Users/n.t.wamsley/Desktop/ttable.arrow", total_psms)
        p = plot()
        plot!(p, total_psms[total_psms[!,:prec_charge].==2,:x0], total_psms[total_psms[!,:prec_charge].==2,:yt], seriestype=:scatter, alpha = 0.1, show = true)
        plot!(p, total_psms[total_psms[!,:prec_charge].==3,:x0], total_psms[total_psms[!,:prec_charge].==3,:yt], seriestype=:scatter, alpha = 0.1, show = true)
        println("size(total_psms) ", size(total_psms))
        println(describe(total_psms))
        # Fit quad model
        window_width = parse(Float64, first(window_widths))
        binned_psms = MergeBins(
            total_psms,
            (-(window_width + 1.0), window_width + 1.0),
            min_bin_size=10,
            min_bin_width=0.1
        )
        #plot( binned_psms[!, :median_x0],
        #binned_psms[!, :median_yt], series_type=:scatter, show = true
        #)
        fitted_model = fitRazoQuadModel(
            window_width,
            binned_psms[!, :median_x0],
            binned_psms[!, :median_x1],
            binned_psms[!, :median_yt],
            λ0=1e-1,
            ϵ1=1e-5,
            ϵ2=1e-4,
            ϵ3=1e-5
        )
        
        results.quad_models[][ms_file_idx] = RazoQuadModel(fitted_model)
        #append!(results.tuning_psms, total_psms)
        
    catch e
        @warn "Quad transmission function fit failed" exception=e
        merge_into_ref_dict!(
            results.quad_models,
            Dict(ms_file_idx => GeneralGaussModel(5.0f0, 0.0f0))
        )
    end
    
    return results
end


function process_quad_results(
    psms::DataFrame,
    precursors::Arrow.Table,
    iso_splines::IsotopeSplineModel
)
    # Sort by scan and precursor index
    sort!(psms, [:scan_idx, :precursor_idx, :iso_idx])
    
    # Add necessary columns
    processed = hcat(psms, addColumns(
        psms[!, :precursor_idx],
        precursors[:mz],
        precursors[:prec_charge],
        precursors[:sulfur_count],
        psms[!, :iso_idx],
        psms[!, :center_mz],
        iso_splines
    ))
    
    # Combine precursor results
    combined = combine(groupby(processed, [:scan_idx, :precursor_idx])) do group
        summarizePrecursor(
            group[!, :iso_idx],
            group[!, :center_mz],
            group[!, :iso_mz],
            group[!, :prec_charge],
            group[!, :weight],
            group[!, :δ]
        )
    end
    
    # Filter and convert types
    filter!(row -> !ismissing(row.yt), combined)
    combined[!, :prec_charge] = UInt8.(combined[!, :prec_charge])
    combined[!, :x0] = Float32.(combined[!, :x0])
    combined[!, :yt] = Float32.(combined[!, :yt])
    combined[!, :x1] = Float32.(combined[!, :x1])
    filter!(row -> row.prec_charge < 4, combined)
    
    return combined
end
function merge_into_ref_dict!(ref_dict::Base.RefValue{Dict{K,V}}, new_dict::Dict{K,V}) where {K,V}
    ref_dict[] = merge(ref_dict[], new_dict)
end

function filter_quad_psms(iso_idx::AbstractVector{UInt8},
    n_matches::AbstractVector{UInt8},
    weight::AbstractVector{Float32}, 
    params::QuadTuningSearchParameters)

    #filter!(row -> row.n_matches >= params.min_quad_tuning_fragments, psms)
    n = length(iso_idx)
    mask = Vector{Bool}(undef, n)
    @inbounds for i in 1:n
        mask[i] = ((iso_idx[i] < 3) && #Only M0 and M1
                (n_matches[i] >= 5) &&#params.min_quad_tuning_fragments) &&  
                (weight[i]>0)) #Non-zero abundance 
    end
    return mask
end

function process_search_results!(results::QuadTuningSearchResults, params::P, search_context::SearchContext, ms_file_idx::Int64) where {P<:QuadTuningSearchParameters}
    # Store quad model in search context
    setQuadTransmissionModel!(search_context, ms_file_idx, results.quad_models[ms_file_idx])
end

function summarize_results!(results::QuadTuningSearchResults, params::P, search_context::SearchContext) where {P<:QuadTuningSearchParameters}
    # Could add summary statistics or plots about quad transmission models if desired
    resetPrecursorArrays!(search_context)
    return nothing
end

function reset_results!(results::QuadTuningSearchResults)
    empty!(results.quad_models)
    empty!(results.tuning_psms)
    empty!(results.window_widths)
end


function performQuadTransmissionSearch(
    spectra::Arrow.Table,
    scan_idx_to_prec_idx::Dictionary{UInt32, Vector{UInt32}},
    search_context::SearchContext,
    params::QuadTuningSearchParameters,
    ms_file_idx::Int64)
    thread_tasks = partition_scans(spectra, Threads.nthreads())
    scan_idxs = Set(keys(scan_idx_to_prec_idx))

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin


            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]
            Hs = getHs(search_data,)
            weights = getTempWeights(search_data)
            precursor_weights = getPrecursorWeights(search_data)
            residuals = getResiduals(search_data)

            for scan_idx in last(thread_task)
                scan_idx ∉ scan_idxs && continue
                
                msn = spectra[:msOrder][scan_idx]
                msn ∉ params.spec_order && continue
                
                ion_idx = selectTransitions!(
                    getIonTemplates(search_data),
                    QuadEstimationTransitionSelection(),
                    PartialPrecCapture(),
                    getFragmentLookupTable(getSpecLib(search_context)),
                    scan_idx_to_prec_idx[scan_idx],
                    getPrecursors(getSpecLib(search_context))[:mz],
                    getPrecursors(getSpecLib(search_context))[:prec_charge],
                    getPrecursors(getSpecLib(search_context))[:sulfur_count],
                    getIsoSplines(search_data),
                    getPrecursorTransmission(search_data),
                    getIsotopes(search_data),
                    (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
                    block_size = 10000
                )
                
                nmatches, nmisses = matchPeaks!(
                    getIonMatches(search_data),
                    getIonMisses(search_data),
                    getIonTemplates(search_data),
                    ion_idx,
                    spectra[:mz_array][scan_idx],
                    spectra[:intensity_array][scan_idx],
                    getMassErrorModel(search_context, ms_file_idx),
                    spectra[:highMz][scan_idx],
                    UInt32(scan_idx),
                    UInt32(ms_file_idx)
                )
                
                sort!(@view(getIonMatches(search_data)[1:nmatches]), by = x->(x.peak_ind, x.prec_id), alg=QuickSort)
                
                if nmatches > 2
            
                    buildDesignMatrix!(
                        Hs, 
                        getIonMatches(search_data), 
                        getIonMisses(search_data), 
                        nmatches, 
                        nmisses, 
                        getIdToCol(search_data)
                    )
                    
                    if getIdToCol(search_data).size > length(weights)
                        new_entries = getIdToCol(search_data).size - length(weights) + 1000
                        resize!(weights, length(weights) + new_entries)
                        resize!(getSpectralScores(search_data), length(getSpectralScores(search_data)) + new_entries)
                        append!(getUnscoredPsms(search_data), [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])
                    end
                    
                    for i in 1:getIdToCol(search_data).size
                        weights[getIdToCol(search_data)[getIdToCol(search_data).keys[i]]] = precursor_weights[getIdToCol(search_data).keys[i]]
                    end
                    
                    initResiduals!(residuals, Hs, weights)

                    solveHuber!(
                        Hs, 
                        residuals, 
                        weights,
                        Float32(100000.0f0), 
                        0.0f0,
                        params.max_iter_newton,
                        params.max_iter_bisection,
                        params.max_iter_outer,
                        params.accuracy_newton,
                        params.accuracy_bisection,
                        10.0,
                        params.max_diff
                    )
                    
                    tuning_results = getTuningResults(search_data)

                    for i in 1:getIdToCol(search_data).size
                        id = getIdToCol(search_data).keys[i]
                        colid = getIdToCol(search_data)[id]
                        precursor_weights[id] = weights[colid]
                        
                        isotope_idx = UInt8(((id - 1) % 3) + 1)
                        pid = UInt32(((id - 1) ÷ 3) + 1)
                        
                        n_matches = sum(Hs.matched[j] for j in Hs.colptr[colid]:(Hs.colptr[colid+1] - 1))
                        
                        push!(tuning_results[:precursor_idx], pid)
                        push!(tuning_results[:weight], weights[colid])
                        push!(tuning_results[:iso_idx], isotope_idx)
                        push!(tuning_results[:scan_idx], scan_idx)
                        push!(tuning_results[:center_mz], spectra[:centerMz][scan_idx])
                        push!(tuning_results[:n_matches], n_matches)
                    end
                end
                
                reset!(getIdToCol(search_data))
                reset!(Hs)
            end
            
            return DataFrame(getTuningResults(search_data))
        end
    end
    
    return vcat(fetch.(tasks)...)
end