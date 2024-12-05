# New SearchDataStructures type for quad tuning
struct QuadTuningDataStructures{T<:AbstractFloat} <: SearchDataStructures
    ion_matches::Vector{FragmentMatch{T}}
    ion_misses::Vector{FragmentMatch{T}}
    id_to_col::ArrayDict{UInt32, UInt16}
    ion_templates::Vector{DetailedFrag{T}}
    iso_splines::IsotopeSplineModel
    unscored_psms::Vector{SimpleUnscoredPSM{T}}
    spectral_scores::Vector{SpectralScoresSimple{Float16}}
    precursor_weights::Vector{T}
    weights::Vector{T}  # For deconvolution
    residuals::Vector{T}  # For deconvolution
end

# Constructor for QuadTuningDataStructures
function initQuadTuningDataStructures(
    iso_splines::IsotopeSplineModel,
    n_precursors::Int64,
    M::Int64)
    
    return QuadTuningDataStructures(
        [FragmentMatch{Float32}() for _ in 1:M],
        [FragmentMatch{Float32}() for _ in 1:M],
        ArrayDict(UInt32, UInt16, n_precursors*3+1),  # *3+1 for isotope variants
        [DetailedFrag{Float32}() for _ in 1:M],
        iso_splines,
        [SimpleUnscoredPSM{Float32}() for _ in 1:5000],
        Vector{SpectralScoresSimple{Float16}}(undef, 5000),
        zeros(Float32, n_precursors*3+1),
        zeros(Float32, 5000),
        zeros(Float32, 5000)
    )
end

# Function to initialize QuadTuningDataStructures for multiple threads
function initQuadTuningSearchContexts(
    iso_splines::IsotopeSplineModel,
    n_precursors::Int64,
    N::Int64,  # number of threads
    M::Int64   # buffer size for matches
)
    return [initQuadTuningDataStructures(iso_splines, n_precursors, M) for _ in 1:N]
end

# Then to modify existing search context:
function convertToQuadTuningContext!(search_context::SearchContext)
    iso_splines = getIsoSplines(search_context.temp_structures[1])  # Get existing iso_splines
    n_precursors = length(getIdToCol(search_context.temp_structures[1]))
    N = length(search_context.temp_structures)
    M = length(getIonMatches(search_context.temp_structures[1]))
    
    # Create new quad tuning structures
    search_context.temp_structures = initQuadTuningSearchContexts(
        iso_splines,
        n_precursors,
        N,
        M
    )
    
    return search_context
end


# Main search types
struct QuadTuningSearch <: TuningMethod end

struct QuadTuningSearchResults <: SearchResults
    quad_models::Base.Ref{Dict{Int64, QuadTransmissionModel}}
    tuning_psms::Base.Ref{DataFrame}
    window_widths::Base.Ref{Set{String}}
end

# And modify the constructor to:
function init_search_results(search_parameters::QuadTuningSearchParameters, search_context::SearchContext, ms_file_idx::Int64)
    return QuadTuningSearchResults(
        Ref(Dict{Int64, QuadTransmissionModel}()),
        Ref(DataFrame()),
        Ref(Set{String}())
    )
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

# Interface implementations
get_parameters(search_type::QuadTuningSearch, params::Any) = QuadTuningSearchParameters(params)



function process_file!(
    results::QuadTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::Arrow.Table) where {P<:QuadTuningSearchParameters}

    convertToQuadTuningContext!(search_context)
    try
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
        
        results.window_widths = window_widths
        
        # Collect PSMs through multiple iterations
        total_psms = DataFrame()
        n_iterations = 0
        
        while n_iterations < 10
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
            
            # Score and filter PSMs
            scorePresearch!(psms)
            getQvalues!(psms[!, :prob], psms[!, :target], psms[!, :q_value])
            filter!(:q_value => x -> x <= 0.01, psms)
            filter!(:target => identity, psms)
            
            # Get best PSMs per precursor
            psms[!, :best_psms] .= false
            for group in groupby(psms, :precursor_idx)
                best_idx = argmax(group.prob)
                group[best_idx, :best_psms] = true
            end
            filter!(row -> row.best_psms, psms)
            
            # Create scan to precursor mapping
            scan_idx_to_prec_idx = getScanToPrecIdx(
                psms[!, :scan_idx],
                psms[!, :precursor_idx],
                spectra[:centerMz],
                spectra[:isolationWidthMz]
            )
            
            # Perform quad transmission search
            quad_psms = performQuadTransmissionSearch(
                spectra,
                scan_idx_to_prec_idx,
                search_context,
                params,
                ms_file_idx
            )
            
            # Filter and process quad search results
            quad_psms = filter_quad_psms(quad_psms, params)
            processed_psms = process_quad_results(
                quad_psms,
                getPrecursors(getSpecLib(search_context)),
                getIsoSplines(search_context)
            )
            
            append!(total_psms, processed_psms)
            
            if nrow(total_psms) > params.min_quad_tuning_psms
                break
            end
            
            n_iterations += 1
        end
        
        if nrow(total_psms) < params.min_quad_tuning_psms
            @warn "Insufficient PSMs for quad model estimation"
            results.quad_models[ms_file_idx] = GeneralGaussModel(5.0f0, 0.0f0)
            return results
        end
        
        # Fit quad model
        window_width = parse(Float64, first(window_widths))
        binned_psms = MergeBins(
            total_psms,
            (-(window_width + 1.0), window_width + 1.0),
            min_bin_size=10,
            min_bin_width=0.1
        )
        
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
        
        results.quad_models[ms_file_idx] = RazoQuadModel(fitted_model)
        append!(results.tuning_psms, total_psms)
        
    catch e
        @warn "Quad transmission function fit failed" exception=e
        results.quad_models[ms_file_idx] = GeneralGaussModel(5.0f0, 0.0f0)
    end
    
    return results
end

function process_search_results!(results::QuadTuningSearchResults, params::P, search_context::SearchContext, ms_file_idx::Int64) where {P<:QuadTuningSearchParameters}
    # Store quad model in search context
    setQuadTransmissionModel!(search_context, ms_file_idx, results.quad_models[ms_file_idx])
end

function summarize_results!(results::QuadTuningSearchResults, params::P, search_context::SearchContext) where {P<:QuadTuningSearchParameters}
    # Could add summary statistics or plots about quad transmission models if desired
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
    ms_file_idx::UInt32
)
    thread_tasks = partition_scans(spectra, Threads.nthreads())
    scan_idxs = Set(keys(scan_idx_to_prec_idx))

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]

            # Initialize working arrays
            Hs = SparseArray(UInt32(5000))
            _weights_ = zeros(Float32, 5000)
            _residuals_ = zeros(Float32, 5000)
            isotopes = zeros(Float32, 5)
            precursor_transmission = zeros(Float32, 5)

            tuning_results = Dict(
                :precursor_idx => UInt32[],
                :scan_idx => UInt32[],
                :weight => Float32[],
                :iso_idx => UInt8[],
                :center_mz => Float32[],
                :n_matches => UInt8[]
            )

            for scan_idx in last(thread_task)
                scan_idx ∉ scan_idxs && continue
                
                # Scan Filtering
                msn = spectra[:msOrder][scan_idx]
                msn ∉ params.spec_order && continue
                
                # Ion Template Selection
                ion_idx = selectTransitions!(
                    getIonTemplates(search_data),
                    QuadEstimationTransitionSelection(),
                    getPrecEstimation(params),
                    getFragmentLookupTable(getSpecLib(search_context)),
                    scan_idx_to_prec_idx[scan_idx],
                    getPrecursors(getSpecLib(search_context))[:mz],
                    getPrecursors(getSpecLib(search_context))[:prec_charge],
                    getPrecursors(getSpecLib(search_context))[:sulfur_count],
                    getIsoSplines(search_data),
                    precursor_transmission,
                    isotopes,
                    (spectra[:lowMz][scan_idx], spectra[:highMz][scan_idx]);
                    block_size = 10000
                )
                
                # Match peaks
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
                    ms_file_idx
                )
                
                sort!(@view(getIonMatches(search_data)[1:nmatches]), by = x->(x.peak_ind, x.prec_id), alg=QuickSort)
                
                # Spectral Deconvolution
                if nmatches > 2
                    buildDesignMatrix!(
                        Hs, 
                        getIonMatches(search_data), 
                        getIonMisses(search_data), 
                        nmatches, 
                        nmisses, 
                        getIdToCol(search_data)
                    )
                    
                    # Resize arrays if needed
                    if getIdToCol(search_data).size > length(_weights_)
                        new_entries = getIdToCol(search_data).size - length(_weights_) + 1000
                        append!(_weights_, zeros(eltype(_weights_), new_entries))
                        append!(getSpectralScores(search_data), Vector{eltype(getSpectralScores(search_data))}(undef, new_entries))
                        append!(getUnscoredPsms(search_data), [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])
                    end
                    
                    # Get hot start weights
                    for i in 1:getIdToCol(search_data).size
                        _weights_[getIdToCol(search_data)[getIdToCol(search_data).keys[i]]] = search_data.precursor_weights[getIdToCol(search_data).keys[i]]
                    end
                    
                    # Deconvolution
                    initResiduals!(_residuals_, Hs, _weights_)
                    solveHuber!(
                        Hs, 
                        _residuals_, 
                        _weights_,
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
                    
                    # Record results
                    for i in 1:getIdToCol(search_data).size
                        id = getIdToCol(search_data).keys[i]
                        colid = getIdToCol(search_data)[id]
                        search_data.precursor_weights[id] = _weights_[colid]
                        
                        isotope_idx = UInt8(((id - 1) % 3) + 1)
                        pid = UInt32(((id - 1) ÷ 3) + 1)
                        
                        n_matches = sum(Hs.matched[j] for j in Hs.colptr[colid]:(Hs.colptr[colid+1] - 1))
                        
                        push!(tuning_results[:precursor_idx], pid)
                        push!(tuning_results[:weight], _weights_[colid])
                        push!(tuning_results[:iso_idx], isotope_idx)
                        push!(tuning_results[:scan_idx], scan_idx)
                        push!(tuning_results[:center_mz], spectra[:centerMz][scan_idx])
                        push!(tuning_results[:n_matches], n_matches)
                    end
                end
                
                # Reset arrays
                reset!(getIdToCol(search_data))
                reset!(Hs)
            end
            
            DataFrame(tuning_results)
        end
    end
    
    return vcat(fetch.(tasks)...)
end
