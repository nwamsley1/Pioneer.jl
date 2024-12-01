"""
Abstract type hierarchy for different spectral library search strategies
"""
abstract type SearchStrategy end
struct StandardLibrarySearch <: SearchStrategy end
struct NceScanningLibrarySearch <: SearchStrategy end
struct MassErrorLibrarySearch <: SearchStrategy end

"""
Common utilities for spectral library searching
"""
module SearchUtils
    using Arrow  # Keep existing dependencies

    function initialize_search_arrays(size::Integer=5000)
        return (
            SparseArray(UInt32(size)),
            zeros(Float32, size),
            zeros(Float32, size),
            zeros(Float32, 5),  # isotopes
            zeros(Float32, 5)   # precursor_transmission
        )
    end

    function ensure_array_capacity!(arrays::Vector, required_size::Integer, increment::Integer=1000)
        for arr in arrays
            if length(arr) < required_size
                append!(arr, [typeof(arr[1])() for _ in 1:increment])
            end
        end
    end

    function check_scan_validity(scan_idx::Integer, mz_array_length::Integer)
        return scan_idx != 0 && scan_idx <= mz_array_length
    end

    getRTWindow(irt::T, irt_tol::U) where {T,U<:AbstractFloat} = Float32(irt - irt_tol), Float32(irt + irt_tol)
end

"""
    library_search(strategy::SearchStrategy, spectra, params; kwargs...) -> DataFrame

Main entry point for spectral library searching with configurable strategies.
"""
function library_search(strategy::SearchStrategy, spectra::Arrow.Table, params::NamedTuple; kwargs...)
    thread_tasks, _ = partitionScansToThreads(
        spectra[:mz_array], spectra[:retentionTime],
        spectra[:centerMz], spectra[:msOrder],
        Threads.nthreads(), 1
    )
    
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra[:msOrder]))
    
    # Get precursors that passed scoring
    precursors_passed = collect_valid_precursors(strategy, spectra, thread_tasks, scan_to_prec_idx, params, kwargs)
    
    # Execute search based on strategy
    results = execute_search(strategy, spectra, thread_tasks, scan_to_prec_idx, precursors_passed, params, kwargs)
    
    return process_results(strategy, results)
end

# Strategy-specific implementations
function execute_search(::StandardLibrarySearch, spectra, thread_tasks, scan_to_prec_idx, precursors_passed, params, kwargs)
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            getPSMS(
                spectra, last(thread_task), kwargs[:precursors], scan_to_prec_idx,
                precursors_passed[thread_id], kwargs[:fragment_lookup_table], 
                kwargs[:rt_to_irt_spline], kwargs[:ms_file_idx], kwargs[:mass_err_model],
                kwargs[:quad_transmission_model], kwargs[:ion_matches][thread_id],
                kwargs[:ion_misses][thread_id], kwargs[:id_to_col][thread_id],
                kwargs[:ion_templates][thread_id], kwargs[:iso_splines],
                kwargs[:scored_psms][thread_id], kwargs[:unscored_psms][thread_id],
                kwargs[:spectral_scores][thread_id], kwargs[:isotope_err_bounds],
                params["min_frag_count"], params["min_spectral_contrast"],
                params["min_log2_matched_ratio"], params["min_topn_of_m"],
                params["max_best_rank"], params["n_frag_isotopes"],
                params["max_frag_rank"], params["abreviate_precursor_calc"],
                kwargs[:irt_tol], Set(2)
            )
        end
    end
    
    return fetch.(tasks)
end

"""
    getPSMS() implementation with improved organization
"""
function getPSMS(spectra::Arrow.Table, thread_task, precursors, scan_to_prec_idx, 
                precursors_passed_scoring, ion_list, args...; kwargs...)
    # Initialize arrays
    Hs, weights, residuals, isotopes, precursor_transmission = SearchUtils.initialize_search_arrays()
    msms_counts = Dict{Int64, Int64}()
    prec_idx, ion_idx, last_val = 0, 0, 0

    for i in thread_task
        # Skip invalid scans
        if !SearchUtils.check_scan_validity(i, length(spectra[:mz_array])) || 
           ismissing(scan_to_prec_idx[i])
            continue
        end

        # Process valid scans
        last_val = process_scan!(
            i, spectra, scan_to_prec_idx, precursors_passed_scoring,
            Hs, weights, residuals, isotopes, precursor_transmission,
            ion_list, args...; kwargs...
        )
    end

    return DataFrame(@view(scored_PSMs[1:last_val]))
end

"""
Core search functionality implementations. Each search type uses the strategy pattern
through SearchStrategy subtypes.
"""

module SearchCore
    using Arrow

    # Strategy-specific search implementations
    function execute_search(::NceScanningLibrarySearch, spectra, thread_tasks, scan_to_prec_idx, 
                          precursors_passed, params, nce_grid; kwargs...)
        all_results = []
        flt = kwargs[:fragment_lookup_table]
        
        for nce in nce_grid
            flt = updateNceModel(flt, PiecewiseNceModel(nce))
            results = run_parallel_search(thread_tasks, spectra, scan_to_prec_idx, 
                                       precursors_passed, flt, params, kwargs)
            
            tasks_out = vcat(results...)
            tasks_out[!, :nce] .= nce
            push!(all_results, tasks_out)
        end
        
        return vcat(all_results...)
    end

    function execute_search(::MassErrorLibrarySearch, spectra, scan_idxs, 
                          precursors_passed, library_fragment_lookup, ms_file_idx, 
                          mass_err_model, search_arrays; kwargs...)
        sorted_indices = sortperm(scan_idxs)
        scan_idxs = scan_idxs[sorted_indices]
        precursors_passed = precursors_passed[sorted_indices]
        
        scan_to_prec_idx = build_scan_precursor_index(scan_idxs, length(spectra[:mz_array]))
        thread_tasks = distribute_work(scan_idxs)
        
        return run_parallel_mass_error_search(thread_tasks, spectra, scan_idxs, 
                                            scan_to_prec_idx, precursors_passed,
                                            library_fragment_lookup, mass_err_model,
                                            search_arrays)
    end

    # Helper functions for chromatogram extraction and processing
    function process_chromatograms!(chromatograms, weights, precs_temp, prec_temp_size, 
                                  IDtoCOL, rt_idx, scan_idx, retention_time)
        for j in 1:prec_temp_size
            rt_idx += 1
            weight = !iszero(IDtoCOL[precs_temp[j]]) ? weights[IDtoCOL[precs_temp[j]]] : zero(Float32)
            
            chromatograms[rt_idx] = ChromObject(
                Float16(retention_time),
                weight,
                scan_idx,
                precs_temp[j]
            )
            
            ensure_capacity!(chromatograms, rt_idx + 1, 500000)
        end
        return rt_idx
    end

    # Core processing functions
    function process_scan!(scan_data, arrays, params; kwargs...)
        nmatches, nmisses = match_peaks!(scan_data, arrays)
        
        if nmatches > 2
            process_matches(nmatches, nmisses, arrays, params)
        end
        
        return nmatches, nmisses
    end
end

"""
High-level search function implementations using the core functionality
"""

function huberTuningSearch(spectra::Arrow.Table, thread_task::Vector{Int64},
                          params::NamedTuple; kwargs...)
    # Initialize arrays and parameters
    Hs, weights, residuals, isotopes, precursor_transmission = 
        SearchUtils.initialize_search_arrays()
    
    tuning_results = initialize_tuning_results()
    rt_params = initialize_rt_params()
    
    # Process each scan
    for scan_idx in thread_task
        process_scan_for_tuning!(
            scan_idx, spectra, rt_params, tuning_results,
            Hs, weights, residuals, isotopes, precursor_transmission,
            params, kwargs
        )
    end
    
    return DataFrame(tuning_results)
end

function QuadTransmissionSearch(spectra::Arrow.Table, thread_task::Vector{Int64},
                              params::NamedTuple; kwargs...)
    # Initialize search arrays and parameters
    search_arrays = SearchUtils.initialize_search_arrays()
    tuning_results = initialize_quad_results()
    
    # Process scans
    for scan_idx in thread_task
        if !is_valid_scan(scan_idx, spectra, kwargs[:scan_idxs], kwargs[:spec_order])
            continue
        end
        
        process_scan_for_quad_transmission!(
            scan_idx, spectra, search_arrays, tuning_results,
            params, kwargs
        )
    end
    
    return DataFrame(tuning_results)
end

"""
Utility functions for array management and data filtering
"""

module FilterUtils
    function filterMatchedIons!(matches::ArrayDict{UInt32, UInt16}, 
                              ion_matches::Vector{FragmentMatch{Float32}},
                              ion_misses::Vector{FragmentMatch{Float32}},
                              match_counts::NamedTuple;
                              min_matched_ions::Int64)
        nmatches_all, nmisses_all = match_counts.matches, match_counts.misses
        
        # Count matches per precursor
        count_matches!(matches, ion_matches, nmatches_all)
        
        # Filter matches and misses
        nmatches = filter_by_count!(ion_matches, matches, min_matched_ions, nmatches_all)
        nmisses = filter_by_count!(ion_misses, matches, min_matched_ions, nmisses_all)
        
        reset!(matches)
        return nmatches_all, nmisses_all, nmatches, nmisses
    end
    
    function collectFragErrs(all_matches::Vector{M}, new_matches::Vector{M}, 
                           nmatches::Int, current_idx::Int) where {M<:MatchIon{Float32}}
        for match_idx in 1:nmatches
            current_idx = add_match!(all_matches, new_matches[match_idx], current_idx)
        end
        return current_idx
    end
end

# Add type information and documentation to existing structs
"""
    LibraryFragmentLookup

Stores precursor and fragment information for library searching.
"""
struct LibraryFragmentLookup
    fragments::Vector{DetailedFrag{Float32}}
    prec_frag_ranges::Vector{UnitRange{UInt32}}
end

