# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.


"""
Implementation of core search functionality and helper methods.
"""

#==========================================================
Core Search Execution
==========================================================#

"""
    execute_search(search_type::SearchMethod, search_context::SearchContext, params::Any)

Execute a mass spec search with the specified search type and parameters.

The execution follows these steps for each MS file:
1. Initialize search results
2. Process file data
3. Update search context with results
4. Reset for next file
5. Summarize overall results
"""
function execute_search(
    search_type::SearchMethod,
    search_context::SearchContext,
    params::PioneerParameters)

    msdr = getMassSpecData(search_context)
    n_files=length(msdr.file_paths)

    search_parameters = get_parameters(search_type, params)
    Random.seed!(1844)
    search_results = init_search_results(search_parameters, search_context)

    n_processed = 0
    n_failed = 0

    for (ms_file_idx, spectra) in ProgressBar(enumerate(msdr))
        # Skip files that have been marked as failed in previous search methods
        if is_file_failed(search_context, ms_file_idx)
            file_name = try
                getMassSpecData(search_context).file_id_to_name[ms_file_idx]
            catch
                "file_$ms_file_idx"
            end
            @debug_l1 "Skipping file $ms_file_idx ($file_name) - marked as failed in previous step"
            continue
        end

        try
            process_file!(search_results, search_parameters, search_context, ms_file_idx, spectra)
            process_search_results!(search_results, search_parameters, search_context, ms_file_idx, spectra)
            n_processed += 1
        catch e
            file_name = try
                getMassSpecData(search_context).file_id_to_name[ms_file_idx]
            catch
                "file_$ms_file_idx"
            end
            bt = catch_backtrace()
            @user_error "File $ms_file_idx ($file_name) failed during $(typeof(search_type)) processing"
            @user_error sprint(showerror, e, bt)
            mark_file_as_failed_if_needed!(search_context, ms_file_idx, e)
            n_failed += 1
            # Continue with next file instead of crashing entire search
        end

        reset_results!(search_results)
    end

    # Log summary of file processing
    if n_failed > 0
        @user_warn "$(typeof(search_type)): $n_processed files processed successfully, $n_failed files failed"
    else
        @debug_l1 "$(typeof(search_type)): $n_processed files processed successfully"
    end

    # Only proceed to summarize if we have some successful files
    if n_processed > 0
        summarize_results!(search_results, search_parameters, search_context)
    else
        @user_warn "$(typeof(search_type)): No files processed successfully - skipping result summarization"
    end

    return nothing#search_results
end

#==========================================================
Required Interface Methods
==========================================================#

"""Templates for required method implementations"""
function get_parameters(search_type::SearchMethod, params::Any)
    error("get_parameters not implemented for search method of type $(typeof(search_type))")
end

function init_search_results(search_parameters::SearchParameters, search_context::SearchContext, ms_file_idx::Int64)
    error("init_search_results not implemented for params of type $(typeof(search_parameters))")
end

function process_file!(results::SearchResults, params::SearchParameters, search_context::SearchContext, ms_file_idx::Int64, spectra::MassSpecData)
    error("process_file! not implemented for params of type $(typeof(params)) and results of type$(typeof(results)) ")
end

function summarize_results!(results::SearchResults, params::SearchParameters, search_context::SearchContext)
    error("summarize_results! not implemented for params of type $(typeof(params)) and results of type$(typeof(results)) ")
end
# ... other interface methods ...

#==========================================================
Helper Functions
==========================================================#

#==========================================================
Error Handling Helpers
==========================================================#

"""
    check_and_skip_failed_file(search_context, ms_file_idx, method_name)

Check if a file should be skipped due to previous failure and log appropriate warning.
Returns true if file should be skipped, false otherwise.
"""
function check_and_skip_failed_file(search_context::SearchContext, ms_file_idx::Int64, method_name::String)
    if is_file_failed(search_context, ms_file_idx)
        file_name = try
            getFileIdToName(getMSData(search_context), ms_file_idx)
        catch
            "file_$ms_file_idx"
        end
        @user_warn "Skipping $method_name for previously failed file: $file_name"
        return true
    end
    return false
end

"""
    handle_search_error!(search_context, ms_file_idx, method_name, error, fallback_function, results)

Handle search method errors gracefully by marking the file as failed and creating fallback results.
"""
function handle_search_error!(search_context::SearchContext, ms_file_idx::Int64, method_name::String, 
                             error::Exception, fallback_function::Function, results)
    file_name = try
        getFileIdToName(getMSData(search_context), ms_file_idx)
    catch
        "file_$ms_file_idx"
    end
    
    reason = "$method_name failed: $(typeof(error))"

    # Mark failed in both tracking systems (SearchContext set + Arrow reference flag)
    markFileFailed!(search_context, ms_file_idx, reason)
    try
        setFailedIndicator!(getMSData(search_context), ms_file_idx, true)
    catch
        # Fallback for legacy contexts
        mark_file_as_failed_if_needed!(search_context, ms_file_idx, reason)
    end

    # Log full error with stacktrace for easier debugging
    bt = catch_backtrace()
    @user_error "$method_name failed for MS data file: $file_name"
    @user_error sprint(showerror, error, bt)
    
    # Create fallback results
    fallback_function(results, ms_file_idx)
end

#==========================================================
Index-Based Failed File Handling Utilities
==========================================================#

"""
    get_valid_indexed_paths(path_array::Vector{String}, search_context::SearchContext)

Get paths for valid files, maintaining index association.
Returns: Vector of (index, path) tuples for files that passed pipeline stages.
"""
function get_valid_indexed_paths(path_array::Vector{String}, search_context::SearchContext)
    valid_indices = get_valid_file_indices(search_context)
    indexed_paths = Tuple{Int64, String}[]
    
    for idx in valid_indices
        if idx <= length(path_array)
            path = path_array[idx]
            if !isempty(path) && isfile(path)
                push!(indexed_paths, (idx, path))
            end
        end
    end
    
    return indexed_paths
end

"""
    get_valid_paths_only(path_array::Vector{String}, search_context::SearchContext)

Get paths for valid files as a simple vector (without index information).
"""
function get_valid_paths_only(path_array::Vector{String}, search_context::SearchContext)
    return [path for (_, path) in get_valid_indexed_paths(path_array, search_context)]
end

"""
    get_valid_file_names_by_indices(search_context::SearchContext)

Get file names for valid files, preserving index order.
Uses the failed file tracking in SearchContext to determine which files are valid.
"""
function get_valid_file_names_by_indices(search_context::SearchContext)
    valid_indices = get_valid_file_indices(search_context)
    all_names = getFileIdToName(getMSData(search_context))
    
    # Return names in index order (indices are already sorted)
    return [all_names[idx] for idx in valid_indices]
end

"""
    partition_scans(ms_table, n_threads)

Partition MS data into chunks for parallel processing.
"""
function partition_scans(ms_table, n_threads; ms_order_select = 2)

    if ms_order_select == 2
    thread_tasks, total_peaks = partitionScansToThreads(
        getMzArrays(ms_table),
        getRetentionTimes(ms_table),
        getCenterMzs(ms_table),
        getMsOrders(ms_table),
        n_threads,
        1
    )
    else
    thread_tasks, total_peaks = partitionScansToThreadsMS1(
        getMzArrays(ms_table),
        getRetentionTimes(ms_table),
        getCenterMzs(ms_table),
        getMsOrders(ms_table),
        n_threads,
        1
    )
    end
    return thread_tasks
end

"""
Specialized partition_scans for IndexedMassSpecData.
This function handles the proper mapping between virtual indices (used by library search)
and actual scan indices (stored in IndexedMassSpecData.scan_indices).
"""
function partition_scans(indexed_data::IndexedMassSpecData, n_threads; ms_order_select = 2)

    # Get properties for the actual scan indices
    actual_scan_indices = indexed_data.scan_indices
    original_data = indexed_data.original_data

    # Extract properties using actual scan indices
    rt_values = Float32[getRetentionTime(original_data, actual_idx) for actual_idx in actual_scan_indices]
    center_mz_values = Union{Missing, Float32}[getCenterMz(original_data, actual_idx) for actual_idx in actual_scan_indices]
    ms_orders = UInt8[getMsOrder(original_data, actual_idx) for actual_idx in actual_scan_indices]
    mz_arrays = [getMzArray(original_data, actual_idx) for actual_idx in actual_scan_indices]

    @debug_l2 "partition_scans(IndexedMassSpecData): Processing $(length(actual_scan_indices)) scans"
    @debug_l2 "partition_scans(IndexedMassSpecData): RT range: $(minimum(rt_values)) - $(maximum(rt_values))"

    # Use the existing partitioning logic but with extracted data
    if ms_order_select == 2
        thread_tasks, total_peaks = partitionScansToThreadsIndexed(
            mz_arrays,
            rt_values,
            center_mz_values,
            ms_orders,
            actual_scan_indices,
            n_threads,
            1
        )
    else
        thread_tasks, total_peaks = partitionScansToThreadsMS1Indexed(
            mz_arrays,
            rt_values,
            center_mz_values,
            ms_orders,
            actual_scan_indices,
            n_threads,
            1
        )
    end

    @debug_l2 "partition_scans(IndexedMassSpecData): Created $(length(thread_tasks)) thread tasks"
    for (i, task) in enumerate(thread_tasks)
        task_range = last(task)
        @debug_l2 "  Thread $i: $(length(task_range)) scans"
    end

    return thread_tasks
end

"""
    initSearchContext(spec_lib, iso_splines, ms_data_reference, n_threads, n_precursors, buffer_size)

Initialize a new search context with simple library search structures.
"""
function initSearchContext(
    spec_lib::SpectralLibrary,
    iso_splines::IsotopeSplineModel,
    ms_data_reference::MassSpecDataReference,
    n_threads::Int64,
    buffer_size::Int64
)
    n_precursors = length(getPrecursors(spec_lib))
    temp_structures = initSimpleSearchContexts(
        iso_splines,
        n_precursors,
        n_threads,
        buffer_size
    )
    
    # Create SearchContext
    search_context = SearchContext(
        spec_lib,
        temp_structures,
        ms_data_reference,
        n_threads,
        n_precursors,
        buffer_size
    )
    
    # Calculate target/decoy statistics from library
    is_decoy_array = getIsDecoy(getPrecursors(spec_lib))
    n_targets = 0
    n_decoys = 0
    
    # Use loop to avoid allocation
    for i in 1:length(is_decoy_array)
        if is_decoy_array[i]
            n_decoys += 1
        else
            n_targets += 1
        end
    end
    
    # Set library statistics
    search_context.n_library_targets = n_targets
    search_context.n_library_decoys = n_decoys
    search_context.library_fdr_scale_factor = Float32(n_targets / max(n_decoys, 1))
    
    @debug_l1 "Library contains $n_targets targets and $n_decoys decoys (FDR scale factor: $(round(search_context.library_fdr_scale_factor, digits=7)))"
    
    return search_context
end

function initSimpleSearchContexts(
    iso_splines::IsotopeSplineModel,
    n_precursors::Int64,
    N::Int64,
    M::Int64)
    [initSimpleSearchContext(iso_splines, n_precursors, M) for _ in 1:N]
end


function initSimpleSearchContext(
    iso_splines::IsotopeSplineModel,
    n_precursors::Int64,
    M::Int64
    )

    SimpleLibrarySearch(
        [FragmentMatch{Float32}() for _ in range(1, M)],
        [FragmentMatch{Float32}() for _ in range(1, M)],
        [FragmentMatch{Float32}() for _ in range(1, M)],
        ArrayDict(UInt32, UInt16, n_precursors),
        Counter(UInt32, UInt8,n_precursors ),
        [DetailedFrag{Float32}() for _ in range(1, M)],
        iso_splines,
        Vector{SimpleScoredPSM{Float32, Float16}}(undef, 5000),
        [SimpleUnscoredPSM{Float32}() for _ in range(1, 5000)],
        Vector{SpectralScoresSimple{Float16}}(undef, 5000),
        Vector{ComplexScoredPSM{Float32, Float16}}(undef, 5000),
        [ComplexUnscoredPSM{Float32}() for _ in range(1, 5000)],
        Vector{SpectralScoresComplex{Float16}}(undef, 5000),
        Vector{Ms1ScoredPSM{Float32, Float16}}(undef, 5000),
        [Ms1UnscoredPSM{Float32}() for _ in range(1, 5000)],
        Vector{SpectralScoresMs1{Float16}}(undef, 5000),
        SparseArray(UInt32(5000)),
        zeros(UInt32, 5000),
        zeros(Float32, n_precursors),
        zeros(Float32, 5000),
        zeros(Float32, 5000),
        zeros(Float32, 5),
        zeros(Float32, 5),
    )
end

#==========================================================
Failed File Management Utilities
==========================================================#

"""
    is_file_failed(search_context, ms_file_idx) -> Bool

Check if a file has been marked as failed in the ArrowTableReference.
"""
function is_file_failed(search_context, ms_file_idx)
    return getFailedIndicator(getMassSpecData(search_context), ms_file_idx)
end

"""
    mark_file_as_failed_if_needed!(search_context, ms_file_idx, reason)

Mark a file as failed in the ArrowTableReference and log the reason.
"""
function mark_file_as_failed_if_needed!(search_context, ms_file_idx, reason)
    setFailedIndicator!(getMassSpecData(search_context), ms_file_idx, true)
    file_name = try
        getMassSpecData(search_context).file_id_to_name[ms_file_idx]
    catch
        "file_$ms_file_idx"
    end
    @user_warn "File $ms_file_idx ($file_name) marked as failed: $reason"
end

"""
    get_valid_file_indices(search_context) -> Vector{Int}

Get indices of files that have not been marked as failed.
"""
function get_valid_file_indices(search_context)
    ms_data = getMassSpecData(search_context)
    indices = [i for i in 1:length(ms_data.file_paths) if !getFailedIndicator(ms_data, i)]
    return sort(indices)  # Ensure consistent ordering across all uses
end

"""
    get_valid_file_paths(search_context, path_accessor_fn) -> Vector{Tuple{Int, String}}

Get (index, path) pairs for valid files using the provided path accessor function.
Returns only non-empty paths from files that haven't failed.
"""
function get_valid_file_paths(search_context, path_accessor_fn)
    valid_indices = get_valid_file_indices(search_context)
    all_paths = path_accessor_fn(getMassSpecData(search_context))
    return [(i, all_paths[i]) for i in valid_indices if i <= length(all_paths) && !isempty(all_paths[i])]
end

"""
    filter_to_valid_files(file_paths::Vector{String}, search_context) -> Vector{String}

Filter a vector of file paths to only include those from valid (non-failed) files.
"""
function filter_to_valid_files(file_paths::Vector{String}, search_context)
    valid_indices = get_valid_file_indices(search_context)
    return [file_paths[i] for i in valid_indices if i <= length(file_paths) && !isempty(file_paths[i])]
end

"""
    get_valid_fold_file_paths(search_context) -> Vector{String}

Get paths to all fold-split second pass PSM files for valid (non-failed) files.

Returns paths to both fold0 and fold1 files for each valid MS file.
Files that don't exist are filtered out (e.g., if a file had no PSMs for one fold).

# Returns
- Vector of paths to fold-specific Arrow files (e.g., ["file1_fold0.arrow", "file1_fold1.arrow", ...])
"""
function get_valid_fold_file_paths(search_context)
    valid_indices = get_valid_file_indices(search_context)
    ms_data = getMassSpecData(search_context)

    fold_paths = String[]
    for idx in valid_indices
        base_path = getSecondPassPsms(ms_data, idx)
        if !isempty(base_path)
            # Add paths for both folds if they exist
            for fold in UInt8[0, 1]
                fold_path = getSecondPassPsmsFold(ms_data, idx, fold)
                if isfile(fold_path)
                    push!(fold_paths, fold_path)
                end
            end
        end
    end
    return fold_paths
end

"""
    safe_process_file!(search_context, ms_file_idx, processing_fn) -> Bool

Safely process a file with error handling. Returns true if successful, false if failed.
Marks file as failed if processing produces no data or throws an error.
"""
function safe_process_file!(search_context, ms_file_idx, processing_fn)
    if is_file_failed(search_context, ms_file_idx)
        return false  # Already failed
    end

    try
        result = processing_fn(ms_file_idx)
        if result === nothing || (isa(result, DataFrame) && nrow(result) == 0)
            mark_file_as_failed_if_needed!(search_context, ms_file_idx, "No valid data produced")
            return false
        end
        return true
    catch e
        mark_file_as_failed_if_needed!(search_context, ms_file_idx, e)
        return false
    end
end


"""
   setDataOutDir!(s::SearchContext, dir::String) -> SearchContext

Set up the output directory structure for a search.

Creates the following directory hierarchy:
- Main output directory (dir)
 - qc_plots/
   - rt_alignment_plots/
   - mass_error_plots/

# Arguments
- `s::SearchContext`: The search context to update
- `dir::String`: Path to main output directory

# Returns
- The modified SearchContext

# Example
setDataOutDir!(search_context, "/path/to/output")
"""
function setDataOutDir!(s::SearchContext, dir::String)
    # Create QC plots directory
    qc_plot_folder = joinpath(dir, "qc_plots")
    !isdir(qc_plot_folder) && mkdir(qc_plot_folder)

    # Create subdirectories for specific plot types 
    rt_alignment_folder = joinpath(qc_plot_folder, "rt_alignment_plots")
    mass_error_folder = joinpath(qc_plot_folder, "mass_error_plots")
    ms1_mass_error_folder = joinpath(qc_plot_folder, "ms1_mass_error_plots")
    !isdir(rt_alignment_folder) && mkdir(rt_alignment_folder)
    !isdir(mass_error_folder) && mkdir(mass_error_folder)
    !isdir(ms1_mass_error_folder) && mkdir(ms1_mass_error_folder)

    # Store paths in search context
    s.qc_plot_folder[] = qc_plot_folder
    s.rt_alignment_plot_folder[] = rt_alignment_folder
    s.mass_err_plot_folder[] = mass_error_folder
    s.ms1_mass_err_plot_folder[] = ms1_mass_error_folder

    temp_data_dir = joinpath(dir, "temp_data")
    # Delete previous temp data if it exists
    try
        rm(temp_data_dir, recursive=true, force=true)
    catch e
        @warn "Could not fully remove previous temp_data" exception=e
    end
    !isdir(temp_data_dir) && mkdir(temp_data_dir)

    s.data_out_dir[] = dir

    return nothing
end
