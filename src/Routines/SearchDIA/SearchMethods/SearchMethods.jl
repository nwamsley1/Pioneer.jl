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
    search_results = init_search_results(search_type, search_parameters, search_context)

    for (ms_file_idx, spectra) in ProgressBar(enumerate(msdr))
        process_file!(search_results, search_parameters, search_context, ms_file_idx, spectra)
        process_search_results!(search_results, search_parameters, search_context, ms_file_idx, spectra)
        reset_results!(search_results)
    end

    summarize_results!(search_results, search_parameters, search_context)

    return nothing#search_results
end

#==========================================================
Required Interface Methods
==========================================================#

"""Templates for required method implementations"""
function get_parameters(search_type::SearchMethod, params::Any)
    error("get_parameters not implemented for search method of type $(typeof(search_type))")
end

# Default: dispatch on (search_type, params) falls back to (params) for backward compatibility
init_search_results(::SearchMethod, search_parameters::SearchParameters, search_context::SearchContext) =
    init_search_results(search_parameters, search_context)

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

"""
    get_all_indexed_paths(path_accessor_fn, search_context) -> Vector{Tuple{Int64, String}}

Iterate every MS file index, returning `(idx, path)` tuples for files whose
`path_accessor_fn(...)` registered a non-empty path. Replaces the old
`get_valid_file_paths` helper, which additionally filtered by a per-file
failure flag — that flag has been removed.
"""
function get_all_indexed_paths(path_accessor_fn, search_context::SearchContext)
    ms_data = getMassSpecData(search_context)
    all_paths = path_accessor_fn(ms_data)
    return [(i, all_paths[i]) for i in 1:length(all_paths) if !isempty(all_paths[i])]
end

"""
    get_all_fold_file_paths(search_context) -> Vector{String}

Collect every existing fold-split second-pass PSM file (fold0/fold1) across
all MS file indices. Files with no on-disk fold output (e.g. the file had no
PSMs in that fold) are skipped.
"""
function get_all_fold_file_paths(search_context::SearchContext)
    ms_data = getMassSpecData(search_context)
    fold_paths = String[]
    for idx in 1:length(getFilePaths(ms_data))
        base_path = getSecondPassPsms(ms_data, idx)
        isempty(base_path) && continue
        for fold in UInt8[0, 1]
            fold_path = getSecondPassPsmsFold(ms_data, idx, fold)
            isfile(fold_path) && push!(fold_paths, fold_path)
        end
    end
    return fold_paths
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
    
    @user_info "Library: $n_targets target precursors, $n_decoys decoys (FDR scale factor: $(round(search_context.library_fdr_scale_factor, digits=4)))"
    
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
        [MassErrSample() for _ in range(1, M)],            # mass_err_samples
        # id_to_col: reset per scan, ~5k active. Hint to 8k = 8192 next pow2,
        # avoids rehash on busy scans.
        SparsePrecMap{UInt16}(sizehint=8192),
        iso_splines,
        [MainUnscoredPSM{Float32}() for _ in range(1, 5000)],
        StructArray{MainSearchScoredPSM{Float32, Float16}}(undef, 5000),  # SoA scored buffer
        Vector{SpectralScoresMainSearch{Float16, Float32}}(undef, 5000),
        # Tuning buffers — slim PSM variants for ParameterTuning/QuadTuning/Integrate paths.
        [TuningUnscoredPSM{Float32}() for _ in range(1, 5000)],
        StructArray{TuningScoredPSM{Float32, Float16}}(undef, 5000),      # SoA scored buffer
        zeros(UInt32, 5000),
        # precursor_weights: cumulative matched-precursors over a run.
        # Hint to 256k — rough upper bound for typical DIA fill; resizes
        # gracefully if exceeded.
        SparsePrecMap{Float32}(sizehint=262144),
        zeros(Float32, 5000),
        zeros(Float32, 5000),
        zeros(Float32, 5000),
        zeros(Float32, 5000),
        zeros(Float32, 5000),
        zeros(Float32, 5),
        zeros(Float32, 5),
        FusedScratch(256),
        SparseArrayFused(UInt32(5000)),
        zeros(Float32, 5000),  # scan_corrected_mz
        zeros(Float32, 5000),  # scan_obs_low
        zeros(Float32, 5000),  # scan_obs_high
    )
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
