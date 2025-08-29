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
    rm(temp_data_dir, recursive=true, force=true)
    !isdir(temp_data_dir) && mkdir(temp_data_dir)

    s.data_out_dir[] = dir

    return nothing
end