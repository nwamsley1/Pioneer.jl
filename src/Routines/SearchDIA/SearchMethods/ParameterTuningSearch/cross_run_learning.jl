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

# Type definitions moved to types.jl

"""
    InitialParameters

Parameters to use for initial search, potentially informed by previous runs.
"""
struct InitialParameters
    bias_estimate::Float32
    bias_search_range::Float32
    initial_tolerance::Float32
    informed_by_history::Bool
end

"""
    store_tuning_results!(history::ParameterHistory, file_idx::Int64, results::TuningResults)

Store the parameter tuning results for a file.
"""
function store_tuning_results!(history::ParameterHistory, file_idx::Int64, results::TuningResults)
    history.file_parameters[file_idx] = results
    
    # Update global statistics if this was successful
    if results.converged
        update_global_statistics!(history)
    end
end

"""
    update_global_statistics!(history::ParameterHistory)

Recompute global statistics from all successful files.
"""
function update_global_statistics!(history::ParameterHistory)
    successful_results = [r for r in values(history.file_parameters) if r.converged]
    
    if isempty(successful_results)
        return
    end
    
    # Extract values
    mass_offsets = [r.mass_offset for r in successful_results]
    tolerances = [mean(r.mass_tolerance) for r in successful_results]
    
    # Compute statistics
    stats = history.global_stats
    stats.median_mass_offset = median(mass_offsets)
    stats.mass_offset_mad = mad(mass_offsets, normalize=false)
    stats.median_tolerance = median(tolerances)
    stats.tolerance_mad = mad(tolerances, normalize=false)
    stats.n_successful_files = length(successful_results)
end

"""
    get_initial_parameters(history::ParameterHistory, params::ParameterTuningSearchParameters)

Get initial parameters for a new file, potentially informed by previous runs.
"""
function get_initial_parameters(
    history::ParameterHistory,
    params::ParameterTuningSearchParameters
)
    # First file or no successful files: use config defaults
    if history.global_stats.n_successful_files == 0
        return InitialParameters(
            0.0f0,  # No bias estimate
            50.0f0,  # Default search range
            getFragTolPpm(params),  # From config
            false  # Not informed by history
        )
    end
    
    stats = history.global_stats
    
    # Start with median values from previous successful runs
    initial_bias = stats.median_mass_offset
    
    # Add safety margin based on observed variability
    bias_uncertainty = 2.0f0 * stats.mass_offset_mad
    tolerance_buffer = 1.5f0 * stats.median_tolerance + 2.0f0 * stats.tolerance_mad
    
    # For bias search, expand range based on observed variability
    bias_search_range = max(50.0f0, 3.0f0 * stats.mass_offset_mad)
    
    # Make sure we don't go below configured minimum
    initial_tolerance = max(tolerance_buffer, getFragTolPpm(params))
    
    return InitialParameters(
        initial_bias,
        bias_search_range,
        initial_tolerance,
        true  # Informed by history
    )
end

"""
    SearchStrategy

Abstract type for different search strategies.
"""
abstract type SearchStrategy end

struct FocusedSearchStrategy <: SearchStrategy end
struct BroadSearchStrategy <: SearchStrategy end

"""
    select_search_strategy(history::ParameterHistory)

Select appropriate search strategy based on parameter variability.
"""
function select_search_strategy(history::ParameterHistory)
    if history.global_stats.n_successful_files < 2
        # Not enough data, use broad search
        return BroadSearchStrategy()
    end
    
    # If all previous files had similar parameters, use focused search
    if history.global_stats.mass_offset_mad < 5.0f0 &&
       history.global_stats.tolerance_mad < 3.0f0
        return FocusedSearchStrategy()
    else
        return BroadSearchStrategy()
    end
end

"""
    generate_cross_run_report(history::ParameterHistory, output_path::String)

Generate a report of cross-run parameter statistics.
"""
function generate_cross_run_report(history::ParameterHistory, output_path::String)
    open(output_path, "w") do io
        println(io, "# Cross-Run Parameter Statistics")
        println(io)
        
        stats = history.global_stats
        println(io, "## Summary")
        println(io, "- Files processed: $(length(history.file_parameters))")
        println(io, "- Successful tuning: $(stats.n_successful_files)")
        println(io)
        
        if stats.n_successful_files > 0
            println(io, "## Mass Calibration")
            println(io, "- Median offset: $(round(stats.median_mass_offset, digits=2)) ppm")
            println(io, "- MAD: $(round(stats.mass_offset_mad, digits=2)) ppm")
            println(io)
            
            println(io, "## Mass Tolerance")
            println(io, "- Median: $(round(stats.median_tolerance, digits=2)) ppm")
            println(io, "- MAD: $(round(stats.tolerance_mad, digits=2)) ppm")
            println(io)
            
            # Recommend strategy
            strategy = select_search_strategy(history)
            if isa(strategy, FocusedSearchStrategy)
                println(io, "## Recommendation")
                println(io, "Parameters are stable across runs. Using focused search strategy.")
            else
                println(io, "## Recommendation")
                println(io, "Parameters vary across runs. Using broad search strategy.")
            end
        end
    end
end

"""
    should_use_informed_parameters(history::ParameterHistory, min_files::Int = 2)

Determine if we have enough history to inform parameter selection.
"""
function should_use_informed_parameters(history::ParameterHistory, min_files::Int = 2)
    return history.global_stats.n_successful_files >= min_files
end