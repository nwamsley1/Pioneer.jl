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
Type definitions for ParameterTuningSearch module.
This file is loaded first to resolve circular dependencies.
"""

# From diagnostics.jl
"""
    ParameterTuningStatus

Simple structure to track parameter tuning outcomes for each file.
"""
struct ParameterTuningStatus
    file_idx::Int64
    file_name::String
    converged::Bool
    used_fallback::Bool
    fallback_reason::String
    n_iterations::Int
    final_psm_count::Int
    final_mass_offset::Float32
    final_mass_tolerance::Tuple{Float32, Float32}
    warnings::Vector{String}
end

"""
    ParameterTuningDiagnostics

Container for parameter tuning diagnostics across all files.
"""
mutable struct ParameterTuningDiagnostics
    file_statuses::Dict{Int64, ParameterTuningStatus}
    n_successful::Int
    n_fallback::Int
    n_failed::Int
    
    function ParameterTuningDiagnostics()
        new(Dict{Int64, ParameterTuningStatus}(), 0, 0, 0)
    end
end

# From cross_run_learning.jl
"""
    TuningResults

Results from parameter tuning for a single file.
"""
struct TuningResults
    mass_offset::Float32
    mass_tolerance::Tuple{Float32, Float32}
    converged::Bool
    psm_count::Int
    iterations::Int
    warnings::Vector{String}
end

"""
    GlobalParameterStats

Statistics computed across all successfully tuned files.
"""
mutable struct GlobalParameterStats
    median_mass_offset::Float32
    mass_offset_mad::Float32
    median_tolerance::Float32
    tolerance_mad::Float32
    n_successful_files::Int
    
    function GlobalParameterStats()
        new(0.0f0, 0.0f0, 0.0f0, 0.0f0, 0)
    end
end

"""
    ParameterHistory

Container for tracking parameter tuning results across files.
"""
mutable struct ParameterHistory
    file_parameters::Dict{Int64, TuningResults}
    global_stats::GlobalParameterStats
    
    function ParameterHistory()
        new(Dict{Int64, TuningResults}(), GlobalParameterStats())
    end
end

# Iteration configuration
"""
    IterationSettings

Configuration for iteration behavior in parameter tuning.
Controls initial mass tolerance, how it scales, and iterations per phase.
Always uses 3 phases: zero bias, positive shift, negative shift.
"""
struct IterationSettings
    init_mass_tol_ppm::Float32               # Initial mass tolerance (moved from fragment_settings)
    mass_tolerance_scale_factor::Float32     # Factor to scale tolerance each iteration
    iterations_per_phase::Int64              # Iterations before phase reset
    
    function IterationSettings(
        init_tol::Float32 = 20.0f0,
        scale_factor::Float32 = 2.0f0,
        iterations_per_phase::Int64 = 3
    )
        @assert init_tol > 0.0f0 "Initial tolerance must be positive"
        @assert scale_factor > 1.0f0 "Scale factor must be greater than 1"
        @assert iterations_per_phase > 0 "Iterations per phase must be positive"
        
        new(init_tol, scale_factor, iterations_per_phase)
    end
end

"""
    IterationState

Tracks state across phases and iterations during parameter tuning.
"""
mutable struct IterationState
    current_phase::Int64
    current_iteration_in_phase::Int64
    total_iterations::Int64
    phase_bias_shifts::Vector{Float32}  # Bias shift applied at start of each phase
    converged::Bool
    
    function IterationState()
        new(1, 0, 0, Float32[], false)
    end
end

# From ParameterTuningSearch.jl
struct ParameterTuningSearch <: TuningMethod end

"""
Results container for parameter tuning search.
Holds mass error models, RT alignment models, associated data, and diagnostics.
"""
struct ParameterTuningSearchResults <: SearchResults 
    mass_err_model::Base.Ref{<:MassErrorModel}
    rt_to_irt_model::Base.Ref{RtConversionModel}
    irt::Vector{Float32}
    rt::Vector{Float32}
    ppm_errs::Vector{Float32}
    qc_plots_folder_path::String
    diagnostics::ParameterTuningDiagnostics
    parameter_history::ParameterHistory
end

"""
Parameters for parameter tuning search.
Configures fragment matching, RT alignment, and general search behavior.
"""
struct ParameterTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # Core parameters from the original struct
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_fraction_transmitted::Float32
    # frag_tol_ppm moved to iteration_settings as init_mass_tol_ppm
    frag_err_quantile::Float32
    min_psms::Int64
    max_q_val::Float32
    max_presearch_iters::Int64
    min_index_search_score::UInt8
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_frags_for_mass_err_estimation::UInt8  # Number of top fragments per precursor for mass error estimation
    n_frag_isotopes::Int64
    intensity_filter_quantile::Float32  # Quantile for filtering fragments by intensity in mass error model
    max_frag_rank::UInt8
    sample_rate::Float32  # Deprecated - kept for compatibility
    topn_peaks::Union{Nothing, Int64}
    initial_scan_count::Int64
    max_parameter_tuning_scans::Int64
    # max_tol_ppm removed - calculated dynamically from iteration settings
    irt_tol::Float32
    spec_order::Set{Int64}
    relative_improvement_threshold::Float32
    spline_degree::Int64
    spline_n_knots::Int64
    spline_fit_outlier_sd::Int64
    irt_tol_sd::Int64
    prec_estimation::P
    iteration_settings::IterationSettings

    function ParameterTuningSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        global_params = params.global_settings
        tuning_params = params.parameter_tuning
        frag_params = tuning_params.fragment_settings
        search_params = tuning_params.search_settings
        rt_params = params.rt_alignment
        
        # Convert isotope error bounds
        isotope_bounds = global_params.isotope_settings.err_bounds_first_pass
        
        # Create precursor estimation type
        prec_estimation = global_params.isotope_settings.partial_capture ? PartialPrecCapture() : FullPrecCapture()
        
        # Extract intensity filter quantile for mass error model fitting (from fragment settings)
        intensity_filter_quantile = hasproperty(frag_params, :intensity_filter_quantile) ? 
            Float32(frag_params.intensity_filter_quantile) : 
            Float32(0.25)  # Default to 0.25 (25th percentile)
        
        # Extract topn_peaks if present
        topn_peaks = hasproperty(search_params, :topn_peaks) ? Int64(search_params.topn_peaks) : nothing
        @info "Using topn_peaks: $(topn_peaks === nothing ? "no topn peaks found" : topn_peaks)"
        # Extract scan count parameters
        initial_scan_count = hasproperty(search_params, :initial_scan_count) ? Int64(search_params.initial_scan_count) : Int64(2500)
        
        # Check for new name first, then fall back to old name for backward compatibility
        max_parameter_tuning_scans = if hasproperty(search_params, :max_parameter_tuning_scans)
            Int64(search_params.max_parameter_tuning_scans)
        elseif hasproperty(search_params, :expanded_scan_count)
            Int64(search_params.expanded_scan_count)
        else
            Int64(8000)  # Updated default from 10000 to 8000
        end
        
        # Extract max fragments for mass error estimation
        max_frags_for_mass_err_estimation = hasproperty(search_params, :max_frags_for_mass_err_estimation) ? 
            UInt8(search_params.max_frags_for_mass_err_estimation) : UInt8(5)
        
        # Extract max q-value for parameter tuning (with fallback to global threshold)
        max_q_value = hasproperty(search_params, :max_q_value) ? 
            Float32(search_params.max_q_value) : 
            Float32(global_params.scoring.q_value_threshold)
        
        # Extract iteration settings (now includes init_mass_tol_ppm)
        # Get initial tolerance from fragment_settings for backward compatibility
        init_tol_from_frag = hasproperty(frag_params, :tol_ppm) ? 
            Float32(frag_params.tol_ppm) : 20.0f0
            
        iteration_settings = if hasproperty(search_params, :iteration_settings)
            iter = search_params.iteration_settings
            IterationSettings(
                hasproperty(iter, :init_mass_tol_ppm) ? 
                    Float32(iter.init_mass_tol_ppm) : init_tol_from_frag,
                hasproperty(iter, :mass_tolerance_scale_factor) ? 
                    Float32(iter.mass_tolerance_scale_factor) : 2.0f0,
                hasproperty(iter, :iterations_per_phase) ? 
                    Int64(iter.iterations_per_phase) : 3
            )
        else
            # Use defaults for backward compatibility
            IterationSettings(init_tol_from_frag, 2.0f0, 3)
        end
        
        # Construct with appropriate type conversions
        new{typeof(prec_estimation)}(
            # Core parameters
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            Float32(global_params.isotope_settings.min_fraction_transmitted),
            # frag_tol_ppm removed - now in iteration_settings
            Float32(search_params.frag_err_quantile),
            Int64(search_params.min_samples),
            max_q_value,  # Use extracted parameter-tuning-specific q-value
            Int64(search_params.max_presearch_iters),
            UInt8(frag_params.min_score),
            Int64(frag_params.min_count),
            Float32(frag_params.min_spectral_contrast),
            Float32(frag_params.min_log2_ratio),
            (Int64(first(frag_params.min_top_n)), Int64(last(frag_params.min_top_n))),
            max_frags_for_mass_err_estimation,  # Extracted from JSON or defaults to 5
            Int64(frag_params.n_isotopes),
            intensity_filter_quantile,  # Fragment intensity filtering threshold
            UInt8(frag_params.max_rank),
            Float32(search_params.sample_rate),  # Deprecated
            topn_peaks,
            initial_scan_count,
            max_parameter_tuning_scans,
            # max_tol_ppm removed - calculated dynamically
            typemax(Float32), # irt_tol default
            Set{Int64}([2]), # spec_order default
            Float32(frag_params.relative_improvement_threshold),
            3,  # spline_degree default
            5,  # spline_n_knots default
            5,  # spline_fit_outlier_sd default
            Int64(rt_params.sigma_tolerance),
            prec_estimation,
            iteration_settings
        )
    end
end

# Accessor functions for ParameterTuningSearchParameters
# Calculate max tolerance dynamically based on iteration settings
function getMaxTolerancePpm(params::ParameterTuningSearchParameters)
    settings = params.iteration_settings
    # Max tolerance is reached at the end of a phase
    return settings.init_mass_tol_ppm * (settings.mass_tolerance_scale_factor ^ settings.iterations_per_phase)
end
getInitialScanCount(params::ParameterTuningSearchParameters) = params.initial_scan_count
getMaxParameterTuningScans(params::ParameterTuningSearchParameters) = params.max_parameter_tuning_scans
getTopNPeaks(params::ParameterTuningSearchParameters) = params.topn_peaks
getMaxFragsForMassErrEstimation(params::ParameterTuningSearchParameters) = params.max_frags_for_mass_err_estimation
getIntensityFilterQuantile(params::ParameterTuningSearchParameters) = params.intensity_filter_quantile
getIterationSettings(params::ParameterTuningSearchParameters) = params.iteration_settings
# Get initial mass tolerance from iteration settings
getFragTolPpm(params::ParameterTuningSearchParameters) = params.iteration_settings.init_mass_tol_ppm

# Override getMaxBestRank for ParameterTuningSearchParameters since it doesn't have max_best_rank field
# This is for PSM filtering in LibrarySearch, not mass error estimation
import Pioneer: getMaxBestRank
getMaxBestRank(params::ParameterTuningSearchParameters) = UInt8(1)  # Default value for PSM filtering