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
    frag_tol_ppm::Float32
    frag_err_quantile::Float32
    min_psms::Int64
    max_q_val::Float32
    max_presearch_iters::Int64
    min_index_search_score::UInt8
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::UInt8
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    sample_rate::Float32  # Deprecated - kept for compatibility
    topn_peaks::Union{Nothing, Int64}
    initial_scan_count::Int64
    expanded_scan_count::Int64
    max_tolerance_ppm::Float32
    irt_tol::Float32
    spec_order::Set{Int64}
    relative_improvement_threshold::Float32
    spline_degree::Int64
    spline_n_knots::Int64
    spline_fit_outlier_sd::Int64
    irt_tol_sd::Int64
    prec_estimation::P

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
        
        # Extract topn_peaks if present
        topn_peaks = hasproperty(search_params, :topn_peaks) ? Int64(search_params.topn_peaks) : nothing
        @info "Using topn_peaks: $(topn_peaks === nothing ? "no topn peaks found" : topn_peaks)"
        # Extract scan count parameters
        initial_scan_count = hasproperty(search_params, :initial_scan_count) ? Int64(search_params.initial_scan_count) : Int64(2500)
        
        expanded_scan_count = hasproperty(search_params, :expanded_scan_count) ? Int64(search_params.expanded_scan_count) : Int64(10000)
        
        # Extract max tolerance parameter
        max_tolerance_ppm = hasproperty(search_params, :max_tolerance_ppm) ? Float32(search_params.max_tolerance_ppm) : Float32(50.0)
        
        # Construct with appropriate type conversions
        new{typeof(prec_estimation)}(
            # Core parameters
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            Float32(global_params.isotope_settings.min_fraction_transmitted),
            Float32(frag_params.tol_ppm),
            Float32(search_params.frag_err_quantile),
            Int64(search_params.min_samples),
            Float32(global_params.scoring.q_value_threshold),
            Int64(search_params.max_presearch_iters),
            UInt8(frag_params.min_score),
            Int64(frag_params.min_count),
            Float32(frag_params.min_spectral_contrast),
            Float32(frag_params.min_log2_ratio),
            (Int64(first(frag_params.min_top_n)), Int64(last(frag_params.min_top_n))),
            UInt8(5), # max_best_rank default - top fragments for mass error estimation
            Int64(frag_params.n_isotopes),
            UInt8(frag_params.max_rank),
            Float32(search_params.sample_rate),  # Deprecated
            topn_peaks,
            initial_scan_count,
            expanded_scan_count,
            max_tolerance_ppm,
            typemax(Float32), # irt_tol default
            Set{Int64}([2]), # spec_order default
            Float32(frag_params.relative_improvement_threshold),
            3,  # spline_degree default
            5,  # spline_n_knots default
            5,  # spline_fit_outlier_sd default
            Int64(rt_params.sigma_tolerance),
            prec_estimation
        )
    end
end

# Accessor functions for ParameterTuningSearchParameters
getMaxTolerancePpm(params::ParameterTuningSearchParameters) = params.max_tolerance_ppm
getInitialScanCount(params::ParameterTuningSearchParameters) = params.initial_scan_count
getExpandedScanCount(params::ParameterTuningSearchParameters) = params.expanded_scan_count
getTopNPeaks(params::ParameterTuningSearchParameters) = params.topn_peaks