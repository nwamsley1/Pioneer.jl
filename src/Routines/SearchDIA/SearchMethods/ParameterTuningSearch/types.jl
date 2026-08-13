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

"""
    IterationState

Tracks state during parameter tuning's scout-then-collect strategy.
"""
mutable struct IterationState
    converged::Bool
    scan_attempt::Int64
    failed_with_exception::Bool
    best_psm_count::Int64
    best_mass_error_model::Union{Nothing, AbstractMassErrorModel}
    best_psms::Union{Nothing, DataFrame}
    best_score::UInt8
    best_scan_count::Int64
    best_fragments::Union{Nothing, Vector{MassErrSample}}
    wide_scout_plot::Any

    function IterationState()
        new(false, 1, false,
            0, nothing, nothing, UInt8(0), 0,
            nothing, nothing)
    end
end

# From ParameterTuningSearch.jl
struct ParameterTuningSearch <: TuningMethod end

"""
Results container for parameter tuning search.
Holds mass error models, RT alignment models, associated data, and diagnostics.
"""
struct ParameterTuningSearchResults <: SearchResults 
    mass_err_model::Base.RefValue{AbstractMassErrorModel}
    rt_to_irt_model::Base.Ref{RtConversionModel}
    irt::Vector{Float32}
    rt::Vector{Float32}
    ppm_errs::Vector{Float32}
    frag_mzs::Vector{Float32}  # Fragment m/z values corresponding to ppm_errs
    rt_plots::Vector{Vector{UInt8}}    # PNG bytes (legacy)
    mass_plots::Vector{Vector{UInt8}}  # PNG bytes (legacy)
    rt_plot_objects::Vector{Any}       # Plot objects accumulated for combined PDF
    mass_plot_objects::Vector{Any}     # Plot objects accumulated for combined PDF
    nce_plot_objects::Vector{Any}      # Plot objects accumulated for combined PDF
    qc_plots_folder_path::String
    diagnostics::ParameterTuningDiagnostics
    parameter_history::ParameterHistory
    current_iteration_state::Base.Ref{Union{Nothing, IterationState}}  # Store iteration state for plot generation
end

"""
Parameters for parameter tuning search.
Configures fragment matching, RT alignment, and general search behavior.
"""
mutable struct ParameterTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_fraction_transmitted::Float32
    frag_err_quantile::Float32
    min_psms::Int64
    max_q_val::Float32
    max_presearch_iters::Int64
    min_index_search_scores::Vector{UInt8}
    current_min_score::UInt8
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_frags_for_mass_err_estimation::UInt8
    n_frag_isotopes::Int64
    intensity_filter_quantile::Float32
    max_frag_rank::UInt8
    irt_tol::Float32
    spec_order::Set{Int64}
    relative_improvement_threshold::Float32
    spline_degree::Int64
    spline_n_knots::Int64
    spline_fit_outlier_sd::Int64
    prec_estimation::P
    lambda_penalty::Float32
    ransac_threshold_psms::Int64
    min_psms_for_spline::Int64

    function ParameterTuningSearchParameters(params::PioneerParameters)
        # Hardcoded — formerly read from global.isotope_settings (now removed
        # from the public schema; both knobs were always at default values).
        prec_estimation = PartialPrecCapture()

        new{typeof(prec_estimation)}(
            (UInt8(1), UInt8(0)),                                          # isotope_err_bounds
            0.25f0,                                                        # min_fraction_transmitted
            Float32(TUNING_FRAG_ERR_QUANTILE),                             # frag_err_quantile
            TUNING_MIN_SAMPLES,                                            # min_psms
            TUNING_MAX_Q_VALUE,                                            # max_q_val
            TUNING_MAX_PRESEARCH_ITERS,                                    # max_presearch_iters
            Vector{UInt8}(collect(TUNING_SCORE_TIERS)),                    # min_index_search_scores
            UInt8(first(TUNING_SCORE_TIERS)),                              # current_min_score
            TUNING_MIN_SPECTRAL_CONTRAST,                                  # min_spectral_contrast
            TUNING_MIN_LOG2_RATIO,                                         # min_log2_matched_ratio
            (Int64(first(TUNING_MIN_TOPN_OF_M)), Int64(last(TUNING_MIN_TOPN_OF_M))), # min_topn_of_m
            TUNING_MAX_FRAGS_FOR_MASS_ERR,                                 # max_frags_for_mass_err_estimation
            TUNING_N_FRAG_ISOTOPES,                                        # n_frag_isotopes
            TUNING_INTENSITY_FILTER_QUANTILE,                              # intensity_filter_quantile
            PARAM_TUNING_MAX_FRAG_RANK,                                    # max_frag_rank
            typemax(Float32),                                              # irt_tol
            Set{Int64}([2]),                                               # spec_order
            TUNING_RELATIVE_IMPROVEMENT_THRESHOLD,                         # relative_improvement_threshold
            3,                                                             # spline_degree
            100,                                                           # spline_n_knots
            5,                                                             # spline_fit_outlier_sd
            prec_estimation,                                               # prec_estimation
            TUNING_LAMBDA_PENALTY,                                         # lambda_penalty
            TUNING_RANSAC_THRESHOLD_PSMS,                                  # ransac_threshold_psms
            TUNING_MIN_PSMS_FOR_SPLINE,                                    # min_psms_for_spline
        )
    end
end

# Accessor functions for ParameterTuningSearchParameters
getMaxFragsForMassErrEstimation(params::ParameterTuningSearchParameters) = params.max_frags_for_mass_err_estimation
getIntensityFilterQuantile(params::ParameterTuningSearchParameters) = params.intensity_filter_quantile

# Optional JSON override accessors

# RT alignment parameter accessors
getLambdaPenalty(params::ParameterTuningSearchParameters) = params.lambda_penalty
getRansacThresholdPsms(params::ParameterTuningSearchParameters) = params.ransac_threshold_psms
getMinPsmsForSpline(params::ParameterTuningSearchParameters) = params.min_psms_for_spline

# RT alignment parameter accessors (used by fit_irt_model in utils.jl)
getRtAlignmentLambdaPenalty(params::ParameterTuningSearchParameters) = params.lambda_penalty
getRtAlignmentRansacThreshold(params::ParameterTuningSearchParameters) = params.ransac_threshold_psms
getRtAlignmentMinPsms(params::ParameterTuningSearchParameters) = params.min_psms_for_spline

# Override getMaxBestRank for ParameterTuningSearchParameters since it doesn't have max_best_rank field
# This is for PSM filtering in LibrarySearch, not mass error estimation
import Pioneer: getMaxBestRank
getMaxBestRank(params::ParameterTuningSearchParameters) = UInt8(1)  # Default value for PSM filtering

# New getters and setters for multi-score support
getMinIndexSearchScores(params::ParameterTuningSearchParameters) = params.min_index_search_scores
getMinIndexSearchScore(params::ParameterTuningSearchParameters) = params.current_min_score
function setCurrentMinScore!(params::ParameterTuningSearchParameters, score::UInt8)
    params.current_min_score = score
end