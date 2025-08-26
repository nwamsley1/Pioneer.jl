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
    BiasSearchStrategy

Configuration for mass bias detection during parameter tuning.
"""
struct BiasSearchStrategy
    initial_tolerance::Float32  # From config (e.g., 20 ppm)
    bias_search_range::Float32  # How far to search for bias (e.g., ±50 ppm)
    bias_search_steps::Vector{Float32}  # Steps to try: [0, -10, 10, -20, 20, ...]
    quick_search_sample_rate::Float32  # Reduced sampling for bias detection (e.g., 0.005)
    min_psms_for_detection::Int64  # Minimum PSMs to consider bias estimate valid
end

"""
    get_default_bias_search_strategy(params::ParameterTuningSearchParameters)

Create a default bias search strategy based on the search parameters.
"""
function get_default_bias_search_strategy(params::ParameterTuningSearchParameters)
    # Create search steps from 0 to ±50 ppm in 10 ppm increments
    search_steps = Float32[0.0]
    for offset in 10.0:10.0:50.0
        push!(search_steps, -Float32(offset))
        push!(search_steps, Float32(offset))
    end
    
    return BiasSearchStrategy(
        getFragTolPpm(params),
        50.0f0,  # Search up to ±50 ppm
        search_steps,
        0.005f0,  # 0.5% sampling for quick search
        100  # Need at least 100 PSMs to trust bias estimate
    )
end

"""
    detect_mass_bias(spectra::MassSpecData, search_context::SearchContext, 
                    params::ParameterTuningSearchParameters, ms_file_idx::Int64,
                    strategy::BiasSearchStrategy)

Detect systematic mass bias by testing multiple offset values.
Returns (best_bias, best_psm_count).
"""
function detect_mass_bias(
    spectra::MassSpecData, 
    search_context::SearchContext, 
    params::ParameterTuningSearchParameters, 
    ms_file_idx::Int64,
    strategy::BiasSearchStrategy = get_default_bias_search_strategy(params)
)
    best_psm_count = 0
    best_bias = 0.0f0
    
    # Create a modified params with reduced sample rate for quick search
    quick_params = create_quick_search_params(params, strategy.quick_search_sample_rate)
    
    for bias_guess in strategy.bias_search_steps
        # Set temporary mass model with bias offset
        temp_mass_model = MassErrorModel(
            bias_guess,
            (strategy.initial_tolerance, strategy.initial_tolerance)
        )
        setMassErrorModel!(search_context, ms_file_idx, temp_mass_model)
        
        # Do limited search to estimate PSM yield
        psm_count = estimate_psm_yield(spectra, search_context, quick_params, ms_file_idx)
        
        @debug "Bias $bias_guess ppm yielded $psm_count PSMs"
        
        if psm_count > best_psm_count
            best_psm_count = psm_count
            best_bias = bias_guess
        end
        
        # Early termination if found good yield
        if psm_count > strategy.min_psms_for_detection * 2
            break
        end
    end
    
    if best_psm_count < strategy.min_psms_for_detection
        @user_warn "Bias detection yielded only $best_psm_count PSMs. May be unreliable."
    end
    
    return best_bias, best_psm_count
end

"""
    estimate_psm_yield(spectra::MassSpecData, search_context::SearchContext,
                      params::ParameterTuningSearchParameters, ms_file_idx::Int64)

Quick PSM count estimation using reduced sampling.
"""
function estimate_psm_yield(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::ParameterTuningSearchParameters,
    ms_file_idx::Int64
)
    # Perform a quick library search with current mass error model
    psms = library_search(spectra, search_context, params, ms_file_idx)
    
    # Return raw count without scoring (faster for estimation)
    return size(psms, 1)
end

"""
    create_quick_search_params(params::ParameterTuningSearchParameters, sample_rate::Float32)

Create modified parameters with reduced sample rate for quick searching.
"""
function create_quick_search_params(params::ParameterTuningSearchParameters, sample_rate::Float32)
    # This is a simplified approach - in practice you might want to create a proper copy
    # For now, we'll just note that the actual implementation would modify the sample_rate field
    # The existing ParameterTuningSearchParameters struct would need a way to override sample_rate
    return params  # Placeholder - actual implementation would create modified copy
end

"""
    adjust_search_window_for_bias(initial_tolerance::Float32, detected_bias::Float32)

Adjust the mass search window based on detected bias to ensure adequate coverage.
Returns (left_tolerance, right_tolerance).
"""
function adjust_search_window_for_bias(initial_tolerance::Float32, detected_bias::Float32)
    # If bias is large relative to tolerance, expand asymmetrically
    if abs(detected_bias) > initial_tolerance * 0.5
        if detected_bias > 0
            # Positive bias: expand right side more
            left_tol = initial_tolerance
            right_tol = initial_tolerance + abs(detected_bias) * 0.5f0
        else
            # Negative bias: expand left side more
            left_tol = initial_tolerance + abs(detected_bias) * 0.5f0
            right_tol = initial_tolerance
        end
    else
        # Small bias: symmetric expansion
        left_tol = right_tol = initial_tolerance
    end
    
    return (left_tol, right_tol)
end

"""
    should_attempt_bias_detection(params::ParameterTuningSearchParameters, 
                                 iteration::Int, previous_failures::Int)

Determine if bias detection should be attempted based on current state.
"""
function should_attempt_bias_detection(
    params::ParameterTuningSearchParameters,
    iteration::Int,
    previous_failures::Int
)
    # Attempt bias detection on first iteration or after repeated failures
    return iteration == 1 || previous_failures >= 2
end