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
    BoundarySamplingCheck

Configuration for validating adequate sampling at tolerance boundaries.
"""
struct BoundarySamplingCheck
    boundary_fraction_threshold::Float32  # e.g., 0.05 (5%)
    boundary_zone_fraction::Float32       # e.g., 0.1 (10% of tolerance)
end

"""
    BoundarySamplingResult

Results from boundary sampling validation.
"""
struct BoundarySamplingResult
    adequate_sampling::Bool
    left_boundary_fraction::Float32
    right_boundary_fraction::Float32
    expansion_needed::Bool
    suggested_expansion_factor::Float32
    diagnostic_message::String
end

"""
    get_default_boundary_checker()

Create default boundary sampling checker with standard thresholds.
"""
function get_default_boundary_checker()
    return BoundarySamplingCheck(
        0.05f0,  # 5% threshold
        0.1f0    # 10% boundary zone
    )
end

"""
    check_boundary_sampling(ppm_errors::Vector{Float32}, 
                          mass_model::MassErrorModel,
                          checker::BoundarySamplingCheck = get_default_boundary_checker())

Validate that we have adequate sampling at the tolerance boundaries.
Returns BoundarySamplingResult.
"""
function check_boundary_sampling(
    ppm_errors::Vector{Float32}, 
    mass_model::MassErrorModel,
    checker::BoundarySamplingCheck = get_default_boundary_checker()
)
    if isempty(ppm_errors)
        return BoundarySamplingResult(
            false, 0.0f0, 0.0f0, true, 2.0f0,
            "No PPM errors to analyze"
        )
    end
    
    # Remove bias to analyze tolerance boundaries
    mass_offset = getMassOffset(mass_model)
    centered_errors = ppm_errors .- mass_offset
    left_tol, right_tol = getLeftTol(mass_model), getRightTol(mass_model)
    
    # Define boundary zones (outer 10% of tolerance range)
    left_boundary = -left_tol * (1 - checker.boundary_zone_fraction)
    right_boundary = right_tol * (1 - checker.boundary_zone_fraction)
    
    # Count matches in boundary zones
    n_left_boundary = count(err -> err < left_boundary, centered_errors)
    n_right_boundary = count(err -> err > right_boundary, centered_errors)
    n_total = length(centered_errors)
    
    # Calculate fractions
    left_fraction = n_left_boundary / n_total
    right_fraction = n_right_boundary / n_total
    
    # Check if we have adequate sampling
    adequate_left = left_fraction < checker.boundary_fraction_threshold
    adequate_right = right_fraction < checker.boundary_fraction_threshold
    
    # Determine expansion factor if needed
    suggested_expansion = 1.0f0
    if !adequate_left || !adequate_right
        # If more than threshold at boundaries, we need to expand
        max_boundary_fraction = max(left_fraction, right_fraction)
        # Scale expansion based on how much we exceed threshold
        suggested_expansion = 1.0f0 + (max_boundary_fraction / checker.boundary_fraction_threshold) * 0.5f0
        suggested_expansion = min(suggested_expansion, 2.0f0)  # Cap at 2x expansion
    end
    
    # Create diagnostic message
    diagnostic_msg = if adequate_left && adequate_right
        "Adequate boundary sampling: left=$(round(left_fraction*100, digits=1))%, right=$(round(right_fraction*100, digits=1))%"
    else
        issues = String[]
        !adequate_left && push!(issues, "left boundary $(round(left_fraction*100, digits=1))% > $(checker.boundary_fraction_threshold*100)%")
        !adequate_right && push!(issues, "right boundary $(round(right_fraction*100, digits=1))% > $(checker.boundary_fraction_threshold*100)%")
        "Inadequate sampling at: " * join(issues, ", ")
    end
    
    return BoundarySamplingResult(
        adequate_left && adequate_right,
        left_fraction,
        right_fraction,
        !adequate_left || !adequate_right,
        suggested_expansion,
        diagnostic_msg
    )
end

"""
    expand_tolerance!(mass_model::MassErrorModel, expansion_factor::Float32)

Expand the mass tolerance by the given factor.
Returns new MassErrorModel.
"""
function expand_tolerance(mass_model::MassErrorModel, expansion_factor::Float32)
    left_tol, right_tol = getLeftTol(mass_model), getRightTol(mass_model)
    
    return MassErrorModel(
        getMassOffset(mass_model),
        (left_tol * expansion_factor, right_tol * expansion_factor)
    )
end

"""
    validate_tolerance_bounds(ppm_errors::Vector{Float32}, 
                            mass_model::MassErrorModel,
                            min_coverage::Float32 = 0.95f0)

Validate that the tolerance bounds cover at least min_coverage fraction of the data.
"""
function validate_tolerance_bounds(
    ppm_errors::Vector{Float32}, 
    mass_model::MassErrorModel,
    min_coverage::Float32 = 0.95f0
)
    if isempty(ppm_errors)
        return false, "No data to validate"
    end
    
    # Remove bias
    mass_offset = getMassOffset(mass_model)
    centered_errors = ppm_errors .- mass_offset
    left_tol, right_tol = getLeftTol(mass_model), getRightTol(mass_model)
    
    # Count how many are within bounds
    n_within = count(err -> -left_tol <= err <= right_tol, centered_errors)
    coverage = n_within / length(centered_errors)
    
    return coverage >= min_coverage, "Coverage: $(round(coverage*100, digits=1))%"
end