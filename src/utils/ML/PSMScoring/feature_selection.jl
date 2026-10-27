# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

#=
Feature selection strategies.
Uses AbstractPSMContainer for data access.
=#

"""
    get_features(strategy::FeatureSelectionStrategy, iteration::Int, total_iterations::Int) -> Vector{Symbol}

Get the features to use for the current iteration.
"""
function get_features(strategy::StaticFeatureSelection, ::Int, ::Int)
    return strategy.features
end

"""
    filter_available_features(features::Vector{Symbol}, psms::AbstractPSMContainer) -> Vector{Symbol}

Filter features to only those available in the container.
"""
function filter_available_features(features::Vector{Symbol}, psms::AbstractPSMContainer)
    return [f for f in features if has_column(psms, f)]
end

# DataFrame convenience wrapper
function filter_available_features(features::Vector{Symbol}, psms::DataFrame)
    return [f for f in features if hasproperty(psms, f)]
end
