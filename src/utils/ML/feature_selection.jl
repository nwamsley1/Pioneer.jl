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

#=
Feature selection strategies.
Uses AbstractPSMContainer for data access.
=#

"""
    get_features(strategy::FeatureSelectionStrategy, iteration::Int, total_iterations::Int) -> Vector{Symbol}

Get the features to use for the current iteration.
"""
function get_features(strategy::IterativeFeatureSelection, iteration::Int, ::Int)
    if iteration < strategy.mbr_start_iteration
        return strategy.base_features
    else
        return vcat(strategy.base_features, strategy.mbr_features)
    end
end

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
