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
PSMContainer abstraction for PSM data storage.
Provides a unified interface for different backing implementations.
=#

#############################################################################
# Abstract Interface
#############################################################################

"""
    AbstractPSMContainer

Abstract type for PSM data containers.
Implementations must provide access to PSM data with standard column operations.

# Required Interface
- `nrows(container)` - Number of PSMs
- `get_column(container, col::Symbol)` - Get column data
- `set_column!(container, col::Symbol, data)` - Set column data
- `has_column(container, col::Symbol)` - Check if column exists
- `get_view(container, indices)` - Get view of subset
- `copy_container(container)` - Create a copy
- `get_cv_folds(container)` - Get unique CV fold values
- `get_fold_indices(container, fold)` - Get indices for a CV fold
- `get_train_indices(container, test_fold)` - Get training indices (all except test_fold)
- `sort_container!(container, cols)` - Sort in place by columns
- `to_dataframe(container)` - Convert to DataFrame for ML training
"""
abstract type AbstractPSMContainer end

# Required columns for PSM scoring
const REQUIRED_PSM_COLUMNS = [
    :target,           # Bool - target/decoy label
    :cv_fold,          # UInt8 - cross-validation fold
    :precursor_idx,    # UInt32 - precursor identifier
    :ms_file_idx,      # UInt32 - MS file index
    :irt_pred,         # Float32 - predicted iRT
    :isotopes_captured # Tuple{Int8,Int8} - isotope info
]

# Optional columns for MBR
const MBR_PSM_COLUMNS = [
    :weights,                  # Vector{Float32} - fragment weights
    :irts,                     # Vector{Float32} - fragment iRTs
    :weight,                   # Float32 - precursor weight
    :log2_intensity_explained  # Float32 - MS1 intensity metric
]

#############################################################################
# DataFrame-backed Implementation
#############################################################################

"""
    DataFramePSMContainer <: AbstractPSMContainer

PSM container backed by a DataFrame.
This is the standard in-memory implementation.
"""
struct DataFramePSMContainer <: AbstractPSMContainer
    data::DataFrame
    validated::Bool

    function DataFramePSMContainer(df::DataFrame, validated::Bool=true)
        if validated
            # Validate required columns
            for col in REQUIRED_PSM_COLUMNS
                if !hasproperty(df, col)
                    error("Missing required column: $col")
                end
            end
        end
        new(df, validated)
    end
end

# Convenience constructor without validation (for internal use)
DataFramePSMContainer(df::DataFrame, ::Val{:unsafe}) = DataFramePSMContainer(df, false)

# Core interface implementation
nrows(container::DataFramePSMContainer) = nrow(container.data)

function get_column(container::DataFramePSMContainer, col::Symbol)
    return container.data[!, col]
end

function set_column!(container::DataFramePSMContainer, col::Symbol, data)
    container.data[!, col] = data
    return nothing
end

function has_column(container::DataFramePSMContainer, col::Symbol)
    return hasproperty(container.data, col)
end

function get_view(container::DataFramePSMContainer, indices::AbstractVector{<:Integer})
    return DataFramePSMContainerView(container, collect(Int, indices))
end

function copy_container(container::DataFramePSMContainer)
    return DataFramePSMContainer(copy(container.data), Val(:unsafe))
end

function get_cv_folds(container::DataFramePSMContainer)
    return unique(container.data[!, :cv_fold])
end

function get_fold_indices(container::DataFramePSMContainer, fold)
    cv_col = container.data[!, :cv_fold]
    return findall(==(fold), cv_col)
end

function get_train_indices(container::DataFramePSMContainer, test_fold)
    cv_col = container.data[!, :cv_fold]
    return findall(!=(test_fold), cv_col)
end

function sort_container!(container::DataFramePSMContainer, cols::Vector{Symbol})
    sort!(container.data, cols)
    return nothing
end

function to_dataframe(container::DataFramePSMContainer)
    return container.data
end

# Access underlying DataFrame (for backward compatibility)
function Base.getproperty(c::DataFramePSMContainer, s::Symbol)
    if s === :data
        return getfield(c, :data)
    elseif s === :validated
        return getfield(c, :validated)
    else
        return get_column(c, s)
    end
end

#############################################################################
# View Implementation (for train/test splits)
#############################################################################

"""
    DataFramePSMContainerView <: AbstractPSMContainer

A view into a DataFramePSMContainer for a subset of rows.
Used for train/test splits without copying data.
"""
struct DataFramePSMContainerView <: AbstractPSMContainer
    parent::DataFramePSMContainer
    indices::Vector{Int}
end

nrows(view::DataFramePSMContainerView) = length(view.indices)

function get_column(view::DataFramePSMContainerView, col::Symbol)
    return @view get_column(view.parent, col)[view.indices]
end

function set_column!(view::DataFramePSMContainerView, col::Symbol, data)
    parent_col = get_column(view.parent, col)
    parent_col[view.indices] = data
    return nothing
end

function has_column(view::DataFramePSMContainerView, col::Symbol)
    return has_column(view.parent, col)
end

function get_view(view::DataFramePSMContainerView, indices::AbstractVector{<:Integer})
    # Compose views
    return DataFramePSMContainerView(view.parent, view.indices[indices])
end

function copy_container(view::DataFramePSMContainerView)
    # Copy creates a new container with just the viewed rows
    return DataFramePSMContainer(copy(view.parent.data[view.indices, :]), Val(:unsafe))
end

function to_dataframe(view::DataFramePSMContainerView)
    return @view view.parent.data[view.indices, :]
end

function get_cv_folds(view::DataFramePSMContainerView)
    return unique(get_column(view, :cv_fold))
end

function get_fold_indices(view::DataFramePSMContainerView, fold)
    cv_col = get_column(view, :cv_fold)
    return findall(==(fold), collect(cv_col))
end

function get_train_indices(view::DataFramePSMContainerView, test_fold)
    cv_col = get_column(view, :cv_fold)
    return findall(!=(test_fold), collect(cv_col))
end

# Convenience: iterate over rows
Base.eachindex(c::AbstractPSMContainer) = 1:nrows(c)

#############################################################################
# Utility Functions
#############################################################################

"""
    select_subset(container::AbstractPSMContainer, mask::BitVector) -> AbstractPSMContainer

Select rows matching the mask.
Returns a copy with only matching rows.
"""
function select_subset(container::AbstractPSMContainer, mask::BitVector)
    indices = findall(mask)
    return copy_container(get_view(container, indices))
end

"""
    get_feature_matrix(container::AbstractPSMContainer, features::Vector{Symbol}) -> Matrix{Float32}

Extract features as a matrix for ML training.
Handles missing values by replacing with zeros.
"""
function get_feature_matrix(container::AbstractPSMContainer, features::Vector{Symbol})
    n = nrows(container)
    m = length(features)
    X = Matrix{Float32}(undef, n, m)
    for (j, feat) in enumerate(features)
        col = get_column(container, feat)
        for i in 1:n
            val = col[i]
            # Handle missing values
            if ismissing(val)
                X[i, j] = 0.0f0
            elseif val isa Bool
                X[i, j] = Float32(val)
            else
                X[i, j] = Float32(val)
            end
        end
    end
    return X
end

"""
    get_labels(container::AbstractPSMContainer) -> Vector{Bool}

Get target labels for ML training.
"""
function get_labels(container::AbstractPSMContainer)
    return collect(Bool, get_column(container, :target))
end
