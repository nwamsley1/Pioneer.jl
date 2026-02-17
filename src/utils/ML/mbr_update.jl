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
MBR feature update strategies.
Extracted from percolatorSortOf.jl lines 839-970.
Uses AbstractPSMContainer for data access.
=#

"""
    initialize_mbr_columns!(psms::AbstractPSMContainer, strategy::MBRUpdateStrategy) -> Nothing

Initialize MBR-related columns in the PSM container.
"""
function initialize_mbr_columns!(psms::AbstractPSMContainer, ::PairBasedMBR)
    n = nrows(psms)
    set_column!(psms, :trace_prob, zeros(Float32, n))
    set_column!(psms, :q_value, zeros(Float64, n))
    set_column!(psms, :MBR_max_pair_prob, zeros(Float32, n))
    set_column!(psms, :MBR_best_irt_diff, zeros(Float32, n))
    set_column!(psms, :MBR_log2_weight_ratio, zeros(Float32, n))
    set_column!(psms, :MBR_log2_explained_ratio, zeros(Float32, n))
    set_column!(psms, :MBR_rv_coefficient, zeros(Float32, n))
    set_column!(psms, :MBR_is_best_decoy, trues(n))
    set_column!(psms, :MBR_num_runs, zeros(Int32, n))
    set_column!(psms, :MBR_transfer_candidate, falses(n))
    set_column!(psms, :MBR_is_missing, falses(n))
    return nothing
end

function initialize_mbr_columns!(psms::AbstractPSMContainer, ::NoMBR)
    n = nrows(psms)
    set_column!(psms, :trace_prob, zeros(Float32, n))
    set_column!(psms, :q_value, zeros(Float64, n))
    return nothing
end

# DataFrame convenience wrapper
function initialize_mbr_columns!(psms::DataFrame, strategy::MBRUpdateStrategy)
    container = DataFramePSMContainer(psms, Val(:unsafe))
    initialize_mbr_columns!(container, strategy)
    return nothing
end

"""
    update_mbr_features!(psms_train, psms_test, iteration, mbr_start_iter, strategy) -> Nothing

Update MBR features based on the strategy.
"""
function update_mbr_features!(psms_train::AbstractPSMContainer, psms_test::AbstractPSMContainer,
                               iteration::Int, mbr_start_iter::Int, strategy::PairBasedMBR)
    if iteration >= mbr_start_iter - 1
        # Compute q-values on test set
        trace_probs = collect(Float32, get_column(psms_test, :trace_prob))
        targets = collect(Bool, get_column(psms_test, :target))
        q_values = Vector{Float64}(undef, length(trace_probs))
        get_qvalues!(trace_probs, targets, q_values)
        set_column!(psms_test, :q_value, q_values)

        # summarize_precursors! needs DataFrame - convert temporarily
        summarize_precursors!(to_dataframe(psms_test), q_cutoff=strategy.q_cutoff)
        summarize_precursors!(to_dataframe(psms_train), q_cutoff=strategy.q_cutoff)
    end
    return nothing
end

function update_mbr_features!(::AbstractPSMContainer, ::AbstractPSMContainer,
                               ::Int, ::Int, ::NoMBR)
    return nothing
end

# DataFrame convenience wrapper
function update_mbr_features!(psms_train::DataFrame, psms_test::DataFrame,
                               iteration::Int, mbr_start_iter::Int, strategy::MBRUpdateStrategy)
    train_container = DataFramePSMContainer(psms_train, Val(:unsafe))
    test_container = DataFramePSMContainer(psms_test, Val(:unsafe))
    update_mbr_features!(train_container, test_container, iteration, mbr_start_iter, strategy)
    return nothing
end

"""
    update_mbr_features_train_only!(psms_train, iteration, mbr_start_iter, strategy) -> Nothing

Update MBR features on training data only (used in Phase 1 of separated training/prediction).
"""
function update_mbr_features_train_only!(psms_train::AbstractPSMContainer,
                                          iteration::Int, mbr_start_iter::Int, strategy::PairBasedMBR)
    if iteration >= mbr_start_iter - 1
        # summarize_precursors! needs DataFrame - convert temporarily
        summarize_precursors!(to_dataframe(psms_train), q_cutoff=strategy.q_cutoff)
    end
    return nothing
end

function update_mbr_features_train_only!(::AbstractPSMContainer, ::Int, ::Int, ::NoMBR)
    return nothing
end

"""
    update_mbr_features_test_only!(psms_test, iteration, mbr_start_iter, strategy) -> Nothing

Update MBR features on test data only (used in Phase 2 of separated training/prediction).
"""
function update_mbr_features_test_only!(psms_test::AbstractPSMContainer,
                                         iteration::Int, mbr_start_iter::Int, strategy::PairBasedMBR)
    if iteration >= mbr_start_iter - 1
        # Compute q-values on test set
        trace_probs = collect(Float32, get_column(psms_test, :trace_prob))
        targets = collect(Bool, get_column(psms_test, :target))
        q_values = Vector{Float64}(undef, length(trace_probs))
        get_qvalues!(trace_probs, targets, q_values)
        set_column!(psms_test, :q_value, q_values)

        # summarize_precursors! needs DataFrame - convert temporarily
        summarize_precursors!(to_dataframe(psms_test), q_cutoff=strategy.q_cutoff)
    end
    return nothing
end

function update_mbr_features_test_only!(::AbstractPSMContainer, ::Int, ::Int, ::NoMBR)
    return nothing
end

# Note: summarize_precursors! stays in percolatorSortOf.jl - operates on DataFrame
# due to complex groupby operations that require DataFrame

#############################################################################
# MBR Utility Functions (used by both in-memory and OOM MBR paths)
#############################################################################

@inline function irt_residual(psms::AbstractDataFrame, idx::Integer)
    return Float32(psms.irt_pred[idx] - psms.irt_obs[idx])
end

function MBR_rv_coefficient(weights_A::AbstractVector{<:Real},
    times_A::AbstractVector{<:Real},
    weights_B::AbstractVector{<:Real},
    times_B::AbstractVector{<:Real})

    # Construct two Nx2 matrices, each row is (weight, time)
    X = hcat(collect(weights_A), collect(times_A))
    Y = hcat(collect(weights_B), collect(times_B))

    # Compute cross-products (Gram matrices)
    Sx = X' * X
    Sy = Y' * Y

    # Numerator: trace(Sx * Sy)
    numerator = tr(Sx * Sy)

    # Denominator: sqrt( trace(Sx*Sx)* trace(Sy*Sy) )
    denominator = sqrt(tr(Sx * Sx) * tr(Sy * Sy))

    # Protect against zero in denominator (e.g. if X or Y is all zeros)
    if denominator == 0
        return 0.0
    end

    return numerator / denominator
end

function pad_equal_length(x::AbstractVector, y::AbstractVector)

    function pad(data, d)
        left  = div(d, 2)
        right = d - left        # put extra on the right if d is odd
        padded_data = vcat(zeros(eltype(data), left), data, zeros(eltype(data), right))
        return padded_data
    end

    d = length(y) - length(x)
    if d == 0
        # No padding needed
        return (x, y)
    elseif d > 0
        # Pad x
        return (pad(x, d), y)
    else
        # Pad y
        return (x, pad(y, abs(d)))
    end
end


function pad_rt_equal_length(x::AbstractVector, y::AbstractVector)

    function pad(data, n, d, rt_step)
        left  = div(d, 2)
        right = d - left        # put extra on the right if d is odd
        padded_data = vcat(zeros(eltype(data), left), data, zeros(eltype(data), right))
        @inbounds @fastmath for i in range(1, left)
            padded_data[i] = padded_data[left+1] - (rt_step * ((left+1) - i))
        end
        @inbounds @fastmath for i in range(1, right)
            padded_data[i+left+n] = padded_data[left+n] + (rt_step * i)
        end
        return padded_data
    end

    nx = length(x)
    ny = length(y)
    d = ny - nx
    rt_step_x = (x[end] - x[1]) / nx
    rt_step_y = (y[end] - y[1]) / ny
    rt_step = nx > 1 ? rt_step_x : rt_step_y

    if d == 0
        # No padding needed
        return (x, y)
    elseif d > 0
        # Pad x
        return (pad(x, nx, d, rt_step), y)
    else
        # Pad y
        return (x, pad(y, ny, abs(d), rt_step))
    end
end
