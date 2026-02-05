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

# Note: summarize_precursors! stays in percolatorSortOf.jl - operates on DataFrame
# due to complex groupby operations that require DataFrame
