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
Training data selection strategies.
Extracted from percolatorSortOf.jl lines 972-1009.
Uses AbstractPSMContainer for data access.
=#

"""
    select_training_data(psms::AbstractPSMContainer, strategy::TrainingDataStrategy, iteration::Int) -> AbstractPSMContainer

Select training data based on the strategy and current iteration.
Returns a copy of the filtered training data.
"""
function select_training_data(psms::AbstractPSMContainer, ::AllDataSelection, ::Int)
    return copy_container(psms)
end

function select_training_data(psms::AbstractPSMContainer, strategy::QValueNegativeMining, iteration::Int)
    if iteration == 1
        return copy_container(psms)
    end

    # Deep copy to avoid modifying original labels
    psms_train = copy_container(psms)

    # Convert worst-scoring targets to negatives using PEP
    trace_probs = collect(Float32, get_column(psms_train, :trace_prob))
    targets = collect(Bool, get_column(psms_train, :target))

    order = sortperm(trace_probs, rev=true)
    sorted_scores = trace_probs[order]
    sorted_targets = targets[order]

    PEPs = Vector{Float32}(undef, length(order))
    get_PEP!(sorted_scores, sorted_targets, PEPs; doSort=false)

    idx_cutoff = findfirst(x -> x >= strategy.min_pep_threshold, PEPs)
    if !isnothing(idx_cutoff)
        worst_idxs = order[idx_cutoff:end]
        targets[worst_idxs] .= false
        set_column!(psms_train, :target, targets)
    end

    # Filter: keep all decoys + targets passing q-value threshold
    q_values = collect(get_column(psms_train, :q_value))
    targets = collect(get_column(psms_train, :target))
    keep_mask = BitVector([(!t) || (t && q <= strategy.max_q_value) for (t, q) in zip(targets, q_values)])

    return select_subset(psms_train, keep_mask)
end

# DataFrame convenience wrapper
function select_training_data(psms::DataFrame, strategy::TrainingDataStrategy, iteration::Int)
    container = DataFramePSMContainer(psms, Val(:unsafe))
    result = select_training_data(container, strategy, iteration)
    return to_dataframe(result)
end
