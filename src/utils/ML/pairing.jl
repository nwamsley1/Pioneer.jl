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
Pairing strategy implementations.
Extracted from percolatorSortOf.jl lines 94-240.
Uses AbstractPSMContainer for data access.
=#

"""
    assign_pairs!(psms::AbstractPSMContainer, strategy::PairingStrategy) -> Nothing

Assign pair IDs to PSMs based on the pairing strategy.
Modifies `psms` in-place by adding/updating :pair_id and :irt_bin_idx columns.
"""
function assign_pairs!(psms::AbstractPSMContainer, strategy::RandomPairing)
    n = nrows(psms)
    set_column!(psms, :pair_id, zeros(UInt32, n))

    # Compute iRT bins
    irt_pred = collect(get_column(psms, :irt_pred))
    irt_bin_idx = getIrtBins(irt_pred, strategy.bin_size)
    set_column!(psms, :irt_bin_idx, irt_bin_idx)

    # Get grouping columns
    cv_fold = collect(get_column(psms, :cv_fold))
    isotopes = collect(get_column(psms, :isotopes_captured))
    precursor_idx = collect(get_column(psms, :precursor_idx))
    pair_ids = collect(get_column(psms, :pair_id))

    # Group by (irt_bin_idx, cv_fold, isotopes_captured)
    # Use a dictionary to track groups
    groups = Dict{Tuple{UInt32, UInt8, Tuple{Int8,Int8}}, Vector{Int}}()
    for i in 1:n
        key = (irt_bin_idx[i], cv_fold[i], isotopes[i])
        if !haskey(groups, key)
            groups[key] = Int[]
        end
        push!(groups[key], i)
    end

    # Assign pair IDs within each group
    last_pair_id = zero(UInt32)
    for indices in values(groups)
        last_pair_id = assign_random_pair_ids_for_indices!(
            pair_ids, precursor_idx, indices, last_pair_id, strategy.seed
        )
    end

    set_column!(psms, :pair_id, pair_ids)
    return nothing
end

function assign_pairs!(psms::AbstractPSMContainer, ::NoPairing)
    n = nrows(psms)
    set_column!(psms, :pair_id, UInt32.(1:n))
    set_column!(psms, :irt_bin_idx, zeros(UInt32, n))
    return nothing
end

# DataFrame version for backward compatibility
function assign_pairs!(psms::DataFrame, strategy::PairingStrategy)
    container = DataFramePSMContainer(psms, Val(:unsafe))
    assign_pairs!(container, strategy)
    return nothing
end

# Helper functions
function getIrtBins(irts::AbstractVector{R}, bin_size::Int = 1000) where {R<:Real}
    sort_idx = sortperm(irts)
    bin_idx, bin_count = zero(UInt32), zero(UInt32)
    bin_idxs = similar(irts, UInt32, length(irts))
    for idx in sort_idx
        bin_count += one(UInt32)
        bin_idxs[idx] = bin_idx
        if bin_count >= bin_size
            bin_idx += one(UInt32)
            bin_count = zero(UInt32)
        end
    end
    return bin_idxs
end

function assign_random_pair_ids_for_indices!(
    pair_ids::Vector{UInt32},
    precursor_idx::AbstractVector{UInt32},
    indices::Vector{Int},
    last_pair_id::UInt32,
    seed::Int
)
    # Get unique precursors in this group
    group_precursors = precursor_idx[indices]
    unique_precursors = unique(group_precursors)
    n_precursors = length(unique_precursors)

    # Shuffle precursors
    perm = randperm(MersenneTwister(seed), n_precursors)
    shuffled_precursors = unique_precursors[perm]

    # Create mapping
    precursor_to_pair = Dict{UInt32, UInt32}()
    n_pairs = n_precursors ÷ 2

    for i in 1:n_pairs
        last_pair_id += one(UInt32)
        precursor_to_pair[shuffled_precursors[2*i - 1]] = last_pair_id
        precursor_to_pair[shuffled_precursors[2*i]] = last_pair_id
    end

    # Handle odd number
    if isodd(n_precursors)
        last_pair_id += one(UInt32)
        precursor_to_pair[shuffled_precursors[end]] = last_pair_id
    end

    # Assign pair IDs for this group
    for idx in indices
        pair_ids[idx] = precursor_to_pair[precursor_idx[idx]]
    end

    return last_pair_id
end
