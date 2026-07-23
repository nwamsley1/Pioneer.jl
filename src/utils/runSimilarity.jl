# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Pioneer.jl is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

const RunPair = Tuple{UInt32, UInt32}

"""
    RunSimilarityAtlas

Experiment-level, IDF-weighted directional run similarity. Rare observed IDs
are stored as shared run-pair weights, while common IDs use a complement
representation that records the fewer missing run pairs.

`observed_ids_by_run` retains only positive-IDF observations for consumers
such as MBR that need both run-pair similarity and precursor-presence queries.
IDs observed in every run have zero IDF and are not retained.
"""
struct RunSimilarityAtlas
    shared_weight::Dict{RunPair, Float32}
    missing_weight::Dict{RunPair, Float32}
    total_weight_by_run::Dict{UInt32, Float32}
    common_weight_by_run::Dict{UInt32, Float32}
    run_ids::BitSet
    observed_ids_by_run::Dict{UInt32, BitSet}
    centrality_by_run::Dict{UInt32, Float32}
end

@inline function _ordered_run_pair(
    left_run::UInt32,
    right_run::UInt32,
)::RunPair
    return left_run < right_run ?
        (left_run, right_run) : (right_run, left_run)
end

"""
    run_similarity(atlas, receiver_run, donor_run)

Return the directional containment of the receiver run's weighted observed-ID
set in the donor run.
"""
@inline function run_similarity(
    atlas::RunSimilarityAtlas,
    receiver_run::UInt32,
    donor_run::UInt32,
)::Float32
    receiver_run in atlas.run_ids || return 0.0f0
    donor_run in atlas.run_ids || return 0.0f0
    denominator = get(atlas.total_weight_by_run, receiver_run, 0.0f0)
    denominator > 0.0f0 || return 0.0f0
    receiver_run == donor_run && return 1.0f0

    pair_key = _ordered_run_pair(receiver_run, donor_run)
    direct_key = (receiver_run, donor_run)
    rare_shared = get(atlas.shared_weight, pair_key, 0.0f0)
    common_shared = get(
        atlas.common_weight_by_run,
        receiver_run,
        0.0f0,
    ) - get(atlas.missing_weight, direct_key, 0.0f0)
    numerator = clamp(rare_shared + common_shared, 0.0f0, denominator)
    return numerator / denominator
end

"""
    run_centrality(atlas, run_idx)

Return a run's mean directional similarity to all other experiment runs.
"""
@inline run_centrality(atlas::RunSimilarityAtlas, run_idx::UInt32)::Float32 =
    get(atlas.centrality_by_run, run_idx, 0.0f0)

"""
    is_observed_in_run(atlas, observed_id, run_idx)

Test whether a positive-IDF ID used to construct the atlas was observed in a
run.
"""
@inline function is_observed_in_run(
    atlas::RunSimilarityAtlas,
    observed_id::Integer,
    run_idx::UInt32,
)::Bool
    observed = get(atlas.observed_ids_by_run, run_idx, nothing)
    observed === nothing && return false
    return Int(observed_id) in observed
end

function _compute_run_centrality(
    shared_weight::Dict{RunPair, Float32},
    missing_weight::Dict{RunPair, Float32},
    total_weight_by_run::Dict{UInt32, Float32},
    common_weight_by_run::Dict{UInt32, Float32},
    run_ids::BitSet,
)::Dict{UInt32, Float32}
    n_peers = max(length(run_ids) - 1, 0)
    n_peers == 0 && return Dict{UInt32, Float32}()

    shared_sum_by_run = Dict{UInt32, Float64}()
    for ((left_run, right_run), weight) in shared_weight
        shared_sum_by_run[left_run] =
            get(shared_sum_by_run, left_run, 0.0) + Float64(weight)
        shared_sum_by_run[right_run] =
            get(shared_sum_by_run, right_run, 0.0) + Float64(weight)
    end

    missing_sum_by_run = Dict{UInt32, Float64}()
    for ((receiver_run, _), weight) in missing_weight
        missing_sum_by_run[receiver_run] =
            get(missing_sum_by_run, receiver_run, 0.0) + Float64(weight)
    end

    centrality_by_run = Dict{UInt32, Float32}()
    for run_idx_int in run_ids
        run_idx = UInt32(run_idx_int)
        denominator = Float64(get(total_weight_by_run, run_idx, 0.0f0))
        if !(denominator > 0.0)
            centrality_by_run[run_idx] = 0.0f0
            continue
        end
        common_weight = Float64(get(
            common_weight_by_run,
            run_idx,
            0.0f0,
        ))
        numerator_sum = get(shared_sum_by_run, run_idx, 0.0) +
            Float64(n_peers) * common_weight -
            get(missing_sum_by_run, run_idx, 0.0)
        centrality_by_run[run_idx] = Float32(clamp(
            numerator_sum / (Float64(n_peers) * denominator),
            0.0,
            1.0,
        ))
    end
    return centrality_by_run
end

"""
    build_run_similarity(observed_ids_by_run)

Build a compact experiment-level run-similarity atlas from the unique IDs
observed in each run. Each ID is weighted by inverse run frequency. IDs
observed in every run have zero IDF and are excluded from retained atlas
evidence.
"""
function build_run_similarity(
    observed_ids_by_run::Dict{UInt32, Vector{UInt32}},
)::RunSimilarityAtlas
    sorted_run_ids = sort!(collect(keys(observed_ids_by_run)))
    run_ids = BitSet(Int(run_idx) for run_idx in sorted_run_ids)
    n_runs = length(sorted_run_ids)

    ids_to_runs = Dict{UInt32, Vector{UInt32}}()
    for run_idx in sorted_run_ids
        for observed_id in observed_ids_by_run[run_idx]
            push!(
                get!(ids_to_runs, observed_id, UInt32[]),
                run_idx,
            )
        end
    end

    informative_ids_by_run = Dict(
        run_idx => BitSet()
        for run_idx in sorted_run_ids
    )
    total_weight_by_run = Dict(run_idx => 0.0 for run_idx in sorted_run_ids)
    common_weight_by_run = Dict(run_idx => 0.0 for run_idx in sorted_run_ids)
    shared_weight = Dict{RunPair, Float64}()
    missing_weight = Dict{RunPair, Float64}()
    missing_runs = UInt32[]

    for (observed_id, observed_runs) in ids_to_runs
        document_frequency = length(observed_runs)
        idf = log(Float64(n_runs + 1) / Float64(document_frequency + 1))
        idf > 0.0 || continue

        for run_idx in observed_runs
            push!(informative_ids_by_run[run_idx], Int(observed_id))
            total_weight_by_run[run_idx] += idf
        end

        shared_updates = document_frequency * (document_frequency - 1) ÷ 2
        missing_updates = document_frequency * (n_runs - document_frequency)
        if shared_updates <= missing_updates
            @inbounds for left_idx in 1:(document_frequency - 1)
                left_run = observed_runs[left_idx]
                for right_idx in (left_idx + 1):document_frequency
                    pair_key = (left_run, observed_runs[right_idx])
                    shared_weight[pair_key] =
                        get(shared_weight, pair_key, 0.0) + idf
                end
            end
        else
            for run_idx in observed_runs
                common_weight_by_run[run_idx] += idf
            end

            empty!(missing_runs)
            posting_idx = 1
            @inbounds for run_idx in sorted_run_ids
                if posting_idx <= document_frequency &&
                   observed_runs[posting_idx] == run_idx
                    posting_idx += 1
                else
                    push!(missing_runs, run_idx)
                end
            end
            @inbounds for receiver_run in observed_runs
                for donor_run in missing_runs
                    direct_key = (receiver_run, donor_run)
                    missing_weight[direct_key] =
                        get(missing_weight, direct_key, 0.0) + idf
                end
            end
        end
    end

    shared_weight_f32 = Dict(
        key => Float32(value)
        for (key, value) in shared_weight
    )
    missing_weight_f32 = Dict(
        key => Float32(value)
        for (key, value) in missing_weight
    )
    total_weight_f32 = Dict(
        key => Float32(value)
        for (key, value) in total_weight_by_run
    )
    common_weight_f32 = Dict(
        key => Float32(value)
        for (key, value) in common_weight_by_run
    )
    centrality_by_run = _compute_run_centrality(
        shared_weight_f32,
        missing_weight_f32,
        total_weight_f32,
        common_weight_f32,
        run_ids,
    )

    return RunSimilarityAtlas(
        shared_weight_f32,
        missing_weight_f32,
        total_weight_f32,
        common_weight_f32,
        run_ids,
        informative_ids_by_run,
        centrality_by_run,
    )
end
