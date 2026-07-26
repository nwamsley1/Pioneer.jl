# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

@inline function _mbr_profile_centroid_dot(
    profile::_MBRRunIntensityProfile,
    centroid::Vector{Float32},
)
    value = 0.0f0
    @inbounds for idx in eachindex(profile.precursor_columns)
        value += profile.values[idx] * centroid[profile.precursor_columns[idx]]
    end
    return value
end

@inline function _mbr_set_dense_centroid!(
    dense_centroid::Vector{Float32},
    centroid::_MBRSparseCentroid,
)
    @inbounds for idx in eachindex(centroid.precursor_columns)
        dense_centroid[centroid.precursor_columns[idx]] = centroid.values[idx]
    end
    return nothing
end

@inline function _mbr_clear_dense_centroid!(
    dense_centroid::Vector{Float32},
    centroid::_MBRSparseCentroid,
)
    @inbounds for precursor_column in centroid.precursor_columns
        dense_centroid[precursor_column] = 0.0f0
    end
    return nothing
end

@inline function _mbr_profile_dot(
    left::_MBRRunIntensityProfile,
    right::_MBRRunIntensityProfile,
)
    value = 0.0f0
    left_idx = 1
    right_idx = 1
    @inbounds while left_idx <= length(left.precursor_columns) &&
                    right_idx <= length(right.precursor_columns)
        left_column = left.precursor_columns[left_idx]
        right_column = right.precursor_columns[right_idx]
        if left_column == right_column
            value += left.values[left_idx] * right.values[right_idx]
            left_idx += 1
            right_idx += 1
        elseif left_column < right_column
            left_idx += 1
        else
            right_idx += 1
        end
    end
    return value
end

function _mbr_run_intensity_profiles(
    postings_by_precursor::Dict{UInt32, Vector{_MBRIntensityPosting}},
    file_ids::Vector{UInt32},
)
    precursor_ids = sort!(collect(keys(postings_by_precursor)))
    length(precursor_ids) <= Int(typemax(Int32)) ||
        error("MBR run clustering exceeded Int32 precursor columns")

    profile_idx_by_file = Dict{UInt32, Int}(
        file_idx => profile_idx
        for (profile_idx, file_idx) in enumerate(file_ids)
    )
    precursor_columns = [Int32[] for _ in file_ids]
    values = [Float32[] for _ in file_ids]

    for (precursor_column, precursor_idx) in enumerate(precursor_ids)
        for posting in postings_by_precursor[precursor_idx]
            profile_idx = get(profile_idx_by_file, posting.file_idx, 0)
            profile_idx == 0 && continue
            value = posting.log2_weight
            isfinite(value) && value != 0.0f0 || continue
            push!(precursor_columns[profile_idx], Int32(precursor_column))
            push!(values[profile_idx], value)
        end
    end

    profiles = _MBRRunIntensityProfile[]
    sizehint!(profiles, length(file_ids))
    for (profile_idx, file_idx) in enumerate(file_ids)
        profile_values = values[profile_idx]
        sum_squares = sum(abs2, profile_values; init = 0.0f0)
        sum_squares > 0.0f0 || continue
        inv_norm = inv(sqrt(sum_squares))
        @inbounds for idx in eachindex(profile_values)
            profile_values[idx] *= inv_norm
        end
        push!(
            profiles,
            _MBRRunIntensityProfile(
                file_idx,
                precursor_columns[profile_idx],
                profile_values,
            ),
        )
    end
    return profiles, length(precursor_ids)
end

function _mbr_sparse_centroid(
    profiles::Vector{_MBRRunIntensityProfile},
    members::Vector{Int},
    centroid_accumulator::Vector{Float64},
    touched_columns::Vector{Int32},
    column_touched::BitVector,
)
    empty!(touched_columns)
    for profile_idx in members
        profile = profiles[profile_idx]
        @inbounds for idx in eachindex(profile.precursor_columns)
            precursor_column = profile.precursor_columns[idx]
            if !column_touched[precursor_column]
                column_touched[precursor_column] = true
                push!(touched_columns, precursor_column)
            end
            centroid_accumulator[precursor_column] +=
                Float64(profile.values[idx])
        end
    end

    sort!(touched_columns)
    norm_squared = 0.0
    @inbounds for precursor_column in touched_columns
        norm_squared += abs2(centroid_accumulator[precursor_column])
    end
    inv_norm = norm_squared > 0.0 ? inv(sqrt(norm_squared)) : 0.0
    precursor_columns = Int32[]
    values = Float32[]
    sizehint!(precursor_columns, length(touched_columns))
    sizehint!(values, length(touched_columns))
    @inbounds for precursor_column in touched_columns
        value = centroid_accumulator[precursor_column]
        if value != 0.0
            push!(precursor_columns, precursor_column)
            push!(values, Float32(value * inv_norm))
        end
        centroid_accumulator[precursor_column] = 0.0
        column_touched[precursor_column] = false
    end
    empty!(touched_columns)
    return _MBRSparseCentroid(precursor_columns, values)
end

function _mbr_spherical_run_cluster(
    profiles::Vector{_MBRRunIntensityProfile},
    members::Vector{Int},
    centroid_accumulator::Vector{Float64},
    touched_columns::Vector{Int32},
    column_touched::BitVector,
    dense_centroid::Vector{Float32},
)
    centroid = _mbr_sparse_centroid(
        profiles,
        members,
        centroid_accumulator,
        touched_columns,
        column_touched,
    )
    _mbr_set_dense_centroid!(dense_centroid, centroid)
    similarity_sum = 0.0
    for profile_idx in members
        similarity_sum += Float64(_mbr_profile_centroid_dot(
            profiles[profile_idx],
            dense_centroid,
        ))
    end
    _mbr_clear_dense_centroid!(dense_centroid, centroid)
    return _MBRSphericalRunCluster(members, centroid, similarity_sum)
end

@inline function _mbr_profile_as_centroid(profile::_MBRRunIntensityProfile)
    return _MBRSparseCentroid(
        copy(profile.precursor_columns),
        copy(profile.values),
    )
end

function _mbr_partition_spherical_run_members(
    profiles::Vector{_MBRRunIntensityProfile},
    parent_members::Vector{Int},
    dense_centroid1::Vector{Float32},
    dense_centroid2::Vector{Float32},
    seed1::Int,
    seed2::Int,
    min_child_size::Int,
)
    members1 = Int[]
    members2 = Int[]
    sizehint!(members1, length(parent_members) >>> 1)
    sizehint!(members2, length(parent_members) >>> 1)
    for profile_idx in parent_members
        if profile_idx == seed1
            push!(members1, profile_idx)
        elseif profile_idx == seed2
            push!(members2, profile_idx)
        elseif _mbr_profile_centroid_dot(
            profiles[profile_idx],
            dense_centroid1,
        ) >= _mbr_profile_centroid_dot(
            profiles[profile_idx],
            dense_centroid2,
        )
            push!(members1, profile_idx)
        else
            push!(members2, profile_idx)
        end
    end
    length(members1) >= min_child_size &&
        length(members2) >= min_child_size || return nothing
    return members1, members2
end

function _mbr_split_spherical_run_cluster(
    profiles::Vector{_MBRRunIntensityProfile},
    parent::_MBRSphericalRunCluster,
    centroid_accumulator::Vector{Float64},
    touched_columns::Vector{Int32},
    column_touched::BitVector,
    dense_centroid1::Vector{Float32},
    dense_centroid2::Vector{Float32};
    split_iterations::Int,
    min_child_size::Int,
)
    length(parent.members) >= 2 * min_child_size || return nothing

    _mbr_set_dense_centroid!(dense_centroid1, parent.centroid)
    seed1 = parent.members[1]
    lowest_parent_similarity = Inf32
    for profile_idx in parent.members
        similarity = _mbr_profile_centroid_dot(
            profiles[profile_idx],
            dense_centroid1,
        )
        if similarity < lowest_parent_similarity
            seed1 = profile_idx
            lowest_parent_similarity = similarity
        end
    end
    _mbr_clear_dense_centroid!(dense_centroid1, parent.centroid)

    seed2 = 0
    lowest_seed_similarity = Inf32
    for profile_idx in parent.members
        profile_idx == seed1 && continue
        similarity = _mbr_profile_dot(
            profiles[seed1],
            profiles[profile_idx],
        )
        if similarity < lowest_seed_similarity
            seed2 = profile_idx
            lowest_seed_similarity = similarity
        end
    end
    seed2 == 0 && return nothing

    centroid1 = _mbr_profile_as_centroid(profiles[seed1])
    centroid2 = _mbr_profile_as_centroid(profiles[seed2])
    child1 = parent
    child2 = parent
    for _ in 1:split_iterations
        _mbr_set_dense_centroid!(dense_centroid1, centroid1)
        _mbr_set_dense_centroid!(dense_centroid2, centroid2)
        partition = _mbr_partition_spherical_run_members(
            profiles,
            parent.members,
            dense_centroid1,
            dense_centroid2,
            seed1,
            seed2,
            min_child_size,
        )
        _mbr_clear_dense_centroid!(dense_centroid1, centroid1)
        _mbr_clear_dense_centroid!(dense_centroid2, centroid2)
        partition === nothing && return nothing
        members1, members2 = partition
        child1 = _mbr_spherical_run_cluster(
            profiles,
            members1,
            centroid_accumulator,
            touched_columns,
            column_touched,
            dense_centroid1,
        )
        child2 = _mbr_spherical_run_cluster(
            profiles,
            members2,
            centroid_accumulator,
            touched_columns,
            column_touched,
            dense_centroid1,
        )
        centroid1 = child1.centroid
        centroid2 = child2.centroid
    end

    gain = Float32(
        (
            child1.similarity_sum +
            child2.similarity_sum -
            parent.similarity_sum
        ) / Float64(length(parent.members)),
    )
    return (child1 = child1, child2 = child2, gain = gain)
end

function _fit_mbr_run_clusters(
    postings_by_precursor::Dict{UInt32, Vector{_MBRIntensityPosting}},
    file_ids::Vector{UInt32};
    min_child_size::Int = MBR_RUN_CLUSTER_MIN_CHILD_SIZE,
    split_iterations::Int = MBR_RUN_CLUSTER_SPLIT_ITERATIONS,
    min_gain::Float32 = MBR_RUN_CLUSTER_MIN_GAIN,
)
    profiles, n_precursor_columns = _mbr_run_intensity_profiles(
        postings_by_precursor,
        file_ids,
    )
    cluster_by_file = Dict{UInt32, UInt32}()
    cluster_sizes = UInt32[]

    if !isempty(profiles)
        centroid_accumulator = zeros(Float64, n_precursor_columns)
        touched_columns = Int32[]
        column_touched = falses(n_precursor_columns)
        dense_centroid1 = zeros(Float32, n_precursor_columns)
        dense_centroid2 = zeros(Float32, n_precursor_columns)
        root = _mbr_spherical_run_cluster(
            profiles,
            collect(eachindex(profiles)),
            centroid_accumulator,
            touched_columns,
            column_touched,
            dense_centroid1,
        )
        pending = _MBRSphericalRunCluster[root]

        while !isempty(pending)
            cluster = pop!(pending)
            split = _mbr_split_spherical_run_cluster(
                profiles,
                cluster,
                centroid_accumulator,
                touched_columns,
                column_touched,
                dense_centroid1,
                dense_centroid2;
                split_iterations = split_iterations,
                min_child_size = min_child_size,
            )
            if split !== nothing && split.gain >= min_gain
                push!(pending, split.child2)
                push!(pending, split.child1)
                continue
            end

            cluster_idx = UInt32(length(cluster_sizes) + 1)
            push!(cluster_sizes, UInt32(length(cluster.members)))
            for profile_idx in cluster.members
                cluster_by_file[profiles[profile_idx].file_idx] = cluster_idx
            end
        end
    end

    for file_idx in file_ids
        haskey(cluster_by_file, file_idx) && continue
        cluster_idx = UInt32(length(cluster_sizes) + 1)
        cluster_by_file[file_idx] = cluster_idx
        push!(cluster_sizes, UInt32(1))
    end
    return (
        cluster_by_file = cluster_by_file,
        cluster_sizes = cluster_sizes,
    )
end

"""
    build_mbr_receiver_run_clusters(file_paths; q_value_threshold)

Cluster runs by their pre-MBR identified-precursor intensity profiles and
retain precursor support counts within each receiver's run cluster.
"""
function build_mbr_receiver_run_clusters(
    file_paths::Vector{String};
    q_value_threshold::Float32,
)
    passed_by_file = Dict{UInt32, BitSet}()
    files_by_precursor = Dict{UInt32, Vector{UInt32}}()
    intensity_postings =
        Dict{UInt32, Vector{_MBRIntensityPosting}}()

    for path in file_paths
        tbl = Arrow.Table(path)
        for column in (:precursor_idx, :ms_file_idx, :qval, :weight)
            hasproperty(tbl, column) ||
                error("MBR run-cluster features require column $column in $path")
        end
        @inbounds for row in eachindex(tbl.precursor_idx)
            file_idx = UInt32(tbl.ms_file_idx[row])
            passed = get!(() -> BitSet(), passed_by_file, file_idx)
            qval = Float32(tbl.qval[row])
            isfinite(qval) && qval <= q_value_threshold || continue
            precursor_idx = UInt32(tbl.precursor_idx[row])
            if !(Int(precursor_idx) in passed)
                push!(passed, Int(precursor_idx))
                push!(
                    get!(() -> UInt32[], files_by_precursor, precursor_idx),
                    file_idx,
                )
            end
            weight = Float32(tbl.weight[row])
            isfinite(weight) && weight > 0.0f0 || continue
            push!(
                get!(
                    () -> _MBRIntensityPosting[],
                    intensity_postings,
                    precursor_idx,
                ),
                _MBRIntensityPosting(file_idx, log2(weight)),
            )
        end
    end

    file_ids = sort!(collect(keys(passed_by_file)))
    isempty(file_ids) && return _MBRReceiverRunClusters()
    fit = _fit_mbr_run_clusters(intensity_postings, file_ids)
    support_by_precursor =
        Dict{UInt32, Dict{UInt32, UInt32}}()
    for (precursor_idx, run_files) in files_by_precursor
        support_by_cluster = Dict{UInt32, UInt32}()
        for file_idx in run_files
            cluster_idx = get(fit.cluster_by_file, file_idx, UInt32(0))
            cluster_idx == UInt32(0) && continue
            support_by_cluster[cluster_idx] =
                get(support_by_cluster, cluster_idx, UInt32(0)) + UInt32(1)
        end
        isempty(support_by_cluster) ||
            (support_by_precursor[precursor_idx] = support_by_cluster)
    end
    return _MBRReceiverRunClusters(
        fit.cluster_by_file,
        fit.cluster_sizes,
        support_by_precursor,
        passed_by_file,
    )
end

@inline function _mbr_receiver_cluster_features(
    clusters::_MBRReceiverRunClusters,
    precursor_idx::UInt32,
    receiver_file::UInt32,
)
    cluster_idx = get(clusters.cluster_by_file, receiver_file, UInt32(0))
    cluster_int = Int(cluster_idx)
    1 <= cluster_int <= length(clusters.cluster_sizes) || return (
        support_count = 0.0f0,
        peer_count = 0.0f0,
        support_fraction = 0.0f0,
    )
    cluster_size = clusters.cluster_sizes[cluster_int]
    cluster_size > UInt32(1) || return (
        support_count = 0.0f0,
        peer_count = 0.0f0,
        support_fraction = 0.0f0,
    )

    peer_count = cluster_size - UInt32(1)
    support_by_cluster = get(
        clusters.support_by_precursor,
        precursor_idx,
        nothing,
    )
    support_count = support_by_cluster === nothing ?
        UInt32(0) :
        get(support_by_cluster, cluster_idx, UInt32(0))
    receiver_passed = get(clusters.passed_by_file, receiver_file, nothing)
    if support_count > UInt32(0) &&
       receiver_passed !== nothing &&
       Int(precursor_idx) in receiver_passed
        support_count -= UInt32(1)
    end
    support_count = min(support_count, peer_count)
    return (
        support_count = Float32(support_count),
        peer_count = Float32(peer_count),
        support_fraction = Float32(support_count) / Float32(peer_count),
    )
end
