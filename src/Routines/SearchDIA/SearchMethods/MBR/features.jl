# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

@inline function _mbr_sqrt_tuple_is_valid(
    spectrum::NTuple{8, Float32},
)
    @inbounds for rank in 1:8
        spectrum[rank] > 0.0f0 && return true
    end
    return false
end

@inline function _mbr_sqrt_tuple(columns, row::Integer)
    return ntuple(rank -> begin
        value = Float32(columns[rank][row])
        isfinite(value) ? max(value, 0.0f0) : 0.0f0
    end, 8)
end

@inline function _mbr_hellinger_from_sqrt(
    receiver::NTuple{8, Float32},
    donor::NTuple{8, Float32},
)
    (_mbr_sqrt_tuple_is_valid(receiver) &&
     _mbr_sqrt_tuple_is_valid(donor)) || return 1.0f0
    bhattacharyya = 0.0f0
    @inbounds for rank in 1:8
        bhattacharyya += receiver[rank] * donor[rank]
    end
    return sqrt(clamp(1.0f0 - bhattacharyya, 0.0f0, 1.0f0))
end

@inline function _mbr_temporal_masked_hellinger(
    receiver_trace::AbstractVector,
    donor::NTuple{8, Float32},
    mask::UInt8,
)
    count_ones(mask) >= 2 || return 1.0f0
    length(receiver_trace) % MBR_TEMPORAL_TRACE_STRIDE == 0 || return 1.0f0

    donor_mass = 0.0f0
    @inbounds for rank in 1:8
        ((mask >> (rank - 1)) & 0x01) == 0x01 || continue
        donor_mass += donor[rank] * donor[rank]
    end
    donor_mass > 0.0f0 || return 1.0f0
    inv_donor_norm = inv(sqrt(donor_mass))

    weighted_bhattacharyya = 0.0f0
    total_weight = 0.0f0
    @inbounds for offset in 1:MBR_TEMPORAL_TRACE_STRIDE:length(receiver_trace)
        weight_value = Float32(receiver_trace[offset])
        weight = isfinite(weight_value) ? max(weight_value, 0.0f0) : 0.0f0
        weight > 0.0f0 || continue
        total_weight += weight

        receiver_mass = 0.0f0
        for rank in 1:8
            ((mask >> (rank - 1)) & 0x01) == 0x01 || continue
            value = Float32(receiver_trace[offset + rank])
            receiver_mass += isfinite(value) ? max(value, 0.0f0) : 0.0f0
        end
        receiver_mass > 0.0f0 || continue
        inv_receiver_norm = inv(sqrt(receiver_mass))

        scan_bhattacharyya = 0.0f0
        for rank in 1:8
            ((mask >> (rank - 1)) & 0x01) == 0x01 || continue
            value = Float32(receiver_trace[offset + rank])
            intensity = isfinite(value) ? max(value, 0.0f0) : 0.0f0
            scan_bhattacharyya +=
                sqrt(intensity) * inv_receiver_norm *
                donor[rank] * inv_donor_norm
        end
        weighted_bhattacharyya += weight * scan_bhattacharyya
    end
    total_weight > 0.0f0 || return 1.0f0
    return sqrt(clamp(
        1.0f0 - weighted_bhattacharyya / total_weight,
        0.0f0,
        1.0f0,
    ))
end

@inline function _mbr_shared_corr_mask(
    receiver_mask::UInt8,
    donor_mask::UInt8,
)
    return receiver_mask & donor_mask
end

@inline function _mbr_run_similarity(
    atlas::Union{Nothing, RunSimilarityAtlas},
    receiver_file::UInt32,
    donor_file::UInt32,
)
    atlas === nothing && return 0.0f0
    return run_similarity(atlas, receiver_file, donor_file)
end

function _mbr_donor_score_floor(
    file_paths::Vector{String};
    donor_q_threshold::Float32 = MBR_DONOR_Q_THRESHOLD,
    require_initial_pass::Bool = false,
    q_value_threshold::Float32 = donor_q_threshold,
)
    scores = Float32[]
    targets = Bool[]
    for path in file_paths
        tbl = Arrow.Table(path)
        hasproperty(tbl, :trace_prob_prepass) ||
            error("MBR donor selection requires :trace_prob_prepass in $path")
        @inbounds for row in eachindex(tbl.trace_prob_prepass)
            if require_initial_pass
                hasproperty(tbl, :qval) && hasproperty(tbl, :global_qval) ||
                    error("Initial-pass donor floor requires q-value columns in $path")
                _mbr_initial_pass(
                    tbl.qval[row],
                    tbl.global_qval[row],
                    q_value_threshold,
                ) || continue
            end
            push!(scores, Float32(tbl.trace_prob_prepass[row]))
            push!(targets, Bool(tbl.target[row]))
        end
    end
    isempty(scores) && return Inf32
    qvalues = similar(scores)
    get_qvalues!(scores, targets, qvalues)
    eligible = targets .& (qvalues .<= donor_q_threshold)
    return any(eligible) ? minimum(scores[eligible]) : Inf32
end

function _mbr_preintegration_donor_files(
    file_paths::Vector{String},
    score_floor::Float32,
)
    donor_files = Dict{UInt32, Vector{Tuple{UInt32, Float32}}}()
    for path in file_paths
        tbl = Arrow.Table(path)
        @inbounds for row in eachindex(tbl.precursor_idx)
            score = Float32(tbl.trace_prob_prepass[row])
            score >= score_floor || continue
            pid = UInt32(tbl.precursor_idx[row])
            file_idx = UInt32(tbl.ms_file_idx[row])
            entries = get!(
                () -> Tuple{UInt32, Float32}[],
                donor_files,
                pid,
            )
            existing = findfirst(entry -> entry[1] == file_idx, entries)
            if existing === nothing
                push!(entries, (file_idx, score))
            elseif score > entries[existing][2]
                entries[existing] = (file_idx, score)
            end
        end
    end
    return donor_files
end

@inline function _mbr_has_cross_run_donor(
    donor_files::Dict{UInt32, Vector{Tuple{UInt32, Float32}}},
    pid::UInt32,
    receiver_file::UInt32,
)
    entries = get(donor_files, pid, nothing)
    entries === nothing && return false
    return any(entry -> entry[1] != receiver_file, entries)
end

function build_mbr_integrated_donor_dict(
    file_paths::Vector{String},
    score_floor::Float32;
    q_value_threshold::Float32,
)
    donor_dict = Dict{UInt32, Vector{_MBRDonorEntry}}()
    for path in file_paths
        tbl = Arrow.Table(path)
        required = (
            :precursor_idx,
            :ms_file_idx,
            :trace_prob_prepass,
            :qval,
            :global_qval,
            :irt_pred,
            MBR_INTEGRATED_WEIGHT_COLUMN,
            MBR_INTEGRATED_LOG2_INTENSITY_EXPLAINED_COLUMN,
            MBR_INTEGRATED_APEX_IRT_COLUMN,
            MBR_INTEGRATED_FRAG_CORR_BITVEC_COLUMN,
            MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_BITVEC_RANK_COLUMN,
            MBR_INTEGRATED_N_SCANS_COLUMN,
        )
        for col in required
            hasproperty(tbl, col) ||
                error("Integrated MBR donor selection requires column $col in $path")
        end
        frag_cols = ntuple(
            rank -> getproperty(tbl, MBR_INTEGRATED_FRAGMENT_SQRT_COLUMNS[rank]),
            8,
        )
        @inbounds for row in eachindex(tbl.precursor_idx)
            qval = Float32(tbl.qval[row])
            global_qval = Float32(tbl.global_qval[row])
            score = Float32(tbl.trace_prob_prepass[row])
            initial_pass =
                isfinite(qval) && qval <= q_value_threshold &&
                isfinite(global_qval) && global_qval <= q_value_threshold
            initial_pass && score >= score_floor || continue

            pid = UInt32(tbl.precursor_idx[row])
            file_idx = UInt32(tbl.ms_file_idx[row])
            integrated_frag_sqrt = _mbr_sqrt_tuple(frag_cols, row)
            weight =
                Float32(getproperty(tbl, MBR_INTEGRATED_WEIGHT_COLUMN)[row])
            irt_obs =
                Float32(getproperty(tbl, MBR_INTEGRATED_APEX_IRT_COLUMN)[row])
            n_scans =
                Float32(getproperty(tbl, MBR_INTEGRATED_N_SCANS_COLUMN)[row])
            irt_pred = Float32(tbl.irt_pred[row])
            donor = _MBRDonorEntry(
                score,
                pid,
                weight,
                Float32(getproperty(
                    tbl,
                    MBR_INTEGRATED_LOG2_INTENSITY_EXPLAINED_COLUMN,
                )[row]),
                irt_pred - irt_obs,
                irt_obs,
                n_scans,
                integrated_frag_sqrt,
                UInt8(getproperty(
                    tbl,
                    MBR_INTEGRATED_FRAG_CORR_BITVEC_COLUMN,
                )[row]),
                UInt16(getproperty(
                    tbl,
                    MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_BITVEC_RANK_COLUMN,
                )[row]),
                file_idx,
            )
            entries = get!(() -> _MBRDonorEntry[], donor_dict, pid)
            existing = findfirst(entry -> entry.ms_file_idx == file_idx, entries)
            if existing === nothing
                push!(entries, donor)
            elseif donor.trace_prob > entries[existing].trace_prob
                entries[existing] = donor
            end
        end
    end
    return donor_dict
end

function _mbr_lod_thresholds(
    file_paths::Vector{String},
    score_floor::Float32;
    q_value_threshold::Float32,
)
    samples_by_file = Dict{UInt32, Vector{Float32}}()
    global_samples = Float32[]
    for path in file_paths
        tbl = Arrow.Table(path)
        weight_column = getproperty(tbl, MBR_INTEGRATED_WEIGHT_COLUMN)
        @inbounds for row in eachindex(tbl.precursor_idx)
            Bool(tbl.target[row]) || continue
            _mbr_initial_pass(
                tbl.qval[row],
                tbl.global_qval[row],
                q_value_threshold,
            ) || continue
            Float32(tbl.trace_prob_prepass[row]) >= score_floor || continue
            weight = Float32(weight_column[row])
            isfinite(weight) && weight > 0.0f0 || continue
            log_weight = log2(weight)
            file_idx = UInt32(tbl.ms_file_idx[row])
            push!(
                get!(() -> Float32[], samples_by_file, file_idx),
                log_weight,
            )
            push!(global_samples, log_weight)
        end
    end

    quantile_index(samples) = clamp(
        ceil(Int, Float64(MBR_LOD_WEIGHT_QUANTILE) * length(samples)),
        1,
        length(samples),
    )
    by_file = Dict{UInt32, Float32}()
    for (file_idx, samples) in samples_by_file
        sort!(samples)
        by_file[file_idx] = samples[quantile_index(samples)]
    end
    global_lod = if isempty(global_samples)
        NaN32
    else
        sort!(global_samples)
        global_samples[quantile_index(global_samples)]
    end
    return (by_file = by_file, global_lod = global_lod)
end

@inline function _mbr_donor_in_file(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pid::UInt32,
    file_idx::UInt32,
)
    entries = get(donor_dict, pid, nothing)
    entries === nothing && return nothing
    @inbounds for donor in entries
        donor.ms_file_idx == file_idx && return donor
    end
    return nothing
end

@inline function _mbr_select_donor(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pid::UInt32,
    receiver_file::UInt32,
    atlas::Union{Nothing, RunSimilarityAtlas},
)
    entries = get(donor_dict, pid, nothing)
    entries === nothing && return nothing
    best = nothing
    best_similarity = -Inf32
    best_score = -Inf32
    @inbounds for donor in entries
        donor.ms_file_idx == receiver_file && continue
        similarity = _mbr_run_similarity(
            atlas,
            receiver_file,
            donor.ms_file_idx,
        )
        if best === nothing ||
           similarity > best_similarity ||
           (similarity == best_similarity && donor.trace_prob > best_score)
            best = donor
            best_similarity = similarity
            best_score = donor.trace_prob
        end
    end
    return best
end

@inline function _mbr_top_scoring_donor(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pid::UInt32,
    receiver_file::UInt32,
)
    entries = get(donor_dict, pid, nothing)
    entries === nothing && return nothing
    best = nothing
    best_score = -Inf32
    @inbounds for donor in entries
        donor.ms_file_idx == receiver_file && continue
        score = donor.trace_prob
        if best === nothing || score > best_score
            best = donor
            best_score = score
        end
    end
    return best
end

@inline function _mbr_worst_alternate_donor(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pid::UInt32,
    receiver_file::UInt32,
    best_donor::_MBRDonorEntry,
)
    entries = get(donor_dict, pid, nothing)
    entries === nothing && return nothing
    worst = nothing
    worst_weight = Inf32
    worst_score = Inf32
    @inbounds for donor in entries
        donor.ms_file_idx == receiver_file && continue
        donor.ms_file_idx == best_donor.ms_file_idx && continue
        weight = isfinite(donor.weight) ? donor.weight : Inf32
        if worst === nothing ||
           weight < worst_weight ||
           (weight == worst_weight && donor.trace_prob < worst_score)
            worst = donor
            worst_weight = weight
            worst_score = donor.trace_prob
        end
    end
    return worst
end

@inline function _mbr_collect_false_donors_from_pool!(
    donors::Vector{Union{Nothing, _MBRDonorEntry}},
    seen_pids::BitSet,
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pool::_MBRIrtPool,
    target_irt::Float32,
    receiver_pid::UInt32,
    receiver_file::UInt32,
    eligibility::_MBRCounterfactualEligibility,
    atlas::Union{Nothing, RunSimilarityAtlas};
    donor_file::Union{Nothing, UInt32} = nothing,
)
    n_found = count(donor -> donor !== nothing, donors)
    @inbounds for pool_idx in _mbr_nearest_pool_order(pool, target_irt)
        pid = pool.pids[pool_idx]
        pid == receiver_pid && continue
        Int(pid) in seen_pids && continue
        _mbr_counterfactual_eligible(eligibility, receiver_file, pid) ||
            continue
        donor = donor_file === nothing ?
            _mbr_select_donor(donor_dict, pid, receiver_file, atlas) :
            _mbr_donor_in_file(donor_dict, pid, donor_file)
        donor === nothing && continue
        n_found += 1
        donors[n_found] = donor
        push!(seen_pids, Int(pid))
        n_found == MBR_N_COUNTERFACTUALS && break
    end
    return n_found
end

function _mbr_false_donors(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pools::_MBRPartnerPools,
    eligibility::_MBRCounterfactualEligibility,
    receiver_pid::UInt32,
    receiver_file::UInt32,
    target_irt::Float32,
    true_donor::_MBRDonorEntry,
    atlas::Union{Nothing, RunSimilarityAtlas},
)
    pid_int = Int(receiver_pid)
    fold = Int(pools.fold_by_pid[pid_int])
    mz_bin = Int(pools.mz_bin_by_pid[pid_int])
    charge = Int(pools.charge_by_pid[pid_int])
    sequence_length = Int(pools.length_by_pid[pid_int])
    donors = Union{Nothing, _MBRDonorEntry}[
        nothing for _ in 1:MBR_N_COUNTERFACTUALS
    ]
    seen_pids = BitSet()
    n_found = 0

    same_file_pool = get(
        pools.file_charge_length_pool,
        (true_donor.ms_file_idx, charge, sequence_length),
        _empty_mbr_irt_pool(),
    )
    n_found = _mbr_collect_false_donors_from_pool!(
        donors,
        seen_pids,
        donor_dict,
        same_file_pool,
        target_irt,
        receiver_pid,
        receiver_file,
        eligibility,
        atlas;
        donor_file = true_donor.ms_file_idx,
    )
    n_found == MBR_N_COUNTERFACTUALS && return donors

    fallback_pools = (
        get(
            pools.fold_mz_charge_length_pool,
            (fold, mz_bin, charge, sequence_length),
            _empty_mbr_irt_pool(),
        ),
        get(
            pools.fold_charge_length_pool,
            (fold, charge, sequence_length),
            _empty_mbr_irt_pool(),
        ),
        get(
            pools.charge_length_pool,
            (charge, sequence_length),
            _empty_mbr_irt_pool(),
        ),
    )
    for pool in fallback_pools
        n_found = _mbr_collect_false_donors_from_pool!(
            donors,
            seen_pids,
            donor_dict,
            pool,
            target_irt,
            receiver_pid,
            receiver_file,
            eligibility,
            atlas,
        )
        n_found == MBR_N_COUNTERFACTUALS && break
    end
    return donors
end

@inline function _mbr_log2_weight_ratio(
    receiver_weight::Float32,
    donor::Union{Nothing, _MBRDonorEntry},
)
    donor !== nothing &&
        isfinite(receiver_weight) && receiver_weight > 0.0f0 &&
        isfinite(donor.weight) && donor.weight > 0.0f0 || return -1.0f0
    return log2(receiver_weight / donor.weight)
end

@inline function _mbr_log2_explained_ratio(
    receiver_explained::Float32,
    donor::Union{Nothing, _MBRDonorEntry},
)
    donor !== nothing &&
        isfinite(receiver_explained) &&
        isfinite(donor.log2_intensity_explained) || return -1.0f0
    return receiver_explained - donor.log2_intensity_explained
end

@inline function _mbr_abs_scan_diff(
    receiver_n_scans::Float32,
    donor::Union{Nothing, _MBRDonorEntry},
)
    donor !== nothing &&
        isfinite(receiver_n_scans) &&
        isfinite(donor.n_scans) || return -1.0f0
    return abs(receiver_n_scans - donor.n_scans)
end

@inline function _mbr_log2_scan_ratio(
    receiver_n_scans::Float32,
    donor::Union{Nothing, _MBRDonorEntry},
)
    donor !== nothing &&
        isfinite(receiver_n_scans) &&
        isfinite(donor.n_scans) || return -1.0f0
    return log2((receiver_n_scans + 1.0f0) / (donor.n_scans + 1.0f0))
end

@inline function _mbr_residual_irt_diff(
    receiver_irt_pred::Float32,
    receiver_irt::Float32,
    donor::Union{Nothing, _MBRDonorEntry},
)
    donor !== nothing &&
        isfinite(receiver_irt_pred) &&
        isfinite(receiver_irt) &&
        isfinite(donor.irt_residual) || return -1.0f0
    return abs((receiver_irt_pred - receiver_irt) - donor.irt_residual)
end

@inline function _mbr_observed_irt_diff(
    receiver_irt::Float32,
    donor::Union{Nothing, _MBRDonorEntry},
)
    donor !== nothing &&
        isfinite(receiver_irt) &&
        isfinite(donor.irt_obs) || return -1.0f0
    return abs(receiver_irt - donor.irt_obs)
end

function _mbr_feature_values(
    receiver_pid::UInt32,
    receiver_weight::Float32,
    receiver_explained::Float32,
    receiver_irt_pred::Float32,
    receiver_irt::Float32,
    receiver_n_scans::Float32,
    receiver_temporal_mean::NTuple{8, Float32},
    receiver_temporal_trace::AbstractVector,
    receiver_corr_mask::UInt8,
    donor::_MBRDonorEntry,
    receiver_file::UInt32,
    atlas::Union{Nothing, RunSimilarityAtlas},
    clusters::_MBRReceiverRunClusters,
    bitvec_rank_table,
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
)
    worst_donor = _mbr_worst_alternate_donor(
        donor_dict,
        donor.precursor_idx,
        receiver_file,
        donor,
    )
    hellinger_donor = _mbr_top_scoring_donor(
        donor_dict,
        donor.precursor_idx,
        receiver_file,
    )
    hellinger_donor === nothing && (hellinger_donor = donor)
    donor_mask = hellinger_donor.frag_corr_bitvec
    shared_mask = _mbr_shared_corr_mask(receiver_corr_mask, donor_mask)
    cluster = _mbr_receiver_cluster_features(
        clusters,
        receiver_pid,
        receiver_file,
    )
    return (
        donor.trace_prob,
        worst_donor === nothing ? -1.0f0 : worst_donor.trace_prob,
        _mbr_run_similarity(atlas, receiver_file, donor.ms_file_idx),
        cluster.support_count,
        cluster.peer_count,
        cluster.support_fraction,
        _mbr_log2_weight_ratio(receiver_weight, donor),
        _mbr_log2_weight_ratio(receiver_weight, worst_donor),
        _mbr_log2_explained_ratio(receiver_explained, donor),
        _mbr_log2_explained_ratio(receiver_explained, worst_donor),
        _mbr_abs_scan_diff(receiver_n_scans, donor),
        _mbr_abs_scan_diff(receiver_n_scans, worst_donor),
        _mbr_log2_scan_ratio(receiver_n_scans, donor),
        _mbr_log2_scan_ratio(receiver_n_scans, worst_donor),
        _mbr_residual_irt_diff(receiver_irt_pred, receiver_irt, donor),
        _mbr_residual_irt_diff(
            receiver_irt_pred,
            receiver_irt,
            worst_donor,
        ),
        _mbr_observed_irt_diff(receiver_irt, donor),
        _mbr_observed_irt_diff(receiver_irt, worst_donor),
        worst_donor === nothing ? 1.0f0 : 0.0f0,
        hellinger_donor.trace_prob,
        _mbr_hellinger_from_sqrt(
            receiver_temporal_mean,
            hellinger_donor.integrated_frag_sqrt,
        ),
        _mbr_temporal_masked_hellinger(
            receiver_temporal_trace,
            hellinger_donor.integrated_frag_sqrt,
            donor_mask,
        ),
        _mbr_temporal_masked_hellinger(
            receiver_temporal_trace,
            hellinger_donor.integrated_frag_sqrt,
            receiver_corr_mask,
        ),
        _mbr_temporal_masked_hellinger(
            receiver_temporal_trace,
            hellinger_donor.integrated_frag_sqrt,
            shared_mask,
        ),
        Float32(hellinger_donor.frag_corr_bitvec_rank),
        Float32(_bitvec_pattern_rank(bitvec_rank_table, shared_mask)),
    )
end

function _mbr_sidecar_columns(n::Int)
    columns = Dict{Symbol, Any}(
        :precursor_idx => zeros(UInt32, n),
        :scan_idx => zeros(UInt32, n),
        :MBR_best_is_missing_true => trues(n),
        :MBR_log2_weight_lod_ratio => fill(-1.0f0, n),
        :MBR_receiver_frag_corr_bitvec_rank => fill(-1.0f0, n),
    )
    for stem in MBR_PAIRED_FEATURE_STEMS
        columns[_mbr_true_feature(stem)] = fill(-1.0f0, n)
    end
    for counterfactual_idx in 1:MBR_N_COUNTERFACTUALS
        columns[_mbr_missing_feature(counterfactual_idx)] = trues(n)
        for stem in MBR_PAIRED_FEATURE_STEMS
            columns[_mbr_false_feature(stem, counterfactual_idx)] =
                fill(-1.0f0, n)
        end
    end
    return columns
end

@inline function _mbr_log2_weight_lod_ratio(
    weight::Float32,
    file_idx::UInt32,
    lod_by_file::Dict{UInt32, Float32},
    global_lod::Float32,
)
    isfinite(weight) && weight > 0.0f0 || return -1.0f0
    lod = get(lod_by_file, file_idx, global_lod)
    isfinite(lod) || return -1.0f0
    return log2(weight) - lod
end

function compute_postintegration_mbr_features!(
    main_path::String,
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pools::_MBRPartnerPools,
    eligibility::_MBRCounterfactualEligibility;
    run_similarity_atlas::Union{Nothing, RunSimilarityAtlas},
    receiver_run_clusters::_MBRReceiverRunClusters,
    lod_log2_weight_by_file::Dict{UInt32, Float32},
    lod_log2_weight_global::Float32,
    bitvec_rank_table = nothing,
)
    main = Arrow.Table(main_path)
    n = length(main.precursor_idx)
    columns = _mbr_sidecar_columns(n)
    columns[:precursor_idx] .= UInt32.(main.precursor_idx)
    columns[:scan_idx] .= UInt32.(main.scan_idx)

    temporal_mean_columns = ntuple(
        rank -> getproperty(
            main,
            MBR_INTEGRATED_TEMPORAL_MEAN_SQRT_COLUMNS[rank],
        ),
        8,
    )
    temporal_trace_column = getproperty(
        main,
        MBR_INTEGRATED_TEMPORAL_TRACE_COLUMN,
    )
    weight_column = getproperty(main, MBR_INTEGRATED_WEIGHT_COLUMN)
    explained_column = getproperty(
        main,
        MBR_INTEGRATED_LOG2_INTENSITY_EXPLAINED_COLUMN,
    )
    irt_column = getproperty(main, MBR_INTEGRATED_APEX_IRT_COLUMN)
    scan_count_column = getproperty(main, MBR_INTEGRATED_N_SCANS_COLUMN)
    corr_mask_column = getproperty(
        main,
        MBR_INTEGRATED_FRAG_CORR_BITVEC_COLUMN,
    )
    corr_rank_column = getproperty(
        main,
        MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_BITVEC_RANK_COLUMN,
    )
    has_irt_pred = hasproperty(main, :irt_pred)

    @inbounds for row in 1:n
        receiver_pid = UInt32(main.precursor_idx[row])
        receiver_file = UInt32(main.ms_file_idx[row])
        true_donor = _mbr_select_donor(
            donor_dict,
            receiver_pid,
            receiver_file,
            run_similarity_atlas,
        )
        true_donor === nothing && continue

        receiver_weight = Float32(weight_column[row])
        receiver_explained = Float32(explained_column[row])
        receiver_irt = Float32(irt_column[row])
        receiver_irt_pred = has_irt_pred ?
            Float32(main.irt_pred[row]) :
            pools.irt_by_pid[Int(receiver_pid)]
        receiver_n_scans = Float32(scan_count_column[row])
        receiver_temporal_mean = _mbr_sqrt_tuple(
            temporal_mean_columns,
            row,
        )
        receiver_temporal_trace = temporal_trace_column[row]
        receiver_corr_mask = UInt8(corr_mask_column[row])
        columns[:MBR_log2_weight_lod_ratio][row] =
            _mbr_log2_weight_lod_ratio(
                receiver_weight,
                receiver_file,
                lod_log2_weight_by_file,
                lod_log2_weight_global,
            )
        columns[:MBR_receiver_frag_corr_bitvec_rank][row] =
            Float32(corr_rank_column[row])

        true_values = _mbr_feature_values(
            receiver_pid,
            receiver_weight,
            receiver_explained,
            receiver_irt_pred,
            receiver_irt,
            receiver_n_scans,
            receiver_temporal_mean,
            receiver_temporal_trace,
            receiver_corr_mask,
            true_donor,
            receiver_file,
            run_similarity_atlas,
            receiver_run_clusters,
            bitvec_rank_table,
            donor_dict,
        )
        columns[:MBR_best_is_missing_true][row] = false
        for (feature_idx, stem) in enumerate(MBR_PAIRED_FEATURE_STEMS)
            columns[_mbr_true_feature(stem)][row] = true_values[feature_idx]
        end

        target_irt = receiver_irt_pred
        false_donors = _mbr_false_donors(
            donor_dict,
            pools,
            eligibility,
            receiver_pid,
            receiver_file,
            target_irt,
            true_donor,
            run_similarity_atlas,
        )
        for counterfactual_idx in 1:MBR_N_COUNTERFACTUALS
            false_donor = false_donors[counterfactual_idx]
            false_donor === nothing && continue
            false_values = _mbr_feature_values(
                receiver_pid,
                receiver_weight,
                receiver_explained,
                receiver_irt_pred,
                receiver_irt,
                receiver_n_scans,
                receiver_temporal_mean,
                receiver_temporal_trace,
                receiver_corr_mask,
                false_donor,
                receiver_file,
                run_similarity_atlas,
                receiver_run_clusters,
                bitvec_rank_table,
                donor_dict,
            )
            columns[_mbr_missing_feature(counterfactual_idx)][row] = false
            for (feature_idx, stem) in enumerate(MBR_PAIRED_FEATURE_STEMS)
                columns[_mbr_false_feature(
                    stem,
                    counterfactual_idx,
                )][row] = false_values[feature_idx]
            end
        end
    end

    sidecar = DataFrame()
    sidecar[!, :precursor_idx] = columns[:precursor_idx]
    sidecar[!, :scan_idx] = columns[:scan_idx]
    sidecar[!, :MBR_best_is_missing_true] =
        columns[:MBR_best_is_missing_true]
    for feature in MBR_SHARED_FEATURES
        sidecar[!, feature] = columns[feature]
    end
    for stem in MBR_PAIRED_FEATURE_STEMS
        sidecar[!, _mbr_true_feature(stem)] =
            columns[_mbr_true_feature(stem)]
    end
    for counterfactual_idx in 1:MBR_N_COUNTERFACTUALS
        missing_col = _mbr_missing_feature(counterfactual_idx)
        sidecar[!, missing_col] = columns[missing_col]
        for stem in MBR_PAIRED_FEATURE_STEMS
            feature = _mbr_false_feature(stem, counterfactual_idx)
            sidecar[!, feature] = columns[feature]
        end
    end
    writeArrow(main_path * MBR_SIDECAR_SUFFIX, sidecar)
    return main_path * MBR_SIDECAR_SUFFIX
end
