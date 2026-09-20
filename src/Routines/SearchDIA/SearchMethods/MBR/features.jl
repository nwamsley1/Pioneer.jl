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
    started = last_progress = time()
    @debug_l1 "MBR donor threshold collection starting: files=$(length(file_paths)) initial_pass=$require_initial_pass"
    scores = Float32[]
    targets = Bool[]
    for (file_idx, path) in enumerate(file_paths)
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
        if time() - last_progress >= 60
            @debug_l1 "MBR donor threshold collection: files=$file_idx/$(length(file_paths)) rows=$(length(scores)) elapsed=$(round(time() - started, digits=2))s"
            last_progress = time()
        end
    end
    @debug_l1 "MBR donor threshold collection complete: rows=$(length(scores)) elapsed=$(round(time() - started, digits=2))s"
    isempty(scores) && return Inf32
    started = time()
    @debug_l1 "MBR donor threshold q-values starting: rows=$(length(scores))"
    qvalues = similar(scores)
    get_qvalues!(scores, targets, qvalues)
    eligible = targets .& (qvalues .<= donor_q_threshold)
    score_floor = any(eligible) ? minimum(scores[eligible]) : Inf32
    @debug_l1 "MBR donor threshold q-values complete: score_floor=$score_floor elapsed=$(round(time() - started, digits=2))s"
    return score_floor
end

function _collect_mbr_donor_files!(
    donor_files::Dict{UInt32, Tuple{UInt32, UInt32}},
    precursor_ids, file_ids, scores, score_floor::Float32,
)
    for (precursor_id, file_id, score) in zip(precursor_ids, file_ids, scores)
        Float32(score) >= score_floor || continue
        pid, file_idx = UInt32(precursor_id), UInt32(file_id)
        donors = get!(donor_files, pid, (file_idx, file_idx))
        # Two distinct runs suffice to find a donor outside any receiver run.
        if donors[1] == donors[2] && file_idx != donors[1]
            donor_files[pid] = (donors[1], file_idx)
        end
    end
    return donor_files
end

function _mbr_preintegration_donor_files(
    file_paths::Vector{String},
    score_floor::Float32,
)
    donor_files = Dict{UInt32, Tuple{UInt32, UInt32}}()
    started = last_progress = time()
    rows_processed = 0
    for (file_idx, path) in enumerate(file_paths)
        tbl = Arrow.Table(path)
        _collect_mbr_donor_files!(
            donor_files, tbl.precursor_idx, tbl.ms_file_idx, tbl.trace_prob_prepass, score_floor,
        )
        rows_processed += length(tbl.precursor_idx)
        if time() - last_progress >= 60
            @debug_l1 "MBR donor indexing: files=$file_idx/$(length(file_paths)) rows=$rows_processed precursors=$(length(donor_files)) elapsed=$(round(time() - started, digits=2))s"
            last_progress = time()
        end
    end
    return donor_files
end

@inline function _mbr_has_cross_run_donor(
    donor_files::Dict{UInt32, Tuple{UInt32, UInt32}},
    pid::UInt32,
    receiver_file::UInt32,
)
    entries = get(donor_files, pid, nothing)
    entries === nothing && return false
    return entries[1] != receiver_file || entries[2] != receiver_file
end

function _collect_mbr_integrated_donors!(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    previous_files::Set{UInt32}, columns, frag_columns,
    score_floor::Float32, q_value_threshold::Float32,
)
    positions = Dict{Tuple{UInt32, UInt32}, Int}()
    current_files = Set{UInt32}()
    @inbounds for row in 1:length(columns.precursor_idx)
        _mbr_initial_pass(columns.qval[row], columns.global_qval[row], q_value_threshold) ||
            continue
        score = Float32(columns.trace_prob[row])
        score >= score_floor || continue
        pid = UInt32(columns.precursor_idx[row])
        file_idx = UInt32(columns.ms_file_idx[row])
        entries = get!(() -> _MBRDonorEntry[], donor_dict, pid)
        key = (pid, file_idx)
        position = get(positions, key, 0)
        if position == 0 && file_idx in previous_files
            existing = findfirst(donor -> donor.ms_file_idx == file_idx, entries)
            position = existing === nothing ? 0 : existing
        end
        positions[key] = position == 0 ? length(entries) + 1 : position
        push!(current_files, file_idx)
        position != 0 && !(score > entries[position].trace_prob) && continue
        irt_obs = Float32(columns.irt_obs[row])
        donor = _MBRDonorEntry(
            score, pid, Float32(columns.weight[row]), Float32(columns.explained[row]),
            Float32(columns.irt_pred[row]) - irt_obs, irt_obs,
            Float32(columns.n_scans[row]), _mbr_sqrt_tuple(frag_columns, row),
            UInt8(columns.frag_mask[row]), UInt16(columns.frag_rank[row]), file_idx,
        )
        if position == 0
            push!(entries, donor)
        else
            entries[position] = donor
        end
    end
    union!(previous_files, current_files)
    return donor_dict
end

@inline function _mbr_lower_weight(left::_MBRDonorEntry, right::_MBRDonorEntry)
    left_weight = isfinite(left.weight) ? left.weight : Inf32
    right_weight = isfinite(right.weight) ? right.weight : Inf32
    return left_weight < right_weight ||
        (left_weight == right_weight && left.trace_prob < right.trace_prob)
end

function _MBRDonorIndex(entries::Dict{UInt32, Vector{_MBRDonorEntry}})
    lookups = Dict{UInt32, _MBRDonorLookup}()
    files = Set{UInt32}()
    for (pid, donors) in entries
        length(donors) <= typemax(UInt32) || error("Too many donor runs for precursor $pid")
        top = (UInt32(0), UInt32(0))
        lowest = (UInt32(0), UInt32(0), UInt32(0))
        for i in eachindex(donors)
            donor = donors[i]
            push!(files, donor.ms_file_idx)
            position = UInt32(i)
            if top[1] == 0 || donor.trace_prob > donors[top[1]].trace_prob
                top = (position, top[1])
            elseif top[2] == 0 || donor.trace_prob > donors[top[2]].trace_prob
                top = (top[1], position)
            end
            if lowest[1] == 0 || _mbr_lower_weight(donor, donors[lowest[1]])
                lowest = (position, lowest[1], lowest[2])
            elseif lowest[2] == 0 || _mbr_lower_weight(donor, donors[lowest[2]])
                lowest = (lowest[1], position, lowest[2])
            elseif lowest[3] == 0 || _mbr_lower_weight(donor, donors[lowest[3]])
                lowest = (lowest[1], lowest[2], position)
            end
        end
        order = issorted(donors; by=donor -> donor.ms_file_idx) ? UInt32[] :
            UInt32.(sortperm(donors; by=donor -> donor.ms_file_idx))
        lookups[pid] = _MBRDonorLookup(order, top, lowest)
    end
    return _MBRDonorIndex(entries, lookups, sort!(collect(files)))
end

function build_mbr_integrated_donor_dict(
    file_paths::Vector{String},
    score_floor::Float32;
    q_value_threshold::Float32,
)
    donor_dict = Dict{UInt32, Vector{_MBRDonorEntry}}()
    previous_files = Set{UInt32}()
    started = last_progress = time()
    rows_processed = 0
    for (file_position, path) in enumerate(file_paths)
        tbl = Arrow.Table(path)
        required = (
            :precursor_idx, :ms_file_idx, :trace_prob_prepass, :qval, :global_qval,
            :irt_pred, MBR_INTEGRATED_WEIGHT_COLUMN,
            MBR_INTEGRATED_LOG2_INTENSITY_EXPLAINED_COLUMN,
            MBR_INTEGRATED_APEX_IRT_COLUMN, MBR_INTEGRATED_FRAG_CORR_BITVEC_COLUMN,
            MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_BITVEC_RANK_COLUMN,
            MBR_INTEGRATED_N_SCANS_COLUMN,
        )
        for col in required
            hasproperty(tbl, col) ||
                error("Integrated MBR donor selection requires column $col in $path")
        end
        columns = (
            precursor_idx=tbl.precursor_idx, ms_file_idx=tbl.ms_file_idx,
            trace_prob=tbl.trace_prob_prepass, qval=tbl.qval, global_qval=tbl.global_qval,
            irt_pred=tbl.irt_pred, weight=getproperty(tbl, MBR_INTEGRATED_WEIGHT_COLUMN),
            explained=getproperty(tbl, MBR_INTEGRATED_LOG2_INTENSITY_EXPLAINED_COLUMN),
            irt_obs=getproperty(tbl, MBR_INTEGRATED_APEX_IRT_COLUMN),
            frag_mask=getproperty(tbl, MBR_INTEGRATED_FRAG_CORR_BITVEC_COLUMN),
            frag_rank=getproperty(tbl, MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_BITVEC_RANK_COLUMN),
            n_scans=getproperty(tbl, MBR_INTEGRATED_N_SCANS_COLUMN),
        )
        frag_columns = ntuple(rank -> getproperty(tbl, MBR_INTEGRATED_FRAGMENT_SQRT_COLUMNS[rank]), 8)
        _collect_mbr_integrated_donors!(
            donor_dict, previous_files, columns, frag_columns, score_floor, q_value_threshold,
        )
        rows_processed += length(tbl.precursor_idx)
        if time() - last_progress >= 60
            @debug_l1 "Post-integration MBR donor collection: files=$file_position/$(length(file_paths)) rows=$rows_processed precursors=$(length(donor_dict)) elapsed=$(round(time() - started, digits=2))s"
            last_progress = time()
        end
    end
    index_started = time()
    @debug_l1 "Post-integration MBR donor lookup construction starting: precursors=$(length(donor_dict))"
    index = _MBRDonorIndex(donor_dict)
    @debug_l1 "Post-integration MBR donor lookup construction complete: elapsed=$(round(time() - index_started, digits=2))s"
    return index
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

@inline function _mbr_donor_position(entries, lookup::_MBRDonorLookup, file_idx::UInt32)
    low, high = 1, length(entries)
    while low <= high
        middle = (low + high) >>> 1
        position = isempty(lookup.file_order) ? middle : Int(lookup.file_order[middle])
        donor_file = entries[position].ms_file_idx
        donor_file == file_idx && return position
        if donor_file < file_idx
            low = middle + 1
        else
            high = middle - 1
        end
    end
    return 0
end

function _mbr_receiver_donors(index::_MBRDonorIndex, receiver_file::UInt32, atlas)
    ranked = Tuple{Float32, UInt32}[
        (_mbr_run_similarity(atlas, receiver_file, file_idx), file_idx)
        for file_idx in index.file_ids if file_idx != receiver_file
    ]
    sort!(ranked; by=first, rev=true)
    equal_similarity = isempty(ranked) || all(entry -> entry[1] == ranked[1][1], ranked)
    finite_similarity = all(entry -> isfinite(entry[1]), ranked)
    return _MBRReceiverDonors(index, receiver_file, ranked, equal_similarity, finite_similarity)
end

_mbr_receiver_donors(donors::Dict{UInt32, Vector{_MBRDonorEntry}}, ::UInt32, atlas) = donors

@inline function _mbr_donor_in_file(index::_MBRDonorIndex, pid::UInt32, file_idx::UInt32)
    entries = get(index.entries, pid, nothing)
    entries === nothing && return nothing
    position = _mbr_donor_position(entries, index.lookups[pid], file_idx)
    return position == 0 ? nothing : entries[position]
end

@inline function _mbr_top_scoring_donor(index::_MBRDonorIndex, pid::UInt32, receiver_file::UInt32)
    entries = get(index.entries, pid, nothing)
    entries === nothing && return nothing
    for position in index.lookups[pid].top_scores
        position == 0 && continue
        donor = entries[position]
        donor.ms_file_idx != receiver_file && return donor
    end
    return nothing
end

@inline function _mbr_worst_alternate_donor(
    index::_MBRDonorIndex, pid::UInt32, receiver_file::UInt32, best_donor::_MBRDonorEntry,
)
    entries = get(index.entries, pid, nothing)
    entries === nothing && return nothing
    for position in index.lookups[pid].lowest_weights
        position == 0 && continue
        donor = entries[position]
        donor.ms_file_idx == receiver_file && continue
        donor.ms_file_idx == best_donor.ms_file_idx && continue
        return donor
    end
    return nothing
end

_mbr_select_donor(index::_MBRDonorIndex, pid::UInt32, receiver_file::UInt32, atlas) =
    _mbr_select_donor(index.entries, pid, receiver_file, atlas)

function _mbr_select_donor(
    donors::_MBRReceiverDonors, pid::UInt32, receiver_file::UInt32, atlas,
)
    index = donors.index
    receiver_file == donors.receiver_file && donors.finite_similarity ||
        return _mbr_select_donor(index.entries, pid, receiver_file, atlas)
    donors.equal_similarity && return _mbr_top_scoring_donor(index, pid, receiver_file)
    entries = get(index.entries, pid, nothing)
    entries === nothing && return nothing
    lookup = index.lookups[pid]
    ranked = donors.ranked_files
    n_donors = length(entries)
    n_donors == 0 && return nothing
    # Sparse precursors can be cheaper to scan than to probe through the run ranking.
    if n_donors <= (length(ranked) ÷ n_donors) * ndigits(n_donors; base=2)
        return _mbr_select_donor(index.entries, pid, receiver_file, atlas)
    end
    rank = 1
    while rank <= length(ranked)
        similarity = ranked[rank][1]
        best_position = 0
        while rank <= length(ranked) && ranked[rank][1] == similarity
            position = _mbr_donor_position(entries, lookup, ranked[rank][2])
            if position != 0 && (best_position == 0 ||
               entries[position].trace_prob > entries[best_position].trace_prob ||
               (entries[position].trace_prob == entries[best_position].trace_prob && position < best_position))
                best_position = position
            end
            rank += 1
        end
        best_position != 0 && return entries[best_position]
    end
    return nothing
end

_mbr_donor_in_file(donors::_MBRReceiverDonors, pid::UInt32, file_idx::UInt32) =
    _mbr_donor_in_file(donors.index, pid, file_idx)
_mbr_top_scoring_donor(donors::_MBRReceiverDonors, pid::UInt32, receiver_file::UInt32) =
    _mbr_top_scoring_donor(donors.index, pid, receiver_file)
_mbr_worst_alternate_donor(
    donors::_MBRReceiverDonors, pid::UInt32, receiver_file::UInt32, best_donor::_MBRDonorEntry,
) = _mbr_worst_alternate_donor(donors.index, pid, receiver_file, best_donor)

@inline function _mbr_pid_already_taken(
    donors::Vector{Union{Nothing, _MBRDonorEntry}},
    n_found::Int,
    pid::UInt32,
)
    @inbounds for k in 1:n_found
        donor = donors[k]
        donor === nothing && continue
        donor.precursor_idx == pid && return true
    end
    return false
end

@inline function _mbr_collect_false_donors_from_pool!(
    donors::Vector{Union{Nothing, _MBRDonorEntry}},
    donor_dict::_MBRDonorCollection,
    pool::_MBRIrtPool,
    target_irt::Float32,
    receiver_pid::UInt32,
    receiver_file::UInt32,
    eligibility::_MBRCounterfactualEligibility,
    atlas::Union{Nothing, RunSimilarityAtlas};
    donor_file::Union{Nothing, UInt32} = nothing,
)
    n_found = count(donor -> donor !== nothing, donors)
    cursor = _mbr_pool_cursor(pool, target_irt)
    @inbounds while true
        pool_idx = _mbr_next_pool_idx!(cursor, pool, target_irt)
        pool_idx == 0 && break
        pid = pool.pids[pool_idx]
        pid == receiver_pid && continue
        # `seen_pids` was a per-row BitSet, but every pid pushed to it was exactly the pid of a
        # donor just written to `donors` — so the set is derivable from `donors[1:n_found]`, a
        # scan of at most MBR_N_COUNTERFACTUALS (3) entries. A BitSet spans min..max of its
        # contents in a Vector{UInt64}, and precursor ids run into the millions, so the old
        # version allocated a large bitmap per row to hold <= 3 scattered values.
        _mbr_pid_already_taken(donors, n_found, pid) && continue
        _mbr_counterfactual_eligible(eligibility, receiver_file, pid) ||
            continue
        donor = donor_file === nothing ?
            _mbr_select_donor(donor_dict, pid, receiver_file, atlas) :
            _mbr_donor_in_file(donor_dict, pid, donor_file)
        donor === nothing && continue
        n_found += 1
        donors[n_found] = donor
        n_found == MBR_N_COUNTERFACTUALS && break
    end
    return n_found
end

function _mbr_false_donors(
    donor_dict::_MBRDonorCollection,
    pools::_MBRPartnerPools,
    eligibility::_MBRCounterfactualEligibility,
    receiver_pid::UInt32,
    receiver_file::UInt32,
    target_irt::Float32,
    true_donor::_MBRDonorEntry,
    atlas::Union{Nothing, RunSimilarityAtlas},
    donors::Vector{Union{Nothing, _MBRDonorEntry}},
)
    pid_int = Int(receiver_pid)
    fold = Int(pools.fold_by_pid[pid_int])
    mz_bin = Int(pools.mz_bin_by_pid[pid_int])
    charge = Int(pools.charge_by_pid[pid_int])
    sequence_length = Int(pools.length_by_pid[pid_int])
    fill!(donors, nothing)      # caller-owned buffer, reused across rows
    n_found = 0

    same_file_pool = get(
        pools.file_charge_length_pool,
        (true_donor.ms_file_idx, charge, sequence_length),
        _empty_mbr_irt_pool(),
    )
    n_found = _mbr_collect_false_donors_from_pool!(
        donors,
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
    donor_dict::_MBRDonorCollection,
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

# Candidate feature columns, with paired features stored in block-major order.
struct _MBRSidecarColumns
    precursor_idx::Vector{UInt32}
    scan_idx::Vector{UInt32}
    shared::Vector{Vector{Float32}}
    missing_flags::Vector{BitVector}
    paired::Vector{Vector{Float32}}
end

function _mbr_sidecar_columns(n::Int)
    return _MBRSidecarColumns(
        zeros(UInt32, n),
        zeros(UInt32, n),
        Vector{Float32}[fill(-1.0f0, n) for _ in eachindex(MBR_SHARED_FEATURES)],
        BitVector[trues(n) for _ in 1:(MBR_N_COUNTERFACTUALS + 1)],
        Vector{Float32}[
            fill(-1.0f0, n) for _ in eachindex(MBR_PAIRED_COLUMN_NAMES)
        ],
    )
end

# Positions of the two shared features within `shared`, resolved once at load time.
const MBR_SHARED_LOD_RATIO_IDX =
    findfirst(==(:MBR_log2_weight_lod_ratio), MBR_SHARED_FEATURES)
const MBR_SHARED_CORR_RANK_IDX =
    findfirst(==(:MBR_receiver_frag_corr_bitvec_rank), MBR_SHARED_FEATURES)

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

function _mbr_candidate_donors(
    qvals, global_qvals, pids, file_indices, donor_dict, contexts,
    atlas, threshold,
)
    rows = Int64[]
    donors = _MBRDonorEntry[]
    @inbounds for row in 1:length(qvals)
        qval = Float32(qvals[row])
        global_qval = Float32(global_qvals[row])
        isfinite(global_qval) && global_qval <= threshold || continue
        isfinite(qval) && qval <= threshold && continue
        receiver_file = UInt32(file_indices[row])
        context = get!(contexts, receiver_file) do
            _mbr_receiver_donors(donor_dict, receiver_file, atlas)
        end
        donor = _mbr_select_donor(context, UInt32(pids[row]), receiver_file, atlas)
        donor === nothing && continue
        push!(rows, row)
        push!(donors, donor)
    end
    return rows, donors
end

function compute_postintegration_mbr_features!(
    main_path::String,
    donor_dict::_MBRDonorCollection,
    pools::_MBRPartnerPools,
    eligibility::_MBRCounterfactualEligibility;
    run_similarity_atlas::Union{Nothing, RunSimilarityAtlas},
    receiver_run_clusters::_MBRReceiverRunClusters,
    lod_log2_weight_by_file::Dict{UInt32, Float32},
    lod_log2_weight_global::Float32,
    bitvec_rank_table = nothing,
    q_value_threshold::Float32,
    stats = nothing,
)
    started = time()
    main = Arrow.Table(main_path)
    n = length(main.precursor_idx)
    contexts = Dict{UInt32, _MBRDonorCollection}()
    rows, true_donors = _mbr_candidate_donors(
        main.qval, main.global_qval, main.precursor_idx, main.ms_file_idx,
        donor_dict, contexts, run_similarity_atlas, q_value_threshold,
    )
    selection_seconds = time() - started
    features_started = time()
    columns = _mbr_sidecar_columns(length(rows))
    columns.precursor_idx .= main.precursor_idx[rows]
    columns.scan_idx .= main.scan_idx[rows]
    fill!(columns.missing_flags[1], false)

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

    # Each file reuses its own counterfactual donor buffer.
    false_donor_buffer = Union{Nothing, _MBRDonorEntry}[
        nothing for _ in 1:MBR_N_COUNTERFACTUALS
    ]

    _rdiag = get(ENV, "PIONEER_MBR_ROW_DIAG", "0") == "1"
    _ra = Base.gc_bytes(); _rt = time()

    @inbounds for (candidate_row, row) in enumerate(rows)
        receiver_pid = UInt32(main.precursor_idx[row])
        receiver_file = UInt32(main.ms_file_idx[row])
        donor_context = contexts[receiver_file]
        true_donor = true_donors[candidate_row]

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
        columns.shared[MBR_SHARED_LOD_RATIO_IDX][candidate_row] =
            _mbr_log2_weight_lod_ratio(
                receiver_weight,
                receiver_file,
                lod_log2_weight_by_file,
                lod_log2_weight_global,
            )
        columns.shared[MBR_SHARED_CORR_RANK_IDX][candidate_row] =
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
            donor_context,
        )
        @inbounds for feature_idx in 1:MBR_N_PAIRED
            columns.paired[feature_idx][candidate_row] = true_values[feature_idx]
        end
        if _rdiag
            MBR_ROW_DIAG[:true_feat_bytes] += Base.gc_bytes() - _ra
            MBR_ROW_DIAG[:true_feat_ms] += round(Int, (time() - _rt) * 1000)
            _ra = Base.gc_bytes(); _rt = time()
        end

        target_irt = receiver_irt_pred
        false_donors = _mbr_false_donors(
            donor_context,
            pools,
            eligibility,
            receiver_pid,
            receiver_file,
            target_irt,
            true_donor,
            run_similarity_atlas,
            false_donor_buffer,
        )
        if _rdiag
            MBR_ROW_DIAG[:false_select_bytes] += Base.gc_bytes() - _ra
            MBR_ROW_DIAG[:false_select_ms] += round(Int, (time() - _rt) * 1000)
            MBR_ROW_DIAG[:n_rows_with_donor] += 1
            _ra = Base.gc_bytes(); _rt = time()
        end
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
                donor_context,
            )
            columns.missing_flags[counterfactual_idx + 1][candidate_row] = false
            offset = counterfactual_idx * MBR_N_PAIRED
            @inbounds for feature_idx in 1:MBR_N_PAIRED
                columns.paired[offset + feature_idx][candidate_row] =
                    false_values[feature_idx]
            end
        end
        if _rdiag
            MBR_ROW_DIAG[:false_feat_bytes] += Base.gc_bytes() - _ra
            MBR_ROW_DIAG[:false_feat_ms] += round(Int, (time() - _rt) * 1000)
            _ra = Base.gc_bytes(); _rt = time()
        end
    end
    if _rdiag
        MBR_ROW_DIAG[:n_files] += 1
        _mbr_row_diag_report()
    end

    sidecar = DataFrame(row_idx = rows)
    sidecar[!, :precursor_idx] = columns.precursor_idx
    sidecar[!, :scan_idx] = columns.scan_idx
    sidecar[!, MBR_MISSING_COLUMN_NAMES[1]] = columns.missing_flags[1]
    for (shared_idx, feature) in enumerate(MBR_SHARED_FEATURES)
        sidecar[!, feature] = columns.shared[shared_idx]
    end
    for feature_idx in 1:MBR_N_PAIRED
        sidecar[!, MBR_PAIRED_COLUMN_NAMES[feature_idx]] =
            columns.paired[feature_idx]
    end
    for counterfactual_idx in 1:MBR_N_COUNTERFACTUALS
        sidecar[!, MBR_MISSING_COLUMN_NAMES[counterfactual_idx + 1]] =
            columns.missing_flags[counterfactual_idx + 1]
        offset = counterfactual_idx * MBR_N_PAIRED
        for feature_idx in 1:MBR_N_PAIRED
            sidecar[!, MBR_PAIRED_COLUMN_NAMES[offset + feature_idx]] =
                columns.paired[offset + feature_idx]
        end
    end
    feature_seconds = time() - features_started
    write_started = time()
    writeArrow(main_path * MBR_SIDECAR_SUFFIX, sidecar)
    if stats !== nothing
        stats[] = (
            rows = n,
            candidates = length(rows),
            selection_seconds = selection_seconds,
            feature_seconds = feature_seconds,
            write_seconds = time() - write_started,
        )
    end
    return main_path * MBR_SIDECAR_SUFFIX
end


# Accumulators for the row-loop diagnostic above. Threaded (parallel_foreach! over files), so these
# counts are approximate under contention — they are for attribution, not exact accounting.
const MBR_ROW_DIAG = Dict{Symbol, Int}(
    :true_feat_bytes => 0,    :true_feat_ms => 0,
    :false_select_bytes => 0, :false_select_ms => 0,
    :false_feat_bytes => 0,   :false_feat_ms => 0,
    :n_rows_with_donor => 0,  :n_files => 0,
)

function _mbr_row_diag_report()
    d = MBR_ROW_DIAG
    gb(x) = round(x / 2^30, digits = 2)
    tot = d[:true_feat_bytes] + d[:false_select_bytes] + d[:false_feat_bytes]
    pct(x) = tot > 0 ? round(100 * x / tot, digits = 1) : 0.0
    @user_info """
    MBR row-loop diagnostic ($(d[:n_files]) file(s), $(d[:n_rows_with_donor]) rows with a donor):
      true-donor featurisation   : $(gb(d[:true_feat_bytes])) GB  $(d[:true_feat_ms]) ms  ($(pct(d[:true_feat_bytes]))%)
      counterfactual SELECTION   : $(gb(d[:false_select_bytes])) GB  $(d[:false_select_ms]) ms  ($(pct(d[:false_select_bytes]))%)
      counterfactual featurisation: $(gb(d[:false_feat_bytes])) GB  $(d[:false_feat_ms]) ms  ($(pct(d[:false_feat_bytes]))%)
      TOTAL row loop             : $(gb(tot)) GB"""
    return nothing
end
