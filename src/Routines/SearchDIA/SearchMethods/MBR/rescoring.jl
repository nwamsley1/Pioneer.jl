# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

"""
    load_postintegration_mbr_candidates(file_paths, q_value_threshold)

Materialise **only the MBR candidate rows** across all files, plus what the caller needs to expand
results back to full length later.

Previously this concatenated every row of every file into one DataFrame and the candidate mask was
applied afterwards. Candidates are only ~10% of rows (76,075 of 742,179 on a 6-file run), so the
frame was ~10x larger than anything downstream used -- and because `parts` and the `vcat` result
coexist, peak was ~2x that again. Extrapolated to 100 files that is ~9 GB to hand ~0.5 GB of
candidates to the model.

Safe because every input to the decision is row-local: `_mbr_candidate_mask_kernel` reads only that
row's q-values and missing flag, the baseline predicate is row-local, and
`_mbr_add_hellinger_contrasts!` compares a row's own four blocks. All cross-run features were already
resolved into each file's sidecar by `compute_postintegration_mbr_features!` before this runs.

Returns `(candidates, masks, n_rows, base_targets, base_decoys)`; `masks[i]` and `n_rows[i]` describe
file `i` so `_write_mbr_recovery_sidecars!` can scatter per-candidate results back over all rows.
"""
function load_postintegration_mbr_candidates(
    file_paths::Vector{String},
    q_value_threshold::Float32,
)
    parts = DataFrame[]
    masks = BitVector[]
    n_rows = Int[]
    base_targets = 0
    base_decoys = 0
    for path in file_paths
        pass1_path = path * PASS1_SIDECAR_SUFFIX
        mbr_path = path * MBR_SIDECAR_SUFFIX
        isfile(pass1_path) || error("Missing MBR Pass-1 sidecar at $pass1_path")
        isfile(mbr_path) || error("Missing MBR feature sidecar at $mbr_path")
        main = Arrow.Table(path)
        pass1 = Arrow.Table(pass1_path)
        mbr = Arrow.Table(mbr_path)
        n = length(main.precursor_idx)
        length(pass1.precursor_idx) == n ||
            error("Pass-1 sidecar row-count mismatch at $pass1_path")
        length(mbr.precursor_idx) == n ||
            error("MBR sidecar row-count mismatch at $mbr_path")
        @inbounds for row in 1:n
            (
                main.precursor_idx[row] == pass1.precursor_idx[row] ==
                    mbr.precursor_idx[row] &&
                main.scan_idx[row] == pass1.scan_idx[row] ==
                    mbr.scan_idx[row]
            ) || error("MBR sidecar misalignment at row $row of $path")
        end

        # Same predicate as the old whole-frame path, evaluated here so only survivors materialise.
        mask = _mbr_candidate_mask_kernel(
            main.global_qval,
            main.qval,
            Tables.getcolumn(mbr, :MBR_best_is_missing_true),
            q_value_threshold,
        )
        # Baseline counts need every row, but they are only two integers. Predicate must match
        # apply_postintegration_mbr_rescoring! exactly: BOTH q-values under threshold.
        @inbounds for row in 1:n
            qv = Float32(main.qval[row])
            gqv = Float32(main.global_qval[row])
            (isfinite(qv) && qv <= q_value_threshold &&
             isfinite(gqv) && gqv <= q_value_threshold) || continue
            Bool(main.target[row]) ? (base_targets += 1) : (base_decoys += 1)
        end
        push!(masks, mask)
        push!(n_rows, n)
        if !any(mask)
            continue
        end

        frame = DataFrame(
            precursor_idx = UInt32.(main.precursor_idx[mask]),
            scan_idx = UInt32.(main.scan_idx[mask]),
            ms_file_idx = UInt32.(main.ms_file_idx[mask]),
            cv_fold = UInt8.(main.cv_fold[mask]),
            target = Bool.(main.target[mask]),
            qval = Float32.(main.qval[mask]),
            global_qval = Float32.(main.global_qval[mask]),
            trace_prob_prepass = Float32.(pass1.trace_prob_prepass[mask]),
            trace_prob_infold = Float32.(pass1.trace_prob_infold[mask]),
        )
        for feature in MBR_RECEIVER_FEATURES
            feature === :trace_prob_infold && continue
            hasproperty(main, feature) || continue
            frame[!, feature] = Tables.getcolumn(main, feature)[mask]
        end
        for feature in Symbol.(Tables.columnnames(mbr))
            feature in (:precursor_idx, :scan_idx) && continue
            frame[!, feature] = Tables.getcolumn(mbr, feature)[mask]
        end
        push!(parts, frame)
    end
    candidates = isempty(parts) ? DataFrame() : vcat(parts...; cols = :union)
    return (
        candidates = candidates,
        masks = masks,
        n_rows = n_rows,
        base_targets = base_targets,
        base_decoys = base_decoys,
    )
end

function load_postintegration_mbr_frame(file_paths::Vector{String})
    parts = DataFrame[]
    for path in file_paths
        pass1_path = path * PASS1_SIDECAR_SUFFIX
        mbr_path = path * MBR_SIDECAR_SUFFIX
        isfile(pass1_path) || error("Missing MBR Pass-1 sidecar at $pass1_path")
        isfile(mbr_path) || error("Missing MBR feature sidecar at $mbr_path")

        main = Arrow.Table(path)
        pass1 = Arrow.Table(pass1_path)
        mbr = Arrow.Table(mbr_path)
        n = length(main.precursor_idx)
        length(pass1.precursor_idx) == n ||
            error("Pass-1 sidecar row-count mismatch at $pass1_path")
        length(mbr.precursor_idx) == n ||
            error("MBR sidecar row-count mismatch at $mbr_path")
        @inbounds for row in 1:n
            (
                main.precursor_idx[row] == pass1.precursor_idx[row] ==
                    mbr.precursor_idx[row] &&
                main.scan_idx[row] == pass1.scan_idx[row] ==
                    mbr.scan_idx[row]
            ) || error("MBR sidecar misalignment at row $row of $path")
        end

        frame = DataFrame(
            precursor_idx = collect(UInt32.(main.precursor_idx)),
            scan_idx = collect(UInt32.(main.scan_idx)),
            ms_file_idx = collect(UInt32.(main.ms_file_idx)),
            cv_fold = collect(UInt8.(main.cv_fold)),
            target = collect(Bool.(main.target)),
            qval = collect(Float32.(main.qval)),
            global_qval = collect(Float32.(main.global_qval)),
            trace_prob_prepass = collect(Float32.(pass1.trace_prob_prepass)),
            trace_prob_infold = collect(Float32.(pass1.trace_prob_infold)),
        )
        for feature in MBR_RECEIVER_FEATURES
            feature === :trace_prob_infold && continue
            hasproperty(main, feature) || continue
            frame[!, feature] = collect(Tables.getcolumn(main, feature))
        end
        for feature in Symbol.(Tables.columnnames(mbr))
            feature in (:precursor_idx, :scan_idx) && continue
            frame[!, feature] = collect(Tables.getcolumn(mbr, feature))
        end
        push!(parts, frame)
    end
    return isempty(parts) ? DataFrame() : vcat(parts...; cols = :union)
end

# Runs over the FULL staged frame, so it is the widest of the MBR row loops. `frame.col` infers as
# AbstractVector and `frame[row, col]` additionally does a Dict{Symbol,Int} lookup per cell, so the
# loop below is moved behind a barrier that receives the columns. Measured at n = 391,370 with a
# 100-column frame: 300.4 ms / 71.7 MB inline versus 0.8 ms / 0.1 MB here, identical output.
function _mbr_candidate_mask(
    frame::DataFrame,
    q_value_threshold::Float32,
)
    return _mbr_candidate_mask_kernel(
        frame.global_qval,
        frame.qval,
        frame[!, :MBR_best_is_missing_true],
        q_value_threshold,
    )
end

function _mbr_candidate_mask_kernel(
    global_qval_col,
    qval_col,
    missing_true_col,
    q_value_threshold::Float32,
)
    n = length(global_qval_col)
    candidates = falses(n)
    @inbounds for row in 1:n
        # Each q-value column used to be fetched twice per row (isfinite, then the comparison).
        global_qval = Float32(global_qval_col[row])
        run_qval = Float32(qval_col[row])
        global_pass = isfinite(global_qval) && global_qval <= q_value_threshold
        run_pass = isfinite(run_qval) && run_qval <= q_value_threshold
        # Left unconditional rather than short-circuited: Bool(missing) throws, so short-circuiting
        # would change behaviour if this column were ever missing-typed (it is a BitVector today).
        true_present = !Bool(missing_true_col[row])
        candidates[row] =
            global_pass && !run_pass && true_present
    end
    return candidates
end

function _mbr_available_feature_sets(frame::DataFrame)
    true_features = Symbol[]
    false_features = [
        Symbol[] for _ in 1:MBR_N_COUNTERFACTUALS
    ]
    for feature in MBR_FTR_FEATURES_TRUE
        hasproperty(frame, feature) || continue
        if feature in MBR_RECEIVER_FEATURES ||
           feature in MBR_SHARED_FEATURES
            all_present = all(
                counterfactual_idx ->
                    hasproperty(frame, feature),
                1:MBR_N_COUNTERFACTUALS,
            )
            all_present || continue
            push!(true_features, feature)
            for counterfactual_idx in 1:MBR_N_COUNTERFACTUALS
                push!(false_features[counterfactual_idx], feature)
            end
            continue
        end

        stem = nothing
        feature_string = String(feature)
        for candidate_stem in MBR_MODEL_PAIRED_FEATURE_STEMS
            if feature_string == candidate_stem * "_true"
                stem = candidate_stem
                break
            end
        end
        stem === nothing && continue
        mapped = Symbol[
            _mbr_false_feature(stem, counterfactual_idx)
            for counterfactual_idx in 1:MBR_N_COUNTERFACTUALS
        ]
        all(feature_name -> hasproperty(frame, feature_name), mapped) ||
            continue
        push!(true_features, feature)
        for counterfactual_idx in 1:MBR_N_COUNTERFACTUALS
            push!(
                false_features[counterfactual_idx],
                mapped[counterfactual_idx],
            )
        end
    end
    isempty(true_features) &&
        error("No complete post-integration MBR feature pairs are available")
    return true_features, false_features
end

@inline function _mbr_block_feature(
    stem::AbstractString,
    counterfactual_idx::Int,
)
    return counterfactual_idx == 0 ?
        _mbr_true_feature(stem) :
        _mbr_false_feature(stem, counterfactual_idx)
end

# Rank and margin of one pairing against the other three, for a single stem. Extracted so the column
# accesses specialise: `frame[row, col]` is a two-argument getindex doing a Dict{Symbol,Int} lookup
# PER CELL, and this runs 4 stems x 4 blocks x ~4 comparisons per row. Measured on the candidate
# frame (n = 40,150): 474.9 ms / 92.7 MB inline versus 3.5 ms / 5.0 MB here, identical output.
#
# `cols` is a Tuple so it arrives concretely typed -- tuples are covariant, so the kernel specialises
# on the real column types even though the caller can only infer `Tuple`.
function _mbr_contrast_rank_margin(cols::Tuple, source_idx::Int, n::Int)
    ranks = fill(-1.0f0, n)
    margins = fill(-1.0f0, n)
    source = cols[source_idx]
    @inbounds for row in 1:n
        source_value = Float32(source[row])
        isfinite(source_value) && source_value >= 0.0f0 || continue
        rank = 1
        best_other = Inf32
        for other_idx in eachindex(cols)
            other_idx == source_idx && continue
            other_value = Float32(cols[other_idx][row])
            isfinite(other_value) && other_value >= 0.0f0 || continue
            other_value < source_value && (rank += 1)
            other_value < best_other && (best_other = other_value)
        end
        ranks[row] = Float32(rank)
        margins[row] = isfinite(best_other) ?
            best_other - source_value :
            -1.0f0
    end
    return ranks, margins
end

function _mbr_add_hellinger_contrasts!(frame::DataFrame)
    n = nrow(frame)
    for base in MBR_HELLINGER_CONTRAST_BASE_STEMS
        source_columns = Symbol[
            _mbr_block_feature(base, counterfactual_idx)
            for counterfactual_idx in 0:MBR_N_COUNTERFACTUALS
        ]
        all(column -> hasproperty(frame, column), source_columns) ||
            continue
        # Fetched once per stem rather than per cell; a Tuple so the kernel specialises (see above).
        cols = Tuple(frame[!, column] for column in source_columns)

        for counterfactual_idx in 0:MBR_N_COUNTERFACTUALS
            rank_column = _mbr_block_feature(
                base * "_rank",
                counterfactual_idx,
            )
            margin_column = _mbr_block_feature(
                base * "_margin",
                counterfactual_idx,
            )
            ranks, margins = _mbr_contrast_rank_margin(
                cols,
                counterfactual_idx + 1,
                n,
            )
            frame[!, rank_column] = ranks
            frame[!, margin_column] = margins
        end
    end
    return frame
end

function _mbr_eval_mask_and_labels(
    scores::Vector{Float32},
    present::BitMatrix,
    top_labels::BitVector,
    top_mask::BitVector,
)
    n_candidates, n_counterfactuals = size(present)
    n_blocks = 1 + n_counterfactuals
    length(scores) == n_blocks * n_candidates ||
        error("MBR score-frame size mismatch")
    eval_mask = falses(length(scores))
    labels = falses(length(scores))
    @inbounds for candidate_idx in 1:n_candidates
        eval_mask[candidate_idx] = top_mask[candidate_idx]
        labels[candidate_idx] = top_labels[candidate_idx]
        best_false_idx = 0
        best_false_score = -Inf32
        for counterfactual_idx in 1:n_counterfactuals
            present[candidate_idx, counterfactual_idx] || continue
            frame_idx =
                counterfactual_idx * n_candidates + candidate_idx
            score = scores[frame_idx]
            rank_score = isfinite(score) ? score : -Inf32
            if best_false_idx == 0 || rank_score > best_false_score
                best_false_idx = frame_idx
                best_false_score = rank_score
            end
        end
        best_false_idx != 0 && (eval_mask[best_false_idx] = true)
    end
    return eval_mask, labels
end

function _mbr_iteration_metrics(
    scores::Vector{Float32},
    present::BitMatrix,
    target_top::BitVector,
    top_labels::BitVector,
    top_mask::BitVector,
)
    n_candidates = length(target_top)
    eval_mask, labels = _mbr_eval_mask_and_labels(
        scores,
        present,
        top_labels,
        top_mask,
    )
    qvalues = fill(Inf32, length(scores))
    peps = fill(Inf32, length(scores))
    eval_indices = findall(eval_mask)
    if !isempty(eval_indices)
        eval_scores = scores[eval_indices]
        eval_labels = labels[eval_indices]
        eval_qvalues = Vector{Float32}(undef, length(eval_indices))
        eval_peps = Vector{Float32}(undef, length(eval_indices))
        get_qvalues!(eval_scores, eval_labels, eval_qvalues)
        get_PEP!(eval_scores, eval_labels, eval_peps)
        qvalues[eval_indices] .= eval_qvalues
        peps[eval_indices] .= eval_peps
    end
    positive_top =
        BitVector(qvalues[1:n_candidates] .<=
                  MBR_SEMISUPERVISED_FTR_THRESHOLD) .&
        target_top
    return (
        qvalues = qvalues,
        peps = peps,
        eval_mask = eval_mask,
        positive_top = positive_top,
        n_positive = count(positive_top),
    )
end

function _mbr_evenly_spaced_sample(
    rows::Vector{Int},
    limit::Int,
)
    limit >= 0 || throw(ArgumentError("MBR training limit must be nonnegative"))
    length(rows) <= limit && return copy(rows)
    limit == 0 && return Int[]
    n = length(rows)
    return Int[
        rows[fld((2 * sample_idx - 1) * n, 2 * limit) + 1]
        for sample_idx in 1:limit
    ]
end

"""
    _mbr_sample_positions(n, limit) -> Vector{Int}

The 1-based positions `_mbr_evenly_spaced_sample` would pick out of a length-`n` list, without
needing the list. Lets the capped branch of `_mbr_training_rows` emit only the rows it keeps rather
than materialising every eligible row index first: selecting 1.875M negatives out of, say, 60M
eligible ones used to allocate a 60M-element `Vector{Int}` (480 MB) purely to index into it.

Positions are strictly increasing, so a single ordered pass over the candidates can match them off
in lockstep.
"""
function _mbr_sample_positions(n::Int, limit::Int)
    limit >= 0 || throw(ArgumentError("MBR training limit must be nonnegative"))
    (n <= limit) && return collect(1:n)
    limit == 0 && return Int[]
    return Int[fld((2 * sample_idx - 1) * n, 2 * limit) + 1 for sample_idx in 1:limit]
end

function _mbr_proportional_caps(
    counts::Vector{Int},
    limit::Int,
)
    limit >= 0 || throw(ArgumentError("MBR training limit must be nonnegative"))
    total = sum(counts)
    total <= limit && return copy(counts)
    total == 0 && return zeros(Int, length(counts))

    caps = Int[fld(count * limit, total) for count in counts]
    remainders = Int[mod(count * limit, total) for count in counts]
    remaining = limit - sum(caps)
    order = sortperm(
        eachindex(counts);
        by = idx -> (-remainders[idx], idx),
    )
    for idx in order
        remaining == 0 && break
        caps[idx] < counts[idx] || continue
        caps[idx] += 1
        remaining -= 1
    end
    return caps
end

function _mbr_training_rows(
    folds::Vector{UInt8},
    positive_top::BitVector,
    receiver_decoy_top::BitVector,
    present::BitMatrix,
    train_fold::UInt8;
    max_positives::Int = MBR_MAX_POSITIVE_TRAIN_PER_FOLD,
    max_negatives::Int = MBR_MAX_NEGATIVE_TRAIN_PER_FOLD,
)
    n_candidates, n_counterfactuals = size(present)
    length(receiver_decoy_top) == n_candidates ||
        error("MBR receiver-decoy mask length mismatch")
    available_positives = 0
    available_receiver_decoys = 0
    negative_counts = zeros(Int, n_counterfactuals)
    @inbounds for candidate_idx in 1:n_candidates
        folds[candidate_idx] == train_fold || continue
        available_positives += positive_top[candidate_idx]
        available_receiver_decoys += receiver_decoy_top[candidate_idx]
        for counterfactual_idx in 1:n_counterfactuals
            present[candidate_idx, counterfactual_idx] || continue
            negative_counts[counterfactual_idx] += 1
        end
    end
    available_negatives =
        available_receiver_decoys + sum(negative_counts)

    if available_positives <= max_positives &&
       available_negatives <= max_negatives
        train_rows = Int[]
        train_labels = Bool[]
        sizehint!(
            train_rows,
            available_positives + available_negatives,
        )
        sizehint!(
            train_labels,
            available_positives + available_negatives,
        )
        @inbounds for candidate_idx in 1:n_candidates
            folds[candidate_idx] == train_fold || continue
            if positive_top[candidate_idx]
                push!(train_rows, candidate_idx)
                push!(train_labels, true)
            elseif receiver_decoy_top[candidate_idx]
                push!(train_rows, candidate_idx)
                push!(train_labels, false)
            end
            for counterfactual_idx in 1:n_counterfactuals
                present[candidate_idx, counterfactual_idx] || continue
                push!(
                    train_rows,
                    counterfactual_idx * n_candidates + candidate_idx,
                )
                push!(train_labels, false)
            end
        end
        return (
            rows = train_rows,
            labels = train_labels,
            available_positives = available_positives,
            used_positives = available_positives,
            available_negatives = available_negatives,
            used_negatives = available_negatives,
            available_receiver_decoys = available_receiver_decoys,
            used_receiver_decoys = available_receiver_decoys,
            available_negatives_by_counterfactual = negative_counts,
            used_negatives_by_counterfactual = copy(negative_counts),
        )
    end

    # Capped branch. Emit only the rows we keep: the previous version pushed every eligible row
    # index into per-class Vector{Int}s and then sampled those down, so hitting the cap allocated
    # 8 bytes per *available* row to select `limit` of them. `_mbr_sample_positions` gives the
    # positions up front, and each class is matched off in one ordered pass. Selection is
    # bit-identical -- same positions, same order.
    all_negative_counts = vcat(
        available_receiver_decoys,
        negative_counts,
    )
    all_negative_caps = _mbr_proportional_caps(
        all_negative_counts,
        max_negatives,
    )
    negative_caps = all_negative_caps[2:end]

    positive_positions = _mbr_sample_positions(available_positives, max_positives)
    receiver_positions = _mbr_sample_positions(available_receiver_decoys, all_negative_caps[1])
    counterfactual_positions = [
        _mbr_sample_positions(negative_counts[k], negative_caps[k])
        for k in 1:n_counterfactuals
    ]

    sampled_positives = Vector{Int}(undef, length(positive_positions))
    sampled_receiver_decoys = Vector{Int}(undef, length(receiver_positions))
    sampled_by_counterfactual = [
        Vector{Int}(undef, length(counterfactual_positions[k]))
        for k in 1:n_counterfactuals
    ]
    # Running ordinal within each class, and how far into that class's position list we are.
    seen_pos = 0; take_pos = 1
    seen_rcv = 0; take_rcv = 1
    seen_cf = zeros(Int, n_counterfactuals)
    take_cf = ones(Int, n_counterfactuals)
    @inbounds for candidate_idx in 1:n_candidates
        folds[candidate_idx] == train_fold || continue
        if positive_top[candidate_idx]
            seen_pos += 1
            if take_pos <= length(positive_positions) &&
               positive_positions[take_pos] == seen_pos
                sampled_positives[take_pos] = candidate_idx
                take_pos += 1
            end
        end
        if receiver_decoy_top[candidate_idx]
            seen_rcv += 1
            if take_rcv <= length(receiver_positions) &&
               receiver_positions[take_rcv] == seen_rcv
                sampled_receiver_decoys[take_rcv] = candidate_idx
                take_rcv += 1
            end
        end
        for counterfactual_idx in 1:n_counterfactuals
            present[candidate_idx, counterfactual_idx] || continue
            seen_cf[counterfactual_idx] += 1
            positions = counterfactual_positions[counterfactual_idx]
            take = take_cf[counterfactual_idx]
            if take <= length(positions) &&
               positions[take] == seen_cf[counterfactual_idx]
                sampled_by_counterfactual[counterfactual_idx][take] =
                    counterfactual_idx * n_candidates + candidate_idx
                take_cf[counterfactual_idx] = take + 1
            end
        end
    end

    sampled_negatives = Int[]
    sizehint!(
        sampled_negatives,
        length(sampled_receiver_decoys) + sum(negative_caps),
    )
    append!(sampled_negatives, sampled_receiver_decoys)
    for counterfactual_idx in 1:n_counterfactuals
        append!(sampled_negatives, sampled_by_counterfactual[counterfactual_idx])
    end

    train_rows = vcat(sampled_positives, sampled_negatives)
    train_labels = vcat(
        trues(length(sampled_positives)),
        falses(length(sampled_negatives)),
    )
    return (
        rows = train_rows,
        labels = train_labels,
        available_positives = available_positives,
        used_positives = length(sampled_positives),
        available_negatives = available_negatives,
        used_negatives = length(sampled_negatives),
        available_receiver_decoys = available_receiver_decoys,
        used_receiver_decoys = length(sampled_receiver_decoys),
        available_negatives_by_counterfactual = negative_counts,
        used_negatives_by_counterfactual = negative_caps,
    )
end

# Row indices into `x` for each CV fold: the block-0 row of every candidate held out in that fold,
# plus each counterfactual block that actually exists for it.
#
# This used to be rebuilt inside _mbr_fit_oof_iteration, i.e. on all 16 calls (2 folds x 8
# semi-supervised iterations), from `folds`, `present` and n_candidates -- all three fixed for the
# whole search, so all 16 builds produced identical vectors. It is now built once. Measured at
# N = 391,370 with 3 counterfactuals: 487.5 ms / 271.2 MB rebuilt versus 11.9 ms / 33.9 MB hoisted.
# The inline version also had no sizehint!, so most of that allocation was geometric regrowth; the
# count-then-fill below removes it.
#
# Only `positive_top` changes between iterations (the semi-supervised label update), and it does not
# enter this calculation -- so the hoist is value-preserving, not an approximation.
function _mbr_test_rows_by_fold(folds::Vector{UInt8}, present::BitMatrix)
    n_candidates, n_counterfactuals = size(present)
    rows_by_fold = Vector{Vector{Int}}()
    for test_fold in UInt8[0, 1]
        n_rows = 0
        @inbounds for candidate_idx in 1:n_candidates
            folds[candidate_idx] == test_fold || continue
            n_rows += 1
            for counterfactual_idx in 1:n_counterfactuals
                present[candidate_idx, counterfactual_idx] && (n_rows += 1)
            end
        end
        rows = Vector{Int}(undef, n_rows)
        position = 0
        @inbounds for candidate_idx in 1:n_candidates
            folds[candidate_idx] == test_fold || continue
            position += 1
            rows[position] = candidate_idx
            for counterfactual_idx in 1:n_counterfactuals
                present[candidate_idx, counterfactual_idx] || continue
                position += 1
                rows[position] =
                    counterfactual_idx * n_candidates + candidate_idx
            end
        end
        push!(rows_by_fold, rows)
    end
    return rows_by_fold
end

function _mbr_fit_oof_iteration(
    x::Matrix{Float32},
    folds::Vector{UInt8},
    positive_top::BitVector,
    receiver_decoy_top::BitVector,
    present::BitMatrix,
    test_rows_by_fold::Vector{Vector{Int}},
)
    n_candidates, n_counterfactuals = size(present)
    n_blocks = 1 + n_counterfactuals
    scores = fill(NaN32, n_blocks * n_candidates)
    last_classifier = nothing

    for (fold_position, test_fold) in enumerate(UInt8[0, 1])
        train_fold = UInt8(1) - test_fold
        test_rows = test_rows_by_fold[fold_position]
        isempty(test_rows) && continue
        training = _mbr_training_rows(
            folds,
            positive_top,
            receiver_decoy_top,
            present,
            train_fold,
        )
        train_rows = training.rows
        train_labels = training.labels
        capped =
            training.used_positives < training.available_positives ||
            training.used_negatives < training.available_negatives
        used_counterfactuals =
            sum(training.used_negatives_by_counterfactual)
        available_counterfactuals =
            sum(training.available_negatives_by_counterfactual)
        @debug_l1 "MBR transfer model OOF fold $(Int(test_fold)): " *
                  "train positives=$(training.used_positives)/" *
                  "$(training.available_positives), " *
                  "receiver decoys=$(training.used_receiver_decoys)/" *
                  "$(training.available_receiver_decoys), " *
                  "counterfactual negatives=$used_counterfactuals/" *
                  "$available_counterfactuals" *
                  (capped ? " (capped)" : "")
        if isempty(train_rows) || length(unique(train_labels)) < 2
            scores[test_rows] .= isempty(train_labels) ?
                0.5f0 :
                (first(train_labels) ? 1.0f0 : 0.0f0)
            continue
        end

        classifier = build_lightgbm_classifier(; SHARED_LGBM_HP...)
        LightGBM.fit!(
            classifier,
            x[train_rows, :],
            _prepare_labels(train_labels);
            verbosity = -1,
        )
        raw = LightGBM.predict(classifier, x[test_rows, :])
        predictions = ndims(raw) == 2 ?
            dropdims(raw; dims = 2) :
            raw
        scores[test_rows] .= Float32.(predictions)
        last_classifier = classifier
    end
    return scores, last_classifier
end

# Function barrier for the scatter below. `candidate_frame[!, col]` infers as AbstractVector, so
# indexing it in a loop IN THE SAME FUNCTION makes every element a dynamic dispatch that boxes its
# result. Measured on a representative 200,000-row x 104-column frame with 12 threads, identical
# output: 806.7 ms / 713.3 MB inline versus 9.7 ms / 79.4 MB through this barrier (82.8x) -- for a
# matrix that is itself only 79.3 MB, i.e. the inline version allocated 9x its own output in boxes.
@inline function _mbr_fill_feature_block!(
    x::Matrix{Float32},
    column,
    offset::Int,
    n_candidates::Int,
    feature_idx::Int,
)
    @inbounds for candidate_idx in 1:n_candidates
        value = column[candidate_idx]
        x[offset + candidate_idx, feature_idx] =
            value === missing ? 0.0f0 : Float32(value)
    end
    return nothing
end

function _mbr_feature_matrix(
    candidate_frame::DataFrame,
    true_features::Vector{Symbol},
    false_features::Vector{Vector{Symbol}},
)
    n_candidates = nrow(candidate_frame)
    n_features = length(true_features)
    n_blocks = 1 + length(false_features)
    x = Matrix{Float32}(
        undef,
        n_blocks * n_candidates,
        n_features,
    )
    Threads.@threads for feature_idx in 1:n_features
        _mbr_fill_feature_block!(
            x,
            candidate_frame[!, true_features[feature_idx]],
            0,
            n_candidates,
            feature_idx,
        )
        for counterfactual_idx in eachindex(false_features)
            _mbr_fill_feature_block!(
                x,
                candidate_frame[
                    !,
                    false_features[counterfactual_idx][feature_idx],
                ],
                counterfactual_idx * n_candidates,
                n_candidates,
                feature_idx,
            )
        end
    end
    return x
end

function _mbr_semisupervised_oof(
    candidate_frame::DataFrame,
    true_features::Vector{Symbol},
    false_features::Vector{Vector{Symbol}},
)
    n_candidates = nrow(candidate_frame)
    x = _mbr_feature_matrix(
        candidate_frame,
        true_features,
        false_features,
    )
    present = falses(n_candidates, MBR_N_COUNTERFACTUALS)
    @inbounds for counterfactual_idx in 1:MBR_N_COUNTERFACTUALS
        missing = candidate_frame[
            !,
            _mbr_missing_feature(counterfactual_idx),
        ]
        for candidate_idx in 1:n_candidates
            present[candidate_idx, counterfactual_idx] =
                !Bool(missing[candidate_idx])
        end
    end
    target_top = BitVector(candidate_frame.target)
    positive_top = copy(target_top)
    receiver_decoy_top = .!target_top
    eval_top_labels = trues(n_candidates)
    eval_top_mask = trues(n_candidates)
    folds = Vector{UInt8}(candidate_frame.cv_fold)
    # Fixed for the whole search: depends only on folds and present, neither of which the
    # semi-supervised loop mutates.
    test_rows_by_fold = _mbr_test_rows_by_fold(folds, present)
    best_state = nothing
    previous_positive = -1

    for iteration in 1:MBR_SEMISUPERVISED_MAX_ITERATIONS
        scores, classifier = _mbr_fit_oof_iteration(
            x,
            folds,
            positive_top,
            receiver_decoy_top,
            present,
            test_rows_by_fold,
        )
        metrics = _mbr_iteration_metrics(
            scores,
            present,
            target_top,
            eval_top_labels,
            eval_top_mask,
        )
        current = (
            iteration = iteration,
            scores = scores,
            metrics = metrics,
            classifier = classifier,
        )
        if best_state === nothing ||
           metrics.n_positive >= best_state.metrics.n_positive
            best_state = current
        end
        @debug_l1 "MBR transfer model iter $iteration: " *
                  "training positives=$(count(positive_top)), " *
                  "receiver decoys=$(count(receiver_decoy_top)); " *
                  "FTR≤$(100 * MBR_SEMISUPERVISED_FTR_THRESHOLD)% " *
                  "target transfers=$(metrics.n_positive)"
        if iteration > 1 && !_scoring_target_gain_sufficient(
            previous_positive,
            metrics.n_positive,
        )
            @debug_l1 "MBR transfer model stopping: target-transfer count " *
                      "did not improve by " *
                      "$(100 * SCORING_SEMISUPERVISED_MIN_TARGET_GAIN)% " *
                      "over $previous_positive; using iteration " *
                      "$(best_state.iteration) with " *
                      "$(best_state.metrics.n_positive) target transfers"
            break
        elseif iteration == MBR_SEMISUPERVISED_MAX_ITERATIONS
            @debug_l1 "MBR transfer model stopping: hit max iterations; " *
                      "using iteration $(best_state.iteration) with " *
                      "$(best_state.metrics.n_positive) target transfers"
            break
        end
        previous_positive = metrics.n_positive
        positive_top = metrics.positive_top
        eval_top_labels = copy(positive_top)
        eval_top_mask = positive_top .| receiver_decoy_top
    end
    return best_state, present
end

function _mbr_top_counterfactual_scores(
    scores::Vector{Float32},
    present::BitMatrix,
)
    n_candidates, n_counterfactuals = size(present)
    top_scores = fill(-Inf32, n_candidates)
    @inbounds for candidate_idx in 1:n_candidates
        for counterfactual_idx in 1:n_counterfactuals
            present[candidate_idx, counterfactual_idx] || continue
            score = scores[
                counterfactual_idx * n_candidates + candidate_idx
            ]
            if isfinite(score) && score > top_scores[candidate_idx]
                top_scores[candidate_idx] = score
            end
        end
    end
    return top_scores
end

function _mbr_combined_error_recovery(
    real_scores::Vector{Float32},
    target_top::BitVector,
    false_scores::Vector{Float32};
    base_targets::Int,
    base_decoys::Int,
    alpha::Float32,
)
    n = length(real_scores)
    combined_qvalues = fill(Inf32, n)
    combined_rates = fill(Inf32, n)
    recovered = falses(n)
    baseline_rate =
        base_targets > 0 ? Float32(base_decoys / base_targets) : Inf32
    if n == 0 || base_targets <= 0 || baseline_rate >= alpha
        return (
            recovered = recovered,
            qvalues = combined_qvalues,
            rates = combined_rates,
            threshold = Inf32,
            base_targets = base_targets,
            base_decoys = base_decoys,
            baseline_error_rate = baseline_rate,
            mbr_targets = 0,
            mbr_decoys = 0,
            false_transfers = 0,
            total_targets = base_targets,
            total_errors = base_decoys,
            combined_error_rate = baseline_rate,
        )
    end

    order = sort(
        findall(isfinite, real_scores);
        by = idx -> real_scores[idx],
        rev = true,
    )
    false_event_scores = Float32[]
    @inbounds for idx in 1:n
        target_top[idx] || continue
        isfinite(real_scores[idx]) && isfinite(false_scores[idx]) || continue
        push!(
            false_event_scores,
            min(real_scores[idx], false_scores[idx]),
        )
    end
    sort!(false_event_scores; rev = true)

    group_ranges = UnitRange{Int}[]
    group_rates = Float32[]
    n_targets = 0
    n_decoys = 0
    false_pointer = 0
    position = 1
    while position <= length(order)
        group_start = position
        threshold = real_scores[order[position]]
        while position <= length(order) &&
              real_scores[order[position]] == threshold
            if target_top[order[position]]
                n_targets += 1
            else
                n_decoys += 1
            end
            position += 1
        end
        while false_pointer < length(false_event_scores) &&
              false_event_scores[false_pointer + 1] >= threshold
            false_pointer += 1
        end
        total_targets = base_targets + n_targets
        total_errors = base_decoys + n_decoys + false_pointer
        push!(group_ranges, group_start:(position - 1))
        push!(
            group_rates,
            total_targets > 0 ?
                Float32(total_errors / total_targets) :
                Inf32,
        )
    end

    running_min = Inf32
    for group_idx in reverse(eachindex(group_rates))
        running_min = min(running_min, group_rates[group_idx])
        for position_idx in group_ranges[group_idx]
            candidate_idx = order[position_idx]
            combined_qvalues[candidate_idx] = running_min
            combined_rates[candidate_idx] = group_rates[group_idx]
        end
    end
    recovered .= combined_qvalues .<= alpha
    threshold = any(recovered) ?
        minimum(real_scores[recovered]) :
        Inf32

    mbr_targets = count(recovered .& target_top)
    mbr_decoys = count(recovered .& .!target_top)
    false_transfers = 0
    if isfinite(threshold)
        @inbounds for idx in 1:n
            target_top[idx] || continue
            if isfinite(real_scores[idx]) &&
               isfinite(false_scores[idx]) &&
               real_scores[idx] >= threshold &&
               false_scores[idx] >= threshold
                false_transfers += 1
            end
        end
    end
    total_targets = base_targets + mbr_targets
    total_errors = base_decoys + mbr_decoys + false_transfers
    combined_rate =
        total_targets > 0 ? Float32(total_errors / total_targets) : Inf32
    return (
        recovered = recovered,
        qvalues = combined_qvalues,
        rates = combined_rates,
        threshold = threshold,
        base_targets = base_targets,
        base_decoys = base_decoys,
        baseline_error_rate = baseline_rate,
        mbr_targets = mbr_targets,
        mbr_decoys = mbr_decoys,
        false_transfers = false_transfers,
        total_targets = total_targets,
        total_errors = total_errors,
        combined_error_rate = combined_rate,
    )
end

"""
    apply_postintegration_mbr_rescoring!(frame; alpha, q_value_threshold)

Train the paired, out-of-fold transfer model and accept recoveries under one
combined precursor error budget. Existing passing rows are never rescored or
removed here.
"""
function apply_postintegration_mbr_rescoring!(
    frame::DataFrame;
    alpha::Float32,
    q_value_threshold::Float32,
    baseline_counts::Union{Nothing, Tuple{Int, Int}} = nothing,
    frame_is_candidates::Bool = false,
)
    n = nrow(frame)
    frame[!, :mbr_recovered] = falses(n)
    frame[!, :MBR_transfer_candidate] = falses(n)
    frame[!, :mbr_target_decoy_prob] = fill(NaN32, n)
    frame[!, :ftr_qval_true] = fill(NaN32, n)
    frame[!, :ftr_pep_true] = fill(NaN32, n)
    frame[!, :mbr_total_error_qval_true] = fill(NaN32, n)
    frame[!, :mbr_total_error_rate_true] = fill(NaN32, n)
    n == 0 && return (
        n_candidates = 0,
        n_recovered = 0,
        base_targets = 0,
        base_decoys = 0,
        baseline_error_rate = 0.0f0,
        mbr_targets = 0,
        mbr_decoys = 0,
        mbr_false_transfers = 0,
        internal_ftr_targets = 0,
        internal_ftr_errors = 0,
        internal_ftr_estimate = NaN32,
        total_targets = 0,
        total_errors = 0,
        combined_error_rate = 0.0f0,
    )

    # Counts over ALL rows. When the caller pre-filtered to candidates it accumulated these
    # per-file during the load, because they cannot be recovered from candidates alone.
    base_targets, base_decoys = if baseline_counts === nothing
        baseline = BitVector(undef, n)
        @inbounds for row in 1:n
            baseline[row] =
                isfinite(Float32(frame.qval[row])) &&
                Float32(frame.qval[row]) <= q_value_threshold &&
                isfinite(Float32(frame.global_qval[row])) &&
                Float32(frame.global_qval[row]) <= q_value_threshold
        end
        (count(baseline .& BitVector(frame.target)),
         count(baseline .& .!BitVector(frame.target)))
    else
        baseline_counts
    end

    # A pre-filtered frame holds candidates only, so re-deriving the mask would be redundant work
    # over a frame that is by construction all-true.
    candidate_mask = frame_is_candidates ? trues(n) :
        _mbr_candidate_mask(frame, q_value_threshold)
    frame[!, :MBR_transfer_candidate] = candidate_mask
    candidate_indices = findall(candidate_mask)
    if isempty(candidate_indices)
        baseline_rate =
            base_targets > 0 ? Float32(base_decoys / base_targets) : 0.0f0
        return (
            n_candidates = 0,
            n_recovered = 0,
            base_targets = base_targets,
            base_decoys = base_decoys,
            baseline_error_rate = baseline_rate,
            mbr_targets = 0,
            mbr_decoys = 0,
            mbr_false_transfers = 0,
            internal_ftr_targets = 0,
            internal_ftr_errors = 0,
            internal_ftr_estimate = NaN32,
            total_targets = base_targets,
            total_errors = base_decoys,
            combined_error_rate = baseline_rate,
        )
    end

    candidates = frame[candidate_indices, :]
    _mbr_add_hellinger_contrasts!(candidates)
    true_features, false_features = _mbr_available_feature_sets(candidates)
    best_state, present = _mbr_semisupervised_oof(
        candidates,
        true_features,
        false_features,
    )
    n_candidates = nrow(candidates)
    real_scores = copy(best_state.scores[1:n_candidates])
    evaluated_top = @view best_state.metrics.eval_mask[1:n_candidates]
    @inbounds for candidate_idx in 1:n_candidates
        evaluated_top[candidate_idx] ||
            (real_scores[candidate_idx] = NaN32)
    end
    false_scores = _mbr_top_counterfactual_scores(
        best_state.scores,
        present,
    )
    combined = _mbr_combined_error_recovery(
        real_scores,
        BitVector(candidates.target),
        false_scores;
        base_targets = base_targets,
        base_decoys = base_decoys,
        alpha = alpha,
    )

    @inbounds for (candidate_position, row) in enumerate(candidate_indices)
        frame[row, :mbr_target_decoy_prob] =
            combined.recovered[candidate_position] ?
            real_scores[candidate_position] :
            NaN32
        frame[row, :mbr_recovered] =
            combined.recovered[candidate_position]
        frame[row, :ftr_qval_true] =
            best_state.metrics.qvalues[candidate_position]
        frame[row, :ftr_pep_true] =
            best_state.metrics.peps[candidate_position]
        frame[row, :mbr_total_error_qval_true] =
            combined.qvalues[candidate_position]
        frame[row, :mbr_total_error_rate_true] =
            combined.rates[candidate_position]
    end

    internal_targets = combined.mbr_targets
    internal_errors = combined.false_transfers
    internal_ftr =
        internal_targets > 0 ?
        Float32(internal_errors / internal_targets) :
        NaN32
    @debug_l1 "MBR combined-error acceptance: baseline=$(base_decoys)/$(base_targets), " *
              "recovered targets=$(combined.mbr_targets), decoys=$(combined.mbr_decoys), " *
              "false transfers=$(combined.false_transfers), " *
              "total=$(combined.total_errors)/$(combined.total_targets)"

    if best_state.classifier !== nothing
        model = LightGBMModel(
            best_state.classifier,
            true_features,
            nothing,
        )
        gains = importance(model)
        if gains !== nothing
            @debug_l1 "MBR transfer-model feature importances (gain):"
            for (feature, gain) in sort(gains; by = x -> -x[2])
                @debug_l1 "  $(rpad(string(feature), 48)) $(round(gain, digits=2))"
            end
        end
    end

    return (
        n_candidates = n_candidates,
        n_recovered = count(combined.recovered),
        base_targets = base_targets,
        base_decoys = base_decoys,
        baseline_error_rate = combined.baseline_error_rate,
        mbr_targets = combined.mbr_targets,
        mbr_decoys = combined.mbr_decoys,
        mbr_false_transfers = combined.false_transfers,
        internal_ftr_targets = internal_targets,
        internal_ftr_errors = internal_errors,
        internal_ftr_estimate = internal_ftr,
        total_targets = combined.total_targets,
        total_errors = combined.total_errors,
        combined_error_rate = combined.combined_error_rate,
    )
end
