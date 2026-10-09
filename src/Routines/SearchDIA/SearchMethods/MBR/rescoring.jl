# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# Candidate features stay in bounded blocks; only labels, IDs and scores span the search.
struct _MBRFeatureStore
    data::DataFrameBlockStore
    schema::DataFrame
end
_MBRFeatureStore(path; budget::Int=64*1024^2) =
    _MBRFeatureStore(DataFrameBlockStore(path, budget), DataFrame())
Base.close(store::_MBRFeatureStore) = close(store.data)

function _mbr_store_features!(store::_MBRFeatureStore, frame::DataFrame)
    columns = unique(vcat(MBR_FTR_FEATURES_TRUE,
        (_mbr_ftr_features_false(k) for k in 1:MBR_N_COUNTERFACTUALS)...))
    filter!(col -> hasproperty(frame, col), columns)
    block = select(frame, columns; copycols=true)
    append!(store.schema, block[1:0, :]; cols=:union)
    _store_dataframe_block!(store.data, block)
    return nothing
end

function _mbr_copy_feature_rows!(dest, src, positions)
    @inbounds for j in axes(src, 2), i in eachindex(positions)
        dest[positions[i], j] = src[i, j]
    end
end

function _mbr_gather_feature_rows(store::_MBRFeatureStore, true_features, false_features,
                                  rows::AbstractVector{Int}, n_candidates::Int)
    dest = Matrix{Float32}(undef, length(rows), length(true_features))
    positions = [Int[] for _ in store.data.row_ends]
    for (position, row) in enumerate(rows)
        1 <= row <= (1 + length(false_features)) * n_candidates || throw(BoundsError(rows, position))
        candidate = mod(row - 1, n_candidates) + 1
        push!(positions[searchsortedfirst(store.data.row_ends, candidate)], position)
    end
    columns = unique(vcat(true_features, false_features...))
    for index in eachindex(positions)
        selected = positions[index]
        isempty(selected) && continue
        first_row = index == 1 ? 1 : store.data.row_ends[index-1] + 1
        n = store.data.row_ends[index] - first_row + 1
        local_rows = [div(rows[p] - 1, n_candidates) * n +
                      mod(rows[p] - 1, n_candidates) + 2 - first_row for p in selected]
        block = _dataframe_block(store.data, index)
        # Match the unioned in-memory table: absent feature values gather as zero.
        if any(col -> !hasproperty(block, col), columns)
            block = copy(block)
            for col in columns
                hasproperty(block, col) || (block[!, col] = zeros(Float32, n))
            end
        end
        gathered = _mbr_gather_feature_rows(block, true_features, false_features, local_rows, n)
        _mbr_copy_feature_rows!(dest, gathered, selected)
    end
    return dest
end

function _mbr_validate_sidecar_rows(
    main_pids, main_scans, pass1_pids, pass1_scans,
    mbr_pids, mbr_scans, row_indices, path,
)
    n = length(main_pids)
    length(pass1_pids) == length(pass1_scans) == n ||
        error("Pass-1 sidecar row-count mismatch at $path")
    @inbounds for row in 1:n
        (main_pids[row] == pass1_pids[row] &&
         main_scans[row] == pass1_scans[row]) ||
            error("Pass-1 sidecar misalignment at row $row of $path")
    end
    length(mbr_pids) == length(mbr_scans) == length(row_indices) ||
        error("MBR sidecar row-count mismatch at $path")
    previous = 0
    for (sidecar_row, row) in enumerate(row_indices)
        row isa Integer && previous < row <= n ||
            error("Invalid MBR sidecar row index $row at $path")
        (main_pids[row] == mbr_pids[sidecar_row] &&
         main_scans[row] == mbr_scans[sidecar_row]) ||
            error("MBR sidecar misalignment at row $row of $path")
        previous = row
    end
    return nothing
end

function _mbr_feature_row_indices(main, pass1, mbr, path)
    rows = hasproperty(mbr, :row_idx) ? mbr.row_idx : (1:length(main.precursor_idx))
    _mbr_validate_sidecar_rows(
        main.precursor_idx, main.scan_idx, pass1.precursor_idx, pass1.scan_idx,
        mbr.precursor_idx, mbr.scan_idx, rows, path,
    )
    return rows
end

function _mbr_candidate_sidecar_mask(
    global_qvals, qvals, missing_true, rows, sparse::Bool, threshold, path,
)
    mask = falses(length(qvals))
    selected = Int[]
    for (sidecar_row, row) in enumerate(rows)
        global_qval = Float32(global_qvals[row])
        qval = Float32(qvals[row])
        candidate = isfinite(global_qval) && global_qval <= threshold &&
            !(isfinite(qval) && qval <= threshold) && !Bool(missing_true[sidecar_row])
        sparse && !candidate && error("Ineligible MBR sidecar candidate at row $row of $path")
        candidate || continue
        mask[row] = true
        push!(selected, sidecar_row)
    end
    return mask, selected
end

function _mbr_baseline_counts(qvals, global_qvals, targets, threshold)
    base_targets = base_decoys = 0
    @inbounds for row in 1:length(qvals)
        qval, global_qval = Float32(qvals[row]), Float32(global_qvals[row])
        (isfinite(qval) && qval <= threshold &&
         isfinite(global_qval) && global_qval <= threshold) || continue
        Bool(targets[row]) ? (base_targets += 1) : (base_decoys += 1)
    end
    return base_targets, base_decoys
end

"""
    load_postintegration_mbr_candidates(file_paths, q_value_threshold)

Load candidate features and retain per-file masks for scattering recovery results
back to their original rows. Accepts indexed candidate sidecars and legacy dense
sidecars. Baseline target/decoy counts are accumulated from all main-table rows.
With `feature_store`, return candidate metadata only and store model features in
blocks of at most `feature_batch_size` rows. The cache budget excludes the current
input block, gathered training/prediction matrices, and Arrow input buffers.
"""
function load_postintegration_mbr_candidates(
    file_paths::Vector{String},
    q_value_threshold::Float32;
    feature_store::Union{Nothing, _MBRFeatureStore} = nothing,
    feature_batch_size::Int = 50_000,
)
    feature_batch_size > 0 || throw(ArgumentError("feature_batch_size must be positive"))
    metadata = DataFrame()
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
        rows = _mbr_feature_row_indices(main, pass1, mbr, path)
        mask, sidecar_rows = _mbr_candidate_sidecar_mask(
            main.global_qval, main.qval, mbr.MBR_best_is_missing_true,
            rows, hasproperty(mbr, :row_idx), q_value_threshold, path,
        )
        targets, decoys = _mbr_baseline_counts(
            main.qval, main.global_qval, main.target, q_value_threshold,
        )
        base_targets += targets
        base_decoys += decoys
        push!(masks, mask)
        push!(n_rows, n)
        isempty(sidecar_rows) && continue
        batches = feature_store === nothing ? (eachindex(sidecar_rows),) :
            Iterators.partition(eachindex(sidecar_rows), feature_batch_size)
        for batch in batches
            selected_sidecar_rows = @view sidecar_rows[batch]
            candidate_rows = rows[selected_sidecar_rows]

            frame = DataFrame(
                precursor_idx = UInt32.(main.precursor_idx[candidate_rows]),
                scan_idx = UInt32.(main.scan_idx[candidate_rows]),
                ms_file_idx = UInt32.(main.ms_file_idx[candidate_rows]),
                cv_fold = UInt8.(main.cv_fold[candidate_rows]),
                target = Bool.(main.target[candidate_rows]),
                qval = Float32.(main.qval[candidate_rows]),
                global_qval = Float32.(main.global_qval[candidate_rows]),
                trace_prob_prepass = Float32.(pass1.trace_prob_prepass[candidate_rows]),
                trace_prob_infold = Float32.(pass1.trace_prob_infold[candidate_rows]),
            )
            for feature in MBR_RECEIVER_FEATURES
                feature === :trace_prob_infold && continue
                hasproperty(main, feature) || continue
                frame[!, feature] = Tables.getcolumn(main, feature)[candidate_rows]
            end
            for feature in Symbol.(Tables.columnnames(mbr))
                feature in (:row_idx, :precursor_idx, :scan_idx) && continue
                frame[!, feature] = Tables.getcolumn(mbr, feature)[selected_sidecar_rows]
            end
            if feature_store === nothing
                push!(parts, frame)
            else
                metadata_columns = [:precursor_idx, :scan_idx, :ms_file_idx, :cv_fold,
                    :target, :qval, :global_qval,
                    (_mbr_missing_feature(k) for k in 1:MBR_N_COUNTERFACTUALS)...]
                _mbr_add_hellinger_contrasts!(frame)
                _mbr_store_features!(feature_store, frame)
                append!(metadata, select(frame, metadata_columns); cols=:union)
            end
        end
    end
    feature_store === nothing || flush(feature_store.data)
    # Recovery-sidecar writing also requires the ID columns when no candidates survive.
    candidates = if feature_store !== nothing && !isempty(metadata)
        metadata
    elseif isempty(parts)
        DataFrame(
            precursor_idx = UInt32[],
            scan_idx = UInt32[],
            ms_file_idx = UInt32[],
            cv_fold = UInt8[],
            target = Bool[],
            qval = Float32[],
            global_qval = Float32[],
            trace_prob_prepass = Float32[],
            trace_prob_infold = Float32[],
        )
    else
        vcat(parts...; cols = :union)
    end
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
        rows = _mbr_feature_row_indices(main, pass1, mbr, path)

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
            feature in (:row_idx, :precursor_idx, :scan_idx) && continue
            values = Tables.getcolumn(mbr, feature)
            expanded = eltype(values) <: Bool ? trues(n) : fill(-1.0f0, n)
            expanded[rows] = values
            frame[!, feature] = expanded
        end
        push!(parts, frame)
    end
    return isempty(parts) ? DataFrame() : vcat(parts...; cols = :union)
end

# Resolve DataFrame columns once so the row loop specializes on their concrete types.
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

function _mbr_available_feature_sets(frame::DataFrame; ion_mobility::Bool = true)
    true_features = Symbol[]
    false_features = [
        Symbol[] for _ in 1:MBR_N_COUNTERFACTUALS
    ]
    for feature in MBR_FTR_FEATURES_TRUE
        hasproperty(frame, feature) || continue
        ion_mobility || !(feature in MBR_ION_MOBILITY_FEATURES) || continue
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

# A concrete column tuple avoids per-cell DataFrame lookups in rank/margin comparisons.
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
        get_score_statistics!(eval_scores, eval_labels, eval_qvalues, eval_peps)
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

Return increasing, evenly spaced positions without allocating all eligible row indices.
A single ordered pass can match these positions to the selected candidates.
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

# Held-out candidate and available counterfactual rows, built once for all iterations.
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
    candidate_frame,
    true_features::Vector{Symbol},
    false_features::Vector{Vector{Symbol}},
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

        phase_started = time()
        @debug_l1 "MBR transfer model OOF fold $(Int(test_fold)) training gather starting: rows=$(length(train_rows)), features=$(length(true_features))"
        train_matrix = _mbr_gather_feature_rows(
            candidate_frame, true_features, false_features, train_rows, n_candidates)
        labels = _prepare_labels(train_labels)
        @debug_l1 "MBR transfer model OOF fold $(Int(test_fold)) training gather complete: elapsed=$(round(time() - phase_started; digits=2))s"

        classifier = build_lightgbm_classifier(; SHARED_LGBM_HP...)
        phase_started = time()
        @debug_l1 "MBR transfer model OOF fold $(Int(test_fold)) fit starting: rows=$(length(train_rows)), threads=$(classifier.num_threads), iterations=$(classifier.num_iterations)"
        LightGBM.fit!(classifier, train_matrix, labels; verbosity=-1)
        _detach_lightgbm_training_data!(classifier)
        train_matrix = nothing
        labels = nothing
        @debug_l1 "MBR transfer model OOF fold $(Int(test_fold)) fit complete: elapsed=$(round(time() - phase_started; digits=2))s"

        prediction_started = time()
        last_progress = prediction_started
        gather_seconds = 0.0
        predict_seconds = 0.0
        @debug_l1 "MBR transfer model OOF fold $(Int(test_fold)) prediction starting: rows=$(length(test_rows)), batch_size=$MBR_PREDICT_ROW_BATCH"
        for batch_start in 1:MBR_PREDICT_ROW_BATCH:length(test_rows)
            batch_stop = min(batch_start + MBR_PREDICT_ROW_BATCH - 1, length(test_rows))
            batch_rows = @view test_rows[batch_start:batch_stop]
            phase_started = time()
            test_matrix = _mbr_gather_feature_rows(
                candidate_frame, true_features, false_features, batch_rows, n_candidates)
            gather_seconds += time() - phase_started
            phase_started = time()
            raw = LightGBM.predict(classifier, test_matrix)
            predict_seconds += time() - phase_started
            test_matrix = nothing
            predictions = ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
            @inbounds scores[batch_rows] .= Float32.(predictions)
            now = time()
            if now - last_progress >= 60
                @debug_l1 "MBR transfer model OOF fold $(Int(test_fold)) prediction progress: rows=$batch_stop/$(length(test_rows)), gather=$(round(gather_seconds; digits=2))s, predict=$(round(predict_seconds; digits=2))s, elapsed=$(round(now - prediction_started; digits=2))s"
                last_progress = now
            end
        end
        @debug_l1 "MBR transfer model OOF fold $(Int(test_fold)) prediction complete: rows=$(length(test_rows)), gather=$(round(gather_seconds; digits=2))s, predict=$(round(predict_seconds; digits=2))s, elapsed=$(round(time() - prediction_started; digits=2))s"
        last_classifier = classifier
    end
    return scores, last_classifier
end

"""
    MBR_PREDICT_ROW_BATCH

Rows per `LightGBM.predict` call in the OOF pass. Bounds the gathered matrix independently of
experiment size; predictions are per-row independent so batching does not change them.
"""
const MBR_PREDICT_ROW_BATCH = 1_000_000

# Typed kernel for the gather below. `candidate_frame[!, col]` infers as AbstractVector, so reading
# it directly in the loop would dispatch per element.
function _mbr_gather_kernel!(
    dest::Matrix{Float32},
    column::AbstractVector,
    feature_idx::Int,
    candidate_rows::Vector{Int},
    dest_rows::Vector{Int},
)
    @inbounds for k in eachindex(candidate_rows, dest_rows)
        value = column[candidate_rows[k]]
        dest[dest_rows[k], feature_idx] =
            value === missing ? 0.0f0 : Float32(value)
    end
    return nothing
end

"""
    _mbr_gather_feature_rows(candidate_frame, true_features, false_features, rows, n_candidates)

Materialise just the requested rows of the block-stacked feature matrix.

Global row `r` identifies block `(r-1) ÷ n_candidates` and candidate
`(r-1) % n_candidates + 1`. Memory scales with the requested rows.

Rows are grouped by block first so each (block, feature) pair reads one concrete column, which is
what lets the kernel specialise.
"""
function _mbr_gather_feature_rows(
    candidate_frame::DataFrame,
    true_features::Vector{Symbol},
    false_features::Vector{Vector{Symbol}},
    rows::AbstractVector{Int},
    n_candidates::Int,
)
    n_features = length(true_features)
    n_blocks = 1 + length(false_features)
    dest = Matrix{Float32}(undef, length(rows), n_features)
    candidate_rows = [Int[] for _ in 1:n_blocks]
    dest_rows = [Int[] for _ in 1:n_blocks]
    @inbounds for (position, row) in enumerate(rows)
        block = div(row - 1, n_candidates)
        (0 <= block < n_blocks) ||
            throw(BoundsError("MBR gather row $row outside $n_blocks blocks of $n_candidates"))
        push!(candidate_rows[block + 1], mod(row - 1, n_candidates) + 1)
        push!(dest_rows[block + 1], position)
    end
    Threads.@threads for feature_idx in 1:n_features
        for block in 0:(n_blocks - 1)
            isempty(candidate_rows[block + 1]) && continue
            column = candidate_frame[
                !,
                block == 0 ? true_features[feature_idx] :
                    false_features[block][feature_idx],
            ]
            _mbr_gather_kernel!(
                dest, column, feature_idx,
                candidate_rows[block + 1], dest_rows[block + 1],
            )
        end
    end
    return dest
end

function _mbr_semisupervised_oof(
    candidate_frame::DataFrame,
    true_features::Vector{Symbol},
    false_features::Vector{Vector{Symbol}};
    feature_source = candidate_frame,
)
    preparation_started = time()
    n_candidates = nrow(candidate_frame)
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
    @debug_l1 "MBR transfer model OOF preparation: rows=$n_candidates, elapsed=$(round(time() - preparation_started; digits=2))s"
    best_state = nothing
    previous_positive = -1

    for iteration in 1:MBR_SEMISUPERVISED_MAX_ITERATIONS
        scores, classifier = _mbr_fit_oof_iteration(
            feature_source,
            true_features,
            false_features,
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
        if metrics.n_positive == 0
            @debug_l1 "MBR transfer model stopping: no confident target transfers remain; using iteration $(best_state.iteration)"
            break
        elseif iteration > 1 && !_scoring_target_gain_sufficient(
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

function _mbr_top_counterfactual_controls(
    scores::Vector{Float32},
    present::BitMatrix,
)
    n_candidates, n_counterfactuals = size(present)
    top_scores = fill(-Inf32, n_candidates)
    top_indices = zeros(UInt8, n_candidates)
    @inbounds for candidate_idx in 1:n_candidates
        for counterfactual_idx in 1:n_counterfactuals
            present[candidate_idx, counterfactual_idx] || continue
            score = scores[
                counterfactual_idx * n_candidates + candidate_idx
            ]
            if isfinite(score) && score > top_scores[candidate_idx]
                top_scores[candidate_idx] = score
                top_indices[candidate_idx] = UInt8(counterfactual_idx)
            end
        end
    end
    return (scores = top_scores, indices = top_indices)
end

function _mbr_top_counterfactual_scores(
    scores::Vector{Float32},
    present::BitMatrix,
)
    return _mbr_top_counterfactual_controls(scores, present).scores
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
    feature_source::Union{Nothing, _MBRFeatureStore} = nothing,
    ion_mobility::Bool = true,
)
    n = nrow(frame)
    frame[!, :mbr_recovered] = falses(n)
    frame[!, :MBR_transfer_candidate] = falses(n)
    frame[!, :mbr_target_decoy_prob] = fill(NaN32, n)
    frame[!, :ftr_qval_true] = fill(NaN32, n)
    frame[!, :ftr_pep_true] = fill(NaN32, n)
    frame[!, :mbr_total_error_qval_true] = fill(NaN32, n)
    frame[!, :mbr_total_error_rate_true] = fill(NaN32, n)
    frame[!, MBR_COUNTERFACTUAL_DECOY_PROB_COLUMN] = fill(NaN32, n)
    frame[!, MBR_COUNTERFACTUAL_DECOY_INDEX_COLUMN] = zeros(UInt8, n)
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

    candidates = frame_is_candidates ? frame : frame[candidate_indices, :]
    preparation_started = time()
    feature_source === nothing || frame_is_candidates ||
        throw(ArgumentError("Stored MBR features require a candidate-only frame"))
    feature_source === nothing && _mbr_add_hellinger_contrasts!(candidates)
    true_features, false_features = _mbr_available_feature_sets(
        feature_source === nothing ? candidates : feature_source.schema; ion_mobility)
    @debug_l1 "MBR transfer model feature preparation: rows=$(nrow(candidates)), features=$(length(true_features)), elapsed=$(round(time() - preparation_started; digits=2))s"
    best_state, present = _mbr_semisupervised_oof(
        candidates,
        true_features,
        false_features;
        feature_source = feature_source === nothing ? candidates : feature_source,
    )
    n_candidates = nrow(candidates)
    real_scores = copy(best_state.scores[1:n_candidates])
    evaluated_top = @view best_state.metrics.eval_mask[1:n_candidates]
    @inbounds for candidate_idx in 1:n_candidates
        evaluated_top[candidate_idx] ||
            (real_scores[candidate_idx] = NaN32)
    end
    counterfactual_controls = _mbr_top_counterfactual_controls(
        best_state.scores,
        present,
    )
    false_scores = counterfactual_controls.scores
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
        if Bool(candidates.target[candidate_position]) &&
           isfinite(combined.threshold) &&
           isfinite(real_scores[candidate_position]) &&
           isfinite(false_scores[candidate_position]) &&
           real_scores[candidate_position] >= combined.threshold &&
           false_scores[candidate_position] >= combined.threshold
            frame[row, MBR_COUNTERFACTUAL_DECOY_PROB_COLUMN] =
                false_scores[candidate_position]
            frame[row, MBR_COUNTERFACTUAL_DECOY_INDEX_COLUMN] =
                counterfactual_controls.indices[candidate_position]
        end
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
