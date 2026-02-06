#=
ScoringWorkspace abstraction for CV fold setup and output array management.
Dispatches on container type to allow future OOM implementations.
=#

#############################################################################
# Abstract Type
#############################################################################

"""
    AbstractScoringWorkspace

Abstract type for scoring workspace management.
Bundles CV fold information, output arrays, and model storage.

Future OOM implementations can specialize this to map folds to file paths
and write output per-file rather than to global arrays.
"""
abstract type AbstractScoringWorkspace end

#############################################################################
# In-Memory Implementation
#############################################################################

"""
    InMemoryScoringWorkspace <: AbstractScoringWorkspace

In-memory workspace that pre-computes fold indices and allocates global output arrays.
"""
struct InMemoryScoringWorkspace <: AbstractScoringWorkspace
    psms::DataFramePSMContainer
    cv_folds::Vector{UInt8}
    fold_indices::Dict{UInt8, Vector{Int}}
    train_indices::Dict{UInt8, Vector{Int}}
    models::Dict{UInt8, Vector{Any}}
    prob_test::Vector{Float32}
    prob_train::Vector{Float32}
    nonMBR_estimates::Vector{Float32}
    MBR_estimates::Vector{Float32}
end

#############################################################################
# Factory Function (dispatch point for future OOM)
#############################################################################

"""
    setup_scoring_workspace(psms::DataFramePSMContainer, config::ScoringConfig)

Create a scoring workspace dispatched on container type.
For in-memory containers, pre-computes fold/train indices and allocates output arrays.
"""
function setup_scoring_workspace(psms::DataFramePSMContainer, config::ScoringConfig)
    cv_folds = get_cv_folds(psms)
    fold_idx = Dict(fold => get_fold_indices(psms, fold) for fold in cv_folds)
    train_idx = Dict(fold => get_train_indices(psms, fold) for fold in cv_folds)
    n = nrows(psms)
    return InMemoryScoringWorkspace(
        psms, cv_folds, fold_idx, train_idx,
        Dict{UInt8, Vector{Any}}(),
        zeros(Float32, n), zeros(Float32, n),
        zeros(Float32, n), zeros(Float32, n)
    )
end

#############################################################################
# Training Data Preparation (dispatch point for future OOM sampling)
#############################################################################

"""
    prepare_training_data!(psms::DataFramePSMContainer, config::ScoringConfig)

Prepare training data from the PSM container (mutates `psms` in place).
In-memory: assigns pairs, sorts for memory locality, and initializes MBR columns.
Future OOM: would handle pairing/sorting/initialization per-file or during sampling.
"""
function prepare_training_data!(psms::DataFramePSMContainer, config::ScoringConfig)
    # Assign target-decoy pairs
    assign_pairs!(psms, config.pairing)
    # Sort for memory locality in MBR feature computation (summarize_precursors!)
    sort_container!(psms, [:pair_id, :isotopes_captured, :precursor_idx, :ms_file_idx])
    # Initialize MBR columns
    initialize_mbr_columns!(psms, config.mbr_update)
    return nothing
end

#############################################################################
# Interface Functions (all overridable for future OOM)
#############################################################################

# CV fold iteration
get_cv_folds(ws::InMemoryScoringWorkspace) = ws.cv_folds

# Training phase access
get_training_view(ws::InMemoryScoringWorkspace, fold::UInt8) = get_view(ws.psms, ws.train_indices[fold])
get_train_indices(ws::InMemoryScoringWorkspace, fold::UInt8) = ws.train_indices[fold]

# Prediction phase access
get_test_view(ws::InMemoryScoringWorkspace, fold::UInt8) = get_view(ws.psms, ws.fold_indices[fold])
get_test_indices(ws::InMemoryScoringWorkspace, fold::UInt8) = ws.fold_indices[fold]

# Output arrays
get_train_output(ws::InMemoryScoringWorkspace) = ws.prob_train
get_test_output(ws::InMemoryScoringWorkspace) = ws.prob_test
get_baseline_output(ws::InMemoryScoringWorkspace) = ws.nonMBR_estimates
get_mbr_output(ws::InMemoryScoringWorkspace) = ws.MBR_estimates

# Model storage
store_fold_models!(ws::InMemoryScoringWorkspace, fold::UInt8, models::Vector{Any}) = (ws.models[fold] = models)
get_fold_models(ws::InMemoryScoringWorkspace, fold::UInt8) = ws.models[fold]
get_all_models(ws::InMemoryScoringWorkspace) = ws.models

# Container access (for finalize_scoring! and other top-level ops)
get_psms(ws::InMemoryScoringWorkspace) = ws.psms

#############################################################################
# Phase-Dispatched Accessors (unify training/prediction fold access)
#############################################################################

get_phase_view(ws::InMemoryScoringWorkspace, ::TrainingPhase, fold::UInt8) = get_training_view(ws, fold)
get_phase_view(ws::InMemoryScoringWorkspace, ::PredictionPhase, fold::UInt8) = get_test_view(ws, fold)

get_phase_indices(ws::InMemoryScoringWorkspace, ::TrainingPhase, fold::UInt8) = ws.train_indices[fold]
get_phase_indices(ws::InMemoryScoringWorkspace, ::PredictionPhase, fold::UInt8) = ws.fold_indices[fold]

get_phase_output(ws::InMemoryScoringWorkspace, ::TrainingPhase) = ws.prob_train
get_phase_output(ws::InMemoryScoringWorkspace, ::PredictionPhase) = ws.prob_test

init_fold_models(ws::InMemoryScoringWorkspace, ::TrainingPhase, ::UInt8, n::Int) = Vector{Any}(undef, n)
init_fold_models(ws::InMemoryScoringWorkspace, ::PredictionPhase, fold::UInt8, ::Int) = ws.models[fold]

commit_fold!(ws::InMemoryScoringWorkspace, ::TrainingPhase, fold::UInt8, models::Vector{Any}) = (ws.models[fold] = models)
commit_fold!(ws::InMemoryScoringWorkspace, ::PredictionPhase, ::UInt8, ::Vector{Any}) = nothing

function store_final_predictions!(ws::InMemoryScoringWorkspace, fold::UInt8, ::PairBasedMBR)
    idx = ws.fold_indices[fold]
    ws.MBR_estimates[idx] = collect(Float32, get_column(get_test_view(ws, fold), :trace_prob))
end

function store_final_predictions!(ws::InMemoryScoringWorkspace, fold::UInt8, ::NoMBR)
    idx = ws.fold_indices[fold]
    ws.prob_test[idx] = collect(Float32, get_column(get_test_view(ws, fold), :trace_prob))
end

#############################################################################
# OOM Training Data Preparation (ArrowFilePSMContainer)
#############################################################################

"""
    prepare_training_data!(container::ArrowFilePSMContainer, config::ScoringConfig)

Prepare training data for the file-backed OOM container. Two phases:

**Phase 1 — Lightweight global pairing**: Read only the 4 columns needed for
pairing (~19 bytes/PSM) from all files, assign globally-unique pair_ids.

**Phase 2 — Per-file sort + sidecar creation**: Load one file at a time,
merge in global pair_ids, sort for MBR locality, write sorted data back,
and create a scores sidecar file with initialized mutable columns.
"""
function prepare_training_data!(container::ArrowFilePSMContainer, config::ScoringConfig)
    file_groups = container.file_groups
    isempty(file_groups) && return nothing

    # ─── Phase 1: Lightweight global pairing ────────────────────────────
    # Read only the 4 columns needed for pairing from each file.
    # This fits in memory even for very large datasets (10M PSMs ≈ 190 MB).
    pairing_dfs = Vector{DataFrame}(undef, length(file_groups))
    file_boundaries = Vector{UnitRange{Int}}(undef, length(file_groups))
    row_offset = 0

    for (i, group) in enumerate(file_groups)
        tbl = Arrow.Table(group.data_path)
        pairing_dfs[i] = DataFrame(
            irt_pred        = collect(tbl.irt_pred),
            cv_fold         = collect(tbl.cv_fold),
            isotopes_captured = collect(tbl.isotopes_captured),
            precursor_idx   = collect(tbl.precursor_idx),
        )
        n = group.n_rows
        file_boundaries[i] = (row_offset + 1):(row_offset + n)
        row_offset += n
    end

    global_df = vcat(pairing_dfs...)
    pairing_dfs = nothing  # allow GC

    # Reuse existing assign_pairs! via a lightweight container
    pairing_container = DataFramePSMContainer(global_df, Val(:unsafe))
    assign_pairs!(pairing_container, config.pairing)

    global_pair_ids = collect(UInt32, global_df[!, :pair_id])
    global_irt_bins = collect(UInt32, global_df[!, :irt_bin_idx])
    global_df = nothing  # allow GC

    # ─── Phase 2: Per-file sort + sidecar creation ──────────────────────
    for (i, group) in enumerate(file_groups)
        rows = file_boundaries[i]

        # Load full data for this one file
        df = DataFrame(Tables.columntable(Arrow.Table(group.data_path)))

        # Merge in globally-computed columns
        df[!, :pair_id] = global_pair_ids[rows]
        df[!, :irt_bin_idx] = global_irt_bins[rows]

        # Sort for MBR memory locality (same order as in-memory path)
        sort!(df, [:pair_id, :isotopes_captured, :precursor_idx, :ms_file_idx])

        # Write sorted main data back (immutable after this point)
        writeArrow(group.data_path, df)

        # Create scores sidecar with initialized mutable columns
        n = nrow(df)
        scores_df = _create_scores_sidecar(n, config.mbr_update)
        writeArrow(group.scores_path, scores_df)

        df = nothing  # allow GC per file
    end

    return nothing
end

"""
    _create_scores_sidecar(n::Int, ::PairBasedMBR) -> DataFrame

Create the mutable scores sidecar DataFrame for PairBasedMBR mode.
Mirrors the columns initialized by `initialize_mbr_columns!(_, ::PairBasedMBR)`.
"""
function _create_scores_sidecar(n::Int, ::PairBasedMBR)
    return DataFrame(
        trace_prob              = zeros(Float32, n),
        q_value                 = zeros(Float64, n),
        nonMBR_trace_prob       = zeros(Float32, n),
        MBR_max_pair_prob       = zeros(Float32, n),
        MBR_best_irt_diff       = zeros(Float32, n),
        MBR_log2_weight_ratio   = zeros(Float32, n),
        MBR_log2_explained_ratio = zeros(Float32, n),
        MBR_rv_coefficient      = zeros(Float32, n),
        MBR_is_best_decoy       = trues(n),
        MBR_num_runs            = zeros(Int32, n),
        MBR_transfer_candidate  = falses(n),
        MBR_is_missing          = falses(n),
    )
end

"""
    _create_scores_sidecar(n::Int, ::NoMBR) -> DataFrame

Create the mutable scores sidecar DataFrame for NoMBR mode.
Mirrors the columns initialized by `initialize_mbr_columns!(_, ::NoMBR)`.
"""
function _create_scores_sidecar(n::Int, ::NoMBR)
    return DataFrame(
        trace_prob = zeros(Float32, n),
        q_value    = zeros(Float64, n),
    )
end

#############################################################################
# ArrowFileScoringWorkspace — OOM wrapper for pair-sampled training
#############################################################################

"""
    ArrowFileScoringWorkspace <: AbstractScoringWorkspace

OOM-aware workspace that samples a representative subset of PSMs (by pair_id)
from an `ArrowFilePSMContainer`, loads the sample into a normal
`InMemoryScoringWorkspace`, and delegates all training/prediction to it.

The `container` field retains a reference to all files on disk for later
per-file prediction (Step 3).
"""
struct ArrowFileScoringWorkspace <: AbstractScoringWorkspace
    container::ArrowFilePSMContainer
    inner::InMemoryScoringWorkspace
end

#############################################################################
# setup_scoring_workspace — ArrowFilePSMContainer dispatch
#############################################################################

"""
    setup_scoring_workspace(psms::ArrowFilePSMContainer, config::ScoringConfig)

Build a scoring workspace for the OOM path. Samples whole pairs within a PSM
budget, loads the sample into an `InMemoryScoringWorkspace`, and wraps it.
"""
function setup_scoring_workspace(psms::ArrowFilePSMContainer, config::ScoringConfig)
    max_psms = min(psms.max_training_psms, nrows(psms))

    # Sample by pair_id (3-pass: count → select → fill)
    sampled = _sample_by_pair_id(psms, max_psms)

    # Create inner InMemoryScoringWorkspace from sampled data
    inner = setup_scoring_workspace(sampled, config)

    return ArrowFileScoringWorkspace(psms, inner)
end

#############################################################################
# _sample_by_pair_id — 3-pass pair-aware sampling
#############################################################################

"""
    _sample_by_pair_id(container::ArrowFilePSMContainer, max_psms::Int) -> DataFramePSMContainer

Sample whole target-decoy pairs from disk files within a PSM budget.

Three passes:
1. Count PSMs per pair_id (reads only the pair_id column).
2. Shuffle pair_ids deterministically, greedily select within budget.
3. Pre-allocate a DataFrame and copy only matching rows from each file.
"""
function _sample_by_pair_id(container::ArrowFilePSMContainer, max_psms::Int)
    # ─── Pass 1: Count PSMs per pair_id ──────────────────────────────────
    pair_counts = Dict{UInt32, Int}()
    for group in container.file_groups
        tbl = Arrow.Table(group.data_path)
        for pid in tbl.pair_id
            pair_counts[pid] = get(pair_counts, pid, 0) + 1
        end
    end

    # ─── Pass 2: Select pair_ids within budget ───────────────────────────
    unique_pids = collect(keys(pair_counts))
    shuffle!(MersenneTwister(1776), unique_pids)

    selected_vec = Vector{UInt32}(undef, length(unique_pids))
    n_selected = 0
    running_count = 0
    for pid in unique_pids
        count = pair_counts[pid]
        if running_count + count > max_psms
            continue  # skip this pair, try smaller ones
        end
        n_selected += 1
        selected_vec[n_selected] = pid
        running_count += count
    end
    resize!(selected_vec, n_selected)

    # BitVector for O(1) lookup — pair_ids are dense (sequential from 1)
    max_pid = maximum(keys(pair_counts))
    is_selected = falses(max_pid)
    for pid in selected_vec
        is_selected[pid] = true
    end

    # ─── Pass 3: Pre-allocate and fill ───────────────────────────────────
    # Read schema from first file group (data + sidecar)
    first_data_tbl = Arrow.Table(container.file_groups[1].data_path)
    first_scores_tbl = Arrow.Table(container.file_groups[1].scores_path)

    # Allocate combined DataFrame with exact row count
    sampled_df = _create_combined_dataframe(first_data_tbl, first_scores_tbl, running_count)

    # Fill from each file
    write_offset = 0
    for group in container.file_groups
        data_tbl = Arrow.Table(group.data_path)
        scores_tbl = Arrow.Table(group.scores_path)

        # Build mask from memory-mapped pair_id column
        pair_id_col = data_tbl.pair_id
        mask = BitVector(undef, length(pair_id_col))
        @inbounds for i in eachindex(pair_id_col)
            mask[i] = is_selected[pair_id_col[i]]
        end
        n_match = sum(mask)
        n_match == 0 && continue

        dest_range = (write_offset + 1):(write_offset + n_match)

        # Copy matching rows from data columns
        for col_name in Tables.columnnames(data_tbl)
            src_col = Tables.getcolumn(data_tbl, col_name)
            sampled_df[!, col_name][dest_range] = src_col[mask]
        end
        # Copy matching rows from sidecar columns
        for col_name in Tables.columnnames(scores_tbl)
            src_col = Tables.getcolumn(scores_tbl, col_name)
            sampled_df[!, col_name][dest_range] = src_col[mask]
        end

        write_offset += n_match
    end

    return DataFramePSMContainer(sampled_df, Val(:unsafe))
end

"""
    _create_combined_dataframe(data_tbl, scores_tbl, n_rows::Int) -> DataFrame

Pre-allocate a DataFrame with columns from both data and sidecar Arrow tables.
"""
function _create_combined_dataframe(data_tbl::Arrow.Table, scores_tbl::Arrow.Table, n_rows::Int)
    df = DataFrame()
    for tbl in (data_tbl, scores_tbl)
        schema = Tables.schema(tbl)
        for (col_name, col_type) in zip(Tables.columnnames(tbl), schema.types)
            if col_type isa Union && Missing <: col_type
                non_missing_type = Base.uniontypes(col_type)
                actual_type = first(t for t in non_missing_type if t !== Missing)
                df[!, col_name] = Vector{Union{Missing, actual_type}}(undef, n_rows)
            else
                df[!, col_name] = Vector{col_type}(undef, n_rows)
            end
        end
    end
    return df
end

#############################################################################
# Delegation — ArrowFileScoringWorkspace → inner InMemoryScoringWorkspace
#############################################################################

# CV fold iteration
get_cv_folds(ws::ArrowFileScoringWorkspace) = get_cv_folds(ws.inner)

# Training/prediction phase access
get_training_view(ws::ArrowFileScoringWorkspace, fold::UInt8) = get_training_view(ws.inner, fold)
get_train_indices(ws::ArrowFileScoringWorkspace, fold::UInt8) = get_train_indices(ws.inner, fold)
get_test_view(ws::ArrowFileScoringWorkspace, fold::UInt8) = get_test_view(ws.inner, fold)
get_test_indices(ws::ArrowFileScoringWorkspace, fold::UInt8) = get_test_indices(ws.inner, fold)

# Output arrays
get_train_output(ws::ArrowFileScoringWorkspace) = get_train_output(ws.inner)
get_test_output(ws::ArrowFileScoringWorkspace) = get_test_output(ws.inner)
get_baseline_output(ws::ArrowFileScoringWorkspace) = get_baseline_output(ws.inner)
get_mbr_output(ws::ArrowFileScoringWorkspace) = get_mbr_output(ws.inner)

# Model storage
store_fold_models!(ws::ArrowFileScoringWorkspace, fold::UInt8, models::Vector{Any}) = store_fold_models!(ws.inner, fold, models)
get_fold_models(ws::ArrowFileScoringWorkspace, fold::UInt8) = get_fold_models(ws.inner, fold)
get_all_models(ws::ArrowFileScoringWorkspace) = get_all_models(ws.inner)

# Container access — returns sampled DataFramePSMContainer
get_psms(ws::ArrowFileScoringWorkspace) = get_psms(ws.inner)

# Phase-dispatched accessors (TrainingPhase delegates to inner sampled workspace;
# PredictionPhase is handled by the specialized process_fold_iterations! method)
get_phase_view(ws::ArrowFileScoringWorkspace, phase::TrainingPhase, fold::UInt8) = get_phase_view(ws.inner, phase, fold)
get_phase_indices(ws::ArrowFileScoringWorkspace, phase::TrainingPhase, fold::UInt8) = get_phase_indices(ws.inner, phase, fold)
get_phase_output(ws::ArrowFileScoringWorkspace, phase::TrainingPhase) = get_phase_output(ws.inner, phase)
init_fold_models(ws::ArrowFileScoringWorkspace, phase::TrainingPhase, fold::UInt8, n::Int) = init_fold_models(ws.inner, phase, fold, n)
commit_fold!(ws::ArrowFileScoringWorkspace, phase::TrainingPhase, fold::UInt8, models::Vector{Any}) = commit_fold!(ws.inner, phase, fold, models)

# store_final_predictions! — no-op for ArrowFileScoringWorkspace
# Predictions are already written to sidecar files during process_fold_iterations!
store_final_predictions!(::ArrowFileScoringWorkspace, ::UInt8, ::PairBasedMBR) = nothing
store_final_predictions!(::ArrowFileScoringWorkspace, ::UInt8, ::NoMBR) = nothing

#############################################################################
# Workspace-level finalize_scoring! — ArrowFileScoringWorkspace
#############################################################################

finalize_scoring!(ws::ArrowFileScoringWorkspace, strategy::PairBasedMBR) =
    _finalize_scoring_arrow!(ws.container, strategy)
finalize_scoring!(ws::ArrowFileScoringWorkspace, strategy::NoMBR) =
    _finalize_scoring_arrow!(ws.container, strategy)

#############################################################################
# OOM Prediction Phase (ArrowFileScoringWorkspace)
#############################################################################

"""
    process_fold_iterations!(::PredictionPhase, workspace::ArrowFileScoringWorkspace, ...)

OOM prediction for a single fold: iterate over all Arrow files, filter to this
fold's PSMs, apply pre-trained models, write predictions to sidecars.

MBR features are computed per-fold across all files between iterations, since
pairs span multiple MS runs.
"""
function process_fold_iterations!(
    ::PredictionPhase,
    workspace::ArrowFileScoringWorkspace,
    fold::UInt8,
    iteration_rounds::Vector{Int},
    config::ScoringConfig,
    mbr_start_iter::Int,
    pbar
)
    models = get_all_models(workspace)
    container = workspace.container

    for itr in eachindex(iteration_rounds)
        # Early exit for non-MBR mode
        if !has_mbr_support(config.mbr_update) && itr > (mbr_start_iter - 1)
            break
        end

        # Per-file prediction for this fold
        for group in container.file_groups
            data_df = DataFrame(Tables.columntable(Arrow.Table(group.data_path)))
            scores_df = DataFrame(Tables.columntable(Arrow.Table(group.scores_path)))

            fold_mask = data_df.cv_fold .== fold
            any(fold_mask) || continue

            # Build temp container with data + scores columns for prediction
            fold_df = data_df[fold_mask, :]
            for col in names(scores_df)
                fold_df[!, Symbol(col)] = scores_df[fold_mask, col]
            end
            temp = DataFramePSMContainer(fold_df, Val(:unsafe))

            preds = predict_scores(models[fold][itr], temp)
            scores_df.trace_prob[fold_mask] .= preds

            # Store baseline at pre-MBR iteration
            if itr == mbr_start_iter - 1
                scores_df.nonMBR_trace_prob[fold_mask] .= preds
            end

            writeArrow(group.scores_path, scores_df)
        end

        # Cross-file MBR for this fold
        if itr >= mbr_start_iter - 1 && has_mbr_support(config.mbr_update)
            _compute_mbr_for_fold!(container, fold, config.mbr_update)
        end

        !isnothing(pbar) && update(pbar)
    end
end

#############################################################################
# Cross-file MBR computation (per-fold)
#############################################################################

"""
    _compute_mbr_for_fold!(container, fold, strategy::PairBasedMBR)

Load MBR-relevant columns for a single fold from ALL files into one DataFrame,
run `summarize_precursors!`, distribute results back per-file.

`summarize_precursors!` groups by `(pair_id, isotopes_captured)` and needs
all PSMs for a pair together. Pairs span files because MBR = match between runs.
Filtering by fold ensures we only update MBR features for the held-out test set
of this fold.
"""
function _compute_mbr_for_fold!(container::ArrowFilePSMContainer, fold::UInt8, strategy::PairBasedMBR)
    mbr_data_cols = [:pair_id, :isotopes_captured, :ms_file_idx, :weight,
                      :log2_intensity_explained, :decoy, :irt_pred, :irt_obs, :weights, :irts]
    mbr_score_cols = [:trace_prob, :q_value, :MBR_max_pair_prob, :MBR_best_irt_diff,
                      :MBR_log2_weight_ratio, :MBR_log2_explained_ratio,
                      :MBR_rv_coefficient, :MBR_is_best_decoy, :MBR_num_runs,
                      :MBR_transfer_candidate, :MBR_is_missing]

    # 1. Collect MBR-relevant columns for this fold from all files
    dfs = DataFrame[]
    # Track per-file fold masks and boundaries for write-back
    file_masks = BitVector[]
    file_offsets = Int[]  # offset into the concatenated fold-filtered DataFrame
    offset = 0
    for group in container.file_groups
        data_tbl = Arrow.Table(group.data_path)
        scores_tbl = Arrow.Table(group.scores_path)

        # Filter to this fold
        cv_folds = collect(data_tbl.cv_fold)
        fold_mask = cv_folds .== fold
        push!(file_masks, fold_mask)
        n_fold = sum(fold_mask)
        n_fold == 0 && (push!(file_offsets, offset); continue)

        df = DataFrame([col => collect(getproperty(data_tbl, col))[fold_mask] for col in mbr_data_cols])

        for col in mbr_score_cols
            df[!, col] = collect(getproperty(scores_tbl, col))[fold_mask]
        end
        push!(dfs, df)
        push!(file_offsets, offset)
        offset += n_fold
    end

    isempty(dfs) && return nothing
    global_df = vcat(dfs...); dfs = nothing

    # 2. Compute q-values (needed by summarize_precursors! for MBR threshold)
    targets = .!global_df.decoy
    get_qvalues!(global_df.trace_prob, targets, global_df.q_value)

    # 3. Run existing summarize_precursors!
    summarize_precursors!(global_df, q_cutoff=strategy.q_cutoff)

    # 4. Write MBR columns back to per-file sidecars (fold rows only)
    for (i, group) in enumerate(container.file_groups)
        n_fold = sum(file_masks[i])
        n_fold == 0 && continue
        rows = (file_offsets[i]+1):(file_offsets[i]+n_fold)

        scores_df = DataFrame(Tables.columntable(Arrow.Table(group.scores_path)))
        for col in mbr_score_cols
            scores_df[file_masks[i], col] = global_df[rows, col]
        end
        writeArrow(group.scores_path, scores_df)
    end
    global_df = nothing
end

#############################################################################
# Finalization — merge final scores back into data files
#############################################################################

"""
    _finalize_scoring_arrow!(container, strategy::PairBasedMBR)

Two-pass finalization for PairBasedMBR:
1. Collect nonMBR baseline + targets from all files, compute global q-values
   to find the probability threshold for transfer candidates
2. Merge final columns (trace_prob, MBR_boosted_trace_prob, MBR_transfer_candidate)
   back into the data Arrow files
"""
function _finalize_scoring_arrow!(container::ArrowFilePSMContainer, strategy::PairBasedMBR)
    # Pass 1: Collect nonMBR baseline + targets from all files
    all_probs = Float32[]
    all_targets = Bool[]
    for group in container.file_groups
        scores_tbl = Arrow.Table(group.scores_path)
        append!(all_probs, collect(scores_tbl.nonMBR_trace_prob))
        data_tbl = Arrow.Table(group.data_path)
        append!(all_targets, collect(data_tbl.target))
    end
    qvals = similar(all_probs)
    get_qvalues!(all_probs, all_targets, qvals)
    pass_mask = qvals .<= strategy.q_cutoff
    prob_thresh = any(pass_mask) ? minimum(all_probs[pass_mask]) : typemax(Float32)

    # Pass 2: Merge final columns into data files
    offset = 0
    for group in container.file_groups
        n = group.n_rows; rows = (offset+1):(offset+n)
        data_df = DataFrame(Tables.columntable(Arrow.Table(group.data_path)))
        scores_df = DataFrame(Tables.columntable(Arrow.Table(group.scores_path)))

        # Transfer candidates: didn't pass on baseline but best pair prob exceeds threshold
        file_pass = pass_mask[rows]
        transfer = .!file_pass .& (scores_df.MBR_max_pair_prob .>= prob_thresh)
        data_df[!, :MBR_transfer_candidate] = transfer

        # Final probability columns
        data_df[!, :trace_prob] = scores_df.nonMBR_trace_prob
        data_df[!, :MBR_boosted_trace_prob] = scores_df.trace_prob

        # Copy MBR feature columns from sidecar to data file
        # These are needed by select_mbr_features() and apply_mbr_filter!()
        for col in (:MBR_max_pair_prob, :MBR_best_irt_diff, :MBR_rv_coefficient,
                    :MBR_log2_weight_ratio, :MBR_log2_explained_ratio,
                    :MBR_is_best_decoy, :MBR_num_runs, :MBR_is_missing)
            data_df[!, col] = scores_df[!, col]
        end

        writeArrow(group.data_path, data_df)
        offset += n
    end
end

function _finalize_scoring_arrow!(container::ArrowFilePSMContainer, ::NoMBR)
    for group in container.file_groups
        data_df = DataFrame(Tables.columntable(Arrow.Table(group.data_path)))
        scores_df = DataFrame(Tables.columntable(Arrow.Table(group.scores_path)))
        data_df[!, :trace_prob] = scores_df.trace_prob
        writeArrow(group.data_path, data_df)
    end
end
