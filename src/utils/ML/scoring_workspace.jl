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
