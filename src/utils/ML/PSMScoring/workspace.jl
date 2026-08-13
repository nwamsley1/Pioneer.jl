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
end

#############################################################################
# Factory Function
#############################################################################

"""
    setup_scoring_workspace(psms::DataFramePSMContainer, config::ScoringConfig)

Create a scoring workspace. Pre-computes fold/train indices and allocates output arrays.
"""
function setup_scoring_workspace(psms::DataFramePSMContainer, config::ScoringConfig)
    cv_folds = get_cv_folds(psms)
    fold_idx = Dict(fold => get_fold_indices(psms, fold) for fold in cv_folds)
    train_idx = Dict(fold => get_train_indices(psms, fold) for fold in cv_folds)
    n = nrows(psms)
    return InMemoryScoringWorkspace(
        psms, cv_folds, fold_idx, train_idx,
        Dict{UInt8, Vector{Any}}(),
        zeros(Float32, n), zeros(Float32, n)
    )
end

#############################################################################
# Training Data Preparation
#############################################################################

"""
    prepare_training_data!(psms::DataFramePSMContainer, config::ScoringConfig)

Prepare training data from the PSM container (mutates `psms` in place).
Sorts for memory locality and initializes output columns.
"""
function prepare_training_data!(psms::DataFramePSMContainer, config::ScoringConfig)
    # Sort for memory locality
    sort_container!(psms, [:isotopes_captured, :precursor_idx, :ms_file_idx])
    # Initialize output columns
    n = nrows(psms)
    set_column!(psms, :trace_prob, zeros(Float32, n))
    set_column!(psms, :q_value, zeros(Float64, n))
    return nothing
end

#############################################################################
# Interface Functions
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

# Model storage
store_fold_models!(ws::InMemoryScoringWorkspace, fold::UInt8, models::Vector{Any}) = (ws.models[fold] = models)
get_fold_models(ws::InMemoryScoringWorkspace, fold::UInt8) = ws.models[fold]
get_all_models(ws::InMemoryScoringWorkspace) = ws.models

# Container access
get_psms(ws::InMemoryScoringWorkspace) = ws.psms

#############################################################################
# Phase-Dispatched Accessors
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

function store_final_predictions!(ws::InMemoryScoringWorkspace, fold::UInt8)
    idx = ws.fold_indices[fold]
    ws.prob_test[idx] = collect(Float32, get_column(get_test_view(ws, fold), :trace_prob))
end

#############################################################################
# OOM Training Data Preparation (ArrowFilePSMContainer)
#############################################################################

"""
    prepare_training_data!(container::ArrowFilePSMContainer, config::ScoringConfig)

Prepare training data for the file-backed OOM container.
Per-file: sort for locality and create scores sidecar file.
"""
function prepare_training_data!(container::ArrowFilePSMContainer, config::ScoringConfig)
    file_groups = container.file_groups
    isempty(file_groups) && return nothing

    for group in file_groups
        df = DataFrame(Tables.columntable(Arrow.Table(group.data_path)))

        # Sort for locality
        fast_df_sort!(df, [:isotopes_captured, :precursor_idx, :ms_file_idx])

        # Write sorted main data back (immutable after this point)
        writeArrow(group.data_path, df)

        # Create scores sidecar with initialized mutable columns
        n = nrow(df)
        scores_df = DataFrame(
            trace_prob = zeros(Float32, n),
            q_value    = zeros(Float64, n),
        )
        writeArrow(group.scores_path, scores_df)

        df = nothing  # allow GC per file
    end

    return nothing
end

#############################################################################
# ArrowFileScoringWorkspace — OOM wrapper for sampled training
#############################################################################

"""
    ArrowFileScoringWorkspace <: AbstractScoringWorkspace

OOM-aware workspace that samples a representative subset of PSMs
from an `ArrowFilePSMContainer`, loads the sample into a normal
`InMemoryScoringWorkspace`, and delegates all training/prediction to it.
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

Build a scoring workspace for the OOM path. Randomly samples PSMs within a
budget, loads the sample into an `InMemoryScoringWorkspace`, and wraps it.
"""
function setup_scoring_workspace(psms::ArrowFilePSMContainer, config::ScoringConfig)
    max_psms = min(psms.max_training_psms, nrows(psms))

    sampled = _sample_psms(psms, max_psms)

    # Create inner InMemoryScoringWorkspace from sampled data
    inner = setup_scoring_workspace(sampled, config)

    return ArrowFileScoringWorkspace(psms, inner)
end

#############################################################################
# _sample_psms — random PSM sampling from disk files
#############################################################################

"""
    _sample_psms(container::ArrowFilePSMContainer, max_psms::Int) -> DataFramePSMContainer

Randomly sample PSMs from disk files within a PSM budget.
Uses reservoir-style selection: compute per-file sample counts proportional
to file size, then randomly select rows from each file.
"""
function _sample_psms(container::ArrowFilePSMContainer, max_psms::Int)
    total = nrows(container)
    sample_frac = min(1.0, max_psms / total)
    rng = MersenneTwister(1776)

    # ─── Pass 1: Determine per-file sample sizes ────────────────────────
    file_sample_counts = Int[]
    running_total = 0
    for group in container.file_groups
        n = min(round(Int, group.n_rows * sample_frac), group.n_rows)
        push!(file_sample_counts, n)
        running_total += n
    end

    # ─── Pass 2: Pre-allocate and fill ───────────────────────────────────
    first_data_tbl = Arrow.Table(container.file_groups[1].data_path)
    first_scores_tbl = Arrow.Table(container.file_groups[1].scores_path)

    sampled_df = _create_combined_dataframe(first_data_tbl, first_scores_tbl, running_total)

    write_offset = 0
    for (gi, group) in enumerate(container.file_groups)
        n_sample = file_sample_counts[gi]
        n_sample == 0 && continue

        data_tbl = Arrow.Table(group.data_path)
        scores_tbl = Arrow.Table(group.scores_path)

        # Random row indices for this file
        all_indices = collect(1:group.n_rows)
        shuffle!(rng, all_indices)
        selected = sort!(all_indices[1:n_sample])

        dest_range = (write_offset + 1):(write_offset + n_sample)

        for col_name in Tables.columnnames(data_tbl)
            src_col = Tables.getcolumn(data_tbl, col_name)
            sampled_df[!, col_name][dest_range] = src_col[selected]
        end
        for col_name in Tables.columnnames(scores_tbl)
            src_col = Tables.getcolumn(scores_tbl, col_name)
            sampled_df[!, col_name][dest_range] = src_col[selected]
        end

        write_offset += n_sample
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

get_cv_folds(ws::ArrowFileScoringWorkspace) = get_cv_folds(ws.inner)
get_training_view(ws::ArrowFileScoringWorkspace, fold::UInt8) = get_training_view(ws.inner, fold)
get_train_indices(ws::ArrowFileScoringWorkspace, fold::UInt8) = get_train_indices(ws.inner, fold)
get_test_view(ws::ArrowFileScoringWorkspace, fold::UInt8) = get_test_view(ws.inner, fold)
get_test_indices(ws::ArrowFileScoringWorkspace, fold::UInt8) = get_test_indices(ws.inner, fold)
get_train_output(ws::ArrowFileScoringWorkspace) = get_train_output(ws.inner)
get_test_output(ws::ArrowFileScoringWorkspace) = get_test_output(ws.inner)
store_fold_models!(ws::ArrowFileScoringWorkspace, fold::UInt8, models::Vector{Any}) = store_fold_models!(ws.inner, fold, models)
get_fold_models(ws::ArrowFileScoringWorkspace, fold::UInt8) = get_fold_models(ws.inner, fold)
get_all_models(ws::ArrowFileScoringWorkspace) = get_all_models(ws.inner)
get_psms(ws::ArrowFileScoringWorkspace) = get_psms(ws.inner)

get_phase_view(ws::ArrowFileScoringWorkspace, phase::TrainingPhase, fold::UInt8) = get_phase_view(ws.inner, phase, fold)
get_phase_indices(ws::ArrowFileScoringWorkspace, phase::TrainingPhase, fold::UInt8) = get_phase_indices(ws.inner, phase, fold)
get_phase_output(ws::ArrowFileScoringWorkspace, phase::TrainingPhase) = get_phase_output(ws.inner, phase)
init_fold_models(ws::ArrowFileScoringWorkspace, phase::TrainingPhase, fold::UInt8, n::Int) = init_fold_models(ws.inner, phase, fold, n)
commit_fold!(ws::ArrowFileScoringWorkspace, phase::TrainingPhase, fold::UInt8, models::Vector{Any}) = commit_fold!(ws.inner, phase, fold, models)

store_final_predictions!(::ArrowFileScoringWorkspace, ::UInt8) = nothing

#############################################################################
# OOM Prediction Phase (ArrowFileScoringWorkspace)
#############################################################################

"""
    process_fold!(::PredictionPhase, workspace::ArrowFileScoringWorkspace, ...)

OOM prediction for a single fold: apply pre-trained model to all Arrow files
in this fold, write predictions to sidecars.
"""
function process_fold!(
    ::PredictionPhase,
    workspace::ArrowFileScoringWorkspace,
    fold::UInt8,
    num_rounds::Int,
    config::ScoringConfig,
    pbar
)
    models = get_all_models(workspace)
    container = workspace.container

    # Per-file prediction for this fold (files are already per-fold)
    for group in get_file_groups_for_fold(container, fold)
        data_df = DataFrame(Arrow.Table(group.data_path))
        scores_df = DataFrame(Tables.columntable(Arrow.Table(group.scores_path)))

        # Build temp container with data + scores columns for prediction
        for col in names(scores_df)
            data_df[!, Symbol(col)] = scores_df[!, col]
        end
        temp = DataFramePSMContainer(data_df, Val(:unsafe))

        preds = predict_scores(models[fold][1], temp)
        scores_df.trace_prob .= preds

        writeArrow(group.scores_path, scores_df)
    end

    !isnothing(pbar) && update(pbar)
end

#############################################################################
# Workspace-level finalize_scoring! — ArrowFileScoringWorkspace
#############################################################################

finalize_scoring!(ws::ArrowFileScoringWorkspace) =
    _finalize_scoring_arrow!(ws.container)

function _finalize_scoring_arrow!(container::ArrowFilePSMContainer)
    for group in container.file_groups
        data_df = DataFrame(Tables.columntable(Arrow.Table(group.data_path)))
        scores_df = DataFrame(Tables.columntable(Arrow.Table(group.scores_path)))
        data_df[!, :trace_prob] = scores_df.trace_prob
        writeArrow(group.data_path, data_df)
    end
end
