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

"""
Fold-specific iRT correction based on current predicted iRT and amino-acid composition.

The refinement model is trained after each precursor rescoring iteration using
only training-fold target precursors that pass the requested FDR threshold. It
predicts an iRT correction term for the held-out fold, which is added to the
current `irt_pred` values before the next rescoring iteration begins.
"""

struct AminoAcidCountIrtModel
    intercept::Float32
    irt_pred_coefficient::Float32
    irt_pred_squared_coefficient::Float32
    coefficients::Dict{String, Float32}
end

mutable struct IrtLinearRefinement{P<:LibraryPrecursors} <: IterationPostProcessStrategy
    precursors::P
    q_value_threshold::Float32
    min_precursors::Int
    qc_plot_dir::Union{Nothing, String}
    qc_run_id::Union{Nothing, UInt32}
    run_labels::Dict{UInt32, String}
    token_cache::Dict{UInt32, Dict{String, Float32}}
end

function IrtLinearRefinement(
    precursors::LibraryPrecursors;
    q_value_threshold::Float32 = 0.01f0,
    min_precursors::Int = 2,
    qc_plot_dir::Union{Nothing, String} = nothing,
    qc_run_id::Union{Nothing, UInt32} = nothing,
    run_labels::Dict{UInt32, String} = Dict{UInt32, String}()
)
    return IrtLinearRefinement(
        precursors,
        q_value_threshold,
        min_precursors,
        qc_plot_dir,
        qc_run_id,
        run_labels,
        Dict{UInt32, Dict{String, Float32}}()
    )
end

function _get_irt_errors(df::AbstractDataFrame)
    if hasproperty(df, :irt_error)
        return Float32.(df.irt_error)
    end
    return abs.(Float32.(df.irt_obs) .- Float32.(df.irt_pred))
end

_sorted_run_ids(run_ids::AbstractVector) = sort!(unique(UInt32.(run_ids)))

@inline function _irt_pred_basis(current_irt_pred::Float32)
    x = Float64(current_irt_pred)
    return x, x * x
end

function _get_run_label(strategy::IrtLinearRefinement, run_id::UInt32)
    return get(strategy.run_labels, run_id, "Run $(run_id)")
end

function _add_error_distribution!(
    p,
    errors::Vector{Float32},
    targets::AbstractVector{Bool},
    subplot_idx::Int;
    title::String
)
    valid_mask = isfinite.(errors)
    target_errors = errors[valid_mask .& targets]
    decoy_errors = errors[valid_mask .& .!targets]

    plot!(
        p,
        title = title,
        xlabel = "Absolute iRT error",
        ylabel = "Density",
        legend = :topright,
        subplot = subplot_idx
    )

    if !isempty(target_errors)
        histogram!(
            p,
            target_errors;
            normalize = :pdf,
            bins = 50,
            alpha = 0.35,
            color = :steelblue,
            linecolor = :steelblue4,
            label = "Targets",
            subplot = subplot_idx
        )
    end

    if !isempty(decoy_errors)
        histogram!(
            p,
            decoy_errors;
            normalize = :pdf,
            bins = 50,
            alpha = 0.35,
            color = :firebrick,
            linecolor = :firebrick4,
            label = "Decoys",
            subplot = subplot_idx
        )
    end

    return nothing
end

function maybe_write_irt_refinement_qc_plot(
    strategy::IrtLinearRefinement,
    run_id::UInt32,
    irt_error_before::Vector{Float32},
    irt_error_after::Vector{Float32},
    targets::AbstractVector{Bool}
)
    isnothing(strategy.qc_plot_dir) && return nothing
    isempty(irt_error_before) && return nothing

    isdir(strategy.qc_plot_dir) || mkpath(strategy.qc_plot_dir)
    plot_path = joinpath(
        strategy.qc_plot_dir,
        "irt_error_refinement_run_$(run_id).png"
    )
    run_label = _get_run_label(strategy, run_id)

    try
        p = plot(
            layout = (1, 2),
            size = (1400, 550),
            titlefont = 12,
            margin = 8Plots.mm
        )
        _add_error_distribution!(
            p,
            irt_error_before,
            targets,
            1;
            title = "Before refinement: $(run_label)"
        )
        _add_error_distribution!(
            p,
            irt_error_after,
            targets,
            2;
            title = "After refinement: $(run_label)"
        )
        savefig(p, plot_path)
    catch e
        @user_warn "Failed to write iRT refinement QC plot" plot_path = plot_path error_type = string(typeof(e))
    end

    return nothing
end

function maybe_write_irt_refinement_qc_plots(
    strategy::IrtLinearRefinement,
    run_ids::AbstractVector,
    irt_error_before::Vector{Float32},
    irt_error_after::Vector{Float32},
    targets::AbstractVector{Bool}
)
    isnothing(strategy.qc_plot_dir) && return nothing

    run_ids_u32 = UInt32.(run_ids)
    run_id = strategy.qc_run_id
    if isnothing(run_id)
        isempty(run_ids_u32) && return nothing
        run_id = first(_sorted_run_ids(run_ids_u32))
    end

    run_mask = run_ids_u32 .== run_id
    any(run_mask) || return nothing

    maybe_write_irt_refinement_qc_plot(
        strategy,
        run_id,
        irt_error_before[run_mask],
        irt_error_after[run_mask],
        targets[run_mask]
    )

    return nothing
end

function precursor_token_counts(
    strategy::IrtLinearRefinement,
    precursor_idx::Integer
)
    pid = UInt32(precursor_idx)
    return get!(strategy.token_cache, pid) do
        sequence = String(getSequence(strategy.precursors)[pid])
        structural_mods = getStructuralMods(strategy.precursors)[pid]
        return _precursor_token_counts(sequence, structural_mods)
    end
end

function _precursor_token_counts(
    sequence::AbstractString,
    structural_mods::Union{Missing, AbstractString}
)
    counts = Dict{String, Float32}()
    position_mods = Dict{Int, Vector{String}}()
    residue_tokens = String[]

    if !ismissing(structural_mods) && !isempty(structural_mods)
        mod_regex = r"\((\d+),([A-Z]|[nc]),([^,\)]+)\)"
        for m in eachmatch(mod_regex, structural_mods)
            position = parse(Int, m.captures[1])
            site = only(m.captures[2])
            mod_name = m.captures[3]

            if site == 'n' || site == 'c'
                token = string(site, "|", mod_name)
                counts[token] = get(counts, token, 0.0f0) + 1.0f0
            elseif 1 <= position <= length(sequence)
                push!(get!(position_mods, position, String[]), mod_name)
            end
        end
    end

    for (position, aa) in enumerate(sequence)
        mods_here = get(position_mods, position, nothing)
        token = if isnothing(mods_here) || isempty(mods_here)
            string(aa)
        else
            string(aa, "|", join(sort(mods_here), "&"))
        end
        push!(residue_tokens, token)
        counts[token] = get(counts, token, 0.0f0) + 1.0f0
    end

    if !isempty(residue_tokens)
        n_term_token = "NTERM|" * first(residue_tokens)
        c_term_token = "CTERM|" * last(residue_tokens)
        counts[n_term_token] = get(counts, n_term_token, 0.0f0) + 1.0f0
        counts[c_term_token] = get(counts, c_term_token, 0.0f0) + 1.0f0
    end

    return counts
end

function fit_irt_refinement_model(
    strategy::IrtLinearRefinement,
    precursor_ids::Vector{UInt32},
    irt_pred_inputs::Vector{Float32},
    irt_corrections::Vector{Float32}
)
    n = length(precursor_ids)
    if n < strategy.min_precursors ||
       n != length(irt_pred_inputs) ||
       n != length(irt_corrections)
        return nothing
    end

    token_to_col = Dict{String, Int}()
    precursor_counts = Vector{Dict{String, Float32}}(undef, n)

    for i in eachindex(precursor_ids)
        counts = precursor_token_counts(strategy, precursor_ids[i])
        precursor_counts[i] = counts
        for token in keys(counts)
            if !haskey(token_to_col, token)
                token_to_col[token] = length(token_to_col) + 1
            end
        end
    end

    if isempty(token_to_col)
        return nothing
    end

    use_quadratic_basis = n >= 3
    irt_basis_cols = use_quadratic_basis ? 2 : 1
    X = zeros(Float64, n, length(token_to_col) + 1 + irt_basis_cols)
    X[:, 1] .= 1.0
    for i in eachindex(precursor_counts)
        linear_term, quadratic_term = _irt_pred_basis(irt_pred_inputs[i])
        X[i, 2] = linear_term
        if use_quadratic_basis
            X[i, 3] = quadratic_term
        end
        for (token, count) in precursor_counts[i]
            X[i, token_to_col[token] + 1 + irt_basis_cols] = Float64(count)
        end
    end

    y = Float64.(irt_corrections)
    β = try
        X \ y
    catch
        return nothing
    end

    coefficients = Dict{String, Float32}()
    for (token, col_idx) in token_to_col
        coefficients[token] = Float32(β[col_idx + 1 + irt_basis_cols])
    end

    return AminoAcidCountIrtModel(
        Float32(β[1]),
        Float32(β[2]),
        use_quadratic_basis ? Float32(β[3]) : 0.0f0,
        coefficients
    )
end

function predict_irt_refinement(
    model::AminoAcidCountIrtModel,
    counts::Dict{String, Float32},
    current_irt_pred::Float32
)
    linear_term, quadratic_term = _irt_pred_basis(current_irt_pred)
    prediction = model.intercept +
                 model.irt_pred_coefficient * Float32(linear_term) +
                 model.irt_pred_squared_coefficient * Float32(quadratic_term)
    for (token, count) in counts
        prediction += get(model.coefficients, token, 0.0f0) * count
    end
    return prediction
end

function predict_irt_refinement(
    strategy::IrtLinearRefinement,
    model::AminoAcidCountIrtModel,
    precursor_idx::Integer,
    current_irt_pred::Float32
)
    return predict_irt_refinement(
        model,
        precursor_token_counts(strategy, precursor_idx),
        current_irt_pred
    )
end

function _passing_precursor_targets(
    precursor_idx::AbstractVector,
    targets::AbstractVector{Bool},
    trace_prob::AbstractVector,
    irt_pred::AbstractVector,
    irt_obs::AbstractVector,
    q_value_threshold::Float32
)
    n = length(trace_prob)
    n == 0 && return UInt32[], Float32[], Float32[]

    q_values = Vector{Float32}(undef, n)
    get_qvalues!(Float32.(trace_prob), collect(Bool, targets), q_values)
    pass_mask = collect(Bool, targets) .& (q_values .<= q_value_threshold)

    sum_pred = Dict{UInt32, Float64}()
    sum_irt = Dict{UInt32, Float64}()
    count_irt = Dict{UInt32, Int}()
    for i in eachindex(pass_mask)
        if pass_mask[i]
            pid = UInt32(precursor_idx[i])
            sum_pred[pid] = get(sum_pred, pid, 0.0) + Float64(irt_pred[i])
            sum_irt[pid] = get(sum_irt, pid, 0.0) + Float64(irt_obs[i])
            count_irt[pid] = get(count_irt, pid, 0) + 1
        end
    end

    precursor_ids = collect(keys(sum_irt))
    sort!(precursor_ids)
    irt_pred_inputs = Float32[sum_pred[pid] / count_irt[pid] for pid in precursor_ids]
    irt_corrections = Float32[
        (sum_irt[pid] / count_irt[pid]) - (sum_pred[pid] / count_irt[pid])
        for pid in precursor_ids
    ]
    return precursor_ids, irt_pred_inputs, irt_corrections
end

function _fit_fold_models(
    df::AbstractDataFrame,
    strategy::IrtLinearRefinement,
    cv_folds::AbstractVector{UInt8}
)
    models = Dict{Tuple{UInt8, UInt32}, Union{Nothing, AminoAcidCountIrtModel}}()
    run_ids = _sorted_run_ids(df.ms_file_idx)
    run_col = UInt32.(df.ms_file_idx)

    for fold in cv_folds
        for run_id in run_ids
            train_mask = (df.cv_fold .!= fold) .& (run_col .== run_id)
            precursor_ids, irt_pred_inputs, irt_corrections = _passing_precursor_targets(
                df.precursor_idx[train_mask],
                df.target[train_mask],
                df.trace_prob[train_mask],
                df.irt_pred[train_mask],
                df.irt_obs[train_mask],
                strategy.q_value_threshold
            )
            models[(fold, run_id)] = fit_irt_refinement_model(
                strategy,
                precursor_ids,
                irt_pred_inputs,
                irt_corrections
            )
        end
    end

    return models
end

function _apply_fold_models_to_dataframe!(
    df::AbstractDataFrame,
    strategy::IrtLinearRefinement,
    models::Dict{Tuple{UInt8, UInt32}, Union{Nothing, AminoAcidCountIrtModel}}
)
    current_pred = Float32.(df.irt_pred)
    refined_pred = copy(current_pred)
    run_col = UInt32.(df.ms_file_idx)

    for ((fold, run_id), model) in models
        isnothing(model) && continue

        fold_mask = (df.cv_fold .== fold) .& (run_col .== run_id)
        if !any(fold_mask)
            continue
        end

        for row_idx in findall(fold_mask)
            correction = predict_irt_refinement(
                strategy,
                model,
                df.precursor_idx[row_idx],
                current_pred[row_idx]
            )
            refined_pred[row_idx] = current_pred[row_idx] + correction
        end
    end

    df[!, :irt_pred] = refined_pred
    if hasproperty(df, :irt_error)
        df[!, :irt_error] = abs.(Float32.(df.irt_obs) .- refined_pred)
    end

    return nothing
end

function apply_iteration_postprocess!(
    strategy::IrtLinearRefinement,
    workspace::InMemoryScoringWorkspace,
    iteration::Int,
    total_iterations::Int
)
    if total_iterations < 1
        return nothing
    end

    df = to_dataframe(get_psms(workspace))
    irt_error_before = _get_irt_errors(df)
    targets = Bool.(df.target)
    run_ids = UInt32.(df.ms_file_idx)
    cv_folds = sort!(collect(get_cv_folds(workspace)))
    models = _fit_fold_models(df, strategy, cv_folds)
    _apply_fold_models_to_dataframe!(df, strategy, models)
    if iteration == total_iterations
        maybe_write_irt_refinement_qc_plots(
            strategy,
            run_ids,
            irt_error_before,
            _get_irt_errors(df),
            targets
        )
    end

    return nothing
end

function _compute_fold_prob_threshold(
    fold_groups::Vector{ArrowFileGroup},
    q_value_threshold::Float32
)
    temp_refs = PSMFileReference[]

    for group in fold_groups
        n = group.n_rows
        n == 0 && continue

        data_tbl = Arrow.Table(group.data_path)
        scores_tbl = Arrow.Table(group.scores_path)
        probs = Vector{Float32}(undef, n)
        targets = Vector{Bool}(undef, n)

        @inbounds for i in 1:n
            probs[i] = Float32(scores_tbl.trace_prob[i])
            targets[i] = Bool(data_tbl.target[i])
        end

        perm = sortperm(probs; rev=true)
        probs .= probs[perm]
        targets .= targets[perm]

        temp_path = tempname() * "_irt_refine_fdr.arrow"
        writeArrow(temp_path, DataFrame(trace_prob = probs, target = targets))
        ref = PSMFileReference(temp_path)
        mark_sorted!(ref, :trace_prob)
        push!(temp_refs, ref)
    end

    if isempty(temp_refs)
        return typemax(Float32)
    end

    merged_path = tempname() * "_irt_refine_fdr_merged.arrow"
    stream_sorted_merge(temp_refs, merged_path, :trace_prob; reverse=true)

    GC.gc(false)
    for ref in temp_refs
        safeRm(file_path(ref), nothing; force=true)
    end

    return _compute_prob_threshold_from_merged(merged_path, q_value_threshold)
end

function _get_group_ms_file_idx(group::ArrowFileGroup)
    group.n_rows == 0 && return nothing
    return UInt32(Arrow.Table(group.data_path).ms_file_idx[1])
end

function _get_group_run_ids(groups::Vector{ArrowFileGroup})
    run_ids = UInt32[]
    for group in groups
        run_id = _get_group_ms_file_idx(group)
        isnothing(run_id) && continue
        push!(run_ids, run_id)
    end
    return _sorted_run_ids(run_ids)
end

function _aggregate_passing_precursors(
    fold_groups::Vector{ArrowFileGroup},
    prob_threshold::Float32
)
    prob_threshold == typemax(Float32) && return UInt32[], Float32[], Float32[]

    sum_pred = Dict{UInt32, Float64}()
    sum_irt = Dict{UInt32, Float64}()
    count_irt = Dict{UInt32, Int}()

    for group in fold_groups
        data_tbl = Arrow.Table(group.data_path)
        scores_tbl = Arrow.Table(group.scores_path)

        @inbounds for i in 1:group.n_rows
            if Bool(data_tbl.target[i]) && Float32(scores_tbl.trace_prob[i]) >= prob_threshold
                pid = UInt32(data_tbl.precursor_idx[i])
                sum_pred[pid] = get(sum_pred, pid, 0.0) + Float64(data_tbl.irt_pred[i])
                sum_irt[pid] = get(sum_irt, pid, 0.0) + Float64(data_tbl.irt_obs[i])
                count_irt[pid] = get(count_irt, pid, 0) + 1
            end
        end
    end

    precursor_ids = collect(keys(sum_irt))
    sort!(precursor_ids)
    irt_pred_inputs = Float32[sum_pred[pid] / count_irt[pid] for pid in precursor_ids]
    irt_corrections = Float32[
        (sum_irt[pid] / count_irt[pid]) - (sum_pred[pid] / count_irt[pid])
        for pid in precursor_ids
    ]
    return precursor_ids, irt_pred_inputs, irt_corrections
end

function _apply_fold_model_to_groups!(
    fold_groups::Vector{ArrowFileGroup},
    strategy::IrtLinearRefinement,
    models::Dict{Tuple{UInt8, UInt32}, Union{Nothing, AminoAcidCountIrtModel}},
    fold::UInt8
)
    irt_error_before = Float32[]
    irt_error_after = Float32[]
    targets = Bool[]
    run_ids = UInt32[]

    for group in fold_groups
        data_df = DataFrame(Tables.columntable(Arrow.Table(group.data_path)))
        nrow(data_df) == 0 && continue

        run_id = UInt32(first(data_df.ms_file_idx))
        append!(irt_error_before, _get_irt_errors(data_df))
        append!(targets, Bool.(data_df.target))
        append!(run_ids, fill(run_id, nrow(data_df)))

        model = get(models, (fold, run_id), nothing)

        if !isnothing(model)
            _apply_fold_models_to_dataframe!(
                data_df,
                strategy,
                Dict((fold, run_id) => model)
            )
            writeArrow(group.data_path, data_df)
        end

        append!(irt_error_after, _get_irt_errors(data_df))
    end

    return irt_error_before, irt_error_after, targets, run_ids
end

function apply_iteration_postprocess!(
    strategy::IrtLinearRefinement,
    workspace::ArrowFileScoringWorkspace,
    iteration::Int,
    total_iterations::Int
)
    if total_iterations < 1
        return nothing
    end

    container = workspace.container
    cv_folds = sort!(collect(get_cv_folds(container)))
    run_ids = _get_group_run_ids(container.file_groups)
    models = Dict{Tuple{UInt8, UInt32}, Union{Nothing, AminoAcidCountIrtModel}}()

    for fold in cv_folds
        train_groups = [group for group in container.file_groups if get_fold_number(group) != fold]
        for run_id in run_ids
            run_train_groups = ArrowFileGroup[
                group for group in train_groups if _get_group_ms_file_idx(group) == run_id
            ]
            prob_threshold = _compute_fold_prob_threshold(run_train_groups, strategy.q_value_threshold)
            precursor_ids, irt_pred_inputs, irt_corrections = _aggregate_passing_precursors(
                run_train_groups,
                prob_threshold
            )
            models[(fold, run_id)] = fit_irt_refinement_model(
                strategy,
                precursor_ids,
                irt_pred_inputs,
                irt_corrections
            )
        end
    end

    irt_error_before = Float32[]
    irt_error_after = Float32[]
    targets = Bool[]
    run_ids_per_row = UInt32[]
    for fold in cv_folds
        test_groups = get_file_groups_for_fold(container, fold)
        fold_before, fold_after, fold_targets, fold_run_ids = _apply_fold_model_to_groups!(
            test_groups,
            strategy,
            models,
            fold
        )
        append!(irt_error_before, fold_before)
        append!(irt_error_after, fold_after)
        append!(targets, fold_targets)
        append!(run_ids_per_row, fold_run_ids)
    end

    sampled_df = to_dataframe(get_psms(workspace))
    _apply_fold_models_to_dataframe!(sampled_df, strategy, models)
    if iteration == total_iterations
        maybe_write_irt_refinement_qc_plots(
            strategy,
            run_ids_per_row,
            irt_error_before,
            irt_error_after,
            targets
        )
    end

    return nothing
end
