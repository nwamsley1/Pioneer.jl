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
    _write_protein_model_fold_label_scatter(df, fold, qc_folder, ss; stage = "iter1", context = "protein_lightgbm", y_col = :precursor_consensus_prefix_shape, y_label = "precursor_consensus_prefix_shape", file_suffix = "prefix_shape")

Write a fold-level `log(pg_score)` scatter against a protein feature
using the semi-supervised training labels from the last fitted iteration.
"""
function _write_protein_model_fold_label_scatter(
    df::AbstractDataFrame,
    fold,
    qc_folder::String,
    ss;
    stage::AbstractString = "iter1",
    context::AbstractString = "protein_lightgbm",
    y_col::Symbol = :precursor_consensus_prefix_shape,
    y_label::AbstractString = "precursor_consensus_prefix_shape",
    file_suffix::AbstractString = "prefix_shape",
    file_idx_to_name::Union{Nothing, AbstractDict{Int64, String}} = nothing
)
    required_cols = (:pg_score, y_col, :target)
    missing_cols = [col for col in required_cols if !hasproperty(df, col)]
    if !isempty(missing_cols)
        @user_warn "Skipping protein model fold label plots due to missing columns" context = context fold = fold missing_columns = missing_cols
        return nothing
    end

    isdir(qc_folder) || mkpath(qc_folder)

    target_positive_mask = ss.positive_mask
    target_negative_mask = ss.mined_negative_mask
    target_confident_positive_mask = (
        hasproperty(ss, :confident_positive_mask) ? ss.confident_positive_mask : target_positive_mask
    ) .& .!target_negative_mask
    target_rescued_positive_mask = target_positive_mask .& .!target_confident_positive_mask
    target_dropped_mask = df.target .& .!target_positive_mask .& .!target_negative_mask
    decoy_mask = .!df.target
    pure_yeast_target_mask = falses(nrow(df))
    has_yeast_named_run = !isnothing(file_idx_to_name) &&
                          any(occursin("Yeast", String(run_name)) for run_name in values(file_idx_to_name))
    if has_yeast_named_run && hasproperty(df, :species) && hasproperty(df, :file_idx) && !isnothing(file_idx_to_name)
        @inbounds for i in eachindex(df.target, df.species, df.file_idx)
            if !df.target[i] || ismissing(df.species[i]) || ismissing(df.file_idx[i])
                continue
            end
            run_name = get(file_idx_to_name, Int64(df.file_idx[i]), "")
            pure_yeast_target_mask[i] = df.target[i] &&
                                        !ismissing(df.species[i]) &&
                                        (strip(String(df.species[i])) == "YEAST") &&
                                        !occursin("Yeast", run_name)
        end
    end
    target_rescued_positive_highlight_mask = target_rescued_positive_mask .& pure_yeast_target_mask
    target_dropped_highlight_mask = target_dropped_mask .& pure_yeast_target_mask
    target_confident_positive_highlight_mask = target_confident_positive_mask .& pure_yeast_target_mask
    target_negative_highlight_mask = target_negative_mask .& pure_yeast_target_mask
    target_rescued_positive_mask = target_rescued_positive_mask .& .!pure_yeast_target_mask
    target_dropped_mask = target_dropped_mask .& .!pure_yeast_target_mask
    target_confident_positive_mask = target_confident_positive_mask .& .!pure_yeast_target_mask
    target_negative_mask = target_negative_mask .& .!pure_yeast_target_mask
    log_pg_score = log.(max.(Float64.(df.pg_score), 1e-6))
    y_values = df[!, y_col]
    target_marker_size = 3.0
    target_negative_marker_size = 3.0
    target_rescued_marker_size = 2.5
    target_dropped_marker_size = 2.5
    decoy_marker_size = 6.0
    plot_size = (1800, 1200)
    plot_dpi = 200

    p_scatter = Plots.plot(
        title = "Protein Model Fold $(fold) $(stage): log(pg_score) vs $(file_suffix)",
        xlabel = "log(pg_score)",
        ylabel = y_label,
        legend = :bottomright,
        background_color = :white,
        foreground_color_legend = nothing,
        size = plot_size,
        dpi = plot_dpi
    )

    if any(target_rescued_positive_mask)
        Plots.scatter!(
            p_scatter,
            log_pg_score[target_rescued_positive_mask],
            y_values[target_rescued_positive_mask];
            label = "Target (+ rescued)",
            color = :dodgerblue,
            alpha = 0.35,
            markersize = target_rescued_marker_size,
            markershape = :diamond,
            markerstrokewidth = 0
        )
    end

    if any(target_dropped_mask)
        Plots.scatter!(
            p_scatter,
            log_pg_score[target_dropped_mask],
            y_values[target_dropped_mask];
            label = "Target (dropped)",
            color = :dodgerblue,
            alpha = 0.35,
            markersize = target_dropped_marker_size,
            markershape = :diamond,
            markerstrokewidth = 0
        )
    end

    if any(target_confident_positive_mask)
        Plots.scatter!(
            p_scatter,
            log_pg_score[target_confident_positive_mask],
            y_values[target_confident_positive_mask];
            label = "Target (+ q/PEP)",
            color = :forestgreen,
            alpha = 0.35,
            markersize = target_marker_size,
            markerstrokewidth = 0
        )
    end

    if any(target_negative_mask)
        Plots.scatter!(
            p_scatter,
            log_pg_score[target_negative_mask],
            y_values[target_negative_mask];
            label = "Target (-)",
            color = :darkorange,
            alpha = 0.30,
            markersize = target_negative_marker_size,
            markerstrokewidth = 0
        )
    end

    if any(decoy_mask)
        Plots.scatter!(
            p_scatter,
            log_pg_score[decoy_mask],
            y_values[decoy_mask];
            label = "Decoy",
            color = :firebrick,
            alpha = 0.20,
            markersize = decoy_marker_size,
            markerstrokewidth = 0
        )
    end

    if any(target_rescued_positive_highlight_mask)
        Plots.scatter!(
            p_scatter,
            log_pg_score[target_rescued_positive_highlight_mask],
            y_values[target_rescued_positive_highlight_mask];
            label = "YEAST Target (+ rescued)",
            color = :dodgerblue,
            alpha = 0.5,
            markersize = 2 * target_rescued_marker_size,
            markershape = :diamond,
            markerstrokewidth = 0.7,
            markerstrokecolor = :black
        )
    end

    if any(target_dropped_highlight_mask)
        Plots.scatter!(
            p_scatter,
            log_pg_score[target_dropped_highlight_mask],
            y_values[target_dropped_highlight_mask];
            label = "YEAST Target (dropped)",
            color = :dodgerblue,
            alpha = 0.5,
            markersize = 2 * target_dropped_marker_size,
            markershape = :diamond,
            markerstrokewidth = 0.7,
            markerstrokecolor = :black
        )
    end

    if any(target_confident_positive_highlight_mask)
        Plots.scatter!(
            p_scatter,
            log_pg_score[target_confident_positive_highlight_mask],
            y_values[target_confident_positive_highlight_mask];
            label = "YEAST Target (+ q/PEP)",
            color = :forestgreen,
            alpha = 0.5,
            markersize = 2 * target_marker_size,
            markerstrokewidth = 0.7,
            markerstrokecolor = :black
        )
    end

    if any(target_negative_highlight_mask)
        Plots.scatter!(
            p_scatter,
            log_pg_score[target_negative_highlight_mask],
            y_values[target_negative_highlight_mask];
            label = "YEAST Target (-)",
            color = :darkorange,
            alpha = 0.5,
            markersize = 2 * target_negative_marker_size,
            markerstrokewidth = 0.7,
            markerstrokecolor = :black
        )
    end

    scatter_path = joinpath(qc_folder, "protein_model_fold_$(fold)_$(stage)_pg_score_vs_$(file_suffix).png")
    savefig(p_scatter, scatter_path)
    @user_info "Wrote protein model fold label scatter context=$(context) fold=$(fold) stage=$(stage) feature=$(y_col) scatter_path=$(scatter_path) n_rows=$(nrow(df))"
    return scatter_path
end

"""
    write_protein_model_fold_label_scatter(df, fold, qc_folder, ss; stage = "iter1", context = "protein_lightgbm")

Write the fold-level `log(pg_score)` vs `precursor_consensus_prefix_shape` scatter
using the semi-supervised training labels from the last fitted iteration.
"""
function write_protein_model_fold_label_scatter(
    df::AbstractDataFrame,
    fold,
    qc_folder::String,
    ss;
    stage::AbstractString = "iter1",
    context::AbstractString = "protein_lightgbm",
    file_idx_to_name::Union{Nothing, AbstractDict{Int64, String}} = nothing
)
    return _write_protein_model_fold_label_scatter(
        df,
        fold,
        qc_folder,
        ss;
        stage = stage,
        context = context,
        y_col = :precursor_consensus_prefix_shape,
        y_label = "precursor_consensus_prefix_shape",
        file_suffix = "prefix_shape",
        file_idx_to_name = file_idx_to_name
    )
end

"""
    write_protein_model_fold_coverage_scatter(df, fold, qc_folder, ss; stage = "iter1", context = "protein_lightgbm")

Write the fold-level `log(pg_score)` vs `coverage_log_ratio` scatter
using the semi-supervised training labels from the last fitted iteration.
"""
function write_protein_model_fold_coverage_scatter(
    df::AbstractDataFrame,
    fold,
    qc_folder::String,
    ss;
    stage::AbstractString = "iter1",
    context::AbstractString = "protein_lightgbm",
    file_idx_to_name::Union{Nothing, AbstractDict{Int64, String}} = nothing
)
    return _write_protein_model_fold_label_scatter(
        df,
        fold,
        qc_folder,
        ss;
        stage = stage,
        context = context,
        y_col = :coverage_log_ratio,
        y_label = "coverage_log_ratio",
        file_suffix = "coverage_log_ratio",
        file_idx_to_name = file_idx_to_name
    )
end

"""
    write_protein_model_fold_peptide_coverage_scatter(df, fold, qc_folder, ss; stage = "iter1", context = "protein_lightgbm")

Write the fold-level `log(pg_score)` vs `peptide_coverage_logit` scatter
using the semi-supervised training labels from the last fitted iteration.
"""
function write_protein_model_fold_peptide_coverage_scatter(
    df::AbstractDataFrame,
    fold,
    qc_folder::String,
    ss;
    stage::AbstractString = "iter1",
    context::AbstractString = "protein_lightgbm",
    file_idx_to_name::Union{Nothing, AbstractDict{Int64, String}} = nothing
)
    return _write_protein_model_fold_label_scatter(
        df,
        fold,
        qc_folder,
        ss;
        stage = stage,
        context = context,
        y_col = :peptide_coverage_logit,
        y_label = "peptide_coverage_logit",
        file_suffix = "peptide_coverage_logit",
        file_idx_to_name = file_idx_to_name
    )
end

function make_protein_iteration_plot_callback(
    train_df::AbstractDataFrame,
    fold,
    qc_folder::String;
    context::AbstractString,
    file_idx_to_name::Union{Nothing, AbstractDict{Int64, String}} = nothing
)
    return function(iteration, plot_state)
        write_protein_model_fold_label_scatter(
            train_df,
            fold,
            qc_folder,
            plot_state;
            stage = "iter_$(iteration)",
            context = context,
            file_idx_to_name = file_idx_to_name
        )
        write_protein_model_fold_coverage_scatter(
            train_df,
            fold,
            qc_folder,
            plot_state;
            stage = "iter_$(iteration)",
            context = context,
            file_idx_to_name = file_idx_to_name
        )
        write_protein_model_fold_peptide_coverage_scatter(
            train_df,
            fold,
            qc_folder,
            plot_state;
            stage = "iter_$(iteration)",
            context = context,
            file_idx_to_name = file_idx_to_name
        )
        return nothing
    end
end
