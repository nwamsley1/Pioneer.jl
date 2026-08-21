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
    load_run_level_protein_training_rows(pg_refs; include_qc_plot_columns = false)

Load only the columns needed for run-level protein-model fitting from
protein-group files. This avoids materializing string-heavy output columns like
`peptide_list` that are not used by the model.
"""
function load_run_level_protein_training_rows(
    pg_refs::Vector{ProteinGroupFileReference};
    include_qc_plot_columns::Bool = false
)
    columns_to_load = Symbol[
        :protein_name,
        :target,
        :n_non_mbr_peptides,
    ]
    append!(columns_to_load, run_level_protein_feature_names())
    if include_qc_plot_columns
        push!(columns_to_load, :species, :file_idx)
    end
    unique!(columns_to_load)

    all_protein_groups = DataFrame()
    for pg_ref in pg_refs
        pg_path = file_path(pg_ref)
        isfile(pg_path) || continue

        tbl = Arrow.Table(pg_path)
        n_rows = length(tbl[:protein_name])
        n_rows == 0 && continue

        available_columns = Set(propertynames(tbl))
        missing_columns = [col for col in columns_to_load if !(col in available_columns)]
        isempty(missing_columns) || error("Protein group file $pg_path is missing required columns: $missing_columns")

        chunk_df = DataFrame()
        for col in columns_to_load
            chunk_df[!, col] = collect(tbl[col])
        end

        if ncol(all_protein_groups) == 0
            all_protein_groups = chunk_df
        else
            append!(all_protein_groups, chunk_df)
        end
    end

    return all_protein_groups
end

"""
    prepare_run_level_protein_training_rows(actual_groups, shadow_groups)

Attach explicit training labels to reported protein groups and append
target-backed counterfactual shadow negatives. Actual MBR-only singleton target
groups remain scoreable but are excluded from the positive seed; treating them
as positives would teach the model the failure mode the shadow controls are
designed to expose.
"""
function prepare_run_level_protein_training_rows(
    actual_groups::DataFrame,
    shadow_groups::DataFrame,
)
    groups = copy(actual_groups)
    nrow(groups) == 0 && return groups

    groups[!, :protein_training_label] = Bool.(groups.target)
    groups[!, :protein_training_eligible] = .!(
        Bool.(groups.target) .&
        Bool.(groups.mbr_only_protein) .&
        (Int.(groups.mbr_recovered_peptides) .== 1)
    )
    groups[!, :protein_shadow_negative] = falses(nrow(groups))

    if nrow(shadow_groups) > 0
        required = Symbol.(names(actual_groups))
        missing_columns = [
            column for column in required
            if !hasproperty(shadow_groups, column)
        ]
        isempty(missing_columns) || error(
            "Counterfactual shadow protein groups are missing columns: " *
            string(missing_columns)
        )
        shadows = select(shadow_groups, required)
        shadows[!, :protein_training_label] = falses(nrow(shadows))
        shadows[!, :protein_training_eligible] = trues(nrow(shadows))
        shadows[!, :protein_shadow_negative] = trues(nrow(shadows))
        append!(groups, shadows; cols = :setequal)
    end

    return groups
end

"""
    perform_run_level_protein_scoring(pg_refs::Vector{ProteinGroupFileReference},
                                      max_in_memory_rows::Int64,
                                      qc_folder::String,
                                      precursors::LibraryPrecursors;
                                      protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}})

Fit and apply the run-level protein model.

# Arguments
- `pg_refs`: Vector of protein group file references
- `max_in_memory_rows`: Maximum number of protein-group rows allowed in memory
- `qc_folder`: Folder for QC plots
- `precursors`: Library precursors
- `protein_to_cv_fold`: Pre-built mapping of proteins to CV folds

Returns `true` only when models were trained and applied to every fold.
"""
function perform_run_level_protein_scoring(
    pg_refs::Vector{ProteinGroupFileReference},
    max_in_memory_rows::Int64,
    qc_folder::String,
    precursors::LibraryPrecursors;
    protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
    file_idx_to_name::Union{Nothing, AbstractDict{Int64, String}} = nothing,
    write_qc_plots::Bool = true,
    train_q_value_threshold::Float32 = 0.01f0,
    min_prefix_shape_neg_threshold_itr::Float32 = -0.20f0,
    min_pep_neg_threshold_itr::Float32 = 0.90f0,
    counterfactual_shadow_groups::DataFrame = DataFrame(),
)
    total_protein_groups = 0
    for ref in pg_refs
        if exists(ref)
            total_protein_groups += row_count(ref)
        end
    end
    total_protein_groups += nrow(counterfactual_shadow_groups)

    max_protein_groups_in_memory_limit = max(max_in_memory_rows, 100_000)

    if total_protein_groups > max_protein_groups_in_memory_limit
        error("Run-level protein scoring does not support out-of-memory fitting. total_protein_groups=$(total_protein_groups) exceeds max_protein_groups_in_memory_limit=$(max_protein_groups_in_memory_limit).")
    end

    actual_protein_groups = load_run_level_protein_training_rows(
        pg_refs;
        include_qc_plot_columns = write_qc_plots
    )
    all_protein_groups = prepare_run_level_protein_training_rows(
        actual_protein_groups,
        counterfactual_shadow_groups,
    )

    eligible = all_protein_groups.protein_training_eligible
    labels = all_protein_groups.protein_training_label
    n_targets = count(eligible .& labels)
    n_decoys = count(eligible .& .!labels)
    n_shadow_negatives = count(all_protein_groups.protein_shadow_negative)
    n_excluded_singletons = count(.!eligible)
    @debug_l1 "Run-level protein training controls: " *
              "counterfactual shadows=$(n_shadow_negatives), " *
              "excluded target MBR-only singleton positives=$(n_excluded_singletons)"
    skip_scoring = !(n_targets > 50 && n_decoys > 50 && nrow(all_protein_groups) > 1000)

    return perform_protein_scoring_multifold(
        all_protein_groups,
        qc_folder,
        pg_refs,
        precursors;
        protein_to_cv_fold = protein_to_cv_fold,
        file_idx_to_name = file_idx_to_name,
        skip_scoring = skip_scoring,
        write_qc_plots = write_qc_plots,
        train_q_value_threshold = train_q_value_threshold,
        min_prefix_shape_neg_threshold_itr = min_prefix_shape_neg_threshold_itr,
        min_pep_neg_threshold_itr = min_pep_neg_threshold_itr
    )
end

function run_protein_scoring!(
    search_context::SearchContext;
    passing_refs::Vector{PSMFileReference},
    protein_ambiguity_candidates::Dict{UInt32, Vector{ProteinKey}} =
        Dict{UInt32, Vector{ProteinKey}}(),
    protein_peptide_opportunities::Dict{
        ProteinKey,
        ProteinPeptideOpportunityCounts
    } = Dict{ProteinKey, ProteinPeptideOpportunityCounts}(),
    max_in_memory_table_mb::Float64,
    q_value_threshold::Float32,
    min_peptides::Int64,
    write_qc_plots::Bool,
    min_pep_neg_threshold_itr::Float32,
    q_value_interpolation_points_per_bin::Int64
)
    isempty(passing_refs) && return nothing

    temp_folder = joinpath(getDataOutDir(search_context), "temp_data")
    passing_proteins_folder = joinpath(temp_folder, "passing_proteins")
    !isdir(passing_proteins_folder) && mkpath(passing_proteins_folder)

    qc_folder = joinpath(dirname(temp_folder), "qc_plots")
    !isdir(qc_folder) && mkpath(qc_folder)

    precursors = getPrecursors(getSpecLib(search_context))
    isempty(protein_peptide_opportunities) &&
        error("Protein scoring requires theoretical unique/shared peptide opportunities")

    file_idx_to_name = Dict{Int64, String}()
    for (file_idx, file_name) in enumerate(getFileIdToName(getMSData(search_context)))
        file_idx_to_name[Int64(file_idx)] = String(file_name)
    end

    (
        pg_refs,
        psm_to_pg_mapping,
        protein_to_cv_fold,
        counterfactual_shadow_groups,
    ) = build_protein_group_tables(
        passing_refs,
        passing_proteins_folder,
        protein_peptide_opportunities,
        precursors = precursors,
        protein_ambiguity_candidates = protein_ambiguity_candidates,
        min_peptides = min_peptides,
        q_value_threshold = q_value_threshold
    )

    indexed_paths = get_all_indexed_paths(getPassingPsms, search_context)
    paired_files = PairedSearchFiles[]
    for (file_idx, psm_path) in indexed_paths
        haskey(psm_to_pg_mapping, psm_path) || continue
        pg_path = psm_to_pg_mapping[psm_path]
        push!(paired_files, PairedSearchFiles(psm_path, pg_path, file_idx))
        setPassingProteins!(getMSData(search_context), file_idx, pg_path)
    end

    isempty(paired_files) && error("No protein groups created during protein inference")

    max_in_memory_rows = estimate_max_rows(max_in_memory_table_mb, file_path(first(pg_refs)))
    @debug_l1 "Memory budget $(max_in_memory_table_mb) MB → max_protein_groups = $max_in_memory_rows"
    protein_model_scoring_succeeded = perform_run_level_protein_scoring(
        pg_refs,
        max_in_memory_rows,
        qc_folder,
        precursors;
        protein_to_cv_fold = protein_to_cv_fold,
        file_idx_to_name = file_idx_to_name,
        write_qc_plots = write_qc_plots,
        train_q_value_threshold = q_value_threshold,
        min_pep_neg_threshold_itr = min_pep_neg_threshold_itr,
        counterfactual_shadow_groups = counterfactual_shadow_groups,
    )

    n_proteins = length(getProteins(getSpecLib(search_context)))
    sorted_pg_scores_path = joinpath(temp_folder, "sorted_pg_scores.arrow")
    spline_result = build_qvalue_spline_from_refs(
        pg_refs,
        :pg_score,
        sorted_pg_scores_path;
        batch_size = 1_000_000,
        compute_pep = true,
        min_pep_points_per_bin = q_value_interpolation_points_per_bin,
        temp_prefix = "preglobal_pg_sidecar",
    )
    spline_result === nothing &&
        error("Cannot build global protein model without run-level scores")
    pg_qval_spline = spline_result.qval_spline
    run_score_floor = _score_floor_for_qvalue(
        pg_qval_spline,
        q_value_threshold,
    )
    precursor_scoring_results = get_results(
        search_context,
        PrecursorScoringSearch,
    )
    run_similarity = if precursor_scoring_results isa PrecursorScoringSearchResults
        precursor_scoring_results.run_similarity[]
    else
        nothing
    end
    n_runs_total = length(getFilePaths(getMSData(search_context)))

    global_pg_score_dict, pg_name_to_global_pg_score =
        build_global_protein_score_dicts(
            pg_refs,
            protein_to_cv_fold,
            n_proteins,
            length(pg_refs),
            protein_model_scoring_succeeded,
            run_score_floor,
            n_experiment_runs = n_runs_total,
            run_similarity = run_similarity,
            q_value_threshold = q_value_threshold,
        )
    search_context.pg_name_to_global_pg_score[] = pg_name_to_global_pg_score

    global_pg_qval_dict = build_protein_global_qval_dict(global_pg_score_dict)
    search_context.global_pg_score_to_qval_dict[] = global_pg_qval_dict

    search_context.pg_score_to_qval[] = pg_qval_spline
    search_context.pg_score_to_pep[] = spline_result.pep_interp

    protein_combined_pipeline = TransformPipeline() |>
        add_dict_column_composite_key(:global_pg_score, [:protein_name, :target, :entrap_id], global_pg_score_dict) |>
        add_dict_column_composite_key(:global_pg_qval, [:protein_name, :target, :entrap_id], global_pg_qval_dict) |>
        add_interpolated_column(:pg_qval, :pg_score, search_context.pg_score_to_qval[]) |>
        add_interpolated_column(:pg_pep, :pg_score, search_context.pg_score_to_pep[]) |>
        filter_by_multiple_thresholds([
            (:global_pg_qval, q_value_threshold),
            (:pg_qval, q_value_threshold)
        ])

    apply_pipeline!(pg_refs, protein_combined_pipeline)

    spline_result = build_qvalue_spline_from_refs(pg_refs, :pg_score, sorted_pg_scores_path;
        batch_size = 1_000_000,
        min_pep_points_per_bin = q_value_interpolation_points_per_bin,
        temp_prefix = "pg_recalc")
    search_context.pg_score_to_qval[] = spline_result.qval_spline

    recalc_pg_pipeline = TransformPipeline() |>
        add_interpolated_column(:pg_qval, :pg_score, search_context.pg_score_to_qval[])

    pg_refs = apply_pipeline_batch(pg_refs, recalc_pg_pipeline, passing_proteins_folder)

    update_psms_with_protein_scores_refs(
        paired_files,
        search_context.pg_name_to_global_pg_score[],
        search_context.pg_score_to_qval[],
        search_context.global_pg_score_to_qval_dict[]
    )

    return nothing
end
