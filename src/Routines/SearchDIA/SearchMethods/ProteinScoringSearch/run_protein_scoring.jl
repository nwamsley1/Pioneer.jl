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
    count_protein_peptides(precursors::LibraryPrecursors)

Count all possible peptides for each protein in the library.

# Arguments
- `precursors`: Library precursors

# Returns
- Dictionary mapping protein keys to sets of peptide sequences
"""
function count_protein_peptides(precursors::LibraryPrecursors)
    protein_to_possible_peptides = Dict{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8}, Set{String}}()

    all_accession_numbers = getAccessionNumbers(precursors)
    all_sequences = getSequence(precursors)
    all_decoys = getIsDecoy(precursors)
    all_entrap_ids = getEntrapmentGroupId(precursors)
    n_precursors = length(all_accession_numbers)

    for i in 1:n_precursors
        protein_names = split(all_accession_numbers[i], ';')
        is_decoy = all_decoys[i]
        entrap_id = all_entrap_ids[i]

        for protein_name in protein_names
            key = (protein_name = String(protein_name), target = !is_decoy, entrap_id = entrap_id)
            if !haskey(protein_to_possible_peptides, key)
                protein_to_possible_peptides[key] = Set{String}()
            end
            push!(protein_to_possible_peptides[key], all_sequences[i])
        end
    end

    return protein_to_possible_peptides
end

"""
    perform_protein_probit_regression(pg_refs::Vector{ProteinGroupFileReference},
                                    max_in_memory_rows::Int64,
                                    qc_folder::String,
                                    precursors::LibraryPrecursors;
                                    protein_to_cv_fold::Union{Nothing, Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}} = nothing)

Perform probit regression on protein groups.

# Arguments
- `pg_refs`: Vector of protein group file references
- `max_in_memory_rows`: Maximum number of protein-group rows allowed in memory
- `qc_folder`: Folder for QC plots
- `precursors`: Library precursors
- `protein_to_cv_fold`: Optional pre-built mapping of proteins to CV folds
"""
function perform_protein_probit_regression(
    pg_refs::Vector{ProteinGroupFileReference},
    max_in_memory_rows::Int64,
    qc_folder::String,
    precursors::LibraryPrecursors;
    protein_to_cv_fold::Union{Nothing, Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}} = nothing,
    file_idx_to_name::Union{Nothing, AbstractDict{Int64, String}} = nothing,
    write_qc_plots::Bool = true,
    log_feature_importance::Bool = true,
    train_q_value_threshold::Float32 = 0.01f0,
    min_prefix_shape_neg_threshold_itr::Float32 = -0.20f0,
    min_pep_neg_threshold_itr::Float32 = 0.90f0
)
    passing_pg_paths = [file_path(ref) for ref in pg_refs]

    total_protein_groups = 0
    for ref in pg_refs
        if exists(ref)
            total_protein_groups += row_count(ref)
        end
    end

    max_protein_groups_in_memory_limit = max(max_in_memory_rows, 100_000)

    if total_protein_groups > max_protein_groups_in_memory_limit
        error("Protein probit out-of-memory processing is not supported. total_protein_groups=$(total_protein_groups) exceeds max_protein_groups_in_memory_limit=$(max_protein_groups_in_memory_limit).")
    end

    all_protein_groups = DataFrame()
    for pg_path in passing_pg_paths
        if isfile(pg_path) && endswith(pg_path, ".arrow")
            append!(all_protein_groups, DataFrame(Tables.columntable(Arrow.Table(pg_path))))
        end
    end

    n_targets = sum(all_protein_groups.target)
    n_decoys = sum(.!all_protein_groups.target)
    skip_scoring = !(n_targets > 50 && n_decoys > 50 && nrow(all_protein_groups) > 1000)

    perform_probit_analysis_multifold(
        all_protein_groups,
        qc_folder,
        pg_refs,
        precursors;
        protein_to_cv_fold = protein_to_cv_fold,
        file_idx_to_name = file_idx_to_name,
        skip_scoring = skip_scoring,
        write_qc_plots = write_qc_plots,
        log_feature_importance = log_feature_importance,
        train_q_value_threshold = train_q_value_threshold,
        min_prefix_shape_neg_threshold_itr = min_prefix_shape_neg_threshold_itr,
        min_pep_neg_threshold_itr = min_pep_neg_threshold_itr
    )
end

function run_protein_scoring!(
    search_context::SearchContext;
    passing_refs::Vector{PSMFileReference},
    max_in_memory_table_mb::Float64,
    q_value_threshold::Float32,
    min_peptides::Int64,
    write_qc_plots::Bool,
    log_feature_importance::Bool,
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
    protein_to_possible_peptides = count_protein_peptides(precursors)

    file_idx_to_name = Dict{Int64, String}()
    for (file_idx, file_name) in enumerate(getFileIdToName(getMSData(search_context)))
        file_idx_to_name[Int64(file_idx)] = String(file_name)
    end

    pg_refs, psm_to_pg_mapping = build_protein_group_tables(
        passing_refs,
        passing_proteins_folder,
        protein_to_possible_peptides,
        min_peptides = min_peptides,
        q_value_threshold = q_value_threshold,
        qc_folder = qc_folder,
        file_idx_to_name = file_idx_to_name
    )

    valid_file_data = get_valid_file_paths(search_context, getPassingPsms)
    paired_files = PairedSearchFiles[]
    for (file_idx, psm_path) in valid_file_data
        haskey(psm_to_pg_mapping, psm_path) || continue
        pg_path = psm_to_pg_mapping[psm_path]
        push!(paired_files, PairedSearchFiles(psm_path, pg_path, file_idx))
        setPassingProteins!(getMSData(search_context), file_idx, pg_path)
    end

    isempty(paired_files) && error("No protein groups created during protein inference")

    psm_paths = [file_path(ref) for ref in passing_refs]
    protein_to_cv_fold = build_protein_cv_fold_mapping(psm_paths, precursors)

    max_in_memory_rows = estimate_max_rows(max_in_memory_table_mb, file_path(first(pg_refs)))
    @user_info "Memory budget $(max_in_memory_table_mb) MB → max_protein_groups = $max_in_memory_rows"
    perform_protein_probit_regression(
        pg_refs,
        max_in_memory_rows,
        qc_folder,
        precursors;
        protein_to_cv_fold = protein_to_cv_fold,
        file_idx_to_name = file_idx_to_name,
        write_qc_plots = write_qc_plots,
        log_feature_importance = log_feature_importance,
        train_q_value_threshold = q_value_threshold,
        min_pep_neg_threshold_itr = min_pep_neg_threshold_itr
    )

    sqrt_n_runs = floor(Int, sqrt(length(pg_refs)))
    n_proteins = length(getProteins(getSpecLib(search_context)))

    global_pg_score_dict, pg_name_to_global_pg_score =
        build_protein_global_score_dicts(pg_refs, sqrt_n_runs, n_proteins)
    search_context.pg_name_to_global_pg_score[] = pg_name_to_global_pg_score

    global_pg_qval_dict = build_protein_global_qval_dict(global_pg_score_dict)
    search_context.global_pg_score_to_qval_dict[] = global_pg_qval_dict

    sorted_pg_scores_path = joinpath(temp_folder, "sorted_pg_scores.arrow")
    spline_result = build_qvalue_spline_from_refs(pg_refs, :pg_score, sorted_pg_scores_path;
        batch_size = 1_000_000, compute_pep = true,
        min_pep_points_per_bin = q_value_interpolation_points_per_bin,
        temp_prefix = "pg_sidecar")
    search_context.pg_score_to_qval[] = spline_result.qval_spline
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

    update_psms_with_probit_scores_refs(
        paired_files,
        search_context.pg_name_to_global_pg_score[],
        search_context.pg_score_to_qval[],
        search_context.global_pg_score_to_qval_dict[]
    )

    return nothing
end
