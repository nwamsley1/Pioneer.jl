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

    # Each precursor's sequence must be added to its protein(s)' peptide set, so rows
    # can't be skipped wholesale — but the (accession_string, target, entrap_id) tuple
    # repeats across the ~6.8M precursors (only ~tens of thousands are distinct). Cache
    # the destination Set vector per distinct tuple: a cache hit skips the split, the
    # per-token String() interning, and the NamedTuple-key dict lookups, leaving only
    # the unavoidable push! of the sequence. Bit-identical (same sets, same members).
    set_cache = Dict{Tuple{String, Bool, UInt8}, Vector{Set{String}}}()
    for i in 1:n_precursors
        is_decoy = all_decoys[i]
        entrap_id = all_entrap_ids[i]
        target = !is_decoy

        target_sets = get!(set_cache, (all_accession_numbers[i], target, entrap_id)) do
            sets = Vector{Set{String}}()
            for protein_name in split(all_accession_numbers[i], ';')
                key = (protein_name = String(protein_name), target = target, entrap_id = entrap_id)
                push!(sets, get!(() -> Set{String}(), protein_to_possible_peptides, key))
            end
            sets
        end

        seq = all_sequences[i]
        for s in target_sets
            push!(s, seq)
        end
    end

    return protein_to_possible_peptides
end

"""
    load_protein_probit_training_rows(pg_refs; include_qc_plot_columns = false)

Load only the columns needed for protein probit fitting from protein-group files.
This avoids materializing string-heavy output columns like `peptide_list` that
are not used by the model.
"""
function load_protein_probit_training_rows(
    pg_refs::Vector{ProteinGroupFileReference};
    include_qc_plot_columns::Bool = false
)
    columns_to_load = Symbol[:protein_name, :target, :n_peptides]
    append!(columns_to_load, protein_probit_feature_names())
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
    perform_protein_probit_regression(pg_refs::Vector{ProteinGroupFileReference},
                                    max_in_memory_rows::Int64,
                                    qc_folder::String,
                                    precursors::LibraryPrecursors;
                                    protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}})

Perform probit regression on protein groups.

# Arguments
- `pg_refs`: Vector of protein group file references
- `max_in_memory_rows`: Maximum number of protein-group rows allowed in memory
- `qc_folder`: Folder for QC plots
- `precursors`: Library precursors
- `protein_to_cv_fold`: Pre-built mapping of proteins to CV folds
"""
function perform_protein_probit_regression(
    pg_refs::Vector{ProteinGroupFileReference},
    max_in_memory_rows::Int64,
    qc_folder::String,
    precursors::LibraryPrecursors;
    protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
    file_idx_to_name::Union{Nothing, AbstractDict{Int64, String}} = nothing,
    write_qc_plots::Bool = true,
    log_feature_importance::Bool = true,
    train_q_value_threshold::Float32 = 0.01f0,
    min_prefix_shape_neg_threshold_itr::Float32 = -0.20f0,
    min_pep_neg_threshold_itr::Float32 = 0.90f0
)
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

    all_protein_groups = load_protein_probit_training_rows(
        pg_refs;
        include_qc_plot_columns = write_qc_plots
    )

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

    pg_refs, psm_to_pg_mapping, protein_to_cv_fold = build_protein_group_tables(
        passing_refs,
        passing_proteins_folder,
        protein_to_possible_peptides,
        precursors = precursors,
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
    # PEP is computed once, AFTER the protein filter (recalc below), on the survivor
    # population — so :pg_pep is conditioned on the same set as :pg_qval. No pre-filter
    # PEP is built here (nothing consumes it before the recalc).
    spline_result = build_qvalue_spline_from_refs(pg_refs, :pg_score, sorted_pg_scores_path;
        batch_size = 1_000_000,
        min_pep_points_per_bin = q_value_interpolation_points_per_bin,
        temp_prefix = "pg_sidecar")
    search_context.pg_score_to_qval[] = spline_result.qval_spline

    # SEQUENTIAL FILTER (proteins): Part 1 filters on (:global_pg_qval ≤ threshold)
    # ONLY, mirroring the PSM V2 filter. The experiment-wide (:pg_qval) filter is
    # deferred to the recalc below, where it is applied to the pg_qval recomputed on
    # the global-passing survivors. (Previously this filtered on both simultaneously.)
    protein_combined_pipeline = TransformPipeline() |>
        add_dict_column_composite_key(:global_pg_score, [:protein_name, :target, :entrap_id], global_pg_score_dict) |>
        add_dict_column_composite_key(:global_pg_qval, [:protein_name, :target, :entrap_id], global_pg_qval_dict) |>
        add_interpolated_column(:pg_qval, :pg_score, search_context.pg_score_to_qval[]) |>
        filter_by_multiple_thresholds([
            (:global_pg_qval, q_value_threshold)
        ])

    apply_pipeline!(pg_refs, protein_combined_pipeline)

    spline_result = build_qvalue_spline_from_refs(pg_refs, :pg_score, sorted_pg_scores_path;
        batch_size = 1_000_000, compute_pep = true,
        min_pep_points_per_bin = q_value_interpolation_points_per_bin,
        temp_prefix = "pg_recalc")
    search_context.pg_score_to_qval[] = spline_result.qval_spline
    search_context.pg_score_to_pep[]  = spline_result.pep_interp

    # SEQUENTIAL FILTER (proteins): recompute experiment-wide pg_qval AND pg_pep on the
    # global-passing survivors, THEN apply the deferred (:pg_qval ≤ threshold) filter.
    recalc_pg_pipeline = TransformPipeline() |>
        add_interpolated_column(:pg_qval, :pg_score, search_context.pg_score_to_qval[]) |>
        add_interpolated_column(:pg_pep,  :pg_score, search_context.pg_score_to_pep[]) |>
        filter_by_multiple_thresholds([(:pg_qval, q_value_threshold)])

    pg_refs = apply_pipeline_batch(pg_refs, recalc_pg_pipeline, passing_proteins_folder)

    update_psms_with_probit_scores_refs(
        paired_files,
        search_context.pg_name_to_global_pg_score[],
        search_context.pg_score_to_qval[],
        search_context.global_pg_score_to_qval_dict[]
    )

    return nothing
end
