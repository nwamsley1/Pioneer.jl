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
    columns_to_load = Symbol[:protein_name, :target, :entrap_id, :n_peptides]
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

# ===== Two-round protein scoring: cross-run features + shadow proteins (GROUP_SCORING=1) =====

# Per (protein, run): global_max + global_top3_logodds of round-1 pg_score over the OTHER runs of the
# same protein (leave-one-out). In-memory — all_protein_groups holds every (protein, run) row.
function add_protein_cross_run_features!(df::DataFrame)
    n = nrow(df)
    gmax = zeros(Float32, n); gt3 = zeros(Float32, n)
    hasent = hasproperty(df, :entrap_id)
    groups = Dict{Tuple{String,Bool,UInt8}, Vector{Int}}()
    @inbounds for i in 1:n
        key = (String(df.protein_name[i]), Bool(df.target[i]), hasent ? UInt8(df.entrap_id[i]) : UInt8(0))
        push!(get!(() -> Int[], groups, key), i)
    end
    pg = df.pg_score
    for (_, idxs) in groups
        m = length(idxs)
        @inbounds for a in 1:m
            i = idxs[a]; mx = 0f0; t1 = 0f0; t2 = 0f0; t3 = 0f0
            for b in 1:m
                b == a && continue
                s = Float32(pg[idxs[b]])
                s > mx && (mx = s)
                if     s > t1; t3 = t2; t2 = t1; t1 = s
                elseif s > t2; t3 = t2; t2 = s
                elseif s > t3; t3 = s end
            end
            gmax[i] = mx
            gt3[i]  = _pos_logit(t1) + _pos_logit(t2) + _pos_logit(t3)
        end
    end
    df[!, :global_max] = gmax
    df[!, :global_top3_logodds] = gt3
    return df
end

# Write the two cross-run columns to the on-disk pg files. all_protein_groups preserves pg_ref row
# order (concatenated per file, skipping missing/empty), so slice the ordered vectors per file the
# same way. Needed so round-2's disk write-back (apply_probit_scores_multifold!) can read them.
function write_protein_cross_run_cols!(pg_refs, gmax::AbstractVector{Float32}, gt3::AbstractVector{Float32})
    off = 0
    for ref in pg_refs
        p = file_path(ref); isfile(p) || continue
        main = DataFrame(Tables.columntable(Arrow.Table(p)))
        nr = nrow(main); nr == 0 && continue
        main[!, :global_max] = collect(gmax[off+1:off+nr])
        main[!, :global_top3_logodds] = collect(gt3[off+1:off+nr])
        writeArrow(p, main)
        off += nr
    end
    return nothing
end

# One opposite-class shadow row per real row (TRAINING ONLY): copy a nearest-pg_score opposite-class
# row's features, then graft the real row's cross-run columns (1:1 marginal → the model can't use the
# cross-run features as a standalone lever). Returns a new DataFrame (real rows then shadows).
function inject_protein_shadows(df::DataFrame, graft_cols::Vector{Symbol})
    n = nrow(df); tgt = collect(Bool.(df.target)); pg = collect(Float32.(df.pg_score))
    tidx = findall(tgt); didx = findall(.!tgt)
    (isempty(tidx) || isempty(didx)) && return df
    tperm = sortperm(view(pg, tidx)); dperm = sortperm(view(pg, didx))
    tsorted = pg[tidx[tperm]]; dsorted = pg[didx[dperm]]
    src = Vector{Int}(undef, n)
    @inbounds for r in 1:n
        if tgt[r]
            pos = _nearest_sorted_idx(dsorted, pg[r]); src[r] = pos == 0 ? r : didx[dperm[pos]]
        else
            pos = _nearest_sorted_idx(tsorted, pg[r]); src[r] = pos == 0 ? r : tidx[tperm[pos]]
        end
    end
    shadows = df[src, :]                    # copy opposite-class rows (their label is the flip)
    for c in graft_cols
        hasproperty(shadows, c) && (shadows[!, c] = df[!, c])   # graft real row's cross-run values
    end
    return vcat(df, shadows)
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

    dump_requested = haskey(ENV, "PROTEIN_GROUP_DUMP")
    all_protein_groups = load_protein_probit_training_rows(
        pg_refs;
        include_qc_plot_columns = write_qc_plots || dump_requested
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

    # PROTEIN_GROUP_DUMP=<dir>: after round-1, write the protein-scoring state (standard features +
    # ROUND-1 pg_score + target + cv_fold + protein key + run/file + species) once, so round-2 /
    # cross-run / shadow / global-model scoring methods can be iterated OFFLINE without re-searching.
    if dump_requested && !skip_scoring
        dump_df = copy(all_protein_groups)
        assign_protein_group_cv_folds!(dump_df, protein_to_cv_fold)
        if file_idx_to_name !== nothing && hasproperty(dump_df, :file_idx)
            dump_df[!, :file_name] = String[get(file_idx_to_name, Int64(fi), "") for fi in dump_df.file_idx]
        end
        dump_dir = ENV["PROTEIN_GROUP_DUMP"]; mkpath(dump_dir)
        dump_path = joinpath(dump_dir, "protein_group_dump.arrow")
        writeArrow(dump_path, dump_df)
        @user_info "protein-group dump written (post round-1): $dump_path ($(nrow(dump_df)) rows, $(ncol(dump_df)) cols)"
    end

    # ROUND 2 (GROUP_SCORING=1): add cross-run protein features from round-1 pg_score, inject shadow
    # proteins for regularization, and re-run the same probit with the two extra features.
    if group_scoring_enabled() && !skip_scoring
        add_protein_cross_run_features!(all_protein_groups)
        write_protein_cross_run_cols!(pg_refs, all_protein_groups.global_max, all_protein_groups.global_top3_logodds)
        df_train = get(ENV, "PROTEIN_SHADOW", "1") != "0" ?
            inject_protein_shadows(all_protein_groups, [:global_max, :global_top3_logodds]) :
            all_protein_groups
        @user_info "protein two-round: round-2 with cross-run features (global_max, global_top3_logodds) + $(nrow(df_train) - nrow(all_protein_groups)) shadow proteins"
        perform_probit_analysis_multifold(
            df_train,
            qc_folder,
            pg_refs,
            precursors;
            protein_to_cv_fold = protein_to_cv_fold,
            file_idx_to_name = file_idx_to_name,
            skip_scoring = skip_scoring,
            write_qc_plots = false,
            log_feature_importance = log_feature_importance,
            train_q_value_threshold = train_q_value_threshold,
            min_prefix_shape_neg_threshold_itr = min_prefix_shape_neg_threshold_itr,
            min_pep_neg_threshold_itr = min_pep_neg_threshold_itr
        )
    end
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

    # GLOBAL_MODEL=1: learned protein-level model on per-run pg_score aggregates (2-fold protein CV),
    # replacing the fixed top-√N log-odds (mirrors the precursor global model).
    global_pg_score_dict, pg_name_to_global_pg_score =
        global_model_enabled() ?
            build_protein_global_model_dict(pg_refs, sqrt_n_runs, n_proteins, protein_to_cv_fold) :
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
