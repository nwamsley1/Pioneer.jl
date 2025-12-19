function identification_metrics(
    precursors_long::DataFrame,
    precursors_wide::DataFrame,
    protein_groups_long::DataFrame,
    protein_groups_wide::DataFrame,
)
    precursor_metrics = Dict(
        "total" => nrow(precursors_long),
        "unique" => nrow(precursors_wide),
    )

    protein_metrics = Dict(
        "total" => nrow(protein_groups_long),
        "unique" => nrow(protein_groups_wide),
    )

    precursor_metrics, protein_metrics
end

function cv_metrics(
    precursors_wide::DataFrame,
    protein_groups_wide::DataFrame,
    quant_col_names,
    protein_quant_col_names,
    cv_groups::Dict{String, Vector{String}},
    dataset_name::AbstractString,
)
    @info "Starting CV metrics" dataset=dataset_name
    if isempty(cv_groups)
        @info "No CV groups found; computing across all runs" dataset=dataset_name
    else
        @info "Using CV groups" dataset=dataset_name groups=cv_groups
    end

    precursor_wide_metrics = compute_wide_metrics(
        precursors_wide,
        quant_col_names;
        table_label = "precursors_wide",
        dataset_name = dataset_name,
    )
    protein_wide_metrics = compute_wide_metrics(
        protein_groups_wide,
        protein_quant_col_names;
        table_label = "protein_groups_wide",
        dataset_name = dataset_name,
    )
    precursor_cv_metrics = compute_cv_metrics(
        precursors_wide, quant_col_names; table_label = "precursors_wide", groups = cv_groups
    )
    protein_cv_metrics = compute_cv_metrics(
        protein_groups_wide, protein_quant_col_names; table_label = "protein_groups_wide", groups = cv_groups
    )

    precursor_wide_metrics, protein_wide_metrics, precursor_cv_metrics, protein_cv_metrics
end

function keap1_metrics(
    precursors_wide::DataFrame,
    protein_groups_wide::DataFrame,
    quant_col_names,
    protein_quant_col_names,
    labels_for_runs::Dict{String, String},
    dataset_name::AbstractString,
)
    @info "Starting KEAP1 metrics" dataset=dataset_name
    precursor_metrics = gene_counts_metrics_by_run(
        precursors_wide,
        quant_col_names,
        labels_for_runs,
        "precursors",
        "KEAP1",
    )
    merge!(
        precursor_metrics,
        gene_counts_metrics_by_run(
            precursors_wide,
            quant_col_names,
            labels_for_runs,
            "precursors",
            "NFE2L2",
        ),
    )

    protein_metrics = gene_counts_metrics_by_run(
        protein_groups_wide,
        protein_quant_col_names,
        labels_for_runs,
        "protein_groups",
        "KEAP1",
    )
    merge!(
        protein_metrics,
        gene_counts_metrics_by_run(
            protein_groups_wide,
            protein_quant_col_names,
            labels_for_runs,
            "protein_groups",
            "NFE2L2",
        ),
    )

    precursor_metrics, protein_metrics
end

function fold_change_metrics_block(
    precursors_wide::DataFrame,
    protein_groups_wide::DataFrame,
    quant_col_names,
    protein_quant_col_names,
    design_entry,
    dataset_name::AbstractString,
)
    @info "Starting fold-change metrics" dataset=dataset_name
    if design_entry === nothing || isempty(design_entry.run_to_condition)
        @warn "No three-proteome design available; skipping fold-change metrics" dataset=dataset_name
        return nothing
    elseif isempty(design_entry.condition_pairs)
        @warn "Three-proteome design missing condition pairs; skipping fold-change metrics" dataset=dataset_name
        return nothing
    end

    precursor_fold_changes = fold_change_metrics_for_table(
        precursors_wide,
        quant_col_names,
        design_entry,
        design_entry.condition_pairs;
        table_label = "precursors",
    )

    protein_fold_changes = fold_change_metrics_for_table(
        protein_groups_wide,
        protein_quant_col_names,
        design_entry,
        design_entry.condition_pairs;
        table_label = "protein_groups",
    )

    if precursor_fold_changes === nothing && protein_fold_changes === nothing
        return nothing
    end

    metrics = Dict{String, Any}()
    precursor_fold_changes !== nothing && (metrics["precursors"] = precursor_fold_changes)
    protein_fold_changes !== nothing && (metrics["protein_groups"] = protein_fold_changes)
    metrics
end

function compute_dataset_metrics(
    dataset_dir::AbstractString,
    dataset_name::AbstractString;
    metric_groups::AbstractVector{<:AbstractString} = DEFAULT_METRIC_GROUPS,
    experimental_design::Dict{String, Any} = Dict{String, Any}(),
    three_proteome_designs = nothing,
    dataset_paths::Dict{String, String} = Dict{String, String}(),
)
    requested_groups = Set(replace.(lowercase.(metric_groups), "-" => "_"))
    need_identification = "identification" in requested_groups
    need_cv = "cv" in requested_groups
    need_keap1 = "keap1" in requested_groups
    need_ftr = "ftr" in requested_groups
    need_runtime = "runtime" in requested_groups
    need_three_proteome = ("fold_change" in requested_groups) || ("three_proteome" in requested_groups)
    need_table_metrics = need_identification || need_cv || need_keap1 || need_ftr || need_three_proteome

    dataset_three_proteome_designs = need_three_proteome ? three_proteome_designs : nothing

    precursor_id_metrics = nothing
    protein_id_metrics = nothing
    keap1_precursor_metrics = nothing
    keap1_protein_metrics = nothing
    ftr_metrics = nothing
    fold_change_metrics = nothing
    runtime_minutes = nothing

    if need_table_metrics
        required_files = if need_identification || need_cv || need_keap1 || need_ftr
            [
                "precursors_long.arrow",
                "precursors_wide.arrow",
                "protein_groups_long.arrow",
                "protein_groups_wide.arrow",
            ]
        else
            ["protein_groups_wide.arrow"]
        end

        missing_files = filter(f -> !isfile(joinpath(dataset_dir, f)), required_files)
        if !isempty(missing_files)
            @warn "Skipping dataset $dataset_name: missing required outputs" missing_files=missing_files
            return nothing
        end

        precursors_long = nothing
        precursors_wide = nothing
        protein_groups_long = nothing

        if need_identification || need_cv || need_keap1 || need_ftr
            precursors_long = read_required_table(joinpath(dataset_dir, "precursors_long.arrow"))
            precursors_wide = read_required_table(joinpath(dataset_dir, "precursors_wide.arrow"))
            protein_groups_long = read_required_table(joinpath(dataset_dir, "protein_groups_long.arrow"))
        end

        protein_groups_wide = read_required_table(joinpath(dataset_dir, "protein_groups_wide.arrow"))

        quant_col_names = if precursors_wide !== nothing
            quant_column_names_from_proteins(precursors_wide)
        else
            quant_column_names_from_proteins(protein_groups_wide)
        end
        protein_quant_col_names = quant_column_names_from_proteins(protein_groups_wide)
        precursor_wide_metrics = nothing
        protein_wide_metrics = nothing
        precursor_cv_metrics = nothing
        protein_cv_metrics = nothing

        if need_identification
            @info "Starting identification metrics" dataset=dataset_name
            precursor_id_metrics, protein_id_metrics = identification_metrics(
                precursors_long,
                precursors_wide,
                protein_groups_long,
                protein_groups_wide,
            )
        end

        if need_cv
            cv_groups = run_groups_for_dataset(
                experimental_design,
                dataset_name;
                three_proteome_designs = dataset_three_proteome_designs,
            )
            precursor_wide_metrics, protein_wide_metrics, precursor_cv_metrics, protein_cv_metrics =
                cv_metrics(
                    precursors_wide,
                    protein_groups_wide,
                    quant_col_names,
                    protein_quant_col_names,
                    cv_groups,
                    dataset_name,
                )
        end

        if need_keap1
            labels_for_runs = experimental_design_for_dataset(experimental_design, dataset_name)
            keap1_precursor_metrics, keap1_protein_metrics = keap1_metrics(
                precursors_wide,
                protein_groups_wide,
                quant_col_names,
                protein_quant_col_names,
                labels_for_runs,
                dataset_name,
            )
        end

        if need_ftr
            @info "Starting false transfer rate metrics" dataset=dataset_name
            ftr_metrics = compute_ftr_metrics(
                dataset_name,
                precursors_wide,
                protein_groups_wide,
                experimental_design,
                dataset_paths,
            )
        end

        if need_three_proteome
            design_entry = three_proteome_design_entry(dataset_three_proteome_designs, dataset_name)
            fold_change_metrics = fold_change_metrics_block(
                precursors_wide,
                protein_groups_wide,
                quant_col_names,
                protein_quant_col_names,
                design_entry,
                dataset_name,
            )
        end
    end

    entrapment_metrics = if "efdr" in requested_groups
        @info "Starting eFDR metrics" dataset=dataset_name
        compute_entrapment_metrics(dataset_dir, dataset_name)
    else
        nothing
    end

    if need_runtime
        report_path = joinpath(dataset_dir, "pioneer_search_report.txt")
        runtime_minutes = runtime_minutes_from_report(report_path)
    end

    metrics = Dict{String, Any}()

    identification_metrics_block = Dict{String, Any}()

    if need_identification
        precursors_identification = Dict{String, Any}()
        protein_identification = Dict{String, Any}()

        precursor_id_metrics !== nothing && merge!(precursors_identification, precursor_id_metrics)
        protein_id_metrics !== nothing && merge!(protein_identification, protein_id_metrics)

        !isempty(precursors_identification) && (identification_metrics_block["precursors"] = precursors_identification)
        !isempty(protein_identification) && (identification_metrics_block["protein_groups"] = protein_identification)
    end

    cv_metrics_block = Dict{String, Any}()
    if need_cv
        precursor_cv_block = Dict{String, Any}()
        protein_cv_block = Dict{String, Any}()

        if precursor_wide_metrics !== nothing
            precursor_cv_block["complete_rows"] = precursor_wide_metrics.complete_rows
            precursor_cv_block["data_completeness"] = precursor_wide_metrics.data_completeness
        end

        if protein_wide_metrics !== nothing
            protein_cv_block["complete_rows"] = protein_wide_metrics.complete_rows
            protein_cv_block["data_completeness"] = protein_wide_metrics.data_completeness
        end

        precursor_cv_metrics !== nothing && (precursor_cv_block["median_cv"] = precursor_cv_metrics.median_cv)
        protein_cv_metrics !== nothing && (protein_cv_block["median_cv"] = protein_cv_metrics.median_cv)

        !isempty(precursor_cv_block) && (cv_metrics_block["precursors"] = precursor_cv_block)
        !isempty(protein_cv_block) && (cv_metrics_block["protein_groups"] = protein_cv_block)
    end

    !isempty(identification_metrics_block) && (metrics["identification"] = identification_metrics_block)
    !isempty(cv_metrics_block) && (metrics["cv"] = cv_metrics_block)
    if need_keap1
        keap1_metrics_block = Dict{String, Any}()
        keap1_precursor_metrics !== nothing && !isempty(keap1_precursor_metrics) &&
            (keap1_metrics_block["precursors"] = keap1_precursor_metrics)
        keap1_protein_metrics !== nothing && !isempty(keap1_protein_metrics) &&
            (keap1_metrics_block["protein_groups"] = keap1_protein_metrics)
        !isempty(keap1_metrics_block) && (metrics["keap1"] = keap1_metrics_block)
    end
    if ftr_metrics !== nothing
        metrics["ftr"] = ftr_metrics
    end

    if fold_change_metrics !== nothing
        metrics["fold_change"] = fold_change_metrics
    end

    if entrapment_metrics !== nothing
        metrics["entrapment"] = entrapment_metrics
    end

    if runtime_minutes !== nothing
        metrics["runtime"] = runtime_minutes
    end

    return metrics
end
