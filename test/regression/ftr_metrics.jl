function ftr_metrics_for_table(
    mbr_df::DataFrame,
    no_mbr_df::DataFrame,
    mbr_quant_cols::AbstractVector{<:Union{Symbol, String}},
    no_mbr_quant_cols::AbstractVector{<:Union{Symbol, String}},
    human_only_runs::AbstractVector{<:AbstractString};
    table_label::AbstractString,
)
    yeast_human_only_mbr = count_yeast_ids(mbr_df, mbr_quant_cols, human_only_runs; table_label = table_label)
    yeast_human_only_no_mbr = count_yeast_ids(no_mbr_df, no_mbr_quant_cols, human_only_runs; table_label = table_label)

    total_ids_human_only_mbr = count_total_ids(mbr_df, mbr_quant_cols, human_only_runs; table_label = table_label)
    total_ids_human_only_no_mbr = count_total_ids(no_mbr_df, no_mbr_quant_cols, human_only_runs; table_label = table_label)

    additional_yeast_in_human_only = max(yeast_human_only_mbr - yeast_human_only_no_mbr, 0)
    additional_ids_in_human_only = max(total_ids_human_only_mbr - total_ids_human_only_no_mbr, 0)
    ftr = additional_ids_in_human_only > 0 ? additional_yeast_in_human_only / additional_ids_in_human_only : 0.0

    return Dict(
        "yeast_ids_human_only_no_mbr" => yeast_human_only_no_mbr,
        "yeast_ids_human_only_mbr" => yeast_human_only_mbr,
        "total_ids_human_only_no_mbr" => total_ids_human_only_no_mbr,
        "total_ids_human_only_mbr" => total_ids_human_only_mbr,
        "additional_yeast_ids_in_human_only" => additional_yeast_in_human_only,
        "additional_ids_in_human_only" => additional_ids_in_human_only,
        "false_transfer_rate" => ftr,
    )
end

function paired_mbr_dataset_paths(
    dataset_name::AbstractString,
    dataset_paths::Dict{String, String},
)
    if occursin("_noMBR_", dataset_name)
        mbr_name = replace(dataset_name, "_noMBR_" => "_MBR_")
        mbr_path = get(dataset_paths, mbr_name, nothing)
        nombr_path = get(dataset_paths, dataset_name, nothing)
        (mbr_path === nothing || nombr_path === nothing) && return nothing
        return (mbr_name, mbr_path, dataset_name, nombr_path)
    elseif occursin("_MBR_", dataset_name)
        nombr_name = replace(dataset_name, "_MBR_" => "_noMBR_")
        mbr_path = get(dataset_paths, dataset_name, nothing)
        nombr_path = get(dataset_paths, nombr_name, nothing)
        (mbr_path === nothing || nombr_path === nothing) && return nothing
        return (dataset_name, mbr_path, nombr_name, nombr_path)
    end

    nothing
end

function human_only_condition_from_design(entry, alt_entry)
    for source in (entry, alt_entry)
        if source isa AbstractDict
            ftr_config = get(source, "false_transfer_rate", nothing)
            if ftr_config isa AbstractDict
                condition = get(ftr_config, "human_only_condition", nothing)
                condition !== nothing && return String(condition)
            end

            condition = get(source, "human_only_condition", nothing)
            condition !== nothing && return String(condition)
        end
    end

    return "human_only"
end

function compute_ftr_metrics(
    dataset_name::AbstractString,
    precursors_wide::DataFrame,
    protein_groups_wide::DataFrame,
    experimental_design::Dict{String, Any},
    dataset_paths::Dict{String, String},
)
    paired_paths = paired_mbr_dataset_paths(dataset_name, dataset_paths)
    paired_paths === nothing && begin
        @warn "FTR metrics require paired MBR and noMBR datasets" dataset=dataset_name
        return nothing
    end

    mbr_name, mbr_path, nombr_name, nombr_path = paired_paths

    precursors_mbr = if dataset_name == mbr_name
        precursors_wide
    else
        arrow_path = joinpath(mbr_path, "precursors_wide.arrow")
        isfile(arrow_path) || begin
            @warn "Missing precursors for MBR dataset; skipping FTR metrics" dataset=mbr_name path=arrow_path
            return nothing
        end
        read_required_table(arrow_path)
    end

    precursors_no_mbr = if dataset_name == nombr_name
        precursors_wide
    else
        arrow_path = joinpath(nombr_path, "precursors_wide.arrow")
        isfile(arrow_path) || begin
            @warn "Missing precursors for noMBR dataset; skipping FTR metrics" dataset=nombr_name path=arrow_path
            return nothing
        end
        read_required_table(arrow_path)
    end

    protein_groups_mbr = if dataset_name == mbr_name
        protein_groups_wide
    else
        arrow_path = joinpath(mbr_path, "protein_groups_wide.arrow")
        isfile(arrow_path) || begin
            @warn "Missing protein groups for MBR dataset; skipping FTR metrics" dataset=mbr_name path=arrow_path
            return nothing
        end
        read_required_table(arrow_path)
    end

    protein_groups_no_mbr = if dataset_name == nombr_name
        protein_groups_wide
    else
        arrow_path = joinpath(nombr_path, "protein_groups_wide.arrow")
        isfile(arrow_path) || begin
            @warn "Missing protein groups for noMBR dataset; skipping FTR metrics" dataset=nombr_name path=arrow_path
            return nothing
        end
        read_required_table(arrow_path)
    end

    mbr_quant_cols = quant_column_names_from_proteins(precursors_mbr)
    nombr_quant_cols = quant_column_names_from_proteins(precursors_no_mbr)

    protein_mbr_quant_cols = quant_column_names_from_proteins(protein_groups_mbr)
    protein_nombr_quant_cols = quant_column_names_from_proteins(protein_groups_no_mbr)

    alt_dataset_name = mbr_name == dataset_name ? nombr_name : mbr_name
    groups = run_groups_for_dataset(experimental_design, dataset_name)
    alt_groups = run_groups_for_dataset(experimental_design, alt_dataset_name)
    entry = experimental_design_entry(experimental_design, dataset_name)
    alt_entry = experimental_design_entry(experimental_design, alt_dataset_name)
    human_only_condition = human_only_condition_from_design(entry, alt_entry)

    primary_groups = isempty(groups) ? alt_groups : groups
    secondary_groups = isempty(groups) ? groups : alt_groups

    human_only_runs = get(primary_groups, human_only_condition, String[])
    isempty(human_only_runs) && (human_only_runs = get(secondary_groups, human_only_condition, String[]))

    if isempty(human_only_runs) && human_only_condition != "human_only"
        human_only_runs = get(primary_groups, "human_only", String[])
        isempty(human_only_runs) && (human_only_runs = get(secondary_groups, "human_only", String[]))
    end

    if isempty(human_only_runs)
        @warn "No human-only runs provided for FTR metrics; skipping" dataset=dataset_name expected_condition=human_only_condition
        return nothing
    end

    precursor_metrics = ftr_metrics_for_table(
        precursors_mbr,
        precursors_no_mbr,
        mbr_quant_cols,
        nombr_quant_cols,
        human_only_runs;
        table_label = "precursors",
    )

    protein_metrics = ftr_metrics_for_table(
        protein_groups_mbr,
        protein_groups_no_mbr,
        protein_mbr_quant_cols,
        protein_nombr_quant_cols,
        human_only_runs;
        table_label = "protein_groups",
    )

    return Dict(
        "precursors" => precursor_metrics,
        "protein_groups" => protein_metrics,
        "mbr_dataset" => mbr_name,
        "nombr_dataset" => nombr_name,
    )
end
