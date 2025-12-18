#!/usr/bin/env julia

using Arrow
using DataFrames
using JSON
using Statistics

include("metrics_helpers.jl")
include("entrapment_metrics.jl")
include("three_proteome_metrics.jl")
using .RegressionMetricsHelpers: condition_columns, count_nonyeast_ids, count_species_ids, count_total_ids, count_yeast_ids, gene_names_column, mean_for_columns, quant_column_names_from_proteins, select_quant_columns, species_column, unique_species_value
using .EntrapmentMetrics: compute_entrapment_metrics
using .ThreeProteomeMetrics: experimental_design_for_dataset, fold_change_metrics_for_table, gene_counts_metrics_by_run, load_experimental_design, load_three_proteome_designs, normalize_metric_label, run_groups_for_dataset, three_proteome_design_entry

const DEFAULT_METRIC_GROUPS = ["identification", "CV", "eFDR", "runtime"]

function read_required_table(path::AbstractString)
    isfile(path) || error("Required file not found: $path")
    DataFrame(Arrow.Table(path))
end

function runtime_minutes_from_report(path::AbstractString)
    if !isfile(path)
        @warn "Runtime report not found; skipping runtime metric" path=path
        return nothing
    end

    for line in eachline(path)
        if startswith(line, "Total Runtime:")
            m = match(r"^Total Runtime:\s*([0-9]+(?:\.[0-9]+)?)\s+minutes", line)
            if m !== nothing
                return parse(Float64, m.captures[1])
            end
        end
    end

    @warn "Total runtime line not found in report" path=path
    nothing
end

function compute_wide_metrics(
    df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}};
    table_label::AbstractString = "wide_table",
    dataset_name::AbstractString = "dataset",
)
    existing_quant_cols = select_quant_columns(df, quant_col_names)
    runs = length(existing_quant_cols)
    available_cols = Symbol.(names(df))

    if runs == 0
        expected = Symbol.(quant_col_names)
        @warn "No quantification columns found in dataset" dataset=dataset_name table_label=table_label expected_quant_cols=expected available_cols=available_cols expected_quant_cols_list=join(expected, ", ") available_cols_list=join(available_cols, ", ")
        return (; runs = 0, complete_rows = 0, data_completeness = 0.0)
    end

    if length(existing_quant_cols) < length(quant_col_names)
        missing_cols = setdiff(Symbol.(quant_col_names), Symbol.(existing_quant_cols))
        expected = Symbol.(quant_col_names)
        available_quant_syms = Symbol.(existing_quant_cols)
        @warn "Missing quantification columns in dataset" dataset=dataset_name table_label=table_label missing_cols=missing_cols expected_quant_cols=expected available_quant_cols=available_quant_syms available_cols=available_cols missing_cols_list=join(missing_cols, ", ") expected_quant_cols_list=join(expected, ", ") available_quant_cols_list=join(available_quant_syms, ", ") available_cols_list=join(available_cols, ", ")
    end

    quant_data = df[:, existing_quant_cols]
    quant_matrix = Matrix(quant_data)

    complete_rows = sum(all(!ismissing, row) for row in eachrow(quant_matrix))
    non_missing_values = count(!ismissing, quant_matrix)
    total_cells = nrow(df) * runs

    (; runs, complete_rows, data_completeness = total_cells > 0 ? non_missing_values / total_cells : 0.0)
end

function compute_cv_metrics(
    df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}};
    table_label::AbstractString = "wide_table",
    groups::Dict{String, Vector{String}} = Dict{String, Vector{String}}(),
)
    function cvs_for_columns(
        columns::AbstractVector{<:Union{Symbol, String}},
        label::AbstractString,
    )
        column_syms = Symbol.(columns)
        runs = length(columns)
        if runs == 0
            return (; runs = 0, rows_evaluated = 0, cvs = Float64[])
        end

        quant_data = df[:, column_syms]
        complete_data = dropmissing(quant_data)
        rows_evaluated = nrow(complete_data)

        first_row_values = rows_evaluated == 0 ? NamedTuple() : NamedTuple(complete_data[1, :])
        @info "Computing CVs for group" table_label=table_label group_label=label quant_columns=column_syms rows_evaluated=rows_evaluated first_row_values=first_row_values
        if rows_evaluated == 0
            return (; runs, rows_evaluated, cvs = Float64[])
        end

        quant_complete = Matrix(complete_data)
        computed_cvs = Float64[]
        for row in eachrow(quant_complete)
            mean_val = mean(row)
            if mean_val != 0
                push!(computed_cvs, std(row) / mean_val)
            end
        end

        (; runs, rows_evaluated, cvs = computed_cvs)
    end

    if isempty(groups)
        existing_quant_cols = select_quant_columns(df, quant_col_names)
        stats = cvs_for_columns(existing_quant_cols, "all_runs")
        median_cv = isempty(stats.cvs) ? 0.0 : median(stats.cvs)
        return (; runs = stats.runs, rows_evaluated = stats.rows_evaluated, median_cv)
    end

    all_runs = Set{Symbol}()
    all_cvs = Float64[]
    total_rows = 0
    for (label, runs) in pairs(groups)
        columns = select_quant_columns(df, runs)
        union!(all_runs, Symbol.(columns))
        stats = cvs_for_columns(columns, label)
        total_rows += stats.rows_evaluated
        append!(all_cvs, stats.cvs)
    end

    median_cv = isempty(all_cvs) ? 0.0 : median(all_cvs)
    (; runs = length(all_runs), rows_evaluated = total_rows, median_cv)
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

    precursors_metrics = nothing
    protein_metrics = nothing
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

        precursors_metrics = Dict{String, Any}()
        protein_metrics = Dict{String, Any}()

        if need_identification
            @info "Starting identification metrics" dataset=dataset_name
            merge!(precursors_metrics, Dict(
                "total" => nrow(precursors_long),
                "unique" => nrow(precursors_wide),
            ))

            merge!(protein_metrics, Dict(
                "total" => nrow(protein_groups_long),
                "unique" => nrow(protein_groups_wide),
            ))
        end

        if need_cv
            @info "Starting CV metrics" dataset=dataset_name
            cv_groups = run_groups_for_dataset(
                experimental_design,
                dataset_name;
                three_proteome_designs = three_proteome_designs,
            )
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
        end

        if precursor_wide_metrics !== nothing
            merge!(precursors_metrics, Dict(
                "runs" => precursor_wide_metrics.runs,
                "complete_rows" => precursor_wide_metrics.complete_rows,
                "data_completeness" => precursor_wide_metrics.data_completeness,
                ))
        end

        if protein_wide_metrics !== nothing
            merge!(protein_metrics, Dict(
                "runs" => protein_wide_metrics.runs,
                "complete_rows" => protein_wide_metrics.complete_rows,
                "data_completeness" => protein_wide_metrics.data_completeness,
                ))
        end

        if precursor_cv_metrics !== nothing
            merge!(precursors_metrics, Dict(
                "cv_runs" => precursor_cv_metrics.runs,
                "median_cv" => precursor_cv_metrics.median_cv,
            ))
        end

        if protein_cv_metrics !== nothing
            merge!(protein_metrics, Dict(
                "cv_runs" => protein_cv_metrics.runs,
                "median_cv" => protein_cv_metrics.median_cv,
            ))
        end

        if need_keap1
            @info "Starting KEAP1 metrics" dataset=dataset_name
            labels_for_runs = experimental_design_for_dataset(experimental_design, dataset_name)
            keap1_precursor_metrics = gene_counts_metrics_by_run(
                precursors_wide,
                quant_col_names,
                labels_for_runs,
                "precursors",
                "KEAP1",
            )
            merge!(
                keap1_precursor_metrics,
                gene_counts_metrics_by_run(
                    precursors_wide,
                    quant_col_names,
                    labels_for_runs,
                    "precursors",
                    "NFE2L2",
                ),
            )

            keap1_protein_metrics = gene_counts_metrics_by_run(
                protein_groups_wide,
                protein_quant_col_names,
                labels_for_runs,
                "protein_groups",
                "KEAP1",
            )
            merge!(
                keap1_protein_metrics,
                gene_counts_metrics_by_run(
                    protein_groups_wide,
                    protein_quant_col_names,
                    labels_for_runs,
                    "protein_groups",
                    "NFE2L2",
                ),
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
            @info "Starting fold-change metrics" dataset=dataset_name
            design_entry = three_proteome_design_entry(three_proteome_designs, dataset_name)
            if design_entry === nothing || isempty(design_entry.run_to_condition)
                @warn "No three-proteome design available; skipping fold-change metrics" dataset=dataset_name
            elseif isempty(design_entry.condition_pairs)
                @warn "Three-proteome design missing condition pairs; skipping fold-change metrics" dataset=dataset_name
            else
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

                if precursor_fold_changes !== nothing || protein_fold_changes !== nothing
                    fold_change_metrics = Dict{String, Any}()
                    precursor_fold_changes !== nothing && (fold_change_metrics["precursors"] = precursor_fold_changes)
                    protein_fold_changes !== nothing && (fold_change_metrics["protein_groups"] = protein_fold_changes)
                end
            end
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

    if precursors_metrics !== nothing
        if keap1_precursor_metrics !== nothing
            merge!(precursors_metrics, keap1_precursor_metrics)
        end

        if !isempty(precursors_metrics)
            metrics["precursors"] = precursors_metrics
        end
    end

    if protein_metrics !== nothing
        if keap1_protein_metrics !== nothing
            merge!(protein_metrics, keap1_protein_metrics)
        end

        if !isempty(protein_metrics)
            metrics["protein_groups"] = protein_metrics
        end
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

function normalize_metric_groups(groups)
    if groups isa AbstractVector
        return [String(g) for g in groups]
    end

    @warn "Invalid metric groups entry; falling back to empty list" groups_type=typeof(groups)
    String[]
end

function metric_preferences(config::Dict, dataset_name::AbstractString)
    entry = get(config, dataset_name, DEFAULT_METRIC_GROUPS)

    if entry isa AbstractDict
        groups = normalize_metric_groups(get(entry, "groups", DEFAULT_METRIC_GROUPS))
        return (; groups)
    end

    groups = normalize_metric_groups(entry)
    (; groups)
end


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

    groups = run_groups_for_dataset(experimental_design, dataset_name)
    alt_groups = run_groups_for_dataset(experimental_design, mbr_name == dataset_name ? nombr_name : mbr_name)
    if isempty(groups)
        groups = alt_groups
    elseif isempty(alt_groups)
        alt_groups = groups
    end

    human_only_runs = get(groups, "human_only", String[])
    
    if isempty(human_only_runs)
        @warn "No human-only runs provided for FTR metrics; skipping" dataset=dataset_name
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

function dataset_dirs_from_root(root::AbstractString)
    !isdir(root) && return String[]
    filter(readdir(root; join=true)) do path
        isdir(path)
    end
end

function resolve_results_root()
    env_root = get(ENV, "PIONEER_RESULTS_DIR", "")
    arg_root = length(ARGS) >= 1 ? ARGS[1] : ""

    if !isempty(env_root)
        return env_root
    elseif !isempty(arg_root)
        return arg_root
    end

    error("Results directory is not specified; set PIONEER_RESULTS_DIR or pass a path argument")
end

function main()
    results_dir = resolve_results_root()
    isdir(results_dir) || error("Results directory does not exist: $results_dir")
    dataset_dirs = dataset_dirs_from_root(results_dir)
    isempty(dataset_dirs) && error("No dataset directories found in $results_dir")

    dataset_paths = Dict{String, String}(basename(path) => path for path in dataset_dirs)

    @info "Using regression results directory" results_dir=results_dir

    metrics_config_path = get(
        ENV,
        "PIONEER_METRICS_CONFIG",
        joinpath(@__DIR__, "..", "..", "pioneer-regression-configs", "metrics_config.json"),
    )
    metric_group_config = if isfile(metrics_config_path)
        JSON.parsefile(metrics_config_path)
    else
        @info "No metrics_config.json found; using default metric groups" metrics_config_path=metrics_config_path
        Dict{String, Any}()
    end

    experimental_design_path = get(
        ENV,
        "PIONEER_EXPERIMENTAL_DESIGN",
        joinpath(@__DIR__, "..", "..", "pioneer-regression-configs", "experimental_designs"),
    )
    experimental_design = load_experimental_design(experimental_design_path)

    three_proteome_designs_path = get(
        ENV,
        "PIONEER_THREE_PROTEOME_DESIGNS",
        joinpath(@__DIR__, "..", "..", "pioneer-regression-configs", "experimental_designs"),
    )
    three_proteome_designs = nothing

    dataset_dirs = filter(dataset_dirs) do path
        dataset_name = basename(path)
        preferences = metric_preferences(metric_group_config, dataset_name)
        if isempty(preferences.groups)
            @info "Skipping dataset without requested metrics" dataset=dataset_name
            return false
        end
        return true
    end

    isempty(dataset_dirs) && error("No dataset directories remain after filtering")

    for dataset_dir in dataset_dirs
        dataset_name = basename(dataset_dir)

        preferences = metric_preferences(metric_group_config, dataset_name)

        metric_groups = preferences.groups
        normalized_groups = Set(replace.(lowercase.(metric_groups), "-" => "_"))
        need_three_proteome = ("fold_change" in normalized_groups) || ("three_proteome" in normalized_groups)
        if need_three_proteome && three_proteome_designs === nothing
            three_proteome_designs = load_three_proteome_designs(three_proteome_designs_path)
        end

        metrics = compute_dataset_metrics(
            dataset_dir,
            dataset_name;
            metric_groups = metric_groups,
            experimental_design = experimental_design,
            three_proteome_designs = three_proteome_designs,
            dataset_paths = dataset_paths,
        )
        metrics === nothing && continue

        output_path = joinpath(dataset_dir, "metrics_$(dataset_name).json")
        open(output_path, "w") do io
            JSON.print(io, metrics)
        end
    end
end

main()
