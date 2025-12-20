#!/usr/bin/env julia

using Arrow
using DataFrames
using JSON
using Statistics

include("metrics_helpers.jl")
include("entrapment_metrics.jl")
include("three_proteome_metrics.jl")
include("ftr_metrics.jl")
include("metrics_pipeline.jl")
using .RegressionMetricsHelpers: condition_columns, count_nonyeast_ids, count_species_ids, count_total_ids, count_yeast_ids, gene_names_column, mean_for_columns, quant_column_names_from_proteins, select_quant_columns, species_column, unique_species_value
using .EntrapmentMetrics: compute_entrapment_metrics
using .ThreeProteomeMetrics: experimental_design_entry, experimental_design_for_dataset, fold_change_metrics_for_table, gene_counts_metrics_by_run, load_experimental_design, load_three_proteome_designs, normalize_metric_label, run_groups_for_dataset, three_proteome_design_entry

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
