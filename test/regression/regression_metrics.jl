#!/usr/bin/env julia

using CSV
using DataFrames
using JSON
using Pkg
using Statistics

const DEFAULT_METRIC_GROUPS = ["identification", "CV", "eFDR"]

function read_required_table(path::AbstractString)
    isfile(path) || error("Required file not found: $path")
    CSV.read(path, DataFrame; delim='\t', missingstring=["", "NA"], ignorerepeated=false)
end

is_numeric_column(col) = begin
    T = nonmissingtype(eltype(col))
    return T <: Number
end

function quant_column_names_from_proteins(df::DataFrame)
    col_names = names(df)
    anchor_idx = findfirst(name -> lowercase(String(name)) == "global_qval", col_names)

    if anchor_idx !== nothing
        quant_start = anchor_idx + 1
        quant_start > ncol(df) && return Symbol[]
        return Symbol.(col_names[quant_start:end])
    end

    numeric_flags = [is_numeric_column(df[:, c]) for c in col_names]
    any(numeric_flags) || return Symbol[]

    last_non_numeric = findlast(!, numeric_flags)
    if last_non_numeric === nothing
        return Symbol.(col_names)
    end

    quant_start = last_non_numeric + 1
    quant_start > ncol(df) && return Symbol[]
    Symbol.(col_names[quant_start:end])
end

function select_quant_columns(df::DataFrame, quant_col_names)
    quant_syms = Symbol.(quant_col_names)
    df_names = names(df)
    df_name_syms = Symbol.(df_names)
    df_names[[i for (i, n) in pairs(df_name_syms) if n in quant_syms]]
end

function compute_wide_metrics(
    df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}};
    table_label::AbstractString = "wide_table",
)
    existing_quant_cols = select_quant_columns(df, quant_col_names)
    runs = length(existing_quant_cols)
    if runs == 0
        return (; runs = 0, complete_rows = 0, data_completeness = 0.0)
    end

    if length(existing_quant_cols) < length(quant_col_names)
        missing_cols = setdiff(Symbol.(quant_col_names), Symbol.(existing_quant_cols))
        @warn "Missing quantification columns in dataset" missing_cols=missing_cols
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
)
    existing_quant_cols = select_quant_columns(df, quant_col_names)
    runs = length(existing_quant_cols)
    if runs == 0
        return (; runs = 0, rows_evaluated = 0, mean_cv = 0.0)
    end

    quant_data = df[:, existing_quant_cols]
    quant_matrix = Matrix(quant_data)
    complete_mask = [all(!ismissing, row) for row in eachrow(quant_matrix)]
    rows_evaluated = count(complete_mask)
    if rows_evaluated == 0
        return (; runs, rows_evaluated, mean_cv = 0.0)
    end

    quant_complete = quant_matrix[complete_mask, :]
    cvs = Float64[]
    for row in eachrow(quant_complete)
        mean_val = mean(row)
        if mean_val != 0
            push!(cvs, std(row) / mean_val)
        end
    end

    mean_cv = isempty(cvs) ? 0.0 : mean(cvs)
    (; runs, rows_evaluated, mean_cv)
end

function load_dataset_config(dataset_dir::AbstractString)
    config_path = joinpath(dataset_dir, "config.json")
    if !isfile(config_path)
        @warn "Skipping entrapment metrics; missing dataset config" dataset_dir=dataset_dir
        return nothing
    end

    try
        return JSON.parsefile(config_path)
    catch err
        @warn "Failed to parse dataset config; skipping entrapment metrics" config_path=config_path error=err
        return nothing
    end
end

function spectral_library_path_from_config(config::Dict, dataset_dir::AbstractString)
    paths_section = get(config, "paths", nothing)
    if paths_section isa AbstractDict
        library_entry = get(paths_section, "library", nothing)
        if library_entry isa AbstractDict
            lib_path = get(library_entry, "value", nothing)
            if lib_path !== nothing
                return isabspath(lib_path) ? lib_path : normpath(joinpath(dataset_dir, lib_path))
            end
        end
    end

    @warn "No spectral library path found in dataset config" dataset_dir=dataset_dir
    nothing
end

function load_entrapment_module(repo_path::AbstractString)
    if !isdir(repo_path)
        @warn "Entrapment analyses repository not found" repo_path=repo_path
        return nothing
    end

    project_path = Base.active_project()
    project_reset = project_path === nothing ? () -> Pkg.activate(; io=devnull) : () -> Pkg.activate(project_path; io=devnull)

    # Temporarily add the repository to LOAD_PATH to make the module discoverable and
    # activate its environment so dependencies are available without polluting the
    # main project environment.
    pushfirst!(Base.LOAD_PATH, repo_path)
    try
        Pkg.activate(repo_path; io=devnull)
        return Base.require(Symbol("EntrapmentAnalyses"))
    catch err
        @warn "Unable to load EntrapmentAnalyses module" repo_path=repo_path error=err
        return nothing
    finally
        deleteat!(Base.LOAD_PATH, findfirst(==(repo_path), Base.LOAD_PATH))
        project_reset()
    end
end

function try_run_both_analyses(run_fn, dataset_dir::AbstractString, lib_path::AbstractString, output_dir::AbstractString)
    attempts = [
        () -> run_fn(dataset_dir, lib_path; output_dir=output_dir),
        () -> run_fn(dataset_dir, lib_path, output_dir),
        () -> run_fn(dataset_dir; spectral_library_path=lib_path, output_dir=output_dir),
        () -> run_fn(dataset_dir, output_dir, lib_path),
        () -> run_fn(dataset_dir, lib_path),
    ]

    for attempt in attempts
        try
            return attempt()
        catch err
            err isa MethodError || err isa UndefKeywordError || begin
                @warn "Entrapment analysis run failed" error=err
                return nothing
            end
        end
    end

    @warn "Unable to call run_both_analyses with available signatures"
    nothing
end

function compute_entrapment_metrics(dataset_dir::AbstractString, dataset_name::AbstractString)
    config = load_dataset_config(dataset_dir)
    config === nothing && return nothing

    lib_path = spectral_library_path_from_config(config, dataset_dir)
    lib_path === nothing && return nothing

    repo_path = get(ENV, "ENTRAPMENT_ANALYSES_PATH", joinpath(@__DIR__, "..", "..", "EntrapmentAnalyses.jl"))
    entrapment_module = load_entrapment_module(repo_path)
    entrapment_module === nothing && return nothing

    if !isdefined(entrapment_module, :run_both_analyses)
        @warn "Entrapment analyses module missing run_both_analyses" repo_path=repo_path
        return nothing
    end

    output_dir = joinpath(dataset_dir, "entrapment_analyses")
    mkpath(output_dir)

    run_fn = getproperty(entrapment_module, :run_both_analyses)
    result = try_run_both_analyses(run_fn, dataset_dir, lib_path, output_dir)

    metrics = Dict(
        "dataset" => dataset_name,
        "spectral_library" => lib_path,
        "output_dir" => output_dir,
    )

    if result !== nothing
        try
            metrics["results"] = JSON.lower(result)
        catch err
            @warn "Unable to serialize entrapment analysis result" error=err
        end
    end

    metrics
end

function compute_dataset_metrics(
    dataset_dir::AbstractString,
    dataset_name::AbstractString;
    metric_groups::AbstractVector{<:AbstractString} = DEFAULT_METRIC_GROUPS,
)
    requested_groups = Set(lowercase.(metric_groups))
    required_files = [
        "precursors_long.tsv",
        "precursors_wide.tsv",
        "protein_groups_long.tsv",
        "protein_groups_wide.tsv",
    ]

    missing_files = filter(f -> !isfile(joinpath(dataset_dir, f)), required_files)
    if !isempty(missing_files)
        @warn "Skipping dataset $dataset_name: missing required outputs" missing_files=missing_files
        return nothing
    end

    precursors_long = read_required_table(joinpath(dataset_dir, "precursors_long.tsv"))
    precursors_wide = read_required_table(joinpath(dataset_dir, "precursors_wide.tsv"))
    protein_groups_long = read_required_table(joinpath(dataset_dir, "protein_groups_long.tsv"))
    protein_groups_wide = read_required_table(joinpath(dataset_dir, "protein_groups_wide.tsv"))

    quant_col_names = nothing
    precursor_wide_metrics = nothing
    protein_wide_metrics = nothing
    precursor_cv_metrics = nothing
    protein_cv_metrics = nothing

    if "cv" in requested_groups
        quant_col_names = quant_column_names_from_proteins(protein_groups_wide)
        precursor_wide_metrics = compute_wide_metrics(
            precursors_wide, quant_col_names; table_label = "precursors_wide"
        )
        protein_wide_metrics = compute_wide_metrics(
            protein_groups_wide, quant_col_names; table_label = "protein_groups_wide"
        )
        precursor_cv_metrics = compute_cv_metrics(
            precursors_wide, quant_col_names; table_label = "precursors_wide"
        )
        protein_cv_metrics = compute_cv_metrics(
            protein_groups_wide, quant_col_names; table_label = "protein_groups_wide"
        )
    end

    precursors_metrics = Dict{String, Any}(
        "total" => nrow(precursors_long),
        "unique" => nrow(precursors_wide),
    )

    protein_metrics = Dict{String, Any}(
        "total" => nrow(protein_groups_long),
        "unique" => nrow(protein_groups_wide),
    )

    entrapment_metrics = nothing

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
            "rows_with_complete_values" => precursor_cv_metrics.rows_evaluated,
            "mean_cv" => precursor_cv_metrics.mean_cv,
        ))
    end

    if protein_cv_metrics !== nothing
        merge!(protein_metrics, Dict(
            "cv_runs" => protein_cv_metrics.runs,
            "rows_with_complete_values" => protein_cv_metrics.rows_evaluated,
            "mean_cv" => protein_cv_metrics.mean_cv,
        ))
    end

    if "efdr" in requested_groups
        entrapment_metrics = compute_entrapment_metrics(dataset_dir, dataset_name)
    end

    # Placeholder for future metric groups (e.g., FTR, entrapment, fold-change, KEAP1)
    # that can reuse the requested_groups set to decide whether to run expensive calculations.

    metrics = Dict(
        "dataset" => dataset_name,
        "precursors" => precursors_metrics,
        "protein_groups" => protein_metrics,
    )

    if entrapment_metrics !== nothing
        metrics["entrapment"] = entrapment_metrics
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

    dataset_dirs = filter(dataset_dirs) do path
        dataset_name = basename(path)
        groups = normalize_metric_groups(get(metric_group_config, dataset_name, DEFAULT_METRIC_GROUPS))
        if isempty(groups)
            @info "Skipping dataset without requested metrics" dataset=dataset_name
            return false
        end
        return true
    end

    isempty(dataset_dirs) && error("No dataset directories remain after filtering")

    for dataset_dir in dataset_dirs
        dataset_name = basename(dataset_dir)

        metric_groups = normalize_metric_groups(
            get(metric_group_config, dataset_name, DEFAULT_METRIC_GROUPS),
        )
        metrics = compute_dataset_metrics(dataset_dir, dataset_name; metric_groups)
        metrics === nothing && continue

        output_path = joinpath(dataset_dir, "metrics_$(dataset_name).json")
        open(output_path, "w") do io
            JSON.print(io, metrics)
        end
    end
end

main()
