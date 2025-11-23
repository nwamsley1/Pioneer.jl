#!/usr/bin/env julia

using Arrow
using CSV
using DataFrames
using JSON
using Pkg
using Statistics
using TOML

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
        lib_path = get(paths_section, "library", nothing)
        if lib_path isa AbstractString
            return isabspath(lib_path) ? lib_path : normpath(joinpath(dataset_dir, lib_path))
        end
    end

    @warn "No spectral library path found in dataset config" dataset_dir=dataset_dir
    nothing
end

function arrow_column_names(path::AbstractString)
    try
        table = Arrow.Table(path; convert=false)
        return Set(Symbol.(propertynames(table)))
    catch err
        @warn "Failed to read Arrow schema for score selection" path=path error=err
        return nothing
    end
end

function precursor_score_pairs(path::AbstractString)
    cols = arrow_column_names(path)
    cols === nothing && return nothing

    required_pairs = [
        (:MBR_boosted_global_prob, :MBR_boosted_global_qval),
        (:MBR_boosted_prec_prob, :MBR_boosted_qval),
    ]

    available_pairs = [pair for pair in required_pairs if all(col -> col in cols, pair)]
    if isempty(available_pairs)
        missing_columns_by_pair = Dict(pair => [col for col in pair if col ∉ cols] for pair in required_pairs)
        @warn "No compatible precursor score/q-value columns found for entrapment analysis" path=path missing_columns_by_pair=missing_columns_by_pair available_columns=join(string.(collect(cols)), ", ")
        return nothing
    end

    available_pairs
end

function protein_score_pairs(path::AbstractString)
    cols = arrow_column_names(path)
    cols === nothing && return nothing

    required_pairs = [
        (:global_pg_score, :global_qval),
        (:pg_score, :qval),
    ]

    available_pairs = [pair for pair in required_pairs if all(col -> col in cols, pair)]
    if isempty(available_pairs)
        missing_columns_by_pair = Dict(pair => [col for col in pair if col ∉ cols] for pair in required_pairs)
        @warn "No compatible protein score/q-value columns found for entrapment analysis" path=path missing_columns_by_pair=missing_columns_by_pair available_columns=join(string.(collect(cols)), ", ")
        return nothing
    end

    available_pairs
end

function entrapment_module_name(repo_path::AbstractString)
    project_path = joinpath(repo_path, "Project.toml")
    if isfile(project_path)
        try
            project = TOML.parsefile(project_path)
            name = get(project, "name", nothing)
            if name isa AbstractString
                return name
            end
        catch err
            @warn "Failed to parse Project.toml for entrapment module name" project_path=project_path error=err
        end
    end

    return "PioneerEntrapment"
end

function load_entrapment_module(repo_path::AbstractString)
    if !isdir(repo_path)
        @warn "Entrapment analyses repository not found" repo_path=repo_path
        return nothing
    end

    original_project = Base.active_project()
    original_load_path = copy(Base.LOAD_PATH)
    temp_env = mktempdir()

    module_name = "EntrapmentAnalyses"
    try
        Pkg.activate(temp_env; io=devnull)
        Pkg.develop(; path=repo_path, io=devnull)
        Pkg.instantiate(; io=devnull)
        pushfirst!(Base.LOAD_PATH, temp_env)
        module_name = entrapment_module_name(repo_path)
        return Base.require(Main, Symbol(module_name))
    catch err
        @warn "Unable to load entrapment module" repo_path=repo_path module_name=module_name error=err
        return nothing
    finally
        empty!(Base.LOAD_PATH)
        append!(Base.LOAD_PATH, original_load_path)
        if original_project === nothing
            Pkg.activate(; io=devnull)
        else
            Pkg.activate(original_project; io=devnull)
        end
    end
end

function compute_entrapment_metrics(dataset_dir::AbstractString, dataset_name::AbstractString)
    config = load_dataset_config(dataset_dir)
    config === nothing && return nothing

    lib_path = spectral_library_path_from_config(config, dataset_dir)
    lib_path === nothing && return nothing

    repo_path = get(ENV, "ENTRAPMENT_ANALYSES_PATH", joinpath(@__DIR__, "..", "..", "PioneerEntrapment.jl"))
    entrapment_module = load_entrapment_module(repo_path)
    entrapment_module === nothing && return nothing

    if !isdefined(entrapment_module, :run_efdr_analysis) ||
       !isdefined(entrapment_module, :run_protein_efdr_analysis)
        @warn "Entrapment analyses module missing required analysis functions" repo_path=repo_path
        return nothing
    end

    precursor_results_path = joinpath(dataset_dir, "precursors_long.arrow")
    protein_results_path = joinpath(dataset_dir, "protein_groups_long.arrow")

    missing_arrows = filter(x -> !isfile(x), [precursor_results_path, protein_results_path])
    if !isempty(missing_arrows)
        @warn "Skipping entrapment metrics; missing required arrow outputs" missing_arrows=missing_arrows
        return nothing
    end

    precursor_pairs = precursor_score_pairs(precursor_results_path)
    protein_pairs = protein_score_pairs(protein_results_path)

    if precursor_pairs === nothing || protein_pairs === nothing
        return nothing
    end

    output_dir = joinpath(dataset_dir, "entrapment_analysis")
    mkpath(output_dir)

    efdr_fn = getproperty(entrapment_module, :run_efdr_analysis)
    protein_fn = getproperty(entrapment_module, :run_protein_efdr_analysis)

    precursor_result = try
        Base.invokelatest(
            efdr_fn,
            precursor_results_path,
            lib_path;
            output_dir = output_dir,
            score_qval_pairs = precursor_pairs,
            plot_formats = [:png],
        )
    catch err
        @warn "Precursor entrapment analysis run failed" error=err
        nothing
    end

    protein_result = try
        Base.invokelatest(
            protein_fn,
            protein_results_path;
            output_dir = output_dir,
            score_qval_pairs = protein_pairs,
            plot_formats = [:png],
        )
    catch err
        @warn "Protein entrapment analysis run failed" error=err
        nothing
    end

    if precursor_result !== nothing
        @info "Precursor entrapment analysis completed"
    end

    if protein_result !== nothing
        @info "Protein entrapment analysis completed"
    end

    function summary_from_arrow(path::AbstractString, qval_col::Symbol, efdr_col::Symbol)
        if !isfile(path)
            @warn "Expected entrapment output missing" path=path
            return nothing
        end

        table = Arrow.Table(path; convert=false)
        df = DataFrame(table)
        if !(qval_col in names(df)) || !(efdr_col in names(df))
            @warn "Missing entrapment summary columns" path=path qval_col=qval_col efdr_col=efdr_col available_columns=join(string.(names(df)), ", ")
            return nothing
        end

        if nrow(df) == 0
            @warn "Entrapment output is empty" path=path
            return nothing
        end

        sort!(df, qval_col)
        qval_value = last(df[:, qval_col])
        efdr_value = last(df[:, efdr_col])

        return Dict(
            "qval" => qval_value,
            "paired_efdr" => efdr_value,
            "difference" => efdr_value - qval_value,
        )
    end

    summaries = Dict{String, Any}()

    precursor_summary = summary_from_arrow(
        joinpath(output_dir, "prec_results_with_efdr.arrow"),
        :MBR_boosted_qval,
        :MBR_boosted_prec_prob_paired_efdr,
    )
    if precursor_summary !== nothing
        summaries["precursors"] = precursor_summary
    end

    global_precursor_summary = summary_from_arrow(
        joinpath(output_dir, "global_results_with_efdr.arrow"),
        :MBR_boosted_global_qval,
        :MBR_boosted_global_prob_paired_edfr,
    )
    if global_precursor_summary !== nothing
        summaries["global_precursors"] = global_precursor_summary
    end

    protein_summary = summary_from_arrow(
        joinpath(output_dir, "protein_results_with_efdr.arrow"),
        :qval,
        :pg_score_paired_efdr,
    )
    if protein_summary !== nothing
        summaries["proteins"] = protein_summary
    end

    global_protein_summary = summary_from_arrow(
        joinpath(output_dir, "global_protein_results_with_efdr.arrow"),
        :global_qval,
        :global_pg_score_paired_efdr,
    )
    if global_protein_summary !== nothing
        summaries["global_proteins"] = global_protein_summary
    end

    metrics = Dict(
        "dataset" => dataset_name,
        "spectral_library" => lib_path,
        "output_dir" => output_dir,
    )

    if !isempty(summaries)
        metrics["summaries"] = summaries
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
