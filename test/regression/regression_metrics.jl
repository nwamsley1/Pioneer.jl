#!/usr/bin/env julia

using CSV
using DataFrames
using JSON

function read_required_table(path::AbstractString)
    isfile(path) || error("Required file not found: $path")
    CSV.read(path, DataFrame; delim='\t', missingstring=["", "NA"], ignorerepeated=true)
end

is_numeric_column(col) = begin
    T = nonmissingtype(eltype(col))
    return T <: Number
end

function quant_columns(df::DataFrame)
    col_names = names(df)
    anchor_idx = findfirst(==("global_qval"), col_names)

    if anchor_idx !== nothing
        quant_start = anchor_idx + 1
        quant_start > ncol(df) && return Int[]
        return collect(quant_start:ncol(df))
    end

    numeric_flags = [is_numeric_column(df[:, c]) for c in col_names]
    any(numeric_flags) || return Int[]

    last_non_numeric = findlast(!, numeric_flags)
    if last_non_numeric === nothing
        return collect(1:ncol(df))
    end

    quant_start = last_non_numeric + 1
    quant_start > ncol(df) && return Int[]
    collect(quant_start:ncol(df))
end

function compute_wide_metrics(df::DataFrame)
    quant_cols = quant_columns(df)
    runs = length(quant_cols)
    if runs == 0
        return (; runs = 0, complete_rows = 0, non_missing_values = 0, data_completeness = nothing)
    end

    quant_data = df[:, quant_cols]
    quant_matrix = Matrix(quant_data)

    complete_rows = sum(all(!ismissing, row) for row in eachrow(quant_matrix))
    non_missing_values = count(!ismissing, quant_matrix)
    total_cells = nrow(df) * runs

    (; runs, complete_rows, non_missing_values, data_completeness = total_cells > 0 ? non_missing_values / total_cells : nothing)
end

function compute_dataset_metrics(dataset_dir::AbstractString, dataset_name::AbstractString)
    required_files = [
        "precursors_long.tsv",
        "precursors_wide.tsv",
        "protein_groups_long.tsv",
        "protein_groups_wide.tsv",
    ]

    missing_files = filter(f -> !isfile(joinpath(dataset_dir, f)), required_files)
    if !isempty(missing_files)
        @warn "Skipping dataset due to missing outputs" dataset=dataset_name missing_files=missing_files
        return nothing
    end

    precursors_long = read_required_table(joinpath(dataset_dir, "precursors_long.tsv"))
    precursors_wide = read_required_table(joinpath(dataset_dir, "precursors_wide.tsv"))
    protein_groups_long = read_required_table(joinpath(dataset_dir, "protein_groups_long.tsv"))
    protein_groups_wide = read_required_table(joinpath(dataset_dir, "protein_groups_wide.tsv"))

    precursor_wide_metrics = compute_wide_metrics(precursors_wide)
    protein_wide_metrics = compute_wide_metrics(protein_groups_wide)

    return Dict(
        "dataset" => dataset_name,
        "precursors" => Dict(
            "total" => nrow(precursors_long),
            "unique" => nrow(precursors_wide),
            "runs" => precursor_wide_metrics.runs,
            "complete_rows" => precursor_wide_metrics.complete_rows,
            "non_missing_values" => precursor_wide_metrics.non_missing_values,
            "data_completeness" => precursor_wide_metrics.data_completeness,
        ),
        "protein_groups" => Dict(
            "total" => nrow(protein_groups_long),
            "unique" => nrow(protein_groups_wide),
            "runs" => protein_wide_metrics.runs,
            "complete_rows" => protein_wide_metrics.complete_rows,
            "non_missing_values" => protein_wide_metrics.non_missing_values,
            "data_completeness" => protein_wide_metrics.data_completeness,
        ),
    )
end

function main()
    results_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(pwd(), "results")
    isdir(results_dir) || error("Results directory does not exist: $results_dir")

    dataset_dirs = filter(name -> isdir(joinpath(results_dir, name)), readdir(results_dir))
    isempty(dataset_dirs) && error("No dataset directories found in $results_dir")

    for dataset_name in dataset_dirs
        dataset_dir = joinpath(results_dir, dataset_name)
        @info "Processing dataset" dataset=dataset_name
        metrics = compute_dataset_metrics(dataset_dir, dataset_name)
        metrics === nothing && continue

        output_path = joinpath(dataset_dir, "metrics_$(dataset_name).json")
        open(output_path, "w") do io
            JSON.print(io, metrics; indent=2)
        end
    end
end

main()