module ThreeProteomeMetrics

using DataFrames
using JSON
using Statistics

include("metrics_helpers.jl")
using .RegressionMetricsHelpers: condition_columns, gene_names_column, mean_for_columns, select_quant_columns, species_column, unique_species_value

function load_experimental_design(path::AbstractString)
    if isdir(path)
        files = filter(f -> endswith(f, ".ED.json") || endswith(f, "_ED.json") || endswith(f, ".json"), readdir(path; join=true))
        if isempty(files)
            @info "No experimental design files found; using run names as labels" experimental_design_dir=path
            return Dict{String, Any}()
        end

        designs = Dict{String, Any}()

        for file in files
            dataset_key = replace(
                replace(replace(basename(file), r"\.ED\.json$" => ""), r"_ED\.json$" => ""),
                r"\.json$" => "",
            )
            parsed = load_experimental_design(file)
            isempty(parsed) && continue
            designs[dataset_key] = parsed
        end

        return designs
    end

    if !isfile(path)
        @info "No experimental design file found; using run names as labels" experimental_design_path=path
        return Dict{String, Any}()
    end

    try
        design = JSON.parsefile(path)
        if design isa AbstractDict
            has_runs = haskey(design, "runs")
            if has_runs && length(design) == 1
                return Dict{String, Any}("runs" => design["runs"])
            elseif has_runs
                return design
            end

            mapped_design = Dict{String, Any}()
            for (k, v) in design
                mapped_design[String(k)] = v
            end
            return mapped_design
        end

        @warn "Experimental design file is not a dictionary; ignoring" experimental_design_path=path
        return Dict{String, Any}()
    catch err
        @warn "Failed to parse experimental design file; ignoring" experimental_design_path=path error=err
        return Dict{String, Any}()
    end
end

normalize_metric_label(label::AbstractString) = replace(strip(label), r"\s+" => "_")

function experimental_design_entry(
    experimental_design::Dict{String, Any},
    dataset_name::AbstractString,
)
    entry = get(experimental_design, dataset_name, nothing)
    if entry === nothing && haskey(experimental_design, "runs")
        entry = experimental_design
    end

    if entry isa AbstractDict
        mapped = Dict{String, Any}()
        for (k, v) in entry
            mapped[String(k)] = v
        end
        return mapped
    end

    Dict{String, Any}()
end

function experimental_design_for_dataset(
    experimental_design::Dict{String, Any},
    dataset_name::AbstractString,
)
    entry = experimental_design_entry(experimental_design, dataset_name)
    runs = get(entry, "runs", nothing)
    if runs isa AbstractDict
        return Dict(String(k) => String(v) for (k, v) in runs)
    end

    Dict{String, String}()
end

function run_groups_for_dataset(
    experimental_design::Dict{String, Any},
    dataset_name::AbstractString;
    three_proteome_designs=nothing,
)
    entry = experimental_design_entry(experimental_design, dataset_name)
    grouping = get(entry, "composition", nothing)
    if grouping isa AbstractDict
        groups = Dict{String, Vector{String}}()
        for (group, runs) in grouping
            if runs isa AbstractVector
                groups[String(group)] = [String(r) for r in runs]
            end
        end

        isempty(groups) || return groups
    end

    runs_mapping = get(entry, "runs", nothing)
    if runs_mapping isa AbstractDict
        groups = Dict{String, Vector{String}}()
        for (run, condition) in runs_mapping
            condition_key = String(condition)
            push!(get!(groups, condition_key, String[]), String(run))
        end

        isempty(groups) || return groups
    end

    if three_proteome_designs !== nothing
        tp_entry = three_proteome_design_entry(three_proteome_designs, dataset_name)
        tp_runs_mapping = tp_entry isa AbstractDict ? get(tp_entry, "runs", nothing) : nothing
        if tp_runs_mapping isa AbstractDict
            groups = Dict{String, Vector{String}}()
            for (run, condition) in tp_runs_mapping
                condition_key = String(condition)
                push!(get!(groups, condition_key, String[]), String(run))
            end
            isempty(groups) || return groups
        end
    end

    Dict{String, Vector{String}}()
end

function gene_counts_by_run(
    wide_df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}},
    gene_term::AbstractString;
    table_label::AbstractString = "table",
)
    gene_col = gene_names_column(wide_df; table_label = table_label)
    gene_col === nothing && return Dict{Symbol, Int}()

    quant_columns = select_quant_columns(wide_df, quant_col_names)
    isempty(quant_columns) && return Dict{Symbol, Int}()

    gene_matches = map(wide_df[:, gene_col]) do val
        val === missing && return false
        occursin(gene_term, String(val))
    end

    matched_rows = findall(identity, gene_matches)
    isempty(matched_rows) && return Dict{Symbol, Int}()

    Dict(col => count(!ismissing, wide_df[matched_rows, col]) for col in quant_columns)
end

function gene_counts_metrics_by_run(
    wide_df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}},
    labels_for_runs::Dict{String, String},
    level::AbstractString,
    gene_term::AbstractString,
)
    counts = gene_counts_by_run(wide_df, quant_col_names, gene_term; table_label = level)
    metrics = Dict{String, Int}()

    for (col, count) in counts
        label = get(labels_for_runs, String(col), String(col))
        metric_name = string(gene_term, "_", level, "_in_", normalize_metric_label(label))
        metrics[metric_name] = count
    end

    metrics
end

function fold_change_metrics_for_table(
    df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}},
    design,
    condition_pairs;
    table_label::AbstractString,
)
    quant_columns = select_quant_columns(df, quant_col_names)
    if isempty(quant_columns)
        @warn "No quantification columns available for fold-change metrics" table=table_label
        return nothing
    end

    condition_to_columns, missing_runs = condition_columns(quant_columns, design.run_to_condition)
    if !isempty(missing_runs)
        @warn "Runs listed in three-proteome design missing from table; skipping those runs" table=table_label missing_runs=missing_runs
    end

    species_col = species_column(df; table_label = table_label)
    species_col === nothing && return nothing

    metrics = Dict{String, Any}()

    for pair in condition_pairs
        numerator_columns = get(condition_to_columns, pair.numerator, Vector{eltype(quant_columns)}())
        denominator_columns = get(condition_to_columns, pair.denominator, Vector{eltype(quant_columns)}())
        if isempty(numerator_columns) || isempty(denominator_columns)
            @warn "Missing runs for condition pair; skipping fold-change computation" table=table_label numerator=pair.numerator denominator=pair.denominator
            continue
        end

        log_ratios = Dict{String, Vector{Float64}}()
        expected_log2_by_species = Dict{String, Float64}()

        for row in eachrow(df)
            species = unique_species_value(row[species_col])
            species === nothing && continue

            expected_ratio = get(pair.expected, species, nothing)
            expected_ratio === nothing && continue

            numerator_mean = mean_for_columns(row, numerator_columns)
            denominator_mean = mean_for_columns(row, denominator_columns)

            if numerator_mean === missing || denominator_mean === missing || denominator_mean == 0
                continue
            end

            observed_ratio = numerator_mean / denominator_mean
            if observed_ratio <= 0 || expected_ratio <= 0
                continue
            end

            observed_log2 = log2(observed_ratio)
            expected_log2 = log2(expected_ratio)

            push!(get!(log_ratios, species, Float64[]), observed_log2)
            expected_log2_by_species[species] = get(expected_log2_by_species, species, expected_log2)
        end

        pair_label = string(
            normalize_metric_label(pair.numerator),
            "_over_",
            normalize_metric_label(pair.denominator),
        )

        pair_metrics = Dict{String, Any}()
        for (species, values) in log_ratios
            median_log2_fc = isempty(values) ? missing : median(values)
            expected_log2 = get(expected_log2_by_species, species, missing)
            deviation = (median_log2_fc === missing || expected_log2 === missing) ? missing : median_log2_fc - expected_log2

            if deviation !== missing
                @info "Fold-change deviation" table=table_label pair=pair_label species=species median_log2_fc=median_log2_fc expected_log2_fc=expected_log2 deviation=deviation entries=length(values)
            end

            pair_metrics[string(lowercase(species), "_median_deviation")] = deviation
            pair_metrics[string(lowercase(species), "_entries")] = length(values)
        end

        isempty(pair_metrics) || (metrics[pair_label] = pair_metrics)
    end

    isempty(metrics) ? nothing : metrics
end

function normalize_three_proteome_design(design::Dict{String, Any})
    run_mapping = Dict{String, String}()
    if haskey(design, "runs") && design["runs"] isa AbstractDict
        for (run, condition) in design["runs"]
            run_mapping[String(run)] = String(condition)
        end
    end

    raw_pairs = if haskey(design, "expected_ratios")
        design["expected_ratios"]
    elseif haskey(design, "condition_pairs")
        design["condition_pairs"]
    else
        Any[]
    end

    condition_pairs =
        Vector{NamedTuple{(:numerator, :denominator, :expected), Tuple{String, String, Dict{String, Float64}}}}()

    for pair in raw_pairs
        pair isa AbstractDict || continue

        numerator = get(pair, "numerator", nothing)
        denominator = get(pair, "denominator", nothing)
        expected = get(pair, "expected", nothing)

        if !(numerator isa AbstractString && denominator isa AbstractString && expected isa AbstractDict)
            continue
        end

        expected_ratios = Dict{String, Float64}()
        for (species, ratio) in expected
            ratio isa Real || continue
            expected_ratios[String(species)] = Float64(ratio)
        end

        push!(condition_pairs, (; numerator = String(numerator), denominator = String(denominator), expected = expected_ratios))
    end

    (; run_mapping, condition_pairs)
end

function load_three_proteome_designs(path::AbstractString)
    if isdir(path)
        files = filter(f -> endswith(f, "_three_proteome.json"), readdir(path; join=true))
        if isempty(files)
            @info "No three-proteome design files found" three_proteome_design_dir=path
            return Dict{String, Any}()
        end

        designs = Dict{String, Any}()
        for file in files
            dataset_key = replace(basename(file), r"_three_proteome\.json$" => "")
            parsed = load_three_proteome_designs(file)
            isempty(parsed) && continue
            designs[dataset_key] = parsed
        end

        return designs
    end

    if !isfile(path)
        @info "No three-proteome design file found" three_proteome_design_path=path
        return Dict{String, Any}()
    end

    try
        design = JSON.parsefile(path)
        design isa AbstractDict || return Dict{String, Any}()
        normalized = Dict{String, Any}()
        for (k, v) in design
            normalized[String(k)] = normalize_three_proteome_design(v)
        end
        return normalized
    catch err
        @warn "Failed to parse three-proteome design file; ignoring" three_proteome_design_path=path error=err
        return Dict{String, Any}()
    end
end

function three_proteome_design_entry(three_proteome_designs, dataset_name::AbstractString)
    if three_proteome_designs === nothing
        return nothing
    end

    entry = get(three_proteome_designs, dataset_name, nothing)
    if entry === nothing
        @warn "No three-proteome design found for dataset" dataset=dataset_name
        return nothing
    end

    if !haskey(entry, :run_mapping) || !haskey(entry, :condition_pairs)
        @warn "Three-proteome design missing required fields" dataset=dataset_name
        return nothing
    end

    (; run_to_condition = entry.run_mapping, condition_pairs = entry.condition_pairs)
end

export load_experimental_design, normalize_metric_label, experimental_design_for_dataset, run_groups_for_dataset
export gene_counts_metrics_by_run, fold_change_metrics_for_table
export load_three_proteome_designs, three_proteome_design_entry

end
