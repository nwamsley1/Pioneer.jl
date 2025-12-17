#!/usr/bin/env julia

using Arrow
using CSV
using DataFrames
using JSON
using Pkg
using Statistics
using TOML

const DEFAULT_METRIC_GROUPS = ["identification", "CV", "eFDR"]

const NON_QUANT_COLUMNS = Set([
    :target,
    :isotopic_mods,
    :prec_mz,
    :global_score,
    :MBR_boosted_global_qval,
    :use_for_protein_quant,
    :precursor_idx,
    :entrapment_group_id,
])

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
        return drop_non_quant_columns(Symbol.(col_names[quant_start:end]))
    end

    numeric_flags = [is_numeric_column(df[:, c]) for c in col_names]
    any(numeric_flags) || return Symbol[]

    last_non_numeric = findlast(!, numeric_flags)
    if last_non_numeric === nothing
        return Symbol.(col_names)
    end

    quant_start = last_non_numeric + 1
    quant_start > ncol(df) && return Symbol[]
    drop_non_quant_columns(Symbol.(col_names[quant_start:end]))
end

function drop_non_quant_columns(col_names::AbstractVector{<:Union{Symbol, String}})
    [name for name in col_names if Symbol(name) ∉ NON_QUANT_COLUMNS]
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
        cols = Symbol.(names(df))
        if !(qval_col in cols) || !(efdr_col in cols)
            @warn "Missing entrapment summary columns" path=path qval_col=qval_col efdr_col=efdr_col available_columns=join(string.(cols), ", ")
            return nothing
        end

        if nrow(df) == 0
            @warn "Entrapment output is empty" path=path
            return nothing
        end

        qvals = collect(df[:, qval_col])
        efdrs = collect(df[:, efdr_col])

        best_index = nothing
        best_qval = nothing
        for i in eachindex(qvals)
            qval = qvals[i]
            efdr = efdrs[i]
            if qval === missing || efdr === missing
                continue
            end

            if best_index === nothing || qval > best_qval
                best_index = i
                best_qval = qval
            end
        end

        if best_index === nothing
            @warn "Entrapment output only contains missing values" path=path qval_col=qval_col efdr_col=efdr_col
            return nothing
        end

        qval_value = qvals[best_index]
        efdr_value = efdrs[best_index]

        if !(qval_value isa Real) || !(efdr_value isa Real)
            @warn "Entrapment summary values are non-numeric" path=path qval_value=qval_value efdr_value=efdr_value
            return nothing
        end

        qval_value = Float64(qval_value)
        efdr_value = Float64(efdr_value)

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
        :MBR_boosted_global_prob_paired_efdr,
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

    metrics = Dict{String, Any}()

    if !isempty(summaries)
        metrics["summaries"] = summaries
    end

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
    need_three_proteome = ("fold_change" in requested_groups) || ("three_proteome" in requested_groups)
    need_tsv_metrics = need_identification || need_cv || need_keap1 || need_ftr || need_three_proteome

    precursors_metrics = nothing
    protein_metrics = nothing
    keap1_precursor_metrics = nothing
    keap1_protein_metrics = nothing
    ftr_metrics = nothing
    fold_change_metrics = nothing

    if need_tsv_metrics
        required_files = if need_identification || need_cv || need_keap1 || need_ftr
            [
                "precursors_long.tsv",
                "precursors_wide.tsv",
                "protein_groups_long.tsv",
                "protein_groups_wide.tsv",
            ]
        else
            ["protein_groups_wide.tsv"]
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
            precursors_long = read_required_table(joinpath(dataset_dir, "precursors_long.tsv"))
            precursors_wide = read_required_table(joinpath(dataset_dir, "precursors_wide.tsv"))
            protein_groups_long = read_required_table(joinpath(dataset_dir, "protein_groups_long.tsv"))
        end

        protein_groups_wide = read_required_table(joinpath(dataset_dir, "protein_groups_wide.tsv"))

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
        numerator = get(pair, "numerator_condition", get(pair, "numerator", nothing))
        denominator = get(pair, "denominator_condition", get(pair, "denominator", nothing))
        numerator === nothing && continue
        denominator === nothing && continue

        expected_raw = get(pair, "species_ratios", get(pair, "expected_ratios", Dict()))
        expected = Dict{String, Float64}()
        if expected_raw isa AbstractDict
            for (species, ratio) in expected_raw
                try
                    expected[uppercase(String(species))] = Float64(ratio)
                catch
                end
            end
        end

        push!(
            condition_pairs,
            (; numerator = String(numerator), denominator = String(denominator), expected),
        )
    end

    (; run_to_condition = run_mapping, condition_pairs)
end

function load_three_proteome_designs(path::AbstractString)
    if isdir(path)
        files = filter(f -> endswith(f, "_ED.json"), readdir(path; join=true))
        isempty(files) && return Dict{String, Any}()

        designs = Dict{String, Any}()
        for file in files
            parsed = load_three_proteome_designs(file)
            if parsed isa NamedTuple
                basename_no_ext = replace(basename(file), r"\.json$" => "")
                if endswith(basename_no_ext, "_ED")
                    designs[basename_no_ext] = parsed
                    designs[basename_no_ext[1:end-3]] = parsed
                else
                    designs[basename_no_ext] = parsed
                end
            elseif parsed isa Dict
                merge!(designs, parsed)
            end
        end

        return designs
    end

    if !isfile(path)
        return Dict{String, Any}()
    end

    try
        parsed = JSON.parsefile(path)
        if parsed isa AbstractDict
            contains_direct_design =
                haskey(parsed, "runs") || haskey(parsed, "expected_ratios") || haskey(parsed, "condition_pairs")
            if contains_direct_design
                return normalize_three_proteome_design(Dict{String, Any}(parsed))
            end

            designs = Dict{String, Any}()
            for (k, v) in parsed
                v isa AbstractDict || continue
                designs[String(k)] = normalize_three_proteome_design(Dict{String, Any}(v))
            end

            return designs
        end
    catch err
        @warn "Failed to parse three-proteome design file; ignoring" design_path=path error=err
    end

    Dict{String, Any}()
end

function three_proteome_design_entry(three_proteome_designs, dataset_name::AbstractString)
    if three_proteome_designs isa NamedTuple
        return three_proteome_designs
    elseif three_proteome_designs isa AbstractDict
        entry = get(three_proteome_designs, dataset_name, nothing)
        if entry === nothing && haskey(three_proteome_designs, "runs")
            entry = three_proteome_designs
        end
        return entry
    end

    nothing
end

function load_experimental_design(path::AbstractString)
    if isdir(path)
        files = filter(f -> endswith(f, ".json"), readdir(path; join=true))
        if isempty(files)
            @info "No experimental design files found; using run names as labels" experimental_design_dir=path
            return Dict{String, Any}()
        end

        designs = Dict{String, Any}()

        for file in files
            dataset_key = replace(replace(basename(file), r"_ED\.json$" => ""), r"\.json$" => "")
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

function normalize_metric_label(label::AbstractString)
    replace(strip(label), r"\s+" => "_")
end

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

function all_runs_from_groups(groupings::Dict{String, Vector{String}})
    collected = String[]
    for runs in values(groupings)
        append!(collected, runs)
    end

    unique(collected)
end

function gene_names_column(df::DataFrame; table_label::AbstractString = "table")
    gene_col_index = findfirst(name -> String(name) == "gene_names", names(df))
    if gene_col_index === nothing
        available_columns = join(string.(names(df)), ", ")
        @warn "Table missing gene_names column; skipping gene-based metrics" table=table_label available_columns=available_columns
        return nothing
    end

    names(df)[gene_col_index]
end

function gene_counts_by_run(
    wide_df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}},
    gene_term::AbstractString,
; table_label::AbstractString = "table")
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

function species_column(df::DataFrame; table_label::AbstractString = "table")
    species_col_index = findfirst(name -> String(name) == "species", names(df))
    if species_col_index === nothing
        available_columns = join(string.(names(df)), ", ")
        @warn "Table missing species column; skipping species-based metrics" table=table_label available_columns=available_columns
        return nothing
    end

    names(df)[species_col_index]
end

function is_yeast_only_species(val)
    val === missing && return false
    parts = split(String(val), ";")
    cleaned = unique(filter(!isempty, uppercase.(strip.(parts))))
    length(cleaned) == 1 && cleaned[1] == "YEAST"
end

function unique_species_value(val)
    val === missing && return nothing
    parts = split(String(val), ";")
    cleaned = unique(filter(!isempty, uppercase.(strip.(parts))))
    length(cleaned) == 1 ? cleaned[1] : nothing
end

function median_for_columns(row::DataFrameRow, cols::AbstractVector)
    values = [row[c] for c in cols if row[c] !== missing]
    isempty(values) && return missing
    median(values)
end

function mean_for_columns(row::DataFrameRow, cols::AbstractVector)
    values = [row[c] for c in cols if row[c] !== missing]
    isempty(values) && return missing
    mean(values)
end

function condition_columns(
    quant_col_names::AbstractVector{<:Union{Symbol, String}},
    run_to_condition::Dict{String, String},
)
    quant_lookup = Dict(String(c) => c for c in quant_col_names)
    condition_to_columns = Dict{String, Vector{eltype(quant_col_names)}}()
    missing_runs = String[]

    for (run, condition) in run_to_condition
        col = get(quant_lookup, run, nothing)
        if col === nothing
            push!(missing_runs, run)
            continue
        end

        cols = get!(condition_to_columns, condition, Vector{eltype(quant_col_names)}())
        push!(cols, col)
    end

    return condition_to_columns, missing_runs
end

function resolve_run_columns(
    df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}},
    requested_runs,
)
    quant_columns = select_quant_columns(df, quant_col_names)
    quant_strings = Set(String.(quant_columns))

    if requested_runs === nothing
        return String.(quant_columns)
    end

    run_names = String.(requested_runs)
    String[x for x in run_names if x in quant_strings]
end

function count_species_ids(
    df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}},
    run_names;
    predicate::Function,
    table_label::AbstractString = "protein_groups",
)
    species_col = species_column(df; table_label = table_label)
    species_col === nothing && return 0

    run_columns = resolve_run_columns(df, quant_col_names, run_names)
    if isempty(run_columns)
        if run_names !== nothing
            available_runs = select_quant_columns(df, quant_col_names)
            @warn "Requested runs not found in table; skipping species-based counts" table=table_label requested_runs=run_names available_runs=available_runs
        end
        return 0
    end

    matches = map(df[:, species_col]) do val
        predicate(val)
    end

    matching_indices = findall(identity, matches)
    isempty(matching_indices) && return 0

    quant_matrix = Matrix(df[matching_indices, run_columns])
    count(!ismissing, quant_matrix)
end

function count_yeast_ids(
    df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}},
    run_names = nothing;
    table_label::AbstractString = "protein_groups",
)
    count_species_ids(df, quant_col_names, run_names; predicate = is_yeast_only_species, table_label = table_label)
end

function count_nonyeast_ids(
    df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}},
    run_names = nothing;
    table_label::AbstractString = "protein_groups",
)
    count_species_ids(df, quant_col_names, run_names; predicate = val -> !is_yeast_only_species(val), table_label = table_label)
end

function count_total_ids(
    df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}},
    run_names = nothing;
    table_label::AbstractString = "protein_groups",
)
    run_columns = resolve_run_columns(df, quant_col_names, run_names)
    if isempty(run_columns)
        if run_names !== nothing
            available_runs = select_quant_columns(df, quant_col_names)
            @warn "Requested runs not found in table; skipping total counts" table=table_label requested_runs=run_names available_runs=available_runs
        end
        return 0
    end

    quant_matrix = Matrix(df[:, run_columns])
    count(!ismissing, quant_matrix)
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
        tsv_path = joinpath(mbr_path, "precursors_wide.tsv")
        isfile(tsv_path) || begin
            @warn "Missing precursors for MBR dataset; skipping FTR metrics" dataset=mbr_name path=tsv_path
            return nothing
        end
        read_required_table(tsv_path)
    end

    precursors_no_mbr = if dataset_name == nombr_name
        precursors_wide
    else
        tsv_path = joinpath(nombr_path, "precursors_wide.tsv")
        isfile(tsv_path) || begin
            @warn "Missing precursors for noMBR dataset; skipping FTR metrics" dataset=nombr_name path=tsv_path
            return nothing
        end
        read_required_table(tsv_path)
    end

    protein_groups_mbr = if dataset_name == mbr_name
        protein_groups_wide
    else
        tsv_path = joinpath(mbr_path, "protein_groups_wide.tsv")
        isfile(tsv_path) || begin
            @warn "Missing protein groups for MBR dataset; skipping FTR metrics" dataset=mbr_name path=tsv_path
            return nothing
        end
        read_required_table(tsv_path)
    end

    protein_groups_no_mbr = if dataset_name == nombr_name
        protein_groups_wide
    else
        tsv_path = joinpath(nombr_path, "protein_groups_wide.tsv")
        isfile(tsv_path) || begin
            @warn "Missing protein groups for noMBR dataset; skipping FTR metrics" dataset=nombr_name path=tsv_path
            return nothing
        end
        read_required_table(tsv_path)
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
    three_proteome_designs = load_three_proteome_designs(three_proteome_designs_path)

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
