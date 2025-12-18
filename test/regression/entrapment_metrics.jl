module EntrapmentMetrics

using Arrow
using DataFrames
using JSON
using Pkg
using TOML

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

export compute_entrapment_metrics

end
