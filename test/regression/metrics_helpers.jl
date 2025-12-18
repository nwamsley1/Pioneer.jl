module RegressionMetricsHelpers

using DataFrames
using Statistics

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

is_numeric_column(col) = begin
    T = nonmissingtype(eltype(col))
    return T <: Number
end

function drop_non_quant_columns(col_names::AbstractVector{<:Union{Symbol, String}})
    [name for name in col_names if Symbol(name) ∉ NON_QUANT_COLUMNS]
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

function select_quant_columns(df::DataFrame, quant_col_names)
    quant_syms = Symbol.(quant_col_names)
    df_names = names(df)
    df_name_syms = Symbol.(df_names)
    df_names[[i for (i, n) in pairs(df_name_syms) if n in quant_syms]]
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

export NON_QUANT_COLUMNS
export is_numeric_column, quant_column_names_from_proteins, drop_non_quant_columns, select_quant_columns
export gene_names_column, species_column
export is_yeast_only_species, unique_species_value
export median_for_columns, mean_for_columns, condition_columns, resolve_run_columns
export count_species_ids, count_yeast_ids, count_nonyeast_ids, count_total_ids

end
