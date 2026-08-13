# Top-level spectral library container types.
#
# A SpectralLibrary holds the partitioned fragment index, precursor metadata,
# protein metadata, and the fragment lookup table (for detailed fragment access
# during deconvolution).

abstract type SpectralLibrary end

# Per-output column-exclusion policy. Lets BuildSpecLib indicate columns that
# should be dropped from final exports (e.g. `entrapment_group_id` when the
# library has no entrapment, `isotopic_mods` when there are no isotope mod
# groups). A policy with empty `excluded_columns` keeps every column.
struct OutputSchemaPolicy
    excluded_columns::Dict{Symbol, Tuple{Vararg{Symbol}}}
end

OutputSchemaPolicy() = OutputSchemaPolicy(Dict{Symbol, Tuple{Vararg{Symbol}}}())
OutputSchemaPolicy(::Nothing) = OutputSchemaPolicy()

function OutputSchemaPolicy(config::AbstractDict)
    excluded_columns = Dict{Symbol, Tuple{Vararg{Symbol}}}()
    precursor_exclusions = Symbol[]
    protein_group_exclusions = Symbol[]

    fasta_digest_params = get(config, "fasta_digest_params", nothing)
    entrapment_r = isnothing(fasta_digest_params) ? nothing : get(fasta_digest_params, "entrapment_r", nothing)
    if !isnothing(entrapment_r) && entrapment_r <= 0
        push!(precursor_exclusions, :entrapment_group_id)
        push!(protein_group_exclusions, :entrap_id)
    end

    isotope_mod_groups = get(config, "isotope_mod_groups", nothing)
    if !isnothing(isotope_mod_groups) && isempty(isotope_mod_groups)
        push!(precursor_exclusions, :isotopic_mods)
    end

    !isempty(precursor_exclusions) && (excluded_columns[:precursors] = Tuple(precursor_exclusions))
    !isempty(protein_group_exclusions) && (excluded_columns[:protein_groups] = Tuple(protein_group_exclusions))

    return OutputSchemaPolicy(excluded_columns)
end

excluded_output_columns(policy::OutputSchemaPolicy, output::Symbol) = get(policy.excluded_columns, output, ())
output_column_enabled(policy::OutputSchemaPolicy, output::Symbol, column::Symbol) = !(column in excluded_output_columns(policy, output))
output_column_enabled(policy::OutputSchemaPolicy, output::Symbol, column::AbstractString) = output_column_enabled(policy, output, Symbol(column))
enabled_output_columns(policy::OutputSchemaPolicy, output::Symbol, columns::AbstractVector) = [column for column in columns if output_column_enabled(policy, output, column)]
function enabled_output_table(policy::OutputSchemaPolicy, output::Symbol, tbl)
    column_names = Symbol.(Tables.columnnames(tbl))
    enabled_names = Tuple(name for name in column_names if output_column_enabled(policy, output, name))
    length(enabled_names) == length(column_names) && return tbl
    columns = Tables.columntable(tbl)
    return NamedTuple{enabled_names}(Tuple(getproperty(columns, name) for name in enabled_names))
end

struct FragmentIndexLibrary{P} <: SpectralLibrary
    presearch_partitioned_index::P
    partitioned_index::P
    precursors::LibraryPrecursors
    proteins::LibraryProteins
    fragment_lookup_table::StandardFragmentLookup
    output_schema_policy::OutputSchemaPolicy
end

struct SplineFragmentIndexLibrary{P} <: SpectralLibrary
    presearch_partitioned_index::P
    partitioned_index::P
    precursors::LibraryPrecursors
    proteins::LibraryProteins
    fragment_lookup_table::SplineFragmentLookup
    output_schema_policy::OutputSchemaPolicy
end

getPresearchPartitionedIndex(sl::SpectralLibrary) = sl.presearch_partitioned_index
getPartitionedIndex(sl::SpectralLibrary) = sl.partitioned_index
getPrecursors(sl::SpectralLibrary) = sl.precursors
getFragmentLookupTable(sl::SpectralLibrary) = sl.fragment_lookup_table
getProteins(sl::SpectralLibrary) = sl.proteins
getOutputSchemaPolicy(sl::SpectralLibrary) = sl.output_schema_policy
