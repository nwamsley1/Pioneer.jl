"""
    koina_sequence(sequence, structural_mods) -> String

Derive the cold prediction-input string from residue annotations such as
`(2,C,Unimod:4)`. Preserve multiple modifications at a residue in input order.
This is not invoked during search; model-specific ion-mobility formatting
continues to belong to the prediction layer.
"""
function koina_sequence(sequence::AbstractString, structural_mods)
    tags = [String[] for _ in 0:length(sequence)]
    annotations = coalesce(structural_mods, "")
    pattern = r"\((\d+),[^,]+,([^()]+)\)"
    parsed = collect(eachmatch(pattern, annotations))
    require(join(m.match for m in parsed) == annotations, "Unsupported modification annotation")
    for annotation in parsed
        position = parse(Int, annotation[1])
        require(0 <= position <= length(sequence), "Modification position is outside the peptide")
        push!(tags[position + 1], "[" * uppercase(annotation[2]) * "]")
    end
    output = IOBuffer()
    print(output, join(tags[1]))
    for (position, residue) in enumerate(sequence)
        print(output, residue, join(tags[position + 1]))
    end
    return String(take!(output))
end

"""
    write_precursors(path, table)

Dictionary-encode string columns and omit `koina_sequence`. If that input column
is supplied, require that it is exactly derivable before removing it. There is
one output schema: no fallback that persists the redundant column.
"""
function write_precursors(path, table)
    names = collect(propertynames(table))
    if :koina_sequence in names
        for i in eachindex(table.sequence)
            require(isequal(koina_sequence(table.sequence[i], table.structural_mods[i]), table.koina_sequence[i]),
                "Koina input must be exactly derivable from sequence and structural modifications")
        end
    end
    filter!(!=(:koina_sequence), names)
    columns = Tuple(begin
        column = getproperty(table, name)
        Base.nonmissingtype(eltype(column)) <: AbstractString ?
            Arrow.DictEncode(collect(column)) : column
    end for name in names)
    Arrow.write(path, NamedTuple{Tuple(names)}(columns);
        metadata = Dict("compact_library_storage" => "1"))
    return nothing
end

"""Map dictionary columns directly. Never expand all precursor strings on load."""
function read_precursors(path)
    table = Arrow.Table(path)
    metadata = Arrow.getmetadata(table)
    require(metadata !== nothing && get(metadata, "compact_library_storage", nothing) == "1",
        "Expected the compact precursor-table schema")
    require(!hasproperty(table, :koina_sequence), "Compact precursor tables omit Koina input strings")
    return table
end
