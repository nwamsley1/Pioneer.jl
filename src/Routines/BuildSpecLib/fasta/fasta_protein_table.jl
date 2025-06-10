module FastaProteinTable

using DataFrames

export build_protein_df

"""
    build_protein_df(entries::Vector{FastaEntry};
                     accession_regex::Union{Regex,Nothing}=nothing,
                     gene_regex::Union{Regex,Nothing}=nothing,
                     protein_regex::Union{Regex,Nothing}=nothing,
                     organism_regex::Union{Regex,Nothing}=nothing)

Create a DataFrame of protein information from FASTA entries.
Each regex is applied to the full FASTA header (identifier and description).
Missing matches yield empty strings. If `accession_regex` is `nothing`,
the full header is used as the accession.
"""
function build_protein_df(entries::Vector{FastaEntry};
                          accession_regex::Union{Regex,Nothing}=nothing,
                          gene_regex::Union{Regex,Nothing}=nothing,
                          protein_regex::Union{Regex,Nothing}=nothing,
                          organism_regex::Union{Regex,Nothing}=nothing)
    n = length(entries)
    gene = Vector{String}(undef, n)
    protein = Vector{String}(undef, n)
    accession = Vector{String}(undef, n)
    organism = Vector{String}(undef, n)
    length_vec = Vector{UInt32}(undef, n)
    sequence = Vector{String}(undef, n)

    for (i, entry) in enumerate(entries)
        header = string(get_id(entry), " ", get_description(entry))

        if accession_regex === nothing
            accession[i] = header
        else
            m = match(accession_regex, header)
            accession[i] = m === nothing ? header : String(first(m.captures))
        end

        if gene_regex === nothing
            gene[i] = ""
        else
            m = match(gene_regex, header)
            gene[i] = m === nothing ? "" : String(first(m.captures))
        end

        if protein_regex === nothing
            protein[i] = ""
        else
            m = match(protein_regex, header)
            protein[i] = m === nothing ? "" : String(first(m.captures))
        end

        if organism_regex === nothing
            organism[i] = ""
        else
            m = match(organism_regex, header)
            organism[i] = m === nothing ? "" : String(first(m.captures))
        end

        sequence[i] = get_sequence(entry)
        length_vec[i] = UInt32(length(sequence[i]))
    end

    return DataFrame(
        gene_name = gene,
        protein_name = protein,
        organism = organism,
        accession = accession,
        length = length_vec,
        sequence = sequence,
    )
end

end # module
