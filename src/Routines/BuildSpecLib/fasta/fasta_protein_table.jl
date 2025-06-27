"""
    build_protein_df(entries::Vector{FastaEntry})

Create a DataFrame of protein information from pre-parsed FASTA entries.
"""
function build_protein_df(entries::Vector{FastaEntry})
    n = length(entries)
    gene = Vector{String}(undef, n)
    protein = Vector{String}(undef, n)
    accession = Vector{String}(undef, n)
    organism = Vector{String}(undef, n)
    length_vec = Vector{UInt32}(undef, n)
    sequence = Vector{String}(undef, n)


    for (i, entry) in enumerate(entries)
        gene[i] = get_gene(entry)
        protein[i] = get_protein(entry)
        accession[i] = get_id(entry)
        organism[i] = get_organism(entry)
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