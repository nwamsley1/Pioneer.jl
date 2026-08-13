# Library protein metadata access.
#
# LibraryProteins wraps an Arrow table of protein-level information
# (gene names, protein names, accession numbers).

struct LibraryProteins
    data::Arrow.Table
end

Base.length(lp::LibraryProteins) = length(lp.data[:accession])

SetProteins(tbl::Arrow.Table) = LibraryProteins(tbl)

getGeneName(lp::LibraryProteins) = lp.data[:gene_name]
getGeneName(lp::LibraryProteins, idx::Integer) = lp.data[:gene_name][idx]
getProteinName(lp::LibraryProteins) = lp.data[:protein_name]
getProteinName(lp::LibraryProteins, idx::Integer) = lp.data[:protein_name][idx]
getAccession(lp::LibraryProteins) = lp.data[:accession]
getAccession(lp::LibraryProteins, idx::Integer) = lp.data[:accession][idx]