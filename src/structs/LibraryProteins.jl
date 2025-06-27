"""
    LibraryProteins
Container for protein-level information loaded from an Arrow table.
"""
struct LibraryProteins
    data::Arrow.Table
end

"Return number of proteins in the library." 
Base.length(lp::LibraryProteins) = length(lp.data[:accession])

"Convert Arrow table to LibraryProteins." 
SetProteins(tbl::Arrow.Table) = LibraryProteins(tbl)

getGeneName(lp::LibraryProteins) = lp.data[:gene_name]
getGeneName(lp::LibraryProteins, idx::Integer) = lp.data[:gene_name][idx]
getProteinName(lp::LibraryProteins) = lp.data[:protein_name]
getProteinName(lp::LibraryProteins, idx::Integer) = lp.data[:protein_name][idx]
getAccession(lp::LibraryProteins) = lp.data[:accession]
getAccession(lp::LibraryProteins, idx::Integer) = lp.data[:accession][idx]
getOrganism(lp::LibraryProteins) = lp.data[:organism]
getOrganism(lp::LibraryProteins, idx::Integer) = lp.data[:organism][idx]
getSequence(lp::LibraryProteins) = lp.data[:sequence]
getSequence(lp::LibraryProteins, idx::Integer) = lp.data[:sequence][idx]
getLength(lp::LibraryProteins) = lp.data[:length]
getLength(lp::LibraryProteins, idx::Integer) = lp.data[:length][idx]