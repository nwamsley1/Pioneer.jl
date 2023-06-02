"""
    digest(sequence::AbstractString; regex::Regex = r"[KR][^P|$]", max_length::Int = 40, min_length::Int = 8, missed_cleavages::Int = 1)

Given an amino acid sequence, `sequence`, and a regular expression that matches enzymatic cleavage sites, finds all peptide cleavage products from the
amino acid sequence. Can set a minimum and maximum length for cleaved peptides. Gets all peptides with 0:N missed cleavages. Default cleavage site
regex is for trypsin. Returns a list of strings that are the peptide cleavage products. 

- `sequence::AbstractString` -- An amino acid sequence
- `regex::Regex` -- Cleavage site regex. 
- `max_length::Int` -- Exclude peptides longer than this
- `min_length::Int` -- Exclude peptides shorter than this
- `missed_cleavages::Int` -- Get peptides with 0 to this many missed cleavages

"""
function digest(sequence::AbstractString; regex::Regex = r"[KR][^P|$]", max_length::Int = 40, min_length::Int = 8, missed_cleavages::Int = 1)
    
    function addPeptide!(peptides::Vector{SubString{String}}, n::Int, sequence::AbstractString, site::Int, previous_sites::Vector{Int64}, min_length::Int, max_length::Int, missed_cleavages::Int)
        for i in 1:min(n, missed_cleavages + 1)
            previous_site = previous_sites[end - i + 1]
            if ((site - previous_site) >= min_length) && ((site - previous_site) <= max_length)
                    @inbounds push!(peptides, @view sequence[previous_site+1:site])
            end
        end

        @inbounds for i in 1:length(previous_sites)-1
            previous_sites[i] = previous_sites[i+1]
        end
        n += 1
        previous_sites[end] = site
        return n
    end

    #Iterator of cleavage sites in `sequence`
    cleavage_sites = eachmatch(regex, sequence, overlap = true)
    #in silico digested peptides from `sequence`
    peptides = Vector{SubString{String}}(undef, 0)
    #Positions in `sequence` of 1:(missed_cleavages + 1) cleavage sites
    previous_sites = zeros(Int64, missed_cleavages + 1)
    previous_sites[1] = 1
    #Number of cleavage sites encountered so far
    n = 1
    
    #Advance through each cleavage_site and see if there are any new peptides
    #that satisfy the max_length, min_length, and missed_cleavages thresholds. 
    for site in cleavage_sites
        n = addPeptide!(peptides, n, sequence, site.offset, previous_sites, min_length, max_length, missed_cleavages)
    end

    #Peptides containint the C-terminus are added here
    n = addPeptide!(peptides, n, sequence, length(sequence), previous_sites, min_length, max_length, missed_cleavages)
    return peptides 
end


using FASTX
using CodecZlib
using Dictionaries
file_path = "/Users/n.t.wamsley/RIS_temp/HAMAD_MAY23/mouse_SIL_List/UP000000589_10090.fasta.gz"
#function parseRule(identifier::String)
#    split(identifier, "|")[2]
#end

struct FastaEntry
    uniprot_id::String
    description::String
    sequence::String
end

getID(fe::FastaEntry) = fe.uniprot_id
getDescription(fe::FastaEntry) = fe.description
getSeq(fe::FastaEntry) = fe.sequence


FastaEntry() = FastaEntry("","","")

"""
    parseFasta(fasta_path::String, parse_identifier::Function = x -> x)::Vector{FastaEntry}

Given a fasta file, loads the identifier, description, and sequence into an in memory representation 
`FastaEntry` for each entry. Returns Vector{FastaEntry}

- `fasta_path::String` -- Path to a ".fasta" or ".fasta.gz"
- `parse_identifier::Function ` -- Function takes a String argument and returns a String 

"""
function parseFasta(fasta_path::String, parse_identifier::Function = x -> split(x,"|")[2])::Vector{FastaEntry}

    function getReader(fasta_path::String)
        if endswith(fasta_path, ".fasta.gz")
            return FASTA.Reader(GzipDecompressorStream(open(fasta_path)))
        elseif endswith(fasta_path, ".fasta")
            return FASTA.Reader(open(fasta_path))
        else
            throw(ErrorException("fasta_path \"$fasta_path\" did not end with `.fasta` or `.fasta.gz`"))
        end
    end

    #I/O for Fasta
    reader = getReader(fasta_path)

    #In memory representation of FastaFile
    fasta = Vector{FastaEntry}()
    @time begin
        for record in reader
            push!(fasta, 
                    FastaEntry(parse_identifier(FASTA.identifier(record)),
                               FASTA.description(record),
                               FASTA.sequence(record))
            )
        end
    end

    return fasta
end

function digestFasta(fasta::Vector{FastaEntry}; regex::Regex = r"[KR][^P|$]", max_length::Int = 40, min_length::Int = 8, missed_cleavages::Int = 1)
    peptides_fasta = Vector{FastaEntry}()
    for entry in fasta
        for peptide in digest(getSeq(entry), regex = regex, max_length = max_length, min_length = min_length, missed_cleavages = missed_cleavages)
            push!(peptides_fasta, FastaEntry(getID(entry), 
                                                "",#getDescription(entry), 
                                                peptide))
        end
    end
    return sort!(peptides_fasta, by = x -> getSeq(x))
end

protien_fasta = parseFasta(file_path)
peptides_fasta = digestFasta(protien_fasta)
struct Protein
    name::String
    pepGroup_ids::Set{UInt32}
end
Protein(name::String) = Protein(name, Set{UInt32}())
getName(p::Protein) = p.name
getPepGroupIDs(p::Protein) = p.pepGroup_ids

struct PeptideGroup
    prot_ids::Set{UInt32}
    pep_ids::Set{UInt32}
    sequence::String
end
PeptideGroup(sequence::String) = PeptideGroup(Set{UInt32}(), Set{UInt32}(), sequence)

getProtIDs(pg::PeptideGroup) = pg.prot_ids
getPepIDs(pg::PeptideGroup) = pg.pep_ids
getSeq(pg::PeptideGroup) = pg.sequence
addProtID!(pg::PeptideGroup, prot_id::UInt32) = push!(pg.prot_ids, prot_id)
addPepID!(pg::PeptideGroup, pep_id::UInt32) = push!(pg.pep_ids, pep_id)

struct Peptide
    sequence::String
    pepGroup_id::UInt32
    prec_ids::Set{UInt32}
end
getSeq(p::Peptide) = p.sequence
getPepGroupID(p::Peptide) = p.pepGroup_id
getPrecIDs(p::Peptide) = p.prec_ids
addPrecID!(p::Peptide, prec_id::UInt32) = push!(p.prec_ids, prec_id)

struct SimplePrecursor
    sequence::String
    charge::UInt8
    isotope::UInt8
    pep_id::UInt32
end

mutable struct PrecursorTable
    id_to_prot::UnorderedDictionary{UInt32, Protein}
    id_to_pepGroup::UnorderedDictionary{UInt32, PeptideGroup}
    id_to_pep::UnorderedDictionary{UInt32, Peptide}
    id_to_prec::UnorderedDictionary{UInt32, SimplePrecursor}
    prot_to_id::UnorderedDictionary{String, UInt32}
    peptides::Set{String}
end

hasProtein(pt::PrecursorTable, prot::String) = haskey(pt.prot_to_id, prot)

function addProtein!(pt::PrecursorTable, prot::Protein, prot_id::UInt32) 
    insert!(pt.prot_to_id, getName(prot), prot_id)
    insert!(pt.id_to_prot, prot_id, prot)
end

getProtein(pt::PrecursorTable, prot_id::UInt32) = pt.id_to_prot[prot_id]
getProtID(pt::PrecursorTable, prot_name::String) = pt.prot_to_id[prot_name]
addPepGroupToProtein!(pt::PrecursorTable, prot_id::UInt32, pepGroup_id::UInt32) = push!(getPepGroupIDs(getProtein(pt, prot_id)), pepGroup_id)

getPepGroup(pt::PrecursorTable, pepGroup_id) = pt.id_to_pepGroup[pepGroup_id]
addPepGroup!(pt::PrecursorTable, pepGroup_id::UInt32, sequence::String) = insert!(pt.id_to_pepGroup, pepGroup_id, PeptideGroup(sequence))
addProtToPepGroup!(pt::PrecursorTable, pepGroup_id::UInt32, prot_id::UInt32) = addProtID!(getPepGroup(pt, pepGroup_id), prot_id)

getPep(pt::PrecursorTable, pep_id) = pt.id_to_pep[pep_id]
addPep!(pt::PrecursorTable, pep_id::UInt32, pepGroup_id::UInt32, sequence::String) = insert!(pt.id_to_pep, pep_id, Peptide(sequence, pepGroup_id, Set{UInt32}()))
addPrecToPep!(pt::PrecursorTable, pep_id::UInt32, prec_id::UInt32) = addPrecID!(getPep(pt, pep_id), prec_id)

PrecursorTable() = PrecursorTable(
    UnorderedDictionary{UInt32, Protein}(),
    UnorderedDictionary{UInt32, PeptideGroup}(),
    UnorderedDictionary{UInt32, Peptide}(),
    UnorderedDictionary{UInt32, SimplePrecursor}(),
    UnorderedDictionary{String, UInt32}(),
    Set{String}()
)

function getProtID(ptable::PrecursorTable, prot_name::String, max_prot_id::UInt32)
    if !hasProtein(ptable, prot_name)
        max_prot_id += UInt32(1)
        prot_id = max_prot_id 
        addProtein!(ptable, Protein(prot_name), prot_id)
    else
        prot_id = getProtID(ptable, prot_name)
    end
    return max_prot_id, prot_id
end

function fixedMods(peptide::String, fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}})
    for mod in fixed_mods
        peptide = replace(peptide, mod[:p]=>mod[:r])
    end
    peptide
end

function buildPrecursorTable!(ptable::PrecursorTable, peptides_fasta::Vector{FastaEntry},  
                                fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
                                var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
                                n::Int)
    max_prot_id, pepGroup_id, pep_id, prec_id = UInt32(1), UInt32(1), UInt32(1), UInt32(1)
    
    #Most recently encoutnered peptide sequence
    #It is assumed that `peptides_fasta` is sorted by the `sequence` field
    #so duplicated peptides are adjacent in the list
    previous_peptide = ""
    for peptide in peptides_fasta
        #Current peptide differs from previous
        if getSeq(peptide) != previous_peptide

            #Apply fixed modifications
            peptide_sequence_fixed_mods = fixedMods(getSeq(peptide), fixed_mods)

            pepGroup_id += UInt32(1)
            addPepGroup!(ptable, pepGroup_id, peptide_sequence_fixed_mods)

            max_prot_id, prot_id = getProtID(ptable, peptide_sequence_fixed_mods, max_prot_id)
            addPepGroupToProtein!(ptable, prot_id, pepGroup_id)
            
            addProtToPepGroup!(ptable, pepGroup_id, prot_id)

            previous_peptide = getSeq(peptide)

            ######
            #Apply Variable modifications
            #######

        else #Duplicated peptide, so add a protein to the peptide group

            #Apply fixed modifications
            peptide_sequence_fixed_mods = fixedMods(getSeq(peptide), fixed_mods)

            max_prot_id, prot_id = getProtID(ptable, peptide_sequence_fixed_mods, max_prot_id)
            addPepGroupToProtein!(ptable, prot_id, pepGroup_id)
            addProtToPepGroup!(ptable, pepGroup_id, prot_id)
        end
    end
end

test_table = PrecursorTable()
fixed_mods = [(p=r"C", r="C[Carb]")]
var_mods = [(p=r"(K$)", r="[Hlys]"), (p=r"(R$)", r="[Harg]")]
@time buildPrecursorTable!(test_table, peptides_fasta, fixed_mods, var_mods, 2)
#=
println(i)

precs = Vector{Precursor{Float64}}()
@time begin
for seq in tryptic
    push!(precs, Precursor(String(seq)))
end

ions = Float64[]
@time begin
for seq in tryptic
    append!(ions, getIonSeries(getResidues(String(seq)), UInt8(1)))
    append!(ions, getIonSeries(getResidues(String(seq)), UInt8(2)))
end
end
test = "MPLSLFRRVLLAVLLLVIIWTLFGPSGLGEELLSLSLASLLPAPASPGPPLALPRLLISNSHACGGSGPPPFLLILVCTAPEHLNQRNAIRASWGAIREARGFRVQTLFLLGKPRRQQLADLSSESAAHRDILQASFQDSYRNLTLKTLSGLNWVNKYCPMARYILKTDDDVYVNVPELVSELIQRGGPSEQWQKGKEAQEETTAIHEEHRGQAVPLLYLGRVHWRVRPTRTPESRHHVSEELWPENWGPFPPYASGTGYVLSISAVQLILKVASRAPPLPLEDVFVGVSARRGGLAPTHCVKLAGATHYPLDRCCYGKFLLTSHKVDPWQMQEAWKLVSGMNGERTAPFCSWLQGFLGTLRCRFIAWFSS"

function digest(sequence, regex = r"[KR][^P]", max_length = 40, min_length = 8)
    [x for x in split(sequence, regex) if ((length(x) >= min_length) & (length(x) <= max_length))]
end

test = "AAAKNNNKNNN"
[first(x) for x in findall(r"[KR][^P]", test)]
=#