using FASTX
using CodecZlib
using Dictionaries
using Dictionaries
using Combinatorics
using Random
using Arrow
using Tables

#=include("src/precursor.jl")
include("src/parseFASTA.jl")
include("src/PrecursorDatabase.jl")
include("src/applyMods.jl")=#


#=

peptides_fasta = digestFasta(parseFasta("/Users/n.t.wamsley/Projects/TEST_DATA/proteomes/UP000000589_10090.fasta.gz"), max_length = 30, min_length = 7)
#peptides_fasta = digestFasta(parseFasta(file_path))
test_table = PrecursorTable()
fixed_mods = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}()
var_mods = [(p=r"(M)", r="[MOx]")]
buildPrecursorTable!(test_table, peptides_fasta, fixed_mods, var_mods, 2)
const charge_facs = Float32[1, 0.9, 0.85, 0.8, 0.75]

function adjustNCE(NCE::T, default_charge::Integer, peptide_charge::Integer) where {T<:AbstractFloat}
    return NCE*(charge_facs[default_charge]/charge_facs[peptide_charge])
end

=#

#buildPrositCSV("/Users/n.t.wamsley/Projects/TEST_DATA/proteomes/UP000000589_10090.fasta.gz", 
#"/Users/n.t.wamsley/Desktop/prosit_mouse_NCE35_correctedNCE_073123.csv", min_length = 7, max_length = 30, nce = 35.0, dynamic_nce = true)
#=buildPrositCSV("/Users/n.t.wamsley/Projects/TEST_DATA/proteomes/uniprotkb_proteome_UP000000589_AND_mode_2023_08_28.fasta.gz",  
                "/Users/n.t.wamsley/Projects/TEST_DATA/proteomes/prositInput/prositMouseIsoformsNCE33_corrected_082823.csv",
                min_length = 7, 
                max_length = 30, 
                nce = 33.0, 
                dynamic_nce = true
                )=#



function buildPrositCSV(fasta::Union{String, PrecursorTable}, f_out::String; min_length::Int = 8, max_length::Int = 30,
                        min_charge::Int = 2, max_charge::Int = 4,
                        fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}} = Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}(), 
                        var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}} = [(p=r"(M)", r="[MOx]")],
                        n_var_mods::Int = 2,
                        nce::Float64 = 30.0,
                        default_charge::Int = 3,
                        dynamic_nce::Bool = true)
    precursor_table = fasta
    if typeof(precursor_table) == String
        precursor_table = PrecursorTable()
        buildPrecursorTable!(precursor_table, fasta,
                            fixed_mods, var_mods, n_var_mods)
    end
    
    open(f_out, "w") do file
        write(file, "accession_number,modified_sequence,collision_energy,precursor_charge,prot_ids,pep_id,decoy\n")
        for (id, pep) in ProgressBar(Base.pairs(precursor_table.id_to_pep))
            sequence = replace(getSeq(pep), r"M\[MOx\]"=>"M(ox)")
            sequence = replace(sequence, r"C\[Carb\]"=>"C")
            unmod_sequence = replace(sequence, r"M(ox)"=>"M")
            prot_id = join([prot_id for prot_id in collect(getProtFromPepID(precursor_table, id))],";")
            accession = join([getName(getProtein(precursor_table, prot_id)) for prot_id in collect(getProtFromPepID(precursor_table , id))],";")
            #Check for illegal amino acid characters
            if (occursin("[H", sequence)) | (occursin("U", sequence)) | (occursin("O", sequence)) |  (occursin("X", sequence)) | occursin("Z", getSeq(pep)) | occursin("B", getSeq(pep))
                continue
            end
            #Enforce length constraints
            #There is a bug here because the sequence length includes "M(ox)". 
            #So the maximum length of a methionine oxidized peptide is actually 30 - 4. 
            if (length(unmod_sequence) >max_length) | (length(unmod_sequence) < min_length)
                continue
            end 
            decoy = isDecoy(pep)
            #if (decoy == false)
            for charge in range(min_charge, max_charge)
                if dynamic_nce
                    NCE = adjustNCE(nce, default_charge, charge)
                    write(file, "$accession,$sequence,$NCE,$charge,$prot_id,$id,$decoy\n")
                else
                    write(file, "$accession,$sequence,$nce,$charge,$prot_id,$id,$decoy\n")
                end
            end
            #end
        end
    end

end

#=open("/Users/n.t.wamsley/Desktop/targets.csv", "w") do file
    write(file, "accession_number,modified_sequence,collision_energy,precursor_charge,prot_ids,pep_id,decoy\n")
    for (id, pep) in ProgressBar(pairs(test_table.id_to_pep))
        sequence = replace(getSeq(pep), r"M\[MOx\]"=>"M(ox)")
        sequence = replace(sequence, r"C\[Carb\]"=>"C")
        prot_id = join([prot_id for prot_id in collect(getProtFromPepID(test_table, id))],";")
        accession = join([getName(getProtein(test_table, prot_id)) for prot_id in collect(getProtFromPepID(test_table, id))],";")
        #Check for illegal amino acid characters
        if (occursin("[H", sequence)) | (occursin("U", sequence)) | (occursin("O", sequence)) |  (occursin("X", sequence)) | occursin("Z", getSeq(pep)) | occursin("B", getSeq(pep))
            continue
        end
        #Enforce length constraints
        if (length(sequence) > 30) | (length(sequence) < 8)
            continue
        end 
        decoy = isDecoy(pep)
        #if (decoy == false)
        for charge in [2, 3, 4]
            NCE = adjustNCE(30.0, 3, charge)
            write(file, "$accession,$sequence,$NCE,$charge,$prot_id,$id,$decoy\n")
        end
        #end
    end
end=#