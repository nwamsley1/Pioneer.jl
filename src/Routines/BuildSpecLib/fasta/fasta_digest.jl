# src/fasta/fasta_digest.jl

"""
Enzymatically digest protein sequences into peptides.

Parameters:
- fasta::Vector{FastaEntry}: FASTA entries to digest
- proteome_id::String: Proteome identifier 
- regex::Regex: Enzyme cleavage pattern
- max_length::Int: Maximum peptide length
- min_length::Int: Minimum peptide length
- missed_cleavages::Int: Maximum missed cleavages allowed

Returns:
- Vector{FastaEntry}: Digested peptide entries
"""
function digest_fasta(fasta::Vector{FastaEntry},
                     proteome_id::String;
                     regex::Regex = r"[KR][^P|$]",
                     max_length::Int = 40,
                     min_length::Int = 8,
                     missed_cleavages::Int = 1)::Vector{FastaEntry}

    function digest_sequence(sequence::AbstractString,
                           regex::Regex,
                           max_length::Int,
                           min_length::Int,
                           missed_cleavages::Int)::Vector{String}  # Return String not SubString
        
        function add_peptide!(peptides::Vector{SubString{String}},
                            n::Int,
                            sequence::AbstractString,
                            site::Int,
                            previous_sites::Vector{Int},
                            min_length::Int,
                            max_length::Int,
                            missed_cleavages::Int)
            
            for i in 1:min(n, missed_cleavages + 1)
                previous_site = previous_sites[end - i + 1]
                if ((site - previous_site) >= min_length) && 
                   ((site - previous_site) <= max_length)
                    # Convert SubString to String when adding to peptides
                    push!(peptides, String(@view sequence[previous_site+1:site]))
                end
            end

            for i in 1:length(previous_sites)-1
                previous_sites[i] = previous_sites[i+1]
            end
            previous_sites[end] = site
            return n + 1
        end

        peptides = Vector{SubString{String}}()
        previous_sites = zeros(Int, missed_cleavages + 1)
        previous_sites[1] = 0
        n = 1

        for site in eachmatch(regex, sequence, overlap = true)
            n = add_peptide!(peptides, n, sequence, site.offset,
                           previous_sites, min_length, max_length,
                           missed_cleavages)
        end

        # Handle C-terminal peptides
        n = add_peptide!(peptides, n, sequence, length(sequence),
                       previous_sites, min_length, max_length,
                       missed_cleavages)

        return peptides
    end
    peptides_fasta = Vector{FastaEntry}()
    base_pep_id = one(UInt32)
    base_prec_id = one(UInt32)
    for entry in fasta
        for peptide in digest_sequence(get_sequence(entry), regex,
                                     max_length, min_length,
                                     missed_cleavages)

            if (occursin("[H", peptide)) | (occursin("U", peptide)) | (occursin("O", peptide)) |  (occursin("X", peptide)) | occursin("Z", peptide) | occursin("B", peptide)
                continue
            end

            push!(peptides_fasta, FastaEntry(
                get_id(entry),
                "", # Skip description to save memory
                proteome_id,
                peptide,  # Now String instead of SubString
                missing, #structural_mods 
                missing, #istopic_mods 
                zero(UInt8),
                base_pep_id,
                base_prec_id,
                zero(UInt8),
                false
            ))
            base_prec_id += one(UInt32)
            base_pep_id += one(UInt32)
        end
    end

    return peptides_fasta
end