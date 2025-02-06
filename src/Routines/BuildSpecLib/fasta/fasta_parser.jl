"""
Parse a FASTA file and return a vector of FastaEntry objects.
    
Parameters:
- fasta_path::String: Path to FASTA file (.fasta or .fasta.gz)
- proteome_id::String: Identifier for the proteome
- parse_identifier::Function: Function to parse entry identifiers

Returns:
- Vector{FastaEntry}: Vector of parsed FASTA entries
"""
function parse_fasta(fasta_path::String, 
                    proteome_id::String,
                    #parse_identifier::Function = x -> split(x,"|")[2]
                    )::Vector{FastaEntry}
    
    function parse_identifier(header::String)
        if occursin("|", x)
            return split(x,"|")[2]
        else
            return first(split(x, " "))
        end
    end
    
    function get_reader(fasta_path::String)
        if endswith(fasta_path, ".fasta.gz")
            return FASTA.Reader(GzipDecompressorStream(open(fasta_path)))
        elseif endswith(fasta_path, ".fasta")
            return FASTA.Reader(open(fasta_path))
        else
            throw(ErrorException("Invalid file extension: must be .fasta or .fasta.gz"))
        end
    end

    reader = get_reader(fasta_path)
    entries = Vector{FastaEntry}()

    for record in reader
        push!(entries, FastaEntry(
            parse_identifier(FASTA.identifier(record)),
            FASTA.description(record),
            proteome_id,
            FASTA.sequence(record),
            zero(UInt32),
            zero(UInt8),
            false
        ))
    end

    return entries
end
