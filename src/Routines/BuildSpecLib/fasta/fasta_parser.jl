"""
    parse_fasta(fasta_path::String, proteome_id::String)::Vector{FastaEntry}

Parse a FASTA file and return a vector of FastaEntry objects.
    
# Parameters
- `fasta_path::String`: Path to FASTA file (.fasta or .fasta.gz)
- `proteome_id::String`: Identifier for the proteome

# Returns
- `Vector{FastaEntry}`: Vector of parsed FASTA entries

# Details
The function:
1. Automatically detects and handles compressed (.fasta.gz) or uncompressed (.fasta) files
2. Parses the headers to extract identifiers
3. Creates FastaEntry objects for each sequence with the specified proteome ID
4. Sets default values for optional fields (modifications, charge, etc.)

The identifier parsing attempts to extract:
- UniProt-style identifiers when headers contain pipe symbols (e.g., ">sp|P12345|GENE_HUMAN")
- The first word before any space in other cases

# Examples
```julia
# Parse a standard FASTA file
human_entries = parse_fasta("human_proteome.fasta", "human")

# Parse a compressed FASTA file
mouse_entries = parse_fasta("mouse_proteome.fasta.gz", "mouse")
```

# Notes
The resulting FastaEntry objects will have default values for fields other than sequence,
identifier, and proteome_id. These can be modified later with functions like add_mods().
"""
function parse_fasta(fasta_path::String, 
                    proteome_id::String,
                    #parse_identifier::Function = x -> split(x,"|")[2]
                    )::Vector{FastaEntry}
    
    function parse_identifier(header::AbstractString)
        if count("|", header) == 2
            return split(header,"|")[2]
        else
            return first(split(header, " "))
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
            one(UInt32),
            missing,
            missing,
            zero(UInt8),
            zero(UInt32),
            zero(UInt32),
            zero(UInt8),
            false
        ))
    end

    return entries
end
