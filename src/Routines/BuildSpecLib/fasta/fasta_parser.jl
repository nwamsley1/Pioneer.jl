# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
   parse_fasta(fasta_path::String, proteome_id::String; accession_regex=nothing,
                 gene_regex=nothing, protein_regex=nothing,
                 organism_regex=nothing)::Vector{FastaEntry}

Parse a FASTA file and return a vector of FastaEntry objects.
    
# Parameters
- `fasta_path::String`: Path to FASTA file (.fasta or .fasta.gz)
- `proteome_id::String`: Identifier for the proteome

# Returns
- `Vector{FastaEntry}`: Vector of parsed FASTA entries

# Details
The function:
1. Automatically detects and handles compressed (.fasta.gz) or uncompressed (.fasta) files
2. Optionally parses headers for accession, gene, protein and organism using the provided regexes
3. Creates FastaEntry objects for each sequence with the specified proteome ID
4. Sets default values for optional fields (modifications, charge, etc.)

# Examples
```julia
# Parse a standard FASTA file
human_entries = parse_fasta("human_proteome.fasta", "human")

# Parse a compressed FASTA file
mouse_entries = parse_fasta("mouse_proteome.fasta.gz", "mouse")
```

# Notes
The resulting `FastaEntry` objects will contain the raw header in the `description`
field and parsed metadata when the corresponding regular expressions are provided.
Fields not parsed default to empty strings. Modifications and charge can be added
later with functions like `add_mods()`.
"""
function parse_fasta(
    fasta_path::String,
    proteome_id::String;
    accession_regex::Union{Regex,Nothing}=nothing,
    gene_regex::Union{Regex,Nothing}=nothing,
    protein_regex::Union{Regex,Nothing}=nothing,
    organism_regex::Union{Regex,Nothing}=nothing,
)::Vector{FastaEntry}
    
    function get_reader(fasta_path::String)
        # Try to open based on extension, but handle mismatched compression
        if endswith(fasta_path, ".gz")
            try
                return FASTA.Reader(GzipDecompressorStream(open(fasta_path)))
            catch e
                if isa(e, ZlibError)
                    # File has .gz extension but isn't actually gzipped
                    return FASTA.Reader(open(fasta_path))
                else
                    rethrow(e)
                end
            end
        elseif endswith(fasta_path, ".fasta") || endswith(fasta_path, ".fa") || endswith(fasta_path, ".fna") || endswith(fasta_path, ".faa")
            # Try as plain text for common FASTA extensions
            return FASTA.Reader(open(fasta_path))
        else
            # Unknown extension - try plain text first, then gzipped
            try
                return FASTA.Reader(open(fasta_path))
            catch
                try
                    return FASTA.Reader(GzipDecompressorStream(open(fasta_path)))
                catch
                    throw(ErrorException("Unable to read file $fasta_path - not a valid FASTA or gzipped FASTA file"))
                end
            end
        end
    end

    reader = get_reader(fasta_path)
    entries = Vector{FastaEntry}()

    for record in reader
        identifier = FASTA.identifier(record)
        description = FASTA.description(record)

        accession = if accession_regex === nothing
            identifier
        else
            m = match(accession_regex, identifier)
            m === nothing ? identifier : String(first(m.captures))
        end

        gene = if gene_regex === nothing
            ""
        else
            m = match(gene_regex, description)
            m === nothing ? "" : String(first(m.captures))
        end

        protein = if protein_regex === nothing
            ""
        else
            m = match(protein_regex, description)
            m === nothing ? "" : String(first(m.captures))
        end

        organism = if organism_regex === nothing
            ""
        else
            m = match(organism_regex, description)
            m === nothing ? "" : String(first(m.captures))
        end

        push!(entries, FastaEntry(
            accession,
            FASTA.description(record),
            gene,
            protein,
            organism,
            proteome_id,
            FASTA.sequence(record),
            one(UInt32),
            missing,
            missing,
            zero(UInt8),
            zero(UInt32),  # base_seq_id (will be assigned later)
            zero(UInt32),  # base_entrap_id (will be assigned later)
            zero(UInt32),  # base_pep_id (will be assigned later)
            zero(UInt32),  # base_prec_id (will be assigned later)
            zero(UInt8),   # entrapment_pair_id
            false,         # is_decoy
        ))
    end

    return entries
end
