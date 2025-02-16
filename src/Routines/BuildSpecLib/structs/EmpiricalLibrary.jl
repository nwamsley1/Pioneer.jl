abstract type EmpiricalLibrary end

"""
    BasicEmpiricalLibrary

A type for storing and accessing empirical spectral library data.

# Fields
- `libdf::DataFrame`: DataFrame containing spectral library data with standardized column names

Essential columns are renamed to internal standard names:
- `PrecursorMz` → `prec_mz`
- `Tr_recalibrated` → `irt`
"""
struct BasicEmpiricalLibrary <: EmpiricalLibrary
    libdf::DataFrame
    
    # Internal constructor to standardize column names
    function BasicEmpiricalLibrary(df::DataFrame)
        # Check if required columns exist
        required_cols = ["PrecursorMz", "Tr_recalibrated", "ModifiedPeptide", "PrecursorCharge"]
        missing_cols = setdiff(required_cols, names(df))
        if !isempty(missing_cols)
            throw(ArgumentError("Missing required columns: $(join(missing_cols, ", "))"))
        end
        
        # Create a copy and rename columns
        new_df = copy(df)
        rename!(new_df, 
            "PrecursorMz" => "prec_mz",
            "Tr_recalibrated" => "irt",
            "PeptideSequence" => "sequence",
            "ModifiedPeptide" => "modified_sequence",
            "PrecursorCharge" => "prec_charge",
            "IonMobility" => "ion_mobility",
            "ProductMz" => "frag_mz",
            "LibraryIntensity" => "library_intensity",
            "FragmentCharge" => "frag_charge",
            "FragmentType" => "frag_type",
            "FragmentSeriesNumber" => "frag_series_number",
            "FragmentLossType" => "frag_loss_type",
            "ProteinGroup" => "protein_group",
            "ProteinName" => "protein_name",
            "Genes" => "genes"
        )
        
        # Create precursor_idx based on unique ModifiedPeptide + PrecursorCharge combinations
        # Create a temporary column combining peptide and charge
        new_df.precursor_key = new_df.modified_sequence .* "_" .* string.(new_df.prec_charge)

        # Create dictionary mapping unique combinations to indices
        unique_precursors = unique(new_df.precursor_key)
        precursor_dict = Dict(key => UInt32(i) for (i, key) in enumerate(unique_precursors))
        
        # Add precursor_idx column
        new_df.precursor_idx = UInt32[precursor_dict[key] for key in new_df.precursor_key]
        
        # Remove temporary key column
        select!(new_df, Not(:precursor_key))
        
        # Reorder columns to put precursor_idx first
        new_df = select(new_df, :precursor_idx, Not(:precursor_idx))

        new(new_df)
    end
end

# Getter methods
"""
    getDF(lib::BasicEmpiricalLibrary)

Get the full DataFrame from the library.
"""
getDF(lib::BasicEmpiricalLibrary) = lib.libdf

"""
    getPrecursorMz(lib::BasicEmpiricalLibrary)

Get the precursor m/z values from the library.
"""
getPrecursorMz(lib::BasicEmpiricalLibrary) = lib.libdf.prec_mz

"""
    getIRT(lib::BasicEmpiricalLibrary)

Get the retention time values from the library.
"""
getIRT(lib::BasicEmpiricalLibrary) = lib.libdf.irt

# Extend Base.sort! for BasicEmpiricalLibrary
"""
    sort!(lib::BasicEmpiricalLibrary, cols)

Sort the library in-place by the specified column(s).

# Arguments
- `lib::BasicEmpiricalLibrary`: The library to sort
- `cols`: Column name(s) or sort specification

Note: This modifies the library in-place.
"""
function Base.sort!(lib::BasicEmpiricalLibrary, cols)
    sort!(lib.libdf, cols)
    return lib
end

"""
    sort(lib::BasicEmpiricalLibrary, cols)

Create a sorted copy of the library.

# Arguments
- `lib::BasicEmpiricalLibrary`: The library to sort
- `cols`: Column name(s) or sort specification

Returns a new BasicEmpiricalLibrary with sorted data.
"""
function Base.sort(lib::BasicEmpiricalLibrary, cols)
    new_df = sort(lib.libdf, cols)
    return BasicEmpiricalLibrary(new_df)
end

# Example usage for creating from CSV:
"""
    BasicEmpiricalLibrary(csv_file::String)

Create a BasicEmpiricalLibrary from a CSV file.

# Arguments
- `csv_file::String`: Path to the CSV file
"""
function BasicEmpiricalLibrary(csv_file::String)
    df = CSV.read(csv_file, DataFrame)
    return BasicEmpiricalLibrary(df)
end

function getProteomeId(sl::EmpiricalLibrary, frag_idx::Integer)
    proteome_ids = sl.libdf[!,:proteome_ids]
    function proteome_id(proteome_ids::Vector{Missing}, frag_idx::Integer)
        return ""
    end
    function proteome_id(proteome_ids::Vector{String}, frag_idx::Integer)
        return proteome_ids[frag_idx]
    end
    return proteome_id(proteome_ids, frag_idx)
end

function getProteinGroupId(sl::EmpiricalLibrary, frag_idx::Integer)
    protein_groups = sl.libdf[!,:protein_group]
    function proteome_id(protein_group::Vector{Missing}, frag_idx::Integer)
        return ""
    end
    function proteome_id(protein_group::Vector{String}, frag_idx::Integer)
        return protein_group[frag_idx]
    end
    return proteome_id(protein_groups, frag_idx)
end

function getCollisionEnergy(sl::EmpiricalLibrary, frag_idx::Integer)
    collision_energies = sl.libdf[!,:collision_energy]
    function _collision_energy(collision_energy::Vector{Missing}, frag_idx::Integer)
        return zero(Float32)
    end
    function _collision_energy(collision_energy::Vector{Float32}, frag_idx::Integer)
        return collision_energy[frag_idx]
    end
    return collision_energies(proteome_idx, frag_idx)
end

function getCollisionEnergy(sl::EmpiricalLibrary, frag_idx::Integer)
    collision_energies = sl.libdf[!,:collision_energy]
    function _collision_energy(collision_energy::Vector{Missing}, frag_idx::Integer)
        return zero(Float32)
    end
    function _collision_energy(collision_energy::Vector{Float32}, frag_idx::Integer)
        return collision_energy[frag_idx]
    end
    return collision_energies(proteome_idx, frag_idx)
end

function getMissedCleavages(sl::EmpiricalLibrary, frag_idx::Integer)
    missed_cleavages = sl.libdf[!,:missed_cleavages]
    function _missed_cleavage(missed_cleavage::Vector{Missing}, frag_idx::Integer)
        return zero(Float32)
    end
    function _missed_cleavage(missed_cleavage::Vector{Float32}, frag_idx::Integer)
        return missed_cleavage[frag_idx]
    end
    return missed_cleavages(proteome_idx, frag_idx)
end

getSequence(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[!,:sequence]::Vector{String}[frag_idx]
getPrecCharge(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[!,:prec_charge]::Vector{UInt8}[frag_idx]
getIsDecoy(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[!,:is_decoy]::Vector{Bool}[frag_idx]
function getEntrapmentGroupIdx(sl::EmpiricalLibrary, frag_idx::Integer)
    sl.libdf[!,:entrapment_group_id]::Vector{UInt8}[frag_idx]
end
getPrecMz(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[!,:prec_mz]::Vector{Float32}[frag_idx]
getSeqLength(sl::EmpiricalLibrary, frag_idx::Integer) = length(getSequence(sl, frag_idx))
getIrt(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[!,:irt]::Vector{Float32}[frag_idx]
getPrecursorSulfurCount(sl::EmpiricalLibrary, frag_idx::Integer) = zero(UInt8) #need to implement #sl.libdf[!,:sulfur_count]::Vector{UInt8}[frag_idx]
getFragMz(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[!,:frag_mz]::Vector{Float32}[frag_idx]
getIonType(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[!,:ion_type]::Vector{String}[frag_idx]
getIsY(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[!,:frag_type]::Vector{Char}[frag_idx] == 'y'
getIsB(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[!,:frag_type]::Vector{Char}[frag_idx] == 'b'
getIsP(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[!,:frag_type]::Vector{Char}[frag_idx] == 'p'
function getAXCZ(sl::EmpiricalLibrary, frag_idx::Integer)
    frag_type = sl.libdf[!,:frag_type]::Vector{Char}[frag_idx]
    if frag_type == 'a'
        return true
    elseif frag_type == 'x'
        return true
    elseif frag_type == 'c'
        return true
    elseif frag_type == 'z'
        return true
    end
end
getNeutralDiff(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[!,:neutral_diff]::Vector{Bool}[frag_idx]
getFragSeriesNumber(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[!,:frag_series_number]::Vector{UInt8}[frag_idx]
getFragCharge(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[!,:frag_charge]::Vector{UInt8}[frag_idx]
isInternal(sl::EmpiricalLibrary, frag_idx::Integer) = false #Eventually should implement this
isImmonium(sl::EmpiricalLibrary, frag_idx::Integer) = false #Eventually should implement this
getInternalInd(sl::EmpiricalLibrary, frag_idx::Integer) = (zero(UInt8), zero(UInt8)) #Eventually should implement this
getFragSulfurCount(sl::EmpiricalLibrary, frag_idx::Integer) = zero(UInt8) #Need to implement 



#getMissedCleavages(sl::EmpiricalLibrary, frag_idx::Integer)