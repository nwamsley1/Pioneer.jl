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
        
        # Add proteome_idx column if missing
        if !hasproperty(new_df, :ProteomeIdentifier)
            new_df[!,:proteome_idx] = zeros(Missing, nrow(new_df))
        else
            rename!(new_df, 
                "ProteomeIdentifier" => "proteome_idx"
            )
        end

        # Add collision_energy  column if missing
        if !hasproperty(new_df, :CollisionEnergy)
            new_df[!,:collision_energy] = zeros(Missing, nrow(new_df))
        else
            rename!(new_df, 
                "CollisionEnergy" => "collision_energy"
            )
        end

        # Add collision_energy  column if missing
        if !hasproperty(new_df, :MissedCleavages)
            new_df[!,:missed_cleavages] = zeros(Missing, nrow(new_df))
        else
            rename!(new_df, 
                "MissedCleavages" => "missed_cleavages"
            )
        end
        
        # Convert columns to the correct types
        new_df[!,:prec_charge] = convert(Vector{UInt8}, new_df[!,:prec_charge])
        new_df[!,:frag_charge] = convert(Vector{UInt8}, new_df[!,:frag_charge])
        new_df[!,:frag_series_number] = convert(Vector{UInt8}, new_df[!,:frag_series_number])
        new_df[!,:prec_mz] = convert(Vector{Float32}, new_df[!,:prec_mz])
        new_df[!,:frag_mz] = convert(Vector{Float32}, new_df[!,:frag_mz])
        new_df[!,:irt] = convert(Vector{Float32}, new_df[!,:irt])
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

##########
#Precursor Getter Methods
##########
function getProteomeId(sl::EmpiricalLibrary, frag_idx::Integer)
    proteome_idx = sl.libdf[!,:proteome_idx]
    function proteome_id(proteome_idx::AbstractVector{Missing}, frag_idx::Integer)
        return ""
    end
    function proteome_id(proteome_idx::AbstractVector{String}, frag_idx::Integer)
        return proteome_idx[frag_idx]
    end
    return proteome_id(proteome_idx, frag_idx)
end

function getProteinGroupId(sl::EmpiricalLibrary, frag_idx::Integer)
    protein_groups = sl.libdf[!,:protein_group]
    function proteome_id(protein_group::AbstractVector{Missing}, frag_idx::Integer)
        return ""
    end
    function proteome_id(protein_group::AbstractVector{String}, frag_idx::Integer)
        return protein_group[frag_idx]
    end
    return proteome_id(protein_groups, frag_idx)
end

function getCollisionEnergy(sl::EmpiricalLibrary, frag_idx::Integer)
    collision_energies = sl.libdf[!,:collision_energy]
    function _collision_energy(collision_energy::AbstractVector{Missing}, frag_idx::Integer)
        return zero(Float32)
    end
    function _collision_energy(collision_energy::AbstractVector{Float32}, frag_idx::Integer)
        return collision_energy[frag_idx]
    end
    return _collision_energy(collision_energies, frag_idx)
end

function getCollisionEnergy(sl::EmpiricalLibrary, frag_idx::Integer)
    collision_energies = sl.libdf[!,:collision_energy]
    function _collision_energy(collision_energy::AbstractVector{Missing}, frag_idx::Integer)
        return zero(Float32)
    end
    function _collision_energy(collision_energy::AbstractVector{Float32}, frag_idx::Integer)
        return collision_energy[frag_idx]
    end
    return collision_energies(proteome_idx, frag_idx)
end

function getMissedCleavages(sl::EmpiricalLibrary, frag_idx::Integer)
    missed_cleavages = sl.libdf[!,:missed_cleavages]
    function _missed_cleavage(missed_cleavage::AbstractVector{Missing}, frag_idx::Integer)
        return zero(Float32)
    end
    function _missed_cleavage(missed_cleavage::AbstractVector{Float32}, frag_idx::Integer)
        return missed_cleavage[frag_idx]
    end
    return _missed_cleavage(missed_cleavages, frag_idx)
end

#Placeholder functions
function parseStructuralMods(sl::EmpiricalLibrary, frag_idx::Integer)
    return ""
end
function parseIsotopicMods(sl::EmpiricalLibrary, frag_idx::Integer)
    return ""
end
function getSulfurCount(sl::EmpiricalLibrary, frag_idx::Integer)
    return zero(UInt8)
end
getSequence(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[frag_idx,:sequence]::String
getPrecCharge(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[frag_idx,:prec_charge]::UInt8
getIsDecoy(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[frag_idx,:is_decoy]::Bool
function getEntrapmentGroupIdx(sl::EmpiricalLibrary, frag_idx::Integer)
    sl.libdf[frag_idx,:entrapment_group_id]::UInt8
end
getPrecMz(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[frag_idx,:prec_mz]::Float32
getSeqLength(sl::EmpiricalLibrary, frag_idx::Integer) = length(getSequence(sl, frag_idx))
getIrt(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[frag_idx,:irt]::Float32
#getMissedCleavages(sl::EmpiricalLibrary, frag_idx::Integer)
##########
#Fragment Getter Methods
##########
getFragMz(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[frag_idx, :frag_mz]::Float32

getFragIntensity(sl::EmpiricalLibrary, frag_idx::Integer) = Float16(sl.libdf[frag_idx, :library_intensity])

function getIonType(sl::EmpiricalLibrary, frag_idx::Integer)
    #frag_type = sl.libdf[frag_idx, :frag_type]::String
    #return UInt16(get(ION_TYPE_DICT, frag_type, 0))
    return zero(UInt16)
end

getIsY(sl::EmpiricalLibrary, frag_idx::Integer) = startswith(sl.libdf[frag_idx, :frag_type], "y")::Bool

getIsB(sl::EmpiricalLibrary, frag_idx::Integer) = startswith(sl.libdf[frag_idx, :frag_type], "b")::Bool

getIsP(sl::EmpiricalLibrary, frag_idx::Integer) = startswith(sl.libdf[frag_idx, :frag_type], "p")::Bool

function getAXCZ(sl::EmpiricalLibrary, frag_idx::Integer)
    ftype = first(sl.libdf[frag_idx, :frag_type])
    return (ftype == 'a' || ftype == 'x' || ftype == 'c' || ftype == 'z')::Bool
end

function getNeutralDiff(sl::EmpiricalLibrary, frag_idx::Integer)
    if !hasproperty(sl.libdf, :frag_loss_type)
        return false
    end
    loss = sl.libdf[frag_idx, :frag_loss_type]
    return (!ismissing(loss) && !isempty(loss))::Bool
end

getFragIndex(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[frag_idx, :frag_series_number]::UInt8

getFragCharge(sl::EmpiricalLibrary, frag_idx::Integer) = sl.libdf[frag_idx, :frag_charge]::UInt8

function isInternal(sl::EmpiricalLibrary, frag_idx::Integer)
    #ftype = sl.libdf[frag_idx, :frag_type]::String
    #return startswith(ftype, "Int")::Bool
    return false
end

function isImmonium(sl::EmpiricalLibrary, frag_idx::Integer)
    #ftype = sl.libdf[frag_idx, :frag_type]::String
    #return startswith(ftype, "Imm")::Bool
    return false
end

function getInternalInd(sl::EmpiricalLibrary, frag_idx::Integer)
    #=
    if !isInternal(sl, frag_idx)
        return (zero(UInt8), zero(UInt8))
    end
    ftype = sl.libdf[frag_idx, :frag_type]::String
    m = match(r"Int/(\d+)-(\d+)", ftype)
    if isnothing(m)
        return (zero(UInt8), zero(UInt8))
    end
    return (parse(UInt8, m[1]), parse(UInt8, m[2]))
    =#
    return (zero(UInt8), zero(UInt8))
end

getFragSulfurCount(sl::EmpiricalLibrary, frag_idx::Integer) = zero(UInt8)  # Placeholder implementation