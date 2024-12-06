module MassSpecDataModule

using Arrow

# Abstract type for mass spectrometry data
abstract type MassSpecData end

# Parameterized struct for mass spectrometry data
struct BasicMassSpecData{T} <: MassSpecData
    data::Arrow.Table
    n::Int  # Number of rows in the table
end

# Constructor method for creating BasicMassSpecData from a file path
function BasicMassSpecData(file_path::String)
    table = Arrow.Table(file_path)
    
    # Extract the type inside the Union{Missing, Float32} (or similar)
    union_type = eltype(eltype(table[:mz_array]))
    float_type = first(filter(t -> t !== Missing, Base.uniontypes(union_type)))
    
    # Determine the length (number of rows)
    n = length(table[:mz_array])  # Use any column, as all have the same length
    
    return BasicMassSpecData{float_type}(table, n)
end

# Define length method
Base.length(ms_data::BasicMassSpecData) = ms_data.n

# Getter methods for individual elements
function mz_array(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:mz_array][scan_idx]
end

function intensity_array(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:intensity_array][scan_idx]
end

function scanHeader(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:scanHeader][scan_idx]
end

function scanNumber(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:scanNumber][scan_idx]
end

function basePeakMz(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:basePeakMz][scan_idx]
end

function basePeakIntensity(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:basePeakIntensity][scan_idx]
end

function packetType(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:packetType][scan_idx]
end

function retentionTime(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:retentionTime][scan_idx]
end

function lowMz(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:lowMz][scan_idx]
end

function highMz(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:highMz][scan_idx]
end

function TIC(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:TIC][scan_idx]
end

function centerMz(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:centerMz][scan_idx]
end

function isolationWidthMz(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:isolationWidthMz][scan_idx]
end

function collisionEnergyField(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:collisionEnergyField][scan_idx]
end

function collisionEnergyEvField(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:collisionEnergyEvField][scan_idx]
end

function msOrder(ms_data::BasicMassSpecData, scan_idx::Integer)
    ms_data.data[:msOrder][scan_idx]
end
