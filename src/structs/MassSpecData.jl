# Abstract type for mass spectrometry data
abstract type MassSpecData end

# Parameterized struct for mass spectrometry data with AbstractFloat constraint
struct BasicMassSpecData{T <: AbstractFloat} <: MassSpecData
    data::Arrow.Table
    n::Int  # Number of rows in the table
end

# Constructor method for creating BasicMassSpecData from a file path
function BasicMassSpecData(file_path::String)
    table = Arrow.Table(file_path)
    
    # Extract the type inside the Union{Missing, Float32/Float64}
    union_type = eltype(eltype(table[:mz_array]))
    float_type = first(filter(t -> t <: AbstractFloat, Base.uniontypes(union_type)))
    
    # Determine the length (number of rows)
    n = length(table[:mz_array])
    
    return BasicMassSpecData{float_type}(table, n)
end

# Define length method
Base.length(ms_data::BasicMassSpecData) = ms_data.n

# Generic type alias for array data that can handle different float types
#=
const MsDataArray{T} = SentinelArrays.ChainedVector{
    SubArray{
        Union{Missing, T}, 
        1, 
        Arrow.Primitive{Union{Missing, T}, Vector{T}}, 
        Tuple{UnitRange{Int64}}, 
        true
    }, 
    Arrow.List{
        SubArray{
            Union{Missing, T}, 
            1, 
            Arrow.Primitive{
                Union{Missing, T}, 
                Vector{T}
            }, 
            Tuple{UnitRange{Int64}}, 
            true
        }, 
        Int32, 
        Arrow.Primitive{
            Union{Missing, T}, 
            Vector{T}
        }
    }
}
    =#

# Generic SubArray type alias
const MsSubArray{T} = SubArray{
    Union{Missing, T}, 
    1, 
    Arrow.Primitive{Union{Missing, T}, Vector{T}}, 
    Tuple{UnitRange{Int64}}, 
    true
}

# Getter methods for individual elements
# MZ array getters
function getMzArrays(ms_data::BasicMassSpecData{T}) where T
    return ms_data.data[:mz_array]
end
function getMzArray(ms_data::BasicMassSpecData{T}, scan_idx::Integer)::MsSubArray{T} where T
    getMzArrays(ms_data)[scan_idx]
end

# Intensity array getters
function getIntensityArrays(ms_data::BasicMassSpecData{T}) where T
    return ms_data.data[:intensity_array]
end
function getIntensityArray(ms_data::BasicMassSpecData{T}, scan_idx::Integer)::MsSubArray{T} where T
    getIntensityArrays(ms_data)[scan_idx]
end

# String data doesn't need type parameter
getScanHeaders(ms_data::BasicMassSpecData)::SentinelArrays.ChainedVector{String, Arrow.List{String, Int32, Vector{UInt8}}} = ms_data.data[:scanHeader]
function getScanHeader(ms_data::BasicMassSpecData, scan_idx::Integer)::String
    getScanHeaders(ms_data)[scan_idx]
end

# Scan number getters
getScanNumbers(ms_data::BasicMassSpecData)::SentinelArrays.ChainedVector{Int32, Arrow.Primitive{Int32, Vector{Int32}}} = ms_data.data[:scanNumber]
function getScanNumber(ms_data::BasicMassSpecData, scan_idx::Integer)::Int32
    getScanNumbers(ms_data)[scan_idx]
end

# Base peak m/z getters
function getBasePeakMzs(ms_data::BasicMassSpecData{T}) where T
    return ms_data.data[:basePeakMz]::SentinelArrays.ChainedVector{T, Arrow.Primitive{T, Vector{T}}}
end
function getBasePeakMz(ms_data::BasicMassSpecData{T}, scan_idx::Integer)::T where T
    getBasePeakMzs(ms_data)[scan_idx]
end

# Base peak intensity getters
function getBasePeakIntensities(ms_data::BasicMassSpecData{T}) where T
    return ms_data.data[:basePeakIntensity]::SentinelArrays.ChainedVector{T, Arrow.Primitive{T, Vector{T}}}
end
function getBasePeakIntensity(ms_data::BasicMassSpecData{T}, scan_idx::Integer)::T where T
    getBasePeakIntensities(ms_data)[scan_idx]
end

# Packet type getters remain Int32
getPacketTypes(ms_data::BasicMassSpecData)::SentinelArrays.ChainedVector{Int32, Arrow.Primitive{Int32, Vector{Int32}}} = ms_data.data[:packetType]
function getPacketType(ms_data::BasicMassSpecData, scan_idx::Integer)::Int32
    getPacketTypes(ms_data)[scan_idx]
end

# Retention time getters
function getRetentionTimes(ms_data::BasicMassSpecData{T}) where T
    return ms_data.data[:retentionTime]::SentinelArrays.ChainedVector{T, Arrow.Primitive{T, Vector{T}}}
end
function getRetentionTime(ms_data::BasicMassSpecData{T}, scan_idx::Integer)::T where T
    getRetentionTimes(ms_data)[scan_idx]
end

# Low m/z getters
function getLowMzs(ms_data::BasicMassSpecData{T}) where T
    return ms_data.data[:lowMz]::SentinelArrays.ChainedVector{T, Arrow.Primitive{T, Vector{T}}}
end
function getLowMz(ms_data::BasicMassSpecData{T}, scan_idx::Integer)::T where T
    getLowMzs(ms_data)[scan_idx]
end

# High m/z getters
function getHighMzs(ms_data::BasicMassSpecData{T}) where T
    return ms_data.data[:highMz]::SentinelArrays.ChainedVector{T, Arrow.Primitive{T, Vector{T}}}
end
function getHighMz(ms_data::BasicMassSpecData{T}, scan_idx::Integer)::T where T
    getHighMzs(ms_data)[scan_idx]
end


# TIC getters
function getTICs(ms_data::BasicMassSpecData{T}) where T
    return ms_data.data[:TIC]::SentinelArrays.ChainedVector{T, Arrow.Primitive{T, Vector{T}}}
end
function getTIC(ms_data::BasicMassSpecData{T}, scan_idx::Integer)::T where T
    getTICs(ms_data)[scan_idx]
end


# Center m/z getters
function getCenterMzs(ms_data::BasicMassSpecData{T}) where T
    return ms_data.data[:centerMz]::SentinelArrays.ChainedVector{Union{Missing, T}, Arrow.Primitive{Union{Missing, T}, Vector{T}}}
end
function getCenterMz(ms_data::BasicMassSpecData{T}, scan_idx::Integer)::Union{Missing, T} where T
    getCenterMzs(ms_data)[scan_idx]
end


# Isolation width getters
function getIsolationWidthMzs(ms_data::BasicMassSpecData{T}) where T
    return ms_data.data[:isolationWidthMz]::SentinelArrays.ChainedVector{Union{Missing, T}, Arrow.Primitive{Union{Missing, T}, Vector{T}}}
end
function getIsolationWidthMz(ms_data::BasicMassSpecData{T}, scan_idx::Integer)::Union{Missing, T} where T
    getIsolationWidthMzs(ms_data)[scan_idx]
end

# Collision energy field getters
function getCollisionEnergyFields(ms_data::BasicMassSpecData{T}) where T
    return ms_data.data[:collisionEnergyField]::SentinelArrays.ChainedVector{Union{Missing, T}, Arrow.Primitive{Union{Missing, T}, Vector{T}}}
end
function getCollisionEnergyField(ms_data::BasicMassSpecData{T}, scan_idx::Integer)::Union{Missing, T} where T
    getCollisionEnergyFields(ms_data)[scan_idx]
end


# Collision energy eV field getters
function getCollisionEnergyEvFields(ms_data::BasicMassSpecData{T}) where T
    return ms_data.data[:collisionEnergyEvField]::SentinelArrays.ChainedVector{Union{Missing, T}, Arrow.Primitive{Union{Missing, T}, Vector{T}}}
end
function getCollisionEnergyEvField(ms_data::BasicMassSpecData{T}, scan_idx::Integer)::Union{Missing, T} where T
    getCollisionEnergyEvFields(ms_data)[scan_idx]
end

# MS order getters remain UInt8
getMsOrders(ms_data::BasicMassSpecData)::SentinelArrays.ChainedVector{UInt8, Arrow.Primitive{UInt8, Vector{UInt8}}} = ms_data.data[:msOrder]
function getMsOrder(ms_data::BasicMassSpecData, scan_idx::Integer)::UInt8
    getMsOrders(ms_data)[scan_idx]
end