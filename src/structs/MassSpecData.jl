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

# Abstract type for mass spectrometry data
abstract type MassSpecData end
abstract type NonIonMobilityData <: MassSpecData end 


# Parameterized struct for mass spectrometry data with AbstractFloat constraint
struct BasicNonIonMobilityMassSpecData{T <: AbstractFloat} <: NonIonMobilityData
    data::Arrow.Table
    n::Int  # Number of rows in the table
end

# Parameterized struct for mass spectrometry data with AbstractFloat constraint
struct BatchNonIonMobilityMassSpecData{T <: AbstractFloat} <: NonIonMobilityData
    data::Arrow.Table
    n::Int  # Number of rows in the table
end
Base.length(ms_data::MassSpecData) = ms_data.n


# Constructor method for creating BasicNonIonMobilityMassSpecData from a file path
function BasicMassSpecData(file_path::String)
    table = Arrow.Table(file_path)
    
    # Extract the type inside the Union{Missing, Float32/Float64}
    union_type = eltype(eltype(table[:mz_array]))
    float_type = first(filter(t -> t <: AbstractFloat, Base.uniontypes(union_type)))
    
    # Determine the length (number of rows)
    n = length(table[:mz_array])
    if typeof(table[:msOrder]) == Arrow.Primitive{UInt8, Vector{UInt8}}
        return BasicNonIonMobilityMassSpecData{float_type}(table, n)
    else
        return BatchNonIonMobilityMassSpecData{float_type}(table, n)
    end
end

#=====================================#
#Basic 
#=====================================#
# Define length method

# Generic SubArray type alias
const BasicMsSubArray{T} = SubArray{
        Union{Missing, T},
        1, 
        Arrow.Primitive{
            Union{Missing, T}, 
            Vector{T}
        }, 
        Tuple{UnitRange{Int64}}, 
        true
    }

# Getter methods for individual elements
# MZ array getters
function getMzArrays(ms_data::BasicNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:mz_array]
end
function getMzArray(ms_data::BasicNonIonMobilityMassSpecData{T}, scan_idx::Integer)::BasicMsSubArray{T} where T
    getMzArrays(ms_data)[scan_idx]
end

# Intensity array getters
function getIntensityArrays(ms_data::BasicNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:intensity_array]
end
function getIntensityArray(ms_data::BasicNonIonMobilityMassSpecData{T}, scan_idx::Integer)::BasicMsSubArray{T} where T
    getIntensityArrays(ms_data)[scan_idx]
end

# String data doesn't need type parameter
getScanHeaders(ms_data::BasicNonIonMobilityMassSpecData):: Arrow.List{String, Int32, Vector{UInt8}} = ms_data.data[:scanHeader]
function getScanHeader(ms_data::BasicNonIonMobilityMassSpecData, scan_idx::Integer)::String
    getScanHeaders(ms_data)[scan_idx]
end

# Scan number getters
getScanNumbers(ms_data::BasicNonIonMobilityMassSpecData)::Arrow.Primitive{Int32, Vector{Int32}} = ms_data.data[:scanNumber]
function getScanNumber(ms_data::BasicNonIonMobilityMassSpecData, scan_idx::Integer)::Int32
    getScanNumbers(ms_data)[scan_idx]
end

# Base peak m/z getters
function getBasePeakMzs(ms_data::BasicNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:basePeakMz]::Arrow.Primitive{T, Vector{T}}
end
function getBasePeakMz(ms_data::BasicNonIonMobilityMassSpecData{T}, scan_idx::Integer)::T where T
    getBasePeakMzs(ms_data)[scan_idx]
end

# Base peak intensity getters
function getBasePeakIntensities(ms_data::BasicNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:basePeakIntensity]::Arrow.Primitive{T, Vector{T}}
end
function getBasePeakIntensity(ms_data::BasicNonIonMobilityMassSpecData{T}, scan_idx::Integer)::T where T
    getBasePeakIntensities(ms_data)[scan_idx]
end

# Packet type getters remain Int32
getPacketTypes(ms_data::BasicNonIonMobilityMassSpecData)::Arrow.Primitive{Int32, Vector{Int32}} = ms_data.data[:packetType]
function getPacketType(ms_data::BasicNonIonMobilityMassSpecData, scan_idx::Integer)::Int32
    getPacketTypes(ms_data)[scan_idx]
end

# Retention time getters
function getRetentionTimes(ms_data::BasicNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:retentionTime]::Arrow.Primitive{T, Vector{T}}
end
function getRetentionTime(ms_data::BasicNonIonMobilityMassSpecData{T}, scan_idx::Integer)::T where T
    getRetentionTimes(ms_data)[scan_idx]
end

# Low m/z getters
function getLowMzs(ms_data::BasicNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:lowMz]::Arrow.Primitive{T, Vector{T}}
end
function getLowMz(ms_data::BasicNonIonMobilityMassSpecData{T}, scan_idx::Integer)::T where T
    getLowMzs(ms_data)[scan_idx]
end

# High m/z getters
function getHighMzs(ms_data::BasicNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:highMz]::Arrow.Primitive{T, Vector{T}}
end
function getHighMz(ms_data::BasicNonIonMobilityMassSpecData{T}, scan_idx::Integer)::T where T
    getHighMzs(ms_data)[scan_idx]
end


# TIC getters
function getTICs(ms_data::BasicNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:TIC]::Arrow.Primitive{T, Vector{T}}
end
function getTIC(ms_data::BasicNonIonMobilityMassSpecData{T}, scan_idx::Integer)::T where T
    getTICs(ms_data)[scan_idx]
end


# Center m/z getters
function getCenterMzs(ms_data::BasicNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:centerMz]::Arrow.Primitive{Union{Missing, T}, Vector{T}}
end
function getCenterMz(ms_data::BasicNonIonMobilityMassSpecData{T}, scan_idx::Integer)::Union{Missing, T} where T
    getCenterMzs(ms_data)[scan_idx]
end


# Isolation width getters
function getIsolationWidthMzs(ms_data::BasicNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:isolationWidthMz]::Arrow.Primitive{Union{Missing, T}, Vector{T}}
end
function getIsolationWidthMz(ms_data::BasicNonIonMobilityMassSpecData{T}, scan_idx::Integer)::Union{Missing, T} where T
    getIsolationWidthMzs(ms_data)[scan_idx]
end

# Collision energy field getters
function getCollisionEnergyFields(ms_data::BasicNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:collisionEnergyField]::Arrow.Primitive{Union{Missing, T}, Vector{T}}
end
function getCollisionEnergyField(ms_data::BasicNonIonMobilityMassSpecData{T}, scan_idx::Integer)::Union{Missing, T} where T
    getCollisionEnergyFields(ms_data)[scan_idx]
end


# Collision energy eV field getters
function getCollisionEnergyEvFields(ms_data::BasicNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:collisionEnergyEvField]::Arrow.Primitive{Union{Missing, T}, Vector{T}}
end
function getCollisionEnergyEvField(ms_data::BasicNonIonMobilityMassSpecData{T}, scan_idx::Integer)::Union{Missing, T} where T
    getCollisionEnergyEvFields(ms_data)[scan_idx]
end

# MS order getters remain UInt8
getMsOrders(ms_data::BasicNonIonMobilityMassSpecData)::Arrow.Primitive{UInt8, Vector{UInt8}} = ms_data.data[:msOrder]
function getMsOrder(ms_data::BasicNonIonMobilityMassSpecData, scan_idx::Integer)::UInt8
    getMsOrders(ms_data)[scan_idx]
end

#=====================================#
#Basic 
#=====================================#
# Generic SubArray type alias
const BatchMsSubArray{T} = SubArray{
    Union{Missing, T}, 
    1, 
    Arrow.Primitive{Union{Missing, T}, Vector{T}}, 
    Tuple{UnitRange{Int64}}, 
    true
}

# Getter methods for individual elements
# MZ array getters
function getMzArrays(ms_data::BatchNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:mz_array]
end
function getMzArray(ms_data::BatchNonIonMobilityMassSpecData{T}, scan_idx::Integer)::BatchMsSubArray{T} where T
    getMzArrays(ms_data)[scan_idx]
end

# Intensity array getters
function getIntensityArrays(ms_data::BatchNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:intensity_array]
end
function getIntensityArray(ms_data::BatchNonIonMobilityMassSpecData{T}, scan_idx::Integer)::BatchMsSubArray{T} where T
    getIntensityArrays(ms_data)[scan_idx]
end

# String data doesn't need type parameter
getScanHeaders(ms_data::BatchNonIonMobilityMassSpecData)::SentinelArrays.ChainedVector{String, Arrow.List{String, Int32, Vector{UInt8}}} = ms_data.data[:scanHeader]
function getScanHeader(ms_data::BatchNonIonMobilityMassSpecData, scan_idx::Integer)::String
    getScanHeaders(ms_data)[scan_idx]
end

# Scan number getters
getScanNumbers(ms_data::BatchNonIonMobilityMassSpecData)::SentinelArrays.ChainedVector{Int32, Arrow.Primitive{Int32, Vector{Int32}}} = ms_data.data[:scanNumber]
function getScanNumber(ms_data::BatchNonIonMobilityMassSpecData, scan_idx::Integer)::Int32
    getScanNumbers(ms_data)[scan_idx]
end

# Base peak m/z getters
function getBasePeakMzs(ms_data::BatchNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:basePeakMz]::SentinelArrays.ChainedVector{T, Arrow.Primitive{T, Vector{T}}}
end
function getBasePeakMz(ms_data::BatchNonIonMobilityMassSpecData{T}, scan_idx::Integer)::T where T
    getBasePeakMzs(ms_data)[scan_idx]
end

# Base peak intensity getters
function getBasePeakIntensities(ms_data::BatchNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:basePeakIntensity]::SentinelArrays.ChainedVector{T, Arrow.Primitive{T, Vector{T}}}
end
function getBasePeakIntensity(ms_data::BatchNonIonMobilityMassSpecData{T}, scan_idx::Integer)::T where T
    getBasePeakIntensities(ms_data)[scan_idx]
end

# Packet type getters remain Int32
getPacketTypes(ms_data::BatchNonIonMobilityMassSpecData)::SentinelArrays.ChainedVector{Int32, Arrow.Primitive{Int32, Vector{Int32}}} = ms_data.data[:packetType]
function getPacketType(ms_data::BatchNonIonMobilityMassSpecData, scan_idx::Integer)::Int32
    getPacketTypes(ms_data)[scan_idx]
end

# Retention time getters
function getRetentionTimes(ms_data::BatchNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:retentionTime]::SentinelArrays.ChainedVector{T, Arrow.Primitive{T, Vector{T}}}
end
function getRetentionTime(ms_data::BatchNonIonMobilityMassSpecData{T}, scan_idx::Integer)::T where T
    getRetentionTimes(ms_data)[scan_idx]
end

# Low m/z getters
function getLowMzs(ms_data::BatchNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:lowMz]::SentinelArrays.ChainedVector{T, Arrow.Primitive{T, Vector{T}}}
end
function getLowMz(ms_data::BatchNonIonMobilityMassSpecData{T}, scan_idx::Integer)::T where T
    getLowMzs(ms_data)[scan_idx]
end

# High m/z getters
function getHighMzs(ms_data::BatchNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:highMz]::SentinelArrays.ChainedVector{T, Arrow.Primitive{T, Vector{T}}}
end
function getHighMz(ms_data::BatchNonIonMobilityMassSpecData{T}, scan_idx::Integer)::T where T
    getHighMzs(ms_data)[scan_idx]
end


# TIC getters
function getTICs(ms_data::BatchNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:TIC]::SentinelArrays.ChainedVector{T, Arrow.Primitive{T, Vector{T}}}
end
function getTIC(ms_data::BatchNonIonMobilityMassSpecData{T}, scan_idx::Integer)::T where T
    getTICs(ms_data)[scan_idx]
end


# Center m/z getters
function getCenterMzs(ms_data::BatchNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:centerMz]::SentinelArrays.ChainedVector{Union{Missing, T}, Arrow.Primitive{Union{Missing, T}, Vector{T}}}
end
function getCenterMz(ms_data::BatchNonIonMobilityMassSpecData{T}, scan_idx::Integer)::Union{Missing, T} where T
    getCenterMzs(ms_data)[scan_idx]
end


# Isolation width getters
function getIsolationWidthMzs(ms_data::BatchNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:isolationWidthMz]::SentinelArrays.ChainedVector{Union{Missing, T}, Arrow.Primitive{Union{Missing, T}, Vector{T}}}
end
function getIsolationWidthMz(ms_data::BatchNonIonMobilityMassSpecData{T}, scan_idx::Integer)::Union{Missing, T} where T
    getIsolationWidthMzs(ms_data)[scan_idx]
end

# Collision energy field getters
function getCollisionEnergyFields(ms_data::BatchNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:collisionEnergyField]::SentinelArrays.ChainedVector{Union{Missing, T}, Arrow.Primitive{Union{Missing, T}, Vector{T}}}
end
function getCollisionEnergyField(ms_data::BatchNonIonMobilityMassSpecData{T}, scan_idx::Integer)::Union{Missing, T} where T
    getCollisionEnergyFields(ms_data)[scan_idx]
end


# Collision energy eV field getters
function getCollisionEnergyEvFields(ms_data::BatchNonIonMobilityMassSpecData{T}) where T
    return ms_data.data[:collisionEnergyEvField]::SentinelArrays.ChainedVector{Union{Missing, T}, Arrow.Primitive{Union{Missing, T}, Vector{T}}}
end
function getCollisionEnergyEvField(ms_data::BatchNonIonMobilityMassSpecData{T}, scan_idx::Integer)::Union{Missing, T} where T
    getCollisionEnergyEvFields(ms_data)[scan_idx]
end

# MS order getters remain UInt8
getMsOrders(ms_data::BatchNonIonMobilityMassSpecData)::SentinelArrays.ChainedVector{UInt8, Arrow.Primitive{UInt8, Vector{UInt8}}} = ms_data.data[:msOrder]
function getMsOrder(ms_data::BatchNonIonMobilityMassSpecData, scan_idx::Integer)::UInt8
    getMsOrders(ms_data)[scan_idx]
end