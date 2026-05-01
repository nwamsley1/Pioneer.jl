# Abstract types and concrete structs for mass spectrometry data access.
#
# Two concrete types wrap Arrow tables from converted mzML files:
# - BasicNonIonMobilityMassSpecData: single-batch Arrow files (Arrow.Primitive columns)
# - BatchNonIonMobilityMassSpecData: multi-batch Arrow files (SentinelArrays.ChainedVector columns)
#
# The constructor BasicMassSpecData() auto-detects which type to use based on
# the msOrder column type in the Arrow table.

abstract type MassSpecData end
abstract type NonIonMobilityData{T <: AbstractFloat} <: MassSpecData end

struct BasicNonIonMobilityMassSpecData{T <: AbstractFloat} <: NonIonMobilityData{T}
    data::Arrow.Table
    n::Int
end

struct BatchNonIonMobilityMassSpecData{T <: AbstractFloat} <: NonIonMobilityData{T}
    data::Arrow.Table
    n::Int
end

Base.length(ms_data::MassSpecData) = ms_data.n

"""
    BasicMassSpecData(file_path::String)

Load an Arrow mass spec file and return the appropriate concrete type.
Auto-detects Basic vs Batch based on the msOrder column layout.
"""
function BasicMassSpecData(file_path::String)
    table = Arrow.Table(file_path)
    n = length(table[:mz_array])

    if n == 0
        @user_warn "Empty MS data file detected: $(basename(file_path)). File contains 0 scans and will be skipped in analysis."
        float_type = Float32
    else
        union_type = eltype(eltype(table[:mz_array]))
        abstract_float_types = filter(t -> t <: AbstractFloat, Base.uniontypes(union_type))
        float_type = isempty(abstract_float_types) ? Float32 : first(abstract_float_types)
    end

    if typeof(table[:msOrder]) == Arrow.Primitive{UInt8, Vector{UInt8}}
        return BasicNonIonMobilityMassSpecData{float_type}(table, n)
    else
        return BatchNonIonMobilityMassSpecData{float_type}(table, n)
    end
end
