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
    cycle_idxs::Union{Nothing, Vector{UInt32}}
    n::Int
end

struct BatchNonIonMobilityMassSpecData{T <: AbstractFloat} <: NonIonMobilityData{T}
    data::Arrow.Table
    cycle_idxs::Union{Nothing, Vector{UInt32}}
    n::Int
end

Base.length(ms_data::MassSpecData) = ms_data.n

const CYCLE_IDX_ISOLATION_WINDOW_TOLERANCE_MZ = 1.0f-3

@inline function _cycle_idx_value(x)
    return ismissing(x) ? 0.0f0 : Float32(x)
end

function compute_cycle_idxs(
    ms_orders::AbstractVector,
    center_mzs::AbstractVector,
    isolation_width_mzs::AbstractVector,
)
    n = length(ms_orders)
    cycle_idxs = Vector{UInt32}(undef, n)
    cycle_idx = UInt32(1)
    has_first_ms2_window = false
    first_center_mz = 0.0f0
    first_isolation_width_mz = 0.0f0

    @inbounds for i in 1:n
        if ms_orders[i] == UInt8(2)
            center_mz = _cycle_idx_value(center_mzs[i])
            isolation_width_mz = _cycle_idx_value(isolation_width_mzs[i])
            if !has_first_ms2_window
                has_first_ms2_window = true
                first_center_mz = center_mz
                first_isolation_width_mz = isolation_width_mz
            elseif abs(center_mz - first_center_mz) <= CYCLE_IDX_ISOLATION_WINDOW_TOLERANCE_MZ &&
                   abs(isolation_width_mz - first_isolation_width_mz) <= CYCLE_IDX_ISOLATION_WINDOW_TOLERANCE_MZ
                cycle_idx += UInt32(1)
            end
        end
        cycle_idxs[i] = cycle_idx
    end

    return cycle_idxs
end

function compute_cycle_idxs(ms_data::MassSpecData)
    return compute_cycle_idxs(getMsOrders(ms_data), getCenterMzs(ms_data), getIsolationWidthMzs(ms_data))
end

function _fallback_cycle_idxs(table::Arrow.Table)
    :cycle_idx in Arrow.names(table) && return nothing
    return compute_cycle_idxs(table[:msOrder], table[:centerMz], table[:isolationWidthMz])
end

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

    cycle_idxs = _fallback_cycle_idxs(table)
    if typeof(table[:msOrder]) == Arrow.Primitive{UInt8, Vector{UInt8}}
        return BasicNonIonMobilityMassSpecData{float_type}(table, cycle_idxs, n)
    else
        return BatchNonIonMobilityMassSpecData{float_type}(table, cycle_idxs, n)
    end
end
