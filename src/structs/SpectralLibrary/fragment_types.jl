# Fragment ion types for spectral libraries.
#
# Type hierarchy:
#   Ion → LibraryIon → LibraryFragmentIon → AltimeterFragment → DetailedFrag / SplineDetailedFrag
#                                         → SimpleFrag (lightweight, used during index building)
#
# DetailedFrag has fixed intensity; SplineDetailedFrag has spline coefficients
# for intensity prediction across collision energies.

abstract type LibraryIon{T<:AbstractFloat} <: Ion{T} end

getPrecCharge(f::LibraryIon)::UInt8 = f.prec_charge

abstract type LibraryFragmentIon{T<:AbstractFloat} <: LibraryIon{T} end

getPrecID(f::LibraryFragmentIon{T}) where {T<:AbstractFloat} = f.prec_id
getPrecMZ(f::LibraryFragmentIon{T}) where {T<:AbstractFloat} = f.prec_mz

"""
    SimpleFrag{T}

Lightweight fragment used during index building. Contains precursor metadata
needed for partitioning (prec_mz, prec_irt) plus a score byte.
"""
struct SimpleFrag{T<:AbstractFloat} <: LibraryFragmentIon{T}
    mz::T
    prec_id::UInt32
    prec_mz::T
    prec_irt::T
    prec_charge::UInt8
    score::UInt8
end
ArrowTypes.arrowname(::Type{SimpleFrag}) = :SimpleFrag
ArrowTypes.JuliaType(::Val{:SimpleFrag}) = SimpleFrag

getScore(pbi::SimpleFrag)::UInt8 = pbi.score
getIRT(pbi::SimpleFrag{T}) where {T<:AbstractFloat} = pbi.prec_irt

# ============================================================================
# AltimeterFragment — detailed fragment with full annotation
# ============================================================================

"""
    AltimeterFragment{T} <: LibraryFragmentIon{T}

Base type for fragments with full ion annotation (type, position, charge, rank).
Subtypes differ in how they store intensity: fixed (DetailedFrag) vs spline (SplineDetailedFrag).
"""
abstract type AltimeterFragment{T} <: LibraryFragmentIon{T} end

getPID(f::AltimeterFragment) = f.prec_id
getMz(f::AltimeterFragment) = f.mz
getIonType(f::AltimeterFragment) = f.ion_type
isY(f::AltimeterFragment) = f.is_y
isB(f::AltimeterFragment) = f.is_b
isP(f::AltimeterFragment) = f.is_p
isIso(f::AltimeterFragment) = f.is_isotope
getFragCharge(f::AltimeterFragment) = f.frag_charge
getIonPosition(f::AltimeterFragment) = f.ion_position
getPrecCharge(f::AltimeterFragment) = f.prec_charge
getRank(f::AltimeterFragment) = f.rank
getSulfurCount(f::AltimeterFragment) = f.sulfur_count

# ============================================================================
# Intensity type dispatch
# ============================================================================

abstract type IntensityDataType end
struct ConstantType <: IntensityDataType end

abstract type NceModel{T<:AbstractFloat} end

struct SplineType{M,T<:AbstractFloat} <: IntensityDataType
    knots::NTuple{M, T}
    nce::T
    degree::Int64
end
getNCE(st::SplineType) = st.nce
getKnots(st::SplineType) = st.knots
getDegree(st::SplineType) = st.degree

# ============================================================================
# CompactFrag — packed transition struct for search pipeline (16 bytes)
# ============================================================================

"""
    CompactFrag{T} <: AltimeterFragment{T}

Lightweight fragment for the search hot path. 16 bytes vs DetailedFrag's 24 bytes.
Drops `ion_type` (not used in search) and packs bools/small integers into two UInt16 fields.

**Bit layout for packed_a (UInt16)**:
    bits 0-3:   flags (is_y=0x01, is_b=0x02, is_p=0x04, is_isotope=0x08)
    bits 4-5:   frag_charge (0-3)
    bits 6-8:   prec_charge (0-7)
    bits 9-11:  sulfur_count (0-7, clamped to 5 in isotope splines)
    bits 12-15: spare

**Bit layout for packed_b (UInt16)**:
    bits 0-5:   ion_position (0-63)
    bits 6-10:  rank (0-31)
    bits 11-15: spare
"""
struct CompactFrag{T<:AbstractFloat} <: AltimeterFragment{T}
    prec_id::UInt32
    mz::T
    intensity::Float16
    packed_a::UInt16   # flags + frag_charge + prec_charge + sulfur_count
    packed_b::UInt16   # ion_position + rank
end

# Bit constants for packed_a flags
const _CF_IS_Y       = UInt16(0x0001)
const _CF_IS_B       = UInt16(0x0002)
const _CF_IS_P       = UInt16(0x0004)
const _CF_IS_ISOTOPE = UInt16(0x0008)

# Getters — all @inline for zero overhead
@inline getPID(f::CompactFrag) = f.prec_id
@inline getMz(f::CompactFrag) = f.mz
@inline getIntensity(f::CompactFrag) = f.intensity
@inline getIntensity(f::CompactFrag, ::ConstantType) = f.intensity

@inline isY(f::CompactFrag)  = (f.packed_a & _CF_IS_Y) != 0
@inline isB(f::CompactFrag)  = (f.packed_a & _CF_IS_B) != 0
@inline isP(f::CompactFrag)  = (f.packed_a & _CF_IS_P) != 0
@inline isIso(f::CompactFrag) = (f.packed_a & _CF_IS_ISOTOPE) != 0
@inline getFragCharge(f::CompactFrag)  = UInt8((f.packed_a >> 4) & 0x03)
@inline getPrecCharge(f::CompactFrag)  = UInt8((f.packed_a >> 6) & 0x07)
@inline getSulfurCount(f::CompactFrag) = UInt8((f.packed_a >> 9) & 0x07)
@inline getIonPosition(f::CompactFrag) = UInt8(f.packed_b & 0x3f)
@inline getRank(f::CompactFrag)        = UInt8((f.packed_b >> 6) & 0x1f)

# ion_type not stored — return 0 (callers should use isY/isB/isP instead)
@inline getIonType(f::CompactFrag) = UInt16(0)

# Constructor from individual fields (used by fillTransitionList!)
@inline function CompactFrag(prec_id::UInt32, mz::T, intensity::Float16,
        is_y::Bool, is_b::Bool, is_p::Bool, is_isotope::Bool,
        frag_charge::UInt8, ion_position::UInt8, prec_charge::UInt8,
        rank::UInt8, sulfur_count::UInt8) where {T<:AbstractFloat}
    packed_a = UInt16(is_y) | (UInt16(is_b) << 1) | (UInt16(is_p) << 2) |
               (UInt16(is_isotope) << 3) |
               (UInt16(min(frag_charge, UInt8(3))) << 4) |
               (UInt16(min(prec_charge, UInt8(7))) << 6) |
               (UInt16(min(sulfur_count, UInt8(7))) << 9)
    packed_b = UInt16(min(ion_position, UInt8(63))) |
               (UInt16(min(rank, UInt8(31))) << 6)
    CompactFrag{T}(prec_id, mz, intensity, packed_a, packed_b)
end

# Default constructor
CompactFrag{T}() where {T<:AbstractFloat} = CompactFrag{T}(
    zero(UInt32), zero(T), zero(Float16), zero(UInt16), zero(UInt16))

ArrowTypes.arrowname(::Type{CompactFrag{Float32}}) = :CompactFrag
ArrowTypes.JuliaType(::Val{:CompactFrag}) = CompactFrag{Float32}

getMZ(f::CompactFrag) = f.mz

# ============================================================================
# SplineCompactFrag — packed spline-interpolated intensity (28 bytes)
# ============================================================================

"""
    SplineCompactFrag{N,T} <: AltimeterFragment{T}

Spline-interpolated fragment with packed metadata. 28 bytes for N=4, Float32
(vs SplineDetailedFrag's 36 bytes). Same packed_a/packed_b layout as CompactFrag.

N is the number of spline coefficients for intensity prediction across NCE values.
"""
struct SplineCompactFrag{N,T<:AbstractFloat} <: AltimeterFragment{T}
    prec_id::UInt32
    mz::T
    intensity::NTuple{N, T}  # spline coefficients
    packed_a::UInt16
    packed_b::UInt16
end

# Getters — same packed layout as CompactFrag
@inline getPID(f::SplineCompactFrag) = f.prec_id
@inline getMz(f::SplineCompactFrag) = f.mz
@inline getMZ(f::SplineCompactFrag) = f.mz
@inline isY(f::SplineCompactFrag)  = (f.packed_a & _CF_IS_Y) != 0
@inline isB(f::SplineCompactFrag)  = (f.packed_a & _CF_IS_B) != 0
@inline isP(f::SplineCompactFrag)  = (f.packed_a & _CF_IS_P) != 0
@inline isIso(f::SplineCompactFrag) = (f.packed_a & _CF_IS_ISOTOPE) != 0
@inline getFragCharge(f::SplineCompactFrag)  = UInt8((f.packed_a >> 4) & 0x03)
@inline getPrecCharge(f::SplineCompactFrag)  = UInt8((f.packed_a >> 6) & 0x07)
@inline getSulfurCount(f::SplineCompactFrag) = UInt8((f.packed_a >> 9) & 0x07)
@inline getIonPosition(f::SplineCompactFrag) = UInt8(f.packed_b & 0x3f)
@inline getRank(f::SplineCompactFrag)        = UInt8((f.packed_b >> 6) & 0x1f)
@inline getIonType(f::SplineCompactFrag) = UInt16(0)

function getIntensity(pf::SplineCompactFrag{N, T}, intensity_type::SplineType{M, T}) where {M, N, T<:AbstractFloat}
    return splevl(getNCE(intensity_type), getKnots(intensity_type), pf.intensity, getDegree(intensity_type))::T
end

SplineCompactFrag{N,T}() where {N,T<:AbstractFloat} = SplineCompactFrag(
    zero(UInt32), zero(T), ntuple(_ -> zero(T), Val(N)), zero(UInt16), zero(UInt16))

# Constructor from individual fields (mirrors CompactFrag's positional-args
# constructor; used by BuildSpecLib's getDetailedFrags spline overload so the
# build path emits SplineCompactFrag directly — no SplineDetailedFrag → Compact
# conversion step).
@inline function SplineCompactFrag(prec_id::UInt32, mz::T, intensity::NTuple{N, T},
        is_y::Bool, is_b::Bool, is_p::Bool, is_isotope::Bool,
        frag_charge::UInt8, ion_position::UInt8, prec_charge::UInt8,
        rank::UInt8, sulfur_count::UInt8) where {N, T<:AbstractFloat}
    packed_a = UInt16(is_y) | (UInt16(is_b) << 1) | (UInt16(is_p) << 2) |
               (UInt16(is_isotope) << 3) |
               (UInt16(min(frag_charge, UInt8(3))) << 4) |
               (UInt16(min(prec_charge, UInt8(7))) << 6) |
               (UInt16(min(sulfur_count, UInt8(7))) << 9)
    packed_b = UInt16(min(ion_position, UInt8(63))) |
               (UInt16(min(rank, UInt8(31))) << 6)
    SplineCompactFrag{N, T}(prec_id, mz, intensity, packed_a, packed_b)
end

ArrowTypes.arrowname(::Type{SplineCompactFrag{4, Float32}}) = :SplineCompactFrag
ArrowTypes.JuliaType(::Val{:SplineCompactFrag}) = SplineCompactFrag{4, Float32}

# ============================================================================
# DetailedFrag — fixed intensity (library build/storage format)
# ============================================================================

"""
    DetailedFrag{T}

Fragment with a single fixed intensity value. Used in standard library lookups
where collision energy is fixed or already applied.
"""
struct DetailedFrag{T<:AbstractFloat} <: AltimeterFragment{T}
    prec_id::UInt32
    mz::T
    intensity::Float16
    ion_type::UInt16
    is_y::Bool
    is_b::Bool
    is_p::Bool
    is_isotope::Bool
    frag_charge::UInt8
    ion_position::UInt8
    prec_charge::UInt8
    rank::UInt8
    sulfur_count::UInt8
end

getIntensity(f::DetailedFrag{T}) where {T<:AbstractFloat} = f.intensity
getIntensity(pf::DetailedFrag, ::ConstantType) = pf.intensity

ArrowTypes.arrowname(::Type{DetailedFrag{Float32}}) = :DetailedFrag
ArrowTypes.JuliaType(::Val{:DetailedFrag}) = DetailedFrag

# Convert DetailedFrag → CompactFrag (must be after DetailedFrag definition)
function CompactFrag(f::DetailedFrag{T}) where {T<:AbstractFloat}
    CompactFrag(f.prec_id, f.mz, f.intensity,
        f.is_y, f.is_b, f.is_p, f.is_isotope,
        f.frag_charge, f.ion_position, f.prec_charge, f.rank, f.sulfur_count)
end

DetailedFrag{T}() where {T<:AbstractFloat} = DetailedFrag(
    zero(UInt32), zero(T), zero(Float16), zero(UInt16),
    false, false, false, false,
    zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8)
)

# ============================================================================
# SplineDetailedFrag — spline-interpolated intensity
# ============================================================================

"""
    SplineDetailedFrag{N,T}

Fragment with spline coefficients for intensity prediction across collision energies.
N is the number of spline coefficients.
"""
struct SplineDetailedFrag{N,T<:AbstractFloat} <: AltimeterFragment{T}
    prec_id::UInt32
    mz::T
    intensity::NTuple{N, T}
    ion_type::UInt16
    is_y::Bool
    is_b::Bool
    is_p::Bool
    is_isotope::Bool
    frag_charge::UInt8
    ion_position::UInt8
    prec_charge::UInt8
    rank::UInt8
    sulfur_count::UInt8
end

ArrowTypes.arrowname(::Type{SplineDetailedFrag{4, Float32}}) = :DetailedFrag
ArrowTypes.JuliaType(::Val{:SplineDetailedFrag}) = DetailedFrag

function getIntensity(pf::SplineDetailedFrag{N, T}, intensity_type::SplineType{M, T}) where {M, N, T<:AbstractFloat}
    return splevl(getNCE(intensity_type), getKnots(intensity_type), pf.intensity, getDegree(intensity_type))::T
end

SplineDetailedFrag{N,T}() where {N,T<:AbstractFloat} = SplineDetailedFrag(
    zero(UInt32), zero(T), NTuple{N, T}(undef), zero(UInt16),
    false, false, false, false,
    zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8)
)

# Convert SplineDetailedFrag → SplineCompactFrag (for loading legacy libraries
# that wrote SplineDetailedFrag to disk; new builds emit SplineCompactFrag
# directly via getDetailedFrags).
function SplineCompactFrag(f::SplineDetailedFrag{N, T}) where {N, T<:AbstractFloat}
    SplineCompactFrag(f.prec_id, f.mz, f.intensity,
        f.is_y, f.is_b, f.is_p, f.is_isotope,
        f.frag_charge, f.ion_position, f.prec_charge, f.rank, f.sulfur_count)
end

# ============================================================================
# I/O
# ============================================================================

"""Load detailed fragments from .jls (preferred) or legacy .jld2 format."""
