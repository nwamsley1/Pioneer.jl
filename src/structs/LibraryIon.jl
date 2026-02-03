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

abstract type LibraryIon{T<:AbstractFloat} <: Ion{T} end

getPrecCharge(f::LibraryIon)::UInt8 = f.prec_charge

abstract type LibraryFragmentIon{T<:AbstractFloat} <: LibraryIon{T} end

getPrecID(f::LibraryFragmentIon{T}) where {T<:AbstractFloat} = f.prec_id
getPrecMZ(f::LibraryFragmentIon{T}) where {T<:AbstractFloat} = f.prec_mz
#getIntensity(f::LibraryFragmentIon{T}) where {T<:AbstractFloat} = f.intensity

#=
function reset!(lf::Vector{L}, last_non_empty::Int64) where {L<:LibraryFragmentIon{Float32}}
    for i in range(1, last_non_empty)
        lf[i] = L()
    end
end
=#
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

#Need this information for each distinct fragment type
#Need this information for each distinct fragment type

"""
   AltimeterFragment{T<:AbstractFloat} <: LibraryFragmentIon{T}

Abstract supertype for fragments used in the Altimeter spectral library system.

# Type Parameters
- `T`: Floating point precision type 

# Required Fields
All concrete subtypes must implement:
- `prec_id::UInt32`: Unique identifier of parent precursor
- `mz::T`: Fragment m/z ratio
- `ion_type::UInt16`: Numeric identifier for ion type
- `is_y::Bool`: True if y-ion
- `is_b::Bool`: True if b-ion 
- `is_p::Bool`: True if precursor ion
- `is_isotope::Bool`: True if isotope peak
- `frag_charge::UInt8`: Fragment charge state
- `ion_position::UInt8`: Position in peptide sequence
- `prec_charge::UInt8`: Precursor charge state
- `rank::UInt8`: Intensity rank among siblings
- `sulfur_count::UInt8`: Number of sulfur atoms

# Interface Methods
Standard getters for required fields:
- `getPID(f::AltimeterFragment)`: Get precursor ID
- `getMz(f::AltimeterFragment)`: Get m/z ratio
- `getIonType(f::AltimeterFragment)`: Get ion type
- `isY(f::AltimeterFragment)`: Check if y-ion
- `isB(f::AltimeterFragment)`: Check if b-ion
- `isP(f::AltimeterFragment)`: Check if precursor ion
- `isIso(f::AltimeterFragment)`: Check if isotope peak
- `getFragCharge(f::AltimeterFragment)`: Get fragment charge
- `getIonPosition(f::AltimeterFragment)`: Get sequence position
- `getPrecCharge(f::AltimeterFragment)`: Get precursor charge
- `getRank(f::AltimeterFragment)`: Get intensity rank
- `sulfurCount(f::AltimeterFragment)`: Get sulfur count

# Concrete Subtypes
- `DetailedFrag`: Fragment with fixed intensity
- `SplineDetailedFrag`: Fragment with intensity modeled via spline coefficients

# Notes
- Base type for spectral library fragment ions with detailed annotations
- Subtypes differ primarily in how they handle intensity information
- Used in fragment lookup systems for spectrum prediction
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

abstract type IntensityDataType end

struct SplineType{M,T<:AbstractFloat} <: IntensityDataType
    knots::NTuple{M, T}
    nce::T
    degree::Int64
end

getNCE(st::SplineType) = st.nce
getKnots(st::SplineType) = st.knots
getDegree(st::SplineType) = st.degree

struct ConstantType <: IntensityDataType end

"""
    DetailedFrag{T<:AbstractFloat} <: AltimeterFragment{T}

A detailed representation of a fragment ion with intensity information.

# Type Parameters
- `T`: Floating point precision type

# Fields
- `prec_id::UInt32`: Unique identifier of parent precursor
- `mz::T`: Fragment m/z ratio
- `intensity::Float16`: Observed or predicted intensity
- `ion_type::UInt16`: Numeric identifier for ion type
- `is_y::Bool`: True if y-ion
- `is_b::Bool`: True if b-ion
- `is_p::Bool`: True if precursor ion
- `is_isotope::Bool`: True if isotope peak
- `frag_charge::UInt8`: Fragment charge state
- `ion_position::UInt8`: Position in peptide sequence
- `prec_charge::UInt8`: Precursor charge state
- `rank::UInt8`: Intensity rank among siblings
- `sulfur_count::UInt8`: Number of sulfur atoms

# Notes
- Unlike SplineDetailedFrag, this type stores a single intensity value rather than spline coefficients
- Used in standard library lookups where collision energy is fixed
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

function getIntensity(
    pf::DetailedFrag{T},
    intensity_type::ConstantType
    ) where {T<:AbstractFloat}
    return pf.intensity
end

getIntensity(f::DetailedFrag{T}) where {T<:AbstractFloat} = f.intensity
ArrowTypes.arrowname(::Type{DetailedFrag{Float32}}) = :DetailedFrag
ArrowTypes.JuliaType(::Val{:DetailedFrag}) = DetailedFrag

#LibraryFragment{T}() where {T<:AbstractFloat} = LibraryFragment(zero(T), zero(UInt8), false, zero(UInt8), zero(UInt8), zero(Float32), zero(UInt8), zero(UInt32), zero(UInt8), zero(UInt8))
DetailedFrag{T}() where {T<:AbstractFloat} = DetailedFrag(

                zero(UInt32), #prec_id

                zero(T), #mz
                zero(Float16), #intensity

                zero(UInt16),
                false,
                false,
                false,
                false,
                
                zero(UInt8), #frag_charge
                zero(UInt8), #ion_position
                zero(UInt8), #prec_charge
                zero(UInt8), #rank
                zero(UInt8)  #sulfur_count
 )
 

"""
    SplineDetailedFrag{N,T<:AbstractFloat} <: AltimeterFragment{T}

Fragment ion with spline coefficients for intensity prediction across collision energies.

# Type Parameters
- `N`: Number of spline coefficients
- `T`: Floating point precision type

# Fields
- `prec_id::UInt32`: Unique identifier of parent precursor
- `mz::T`: Fragment m/z ratio
- `intensity::NTuple{N,T}`: Spline coefficients for intensity prediction
- `ion_type::UInt16`: Numeric identifier for ion type
- `is_y::Bool`: True if y-ion
- `is_b::Bool`: True if b-ion
- `is_p::Bool`: True if precursor ion
- `is_isotope::Bool`: True if isotope peak
- `frag_charge::UInt8`: Fragment charge state
- `ion_position::UInt8`: Position in peptide sequence
- `prec_charge::UInt8`: Precursor charge state
- `rank::UInt8`: Intensity rank among siblings
- `sulfur_count::UInt8`: Number of sulfur atoms
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

#=
function load_detailed_frags(filename::String)
    # Define base abstract type if it doesn't exist
    if !isdefined(Main, :AltimeterFragment)
        @eval Main abstract type AltimeterFragment{T<:AbstractFloat} end
    end

    # Define SplineDetailedFrag type if it doesn't exist in workspace
    if !isdefined(Main, :SplineDetailedFrag)
        @eval Main struct Pioneer.SplineDetailedFrag{N,T<:AbstractFloat} <: AltimeterFragment{T}
            prec_id::UInt32
            mz::T
            intensity::NTuple{N,T}
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
    end
    
    # Load data directly - no conversion needed since we've defined the type
    jldopen(filename, "r") do file
        read(file, "data")
    end
end
=#

function load_detailed_frags(filename::String)
    # Determine paths for new (.jls) and legacy (.jld2) formats
    jls_path = if endswith(filename, ".jld2")
        replace(filename, ".jld2" => ".jls")
    elseif endswith(filename, ".jls")
        filename
    else
        filename * ".jls"
    end

    jld2_path = if endswith(filename, ".jls")
        replace(filename, ".jls" => ".jld2")
    elseif endswith(filename, ".jld2")
        filename
    else
        filename * ".jld2"
    end

    # Try new .jls format first
    if isfile(jls_path)
        return deserialize_from_jls(jls_path)
    end

    # Fall back to legacy .jld2 format with warning
    if isfile(jld2_path)
        @warn "Loading legacy JLD2 format from $jld2_path. Consider rebuilding library for better performance."

        # Define base abstract type if it doesn't exist (for JLD2 type reconstruction)
        if !isdefined(Main, :AltimeterFragment)
            @eval Main abstract type AltimeterFragment{T<:AbstractFloat} end
        end

        return jldopen(jld2_path, "r") do file
            read(file, "data")
        end
    end

    error("Fragment file not found: $jls_path or $jld2_path")
end

ArrowTypes.arrowname(::Type{SplineDetailedFrag{4, Float32}}) = :DetailedFrag
ArrowTypes.JuliaType(::Val{:SplineDetailedFrag}) = DetailedFrag
#This needs to be eddited for the splines
function getIntensity(
                    pf::SplineDetailedFrag{N, T},
                    intensity_type::SplineType{M, T}
                    ) where {M, N, T<:AbstractFloat}
   return splevl(getNCE(intensity_type), getKnots(intensity_type), pf.intensity, getDegree(intensity_type))::T
end

#LibraryFragment{T}() where {T<:AbstractFloat} = LibraryFragment(zero(T), zero(UInt8), false, zero(UInt8), zero(UInt8), zero(Float32), zero(UInt8), zero(UInt32), zero(UInt8), zero(UInt8))
SplineDetailedFrag{N,T}() where {N,T<:AbstractFloat} = SplineDetailedFrag(

                zero(UInt32), #prec_id

                zero(T), #mz
                NTuple{N, T}(undef), #intensity

                zero(UInt16),
                false,
                false,
                false,
                false,
                
                zero(UInt8), #frag_charge
                zero(UInt8), #ion_position
                zero(UInt8), #prec_charge
                zero(UInt8), #rank
                zero(UInt8)  #sulfur_count
 )

# Generic conversion for any AltimeterFragment to DetailedFrag
function convert_to_detailed(frag::DetailedFrag{T}, ::ConstantType) where {T <: AbstractFloat}
    DetailedFrag(
        getPID(frag),
        getMz(frag),
        Float16(frag.intensity),  # Default intensity
        getIonType(frag),
        isY(frag),
        isB(frag),
        isP(frag),
        isIso(frag),
        getFragCharge(frag),
        getIonPosition(frag),
        getPrecCharge(frag),
        getRank(frag),
        getSulfurCount(frag)
    )
end

# Specialized version for SplineDetailedFrag that handles intensity calculation
function convert_to_detailed(
    frag::SplineDetailedFrag{N,T}, 
    spline_data::SplineType{M, T},
) where {N,M,T}
    DetailedFrag(
        getPID(frag),
        getMz(frag),
        Float16(getIntensity(frag, spline_data)),
        getIonType(frag),
        isY(frag),
        isB(frag),
        isP(frag),
        isIso(frag),
        getFragCharge(frag),
        getIonPosition(frag),
        getPrecCharge(frag),
        getRank(frag),
        getSulfurCount(frag)
    )
end
abstract type LibraryFragmentLookup end

"""
    StandardFragmentLookup{T<:AbstractFloat} <: LibraryFragmentLookup

A basic container for fragment ion data with fixed intensities.

# Type Parameters
- `T`: Floating point precision type

# Fields
- `frags::Vector{DetailedFrag{T}}`: Vector of fragment ions
- `prec_frag_ranges::Vector{UInt64}`: Index ranges for fragments belonging to each precursor

# Methods
- `getFrag(lfp::StandardFragmentLookup, prec_idx::Integer)`: Get fragment at index
- `getFragments(lfp::StandardFragmentLookup)`: Get full fragment vector
- `getPrecFragRange(lfp::StandardFragmentLookup, prec_idx::Integer)`: Get range of indices for precursor

# Notes
- Alternative to SplineFragmentLookup
- Fragment intensities are stored directly rather than as spline coefficients
"""
struct StandardFragmentLookup{T<:AbstractFloat} <: LibraryFragmentLookup
    frags::Vector{DetailedFrag{T}}
    prec_frag_ranges::Vector{UInt64}
end
getFrag(lfp::StandardFragmentLookup{<:AbstractFloat}, prec_idx::Integer) = lfp.frags[prec_idx]
getFragments(lfp::StandardFragmentLookup{<:AbstractFloat}) = lfp.frags
getPrecFragRange(lfp::StandardFragmentLookup, prec_idx::Integer)::UnitRange{UInt64} = range(lfp.prec_frag_ranges[prec_idx], lfp.prec_frag_ranges[prec_idx+1]-one(UInt64))

"""
   abstract type NceModel{T<:AbstractFloat} end

Abstract type for normalized collision energy prediction models.

# Type Parameters
- `T`: Floating point precision type used by the model

# Interface Requirements
All subtypes must implement:
- Call operator: `(model::ConcreteModel)(x::AbstractFloat, charge::Integer)`
- Vector call operator: `(model::ConcreteModel)(x::AbstractVector, charge::AbstractVector)`
"""
abstract type NceModel{T<:AbstractFloat} end

function getSplineData(lfp::StandardFragmentLookup, prec_charge::UInt8, prec_mz::T) where {T<:AbstractFloat}
    return ConstantType()
end

function getSplineData(lfp::StandardFragmentLookup)
    return ConstantType()
end

# Overloads that accept nce_model for API consistency (unused by StandardFragmentLookup)
function getSplineData(
    lfp::StandardFragmentLookup,
    nce_model::NceModel{T},
    prec_charge::UInt8,
    prec_mz::T
) where {T<:AbstractFloat}
    return ConstantType()
end

function getSplineData(
    lfp::StandardFragmentLookup,
    nce_model::NceModel{T}
) where {T<:AbstractFloat}
    return ConstantType()
end



"""
   PiecewiseNceModel{T<:AbstractFloat} <: NceModel{T}

A piecewise model for normalized collision energy prediction that includes charge dependence.

# Type Parameters
- `T`: Floating point precision type

# Fields
- `breakpoint::T`: x-value where model transitions from linear to constant
- `left_slope::T`: Slope of linear component for x ≤ breakpoint
- `left_intercept::T`: Intercept of linear component for x ≤ breakpoint
- `right_value::T`: Constant value for x > breakpoint
- `charge_slope::T`: Linear charge dependence coefficient

# Model
For a given mass-to-charge ratio x and charge state z:
- When x ≤ breakpoint: f(x,z) = left_slope * x + left_intercept + charge_slope * z
- When x > breakpoint: f(x,z) = right_value + charge_slope * z

# Construction
   PiecewiseNceModel(x::T) where {T<:AbstractFloat}

Creates a PiecewiseNceModel that will return the input value x when evaluated at any point. 
This is achieved by setting:
- breakpoint = 0
- left_slope = 1
- left_intercept = 0
- right_value = x
- charge_slope = 0

This ensures f(m,z) = x for any mass m and charge z.

# Methods
   (f::PiecewiseNceModel)(x::AbstractFloat, charge::Integer)
   (f::PiecewiseNceModel)(x::AbstractVector, charge::AbstractVector)
   fit_nce_model(pwlm::PiecewiseNceModel, x, y, charge, breakpoint)

See also: [`fit_nce_model`](@ref), [`NceModel`](@ref)
"""
struct PiecewiseNceModel{T<:AbstractFloat} <: NceModel{T}
    breakpoint::T
    left_slope::T
    left_intercept::T
    right_value::T
    charge_slope::T
 end



if !@isdefined(SplineFragmentLookup)
"""
    SplineFragmentLookup{N,M,C,T<:AbstractFloat} <: LibraryFragmentLookup

A container for fragment ion data that uses spline interpolation for intensity predictions.

# Type Parameters
- `N`: Number of spline coefficients per fragment
- `M`: Number of spline knots
- `C`: Number of NCE models (typically one per charge state)
- `T`: Floating point precision type

# Fields
- `frags::Vector{SplineDetailedFrag{N,T}}`: Vector of fragment ions with spline coefficients
- `prec_frag_ranges::Vector{UInt64}`: Index ranges for fragments belonging to each precursor
- `knots::NTuple{M, T}`: Knot points for spline interpolation
- `nce_models::NceModelContainer{C,T}`: Models for predicting collision energy by charge state
"""
struct SplineFragmentLookup{N,M,T<:AbstractFloat} <: LibraryFragmentLookup
    frags::Vector{SplineDetailedFrag{N,T}}
    prec_frag_ranges::Vector{UInt64}
    knots::NTuple{M, T}
    degree::Int64
end
end
getDegree(lfp::SplineFragmentLookup) = lfp.degree
#getKnots(lfp::SplineFragmentLookup) = lfp.knots
getFrag(lfp::SplineFragmentLookup, prec_idx::Integer) = lfp.frags[prec_idx]
getFragments(lfp::SplineFragmentLookup) = lfp.frags
getPrecFragRange(lfp::SplineFragmentLookup, prec_idx::Integer)::UnitRange{UInt64} = range(lfp.prec_frag_ranges[prec_idx], lfp.prec_frag_ranges[prec_idx+1]-one(UInt64))


# Version with precursor-specific NCE
function getSplineData(
    lfp::SplineFragmentLookup{N,M,T},
    nce_model::NceModel{T},
    prec_charge::UInt8,
    prec_mz::T
) where {N,M,T<:AbstractFloat}
    return SplineType(
        getKnots(lfp),
        nce_model(prec_mz, prec_charge),  # Direct call to model
        getDegree(lfp)
    )
end

# Version with default NCE (for cases without specific precursor)
function getSplineData(
    lfp::SplineFragmentLookup{N,M,T},
    nce_model::NceModel{T}
) where {N,M,T<:AbstractFloat}
    return SplineType(
        getKnots(lfp),
        nce_model(),  # Call with no arguments
        getDegree(lfp)
    )
end

getKnots(lfp::SplineFragmentLookup) = lfp.knots

# updateNceModel removed - NCE model is now passed as a parameter, not stored in the lookup
"""
    PrecursorBinItem{T<:AbstractFloat}

Item in a precursor bin. Minimal information required to know if the fragment has a correct precursor mass, and if so, 
what the precursor ID is.  

### Fields

- prec_id::UInt32 -- Unique identifier for the precursor
- prec_mz::T -- m/z of the precursor

### Examples

### GetterMethods

- getPrecID(pbi::PrecursorBinItem) = pbi.prec_id
- getPrecMZ(pbi::PrecursorBinItem) = pbi.prec_mz
"""
struct PrecursorBinFragment{T<:AbstractFloat} <: LibraryFragmentIon{T}
    prec_id::UInt32
    prec_mz::T #Only need to tell if the peptide is in the quad isolation window
    score::UInt8 
    charge::UInt8
end

getScore(pbi::PrecursorBinFragment{T}) where {T<:AbstractFloat} = pbi.score


abstract type LibraryPrecursors end

struct StandardLibraryPrecursors <: LibraryPrecursors
    data::Arrow.Table
    n::Int64
    accession_numbers_to_pid::Dictionary{String, UInt32}
    pid_to_cv_fold::Vector{UInt8}
    function StandardLibraryPrecursors(precursor_table::Arrow.Table)
        try
            #precursor_table = Arrow.Table(precursor_table_path)
            n = length(precursor_table[:sequence])
            accession_numbers = precursor_table[:accession_numbers]::Arrow.List{String, Int32, Vector{UInt8}}

            #Maps each unique accession number group to an UInt32 
            unique_proteins = unique(accession_numbers);
            accession_number_to_pgid = Dictionary(
                unique_proteins, range(one(UInt32), UInt32(length(unique_proteins)))
            );

            #precursor idxs to cross-validation folds
            #all precursors corresponding to a given protein-group end up in the same cross validation fold
            pg_to_cv_fold = Dictionary{String, UInt8}()
            cv_folds = UInt8[0, 1]
            Random.seed!(1776)
            for pg in unique_proteins
                insert!(pg_to_cv_fold, pg, rand(cv_folds))
            end
            pid_to_cv_fold = Vector{UInt8}(undef, n)
            for pid in range(1, n)
                pid_to_cv_fold[pid] = pg_to_cv_fold[accession_numbers[pid]]
            end
            if length(keys(accession_number_to_pgid)) <= 1
                @user_warn "Library did not include protein accession numbers. Seeting cross-validation folds based on precursor_idx"
                for pid in range(1, n)
                    pid_to_cv_fold[pid] = rand(cv_folds)
                end
            end
            new(
                precursor_table, n, accession_number_to_pgid, pid_to_cv_fold
            )
        catch e
            @user_warn "Failed to load precursor_table_path: $precursor_table_path"
            throw(e)
        end
    end
end

#=
struct PlexedLibraryPrecursors <: LibraryPrecursors
    data::Arrow.Table
    n::Int64
    accession_numbers_to_pid::Dictionary{String, UInt32}
    pid_to_cv_fold::Vector{UInt8}
    function PlexedLibraryPrecursors(precursor_table::Arrow.Table)
        try
            #precursor_table = Arrow.Table(precursor_table_path)
            n = length(precursor_table[:sequence])
            accession_numbers = precursor_table[:accession_numbers]::Arrow.List{String, Int32, Vector{UInt8}}

            #Maps each unique accession number group to an UInt32 
            unique_proteins = unique(accession_numbers);
            accession_number_to_pgid = Dictionary(
                unique_proteins, range(one(UInt32), UInt32(length(unique_proteins)))
            );

            #precursor idxs to cross-validation folds
            #all precursors corresponding to a given protein-group end up in the same cross validation fold
            pg_to_cv_fold = Dictionary{String, UInt8}()
            cv_folds = UInt8[0, 1]
            Random.seed!(1776)
            for pg in unique_proteins
                insert!(pg_to_cv_fold, pg, rand(cv_folds))
            end
            pid_to_cv_fold = Vector{UInt8}(undef, n)
            for pid in range(1, n)
                pid_to_cv_fold[pid] = pg_to_cv_fold[accession_numbers[pid]]
            end
            if length(keys(accession_number_to_pgid)) <= 1
                @user_warn "Library did not include protein accession numbers. Seeting cross-validation folds based on precursor_idx"
                for pid in range(1, n)
                    pid_to_cv_fold[pid] = rand(cv_folds)
                end
            end
            new(
                precursor_table, n, accession_number_to_pgid, pid_to_cv_fold
            )
        catch e
            @user_warn "Failed to load precursor_table_path: $precursor_table_path"
            throw(e)
        end
    end
end
=#
function SetPrecursors(precursor_table::Arrow.Table)
    #if "plex" in names(DataFrame(precursor_table))
    #    return PlexedLibraryPrecursors(precursor_table)
    #else
        return StandardLibraryPrecursors(precursor_table)
    #end
end
# Define length method
import Base.length
Base.length(ms_data::LibraryPrecursors) = ms_data.n
getProteinGroupId(lp::LibraryPrecursors, accession_numbers::String)::UInt32 = lp.accession_numbers_to_pid[accession_numbers]
getCvFold(lp::LibraryPrecursors, precursor_idx::I) where {I<:Integer} = lp.pid_to_cv_fold[precursor_idx]
#getProteomeIdentifiers(lp::LibraryPrecursors)::Arrow.List{String, Int32, Vector{UInt8}} = lp.data[:proteome_identifiers]
getProteomeIdentifiers(lp::LibraryPrecursors)::Arrow.List{S,Int32,Array{UInt8,1}} where {S<:AbstractString} = lp.data[:proteome_identifiers]
getAccessionNumbers(lp::LibraryPrecursors)::Arrow.List{String, Int32, Vector{UInt8}} = lp.data[:accession_numbers]
#getSequence(lp::LibraryPrecursors)::Arrow.List{String, Int32, Vector{UInt8}} = lp.data[:sequence]
getSequence(lp::LibraryPrecursors)::Arrow.List{S,Int32,Array{UInt8,1}} where {S<:AbstractString} = lp.data[:sequence]
getStartIdx(lp::LibraryPrecursors)::Arrow.Primitive{UInt32, Vector{UInt32}} = lp.data[:start_idx]
getStartIdx(lp::LibraryPrecursors, idx::I) where {I<:Integer} = lp.data[:start_idx][idx]
getStructuralMods(lp::LibraryPrecursors)::Arrow.List{Union{Missing, String}, Int32, Vector{UInt8}} = lp.data[:structural_mods]
getCharge(lp::LibraryPrecursors)::Arrow.Primitive{UInt8, Vector{UInt8}}  = lp.data[:prec_charge]
getCollisionEnergy(lp::LibraryPrecursors)::Arrow.Primitive{Float32, Vector{Float32}} = lp.data[:collision_energy]
getIsDecoy(lp::LibraryPrecursors)::Arrow.BoolVector{Bool}  = lp.data[:is_decoy]
getEntrapmentGroupId(lp::LibraryPrecursors)::Arrow.Primitive{UInt8, Vector{UInt8}} = lp.data[:entrapment_group_id]
getBasePepId(lp::LibraryPrecursors)::Arrow.Primitive{UInt32, Vector{UInt32}} = lp.data[:base_pep_id]
getMz(lp::LibraryPrecursors)::Arrow.Primitive{Float32, Vector{Float32}}  = lp.data[:mz]
getLength(lp::LibraryPrecursors)::Arrow.Primitive{UInt8, Vector{UInt8}}  = lp.data[:length]
getMissedCleavages(lp::LibraryPrecursors)::Arrow.Primitive{UInt8, Vector{UInt8}} = lp.data[:missed_cleavages]
getNumEnzymaticTermini(lp::LibraryPrecursors)::Arrow.Primitive{UInt8, Vector{UInt8}} = lp.data[:num_enzymatic_termini]
getIrt(lp::LibraryPrecursors)::Arrow.Primitive{Float32, Vector{Float32}} = lp.data[:irt]
getSulfurCount(lp::LibraryPrecursors)::Arrow.Primitive{UInt8, Vector{UInt8}} = lp.data[:sulfur_count]
getIsotopicMods(lp::LibraryPrecursors)::Arrow.List{Union{Missing, String}, Int32, Vector{UInt8}} = lp.data[:isotopic_mods]
getPartnerPrecursorIdx(lp::LibraryPrecursors)::Arrow.Primitive{Union{Missing, I}, Vector{I}} where {I<:Integer} = lp.data[:partner_precursor_idx]
# Define two different versions of the function
function extract_pair_idx(pair_idx_column::Arrow.Primitive{Union{Missing, UInt32}, Array{UInt32,1}}, idx)
    value = pair_idx_column[idx]
    return ismissing(value) ? zero(UInt32) : convert(UInt32, value)
end

function extract_pair_idx(pair_idx_column::Arrow.Primitive{UInt32, Vector{UInt32}}, idx)
    return pair_idx_column[idx]
end

# Add a generic fallback for other types
function extract_pair_idx(pair_idx_column, idx)
    value = pair_idx_column[idx]
    if isa(value, UInt32)
        return value
    elseif ismissing(value)
        return zero(UInt32)
    else
        return convert(UInt32, value)
    end
end
getPairIdx(lp::LibraryPrecursors) = lp.data[:pair_id]
#getPlex(lp::PlexedLibraryPrecursors)::Arrow.Primitive{I, Vector{I}} where {I<:Integer} = lp.data[:plex]

