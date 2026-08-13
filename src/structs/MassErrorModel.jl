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

abstract type AbstractMassErrorModel end

struct SimpleMassErrorModel{T<:AbstractFloat} <: AbstractMassErrorModel
    mass_offset::T
    mass_tolerance::Tuple{T, T}
end

const MassErrorModel = SimpleMassErrorModel

getRightTol(mem::SimpleMassErrorModel) = last(mem.mass_tolerance)
getLeftTol(mem::SimpleMassErrorModel) = first(mem.mass_tolerance)
getMassOffset(mem::SimpleMassErrorModel) = last(mem.mass_offset)

function getMassCorrection(mem::SimpleMassErrorModel)
    return mem.mass_offset
end

function getCorrectedMz(mem::SimpleMassErrorModel, mz::Float32)
    return Float32(mz - getMassOffset(mem)*(mz/1f6))
end

# Bounds for the theoretical mass
function getMzBoundsReverse(mem::SimpleMassErrorModel, mass::Float32)
    ppm = mass/1f6
    r_tol = getRightTol(mem)*ppm
    l_tol = getLeftTol(mem)*ppm
    return Float32(mass - r_tol), Float32(mass + l_tol)
end

# Bounds for the empirical mass
function getMzBounds(mem::SimpleMassErrorModel, mass::Float32)
    ppm = mass/1f6
    r_tol = getRightTol(mem)*ppm
    l_tol = getLeftTol(mem)*ppm
    return Float32(mass - l_tol), Float32(mass + r_tol)
end

# 3-arg forwarding: SimpleMassErrorModel ignores intensity
getCorrectedMz(mem::SimpleMassErrorModel, mz::Float32, ::Float32) = getCorrectedMz(mem, mz)
getMzBoundsReverse(mem::SimpleMassErrorModel, mass::Float32, ::Float32) = getMzBoundsReverse(mem, mass)
getCorrectedMz(mem::SimpleMassErrorModel, mz::Float32, intensity::Float32, ::Float32) =
    getCorrectedMz(mem, mz, intensity)
getMzBoundsReverse(mem::SimpleMassErrorModel, mass::Float32, log2I::Float32, ::Float32) =
    getMzBoundsReverse(mem, mass, log2I)

@inline function getCorrectedMzAndBounds(mem::SimpleMassErrorModel, mz::Float32, ::Float32)
    corrected = getCorrectedMz(mem, mz)
    low, high = getMzBoundsReverse(mem, corrected)
    return corrected, low, high
end

@inline function getCorrectedMzAndBounds(mem::SimpleMassErrorModel, mz::Float32, intensity::Float32, ::Float32)
    return getCorrectedMzAndBounds(mem, mz, intensity)
end

#==========================================================
Linear Da mass error model: bias_da(mz) = a + b × mz, constant Da tolerance.
Used during parameter tuning after the wide scout discovers the m/z-dependent bias.
==========================================================#

struct LinearDaMassErrorModel{T<:AbstractFloat} <: AbstractMassErrorModel
    intercept::T      # Da bias at mz=0
    slope::T          # Da bias per Da of m/z
    tolerance_da::T   # symmetric half-width in Da
end

function getCorrectedMz(m::LinearDaMassErrorModel, mz::Float32)
    return Float32(mz - (m.intercept + m.slope * mz))
end

function getMzBoundsReverse(m::LinearDaMassErrorModel, mass::Float32)
    return Float32(mass - m.tolerance_da), Float32(mass + m.tolerance_da)
end

function getMzBounds(m::LinearDaMassErrorModel, mass::Float32)
    return Float32(mass - m.tolerance_da), Float32(mass + m.tolerance_da)
end

getRightTol(m::LinearDaMassErrorModel) = m.tolerance_da
getLeftTol(m::LinearDaMassErrorModel) = m.tolerance_da
getMassOffset(m::LinearDaMassErrorModel) = m.intercept
getMassCorrection(m::LinearDaMassErrorModel) = m.intercept

# 3-arg forwarding: ignores intensity
getCorrectedMz(m::LinearDaMassErrorModel, mz::Float32, ::Float32) = getCorrectedMz(m, mz)
getMzBoundsReverse(m::LinearDaMassErrorModel, mass::Float32, ::Float32) = getMzBoundsReverse(m, mass)
getCorrectedMz(m::LinearDaMassErrorModel, mz::Float32, intensity::Float32, ::Float32) =
    getCorrectedMz(m, mz, intensity)
getMzBoundsReverse(m::LinearDaMassErrorModel, mass::Float32, log2I::Float32, ::Float32) =
    getMzBoundsReverse(m, mass, log2I)

@inline function getCorrectedMzAndBounds(m::LinearDaMassErrorModel, mz::Float32, ::Float32)
    corrected = getCorrectedMz(m, mz)
    low, high = getMzBoundsReverse(m, corrected)
    return corrected, low, high
end

@inline function getCorrectedMzAndBounds(m::LinearDaMassErrorModel, mz::Float32, intensity::Float32, ::Float32)
    return getCorrectedMzAndBounds(m, mz, intensity)
end

# Scoring fallbacks
laplace_log_density(::LinearDaMassErrorModel, ::Float32, ::Float32, ::Float32) = Float32(0)
get_default_top3_ll(::LinearDaMassErrorModel) = Float32(0)

#==========================================================
Linear Da bias + ppm tolerance: bias_da(mz) = a + b × mz, then ±tol_ppm on corrected.
Corrects m/z-dependent bias from the wide scout, then uses the initial ppm
tolerance for the collection window (same width as SimpleMassErrorModel).
==========================================================#

struct LinearBiasPpmTolMassErrorModel{T<:AbstractFloat} <: AbstractMassErrorModel
    intercept::T        # Da bias at mz=0
    slope::T            # Da bias per Da of m/z
    tolerance_ppm::T    # symmetric ppm half-width on corrected mass
end

function getCorrectedMz(m::LinearBiasPpmTolMassErrorModel, mz::Float32)
    return Float32(mz - (m.intercept + m.slope * mz))
end

function getMzBoundsReverse(m::LinearBiasPpmTolMassErrorModel, mass::Float32)
    ppm = mass / 1f6
    hw = m.tolerance_ppm * ppm
    return Float32(mass - hw), Float32(mass + hw)
end

function getMzBounds(m::LinearBiasPpmTolMassErrorModel, mass::Float32)
    ppm = mass / 1f6
    hw = m.tolerance_ppm * ppm
    return Float32(mass - hw), Float32(mass + hw)
end

getRightTol(m::LinearBiasPpmTolMassErrorModel) = m.tolerance_ppm
getLeftTol(m::LinearBiasPpmTolMassErrorModel) = m.tolerance_ppm
getMassOffset(m::LinearBiasPpmTolMassErrorModel) = m.intercept
getMassCorrection(m::LinearBiasPpmTolMassErrorModel) = m.intercept

getCorrectedMz(m::LinearBiasPpmTolMassErrorModel, mz::Float32, ::Float32) = getCorrectedMz(m, mz)
getMzBoundsReverse(m::LinearBiasPpmTolMassErrorModel, mass::Float32, ::Float32) = getMzBoundsReverse(m, mass)
getCorrectedMz(m::LinearBiasPpmTolMassErrorModel, mz::Float32, intensity::Float32, ::Float32) =
    getCorrectedMz(m, mz, intensity)
getMzBoundsReverse(m::LinearBiasPpmTolMassErrorModel, mass::Float32, log2I::Float32, ::Float32) =
    getMzBoundsReverse(m, mass, log2I)

@inline function getCorrectedMzAndBounds(m::LinearBiasPpmTolMassErrorModel, mz::Float32, ::Float32)
    corrected = getCorrectedMz(m, mz)
    low, high = getMzBoundsReverse(m, corrected)
    return corrected, low, high
end

@inline function getCorrectedMzAndBounds(m::LinearBiasPpmTolMassErrorModel, mz::Float32, intensity::Float32, ::Float32)
    return getCorrectedMzAndBounds(m, mz, intensity)
end

laplace_log_density(::LinearBiasPpmTolMassErrorModel, ::Float32, ::Float32, ::Float32) = Float32(0)
get_default_top3_ll(::LinearBiasPpmTolMassErrorModel) = Float32(0)
