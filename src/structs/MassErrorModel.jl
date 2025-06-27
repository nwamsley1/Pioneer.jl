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


struct MassErrorModel{T<:AbstractFloat}
    mass_offset::T#UniformSpline{N, T}
    mass_tolerance::Tuple{T, T}#UniformSpline{N, T}
end

getRightTol(mem::MassErrorModel) = last(mem.mass_tolerance)
getLeftTol(mem::MassErrorModel) = first(mem.mass_tolerance)
getMassOffset(mem::MassErrorModel) = last(mem.mass_offset)

function getMassCorrection(mem::MassErrorModel)
    return mem.mass_offset
end

function getLocation(mem::MassErrorModel)
    return mem.location
end

#Adjust a theoretical mass to be
#=
"""
    (mem::MassErrorModel)(mass::Float32) -> Tuple{Float32, Float32}

Calculate mass tolerance window based on a mass error model.

Given an input observed m/z ratio, this function calculates the lower and upper bounds of a 
matching theoretical m/z ratio tolerance window.

# Arguments
- `mass::Float32`: The target m/z value to calculate tolerances for

# Returns
- `Tuple{Float32, Float32}`: A tuple containing `(lower_bound, upper_bound)` where:
  - `lower_bound`: The mass minus the left tolerance
  - `upper_bound`: The mass plus the right tolerance

# Details
The function:
1. Converts the mass to ppm scale (dividing by 1e6)
2. Applies a mass offset correction
3. Calculates asymmetric tolerance windows using the model's left and right tolerance parameters
4. Returns the lower and upper bounds as Float32 values

# Examples
```julia
mem = MassErrorModel(0.3, 10.0, 10.0)  # offset, left_tol, right_tol
lower, upper = mem(1000.0f0)  # Calculate window for mass 1000 Da
```
"""
function (mem::MassErrorModel)(mass::Float32)
    ppm_norm = Float32(1e6)
    ppm = mass/ppm_norm
    mass -= getMassOffset(mem)*ppm
    r_tol = getRightTol(mem)*ppm
    l_tol = getLeftTol(mem)*ppm
    return Float32(mass - l_tol), Float32(mass + r_tol)
end
=#
#Correct an empeirical mass. 
function getCorrectedMz(mem::MassErrorModel, mz::Float32)
    return Float32(mz - getMassOffset(mem)*(mz/1e6))
end

"""
   getMzBounds(mem::MassErrorModel, mass::Float32)

Given a theoretical mass and a `MassErrorModel`, gets the minimum and maximum bounds for an expected empirical mass.
Critically, assumes the empirical mass has already been corrected using mass offset. See `getCorrectedMz`. 

### Input

- `mem::MassErrorModel`: -- Model for the mass error of an ion 
- `mass::Float32` -- A theoretical mass/mz to get boundaries for 
### Output
Tuple{Float32, Float32}
A tuple with the lower and upper boundary respectively. 
### Notes

- Suppose the mass error was 3 ppm. And the tolerance was 10 ppm and 5ppm on the left and right hand sides respectively. 
A theoretical mass of 1000000.0f0 m/z, would have a tolerance of (999990.0f0, 1000005.0f0). A theoretical mass falling 
### Algorithm 

### Examples 

"""
#Bounds for the theoretical mass 
function getMzBoundsReverse(mem::MassErrorModel, mass::Float32)
    ppm = mass/(1e6)
    r_tol = getRightTol(mem)*ppm
    l_tol = getLeftTol(mem)*ppm
    return Float32(mass - r_tol), Float32(mass + l_tol)
end

#Bounds for the empirical mass 
function getMzBounds(mem::MassErrorModel, mass::Float32)
    ppm = mass/(1e6)
    r_tol = getRightTol(mem)*ppm
    l_tol = getLeftTol(mem)*ppm
    return Float32(mass - l_tol), Float32(mass + r_tol)
end