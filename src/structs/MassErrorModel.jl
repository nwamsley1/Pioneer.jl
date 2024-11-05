
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
function (mem::MassErrorModel)(mass::Float32)
    ppm_norm = Float32(1e6)
    ppm = mass/ppm_norm
    mass += getMassOffset(mem)*ppm
    r_tol = getRightTol(mem)*ppm
    l_tol = getLeftTol(mem)*ppm
    return Float32(mass - l_tol), Float32(mass + r_tol)
end

#Correct an empeirical mass. 
function getCorrectedMz(mem::MassErrorModel, mz::Float32)
    return  Float32(mz - getMassOffset(mem)*(mz/1e6))
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
function getMzBounds(mem::MassErrorModel, mass::Float32)
    ppm = mass/(1e6)
    r_tol = getRightTol(mem)*ppm
    l_tol = getLeftTol(mem)*ppm
    return Float32(mass - l_tol), Float32(mass + r_tol)
end
