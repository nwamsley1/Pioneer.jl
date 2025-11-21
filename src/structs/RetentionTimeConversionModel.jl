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

#Uniform Cubic Smoothing Spline
#Given number of evenly spaced control points and data (time, value)
#solve coefficients for B-spline basis. Then build a speedy implementation. 
struct UniformSpline{N, T<:AbstractFloat} 
    coeffs::SVector{N, T}
    degree::Int64
    first::T
    last::T
    bin_width::T
end

abstract type RtConversionModel end
getModel(rt::RtConversionModel) = rt.model

struct SplineRtConversionModel <: RtConversionModel
    model::UniformSpline
end
(s::SplineRtConversionModel)(x::AbstractFloat) = s.model(x)

struct LinearRtConversionModel <: RtConversionModel
    slope::Float32
    intercept::Float32
end

function (m::LinearRtConversionModel)(rt::AbstractFloat)
    return m.intercept + m.slope * rt
end

# Override getModel for LinearRtConversionModel since it doesn't wrap another model
getModel(m::LinearRtConversionModel) = m

struct IdentityModel <: RtConversionModel
    model::Function
    function IdentityModel()
        new(x::Float32 -> x::Float32)
    end
end

(i::IdentityModel)(x::AbstractFloat) = i.model(x)

RtConversionModel() = IdentityModel()

"""
    IrtRefinementModel

File-specific model to refine library iRT predictions using amino acid composition.

# Fields
- `use_refinement::Bool`: Whether refinement improves validation MAE
- `aa_coefficients::Dict{Char, Float32}`: Per-AA weights (20 standard AAs)
- `intercept::Float32`: Model intercept
- `irt_coefficient::Float32`: Weight for library_irt feature
- `mae_original::Float32`: Validation MAE without refinement
- `mae_refined::Float32`: Validation MAE with refinement
- `r2_train::Float32`: Training R²
- `r2_val::Float32`: Validation R²

# Callable Interface
Model is callable: `refined_irt = model(sequence::String, library_irt::Float32)`

# Algorithm
Predicts error = library_irt - observed_irt, then:
refined_irt = library_irt - predicted_error

# Example
```julia
model = IrtRefinementModel(true, aa_weights, 0.5f0, 0.1f0, ...)
refined = model("PEPTIDE", 50.0f0)  # Returns refined iRT
```
"""
struct IrtRefinementModel
    use_refinement::Bool
    aa_coefficients::Dict{Char, Float32}
    intercept::Float32
    irt_coefficient::Float32
    mae_original::Float32
    mae_refined::Float32
    r2_train::Float32
    r2_val::Float32
end

"""
    (model::IrtRefinementModel)(sequence::String, library_irt::Float32) -> Float32

Apply iRT refinement to a sequence. Zero-allocation via Dict lookup.
"""
function (model::IrtRefinementModel)(sequence::String, library_irt::Float32)::Float32
    if !model.use_refinement
        return library_irt
    end

    # Calculate predicted error
    error_pred = model.intercept + model.irt_coefficient * library_irt

    # Add AA contributions
    for aa in sequence
        if haskey(model.aa_coefficients, aa)
            error_pred += model.aa_coefficients[aa]
        end
    end

    # Return refined iRT (subtract predicted error)
    return library_irt - error_pred
end

# 20 standard amino acids for model training
const STANDARD_AAS = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
                       'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']