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

function PiecewiseNceModel(x::T) where {T<:AbstractFloat}
    PiecewiseNceModel(zero(T), one(T), zero(T), x, zero(T))
end

"""
   fit_nce_model(pwlm::PiecewiseNceModel{T}, x::AbstractVector, y::AbstractVector, 
                charge::AbstractVector, breakpoint::Real) where {T<:AbstractFloat}

Fit a piecewise model with charge dependence to normalized collision energy data, ensuring continuity at the breakpoint.

# Arguments
- `pwlm::PiecewiseNceModel{T}`: Template for the model type to fit
- `x::AbstractVector`: Mass-to-charge ratio values
- `y::AbstractVector`: Observed collision energy values
- `charge::AbstractVector`: Charge states
- `breakpoint::Real`: Point where model transitions from linear to constant behavior

# Returns
- `PiecewiseNceModel{T}`: Fitted model with parameters:
   - left_slope: Slope of linear region (x ≤ breakpoint)
   - left_intercept: Intercept of linear region
   - right_value: Constant value for x > breakpoint (= left_slope * breakpoint + left_intercept)
   - charge_slope: Linear charge dependence coefficient

# Model
For mass-to-charge ratio x and charge state z:
- When x ≤ breakpoint: f(x,z) = left_slope * x + left_intercept + charge_slope * z
- When x > breakpoint: f(x,z) = right_value + charge_slope * z

The model is fit by minimizing the sum of squared residuals while maintaining continuity 
at the breakpoint through the constraint: right_value = left_slope * breakpoint + left_intercept.

# Notes
- Initial parameter estimates are obtained using linear regression on the left region
- Optimization is performed using the LBFGS algorithm
- Charge dependence is linear and consistent across both regions
- Continuity at breakpoint is enforced by construction

See also: [`PiecewiseNceModel`](@ref)
"""
function fit_nce_model(
    pwlm::PiecewiseNceModel{T},
    x::AbstractVector,
    y::AbstractVector,
    charge::AbstractVector,
    breakpoint::Real
) where {T<:AbstractFloat}
    # Define the objective function to minimize sum of squared residuals
    function objective(params)
        # params[1] = left_slope
        # params[2] = left_intercept
        # params[3] = charge_slope
        
        # Calculate the right_value from continuity constraint
        right_value = params[1] * breakpoint + params[2]
        
        # Calculate residuals for left side
        left_mask = x .<= breakpoint
        y_pred_left = params[1] .* x[left_mask] .+
                     params[3] .* charge[left_mask] .+
                     params[2]
        residuals_left = y_pred_left .- y[left_mask]
        
        # Calculate residuals for right side
        right_mask = x .> breakpoint
        y_pred_right = fill(right_value, sum(right_mask)) .+
                     params[3] .* charge[right_mask]
        residuals_right = y_pred_right .- y[right_mask]
        
        # Total sum of squared residuals
        return sum(residuals_left.^2) + sum(residuals_right.^2)
    end
    
    # Initial guess using simple linear regression
    left_mask = x .<= breakpoint
    X_left = [ones(sum(left_mask)) x[left_mask] charge[left_mask]]
    initial_coef = X_left \ y[left_mask]
    
    initial_guess = [
        initial_coef[2],  # left_slope
        initial_coef[1],  # left_intercept
        initial_coef[3]   # charge_slope
    ]
    
    # Optimize
    result = optimize(objective, initial_guess, LBFGS())
    optimal_params = Optim.minimizer(result)
    
    # Calculate right_value using continuity constraint
    right_value = optimal_params[1] * breakpoint + optimal_params[2]
    
    return PiecewiseNceModel(
        T(float(breakpoint)),
        T(optimal_params[1]),  # left_slope
        T(optimal_params[2]),  # left_intercept
        T(right_value),        # right_value from continuity
        T(optimal_params[3])   # charge_slope
    )
end
 
function (f::PiecewiseNceModel)()
    return f.right_value
end

function (f::PiecewiseNceModel)(x::AbstractFloat, charge::Integer)
base = if x <= f.breakpoint
    f.left_slope * x + f.left_intercept
else
    f.right_value
end
return base + f.charge_slope * charge
end
 
 # Add method for vectors
 function (f::PiecewiseNceModel)(x::AbstractVector, charge::AbstractVector)
    return map((x,c) -> f(x,c), x, charge)
 end

# `medians` is a Vector, not an NTuple{N,T}, deliberately. It used to carry the bin count as a type
# parameter, which specialised every method touching an NCE model per dataset -- the count depends on
# how many m/z bins had enough PSMs, so a single trace can produce both {1} and {2}.
#
# Nothing was gained by it: the access below is a single runtime-indexed load (`medians[offset+idx-1]`),
# which a tuple cannot do without spilling to the stack; `offsets` and `n_bins` already carry the
# structure at runtime; the value is built as a Vector and converted at the end; there is one model
# per file, fetched once per search method rather than per precursor; and it is stored in a
# `Dict{Int64, NceModel}` whose abstract value type means the call site dispatches dynamically anyway.
struct BinnedMedianNceModel{T<:AbstractFloat} <: NceModel{T}
    medians::Vector{T}
    offsets::NTuple{6, UInt8}
    n_bins::NTuple{6, UInt8}
    mz_min::NTuple{6, T}
    bin_width::NTuple{6, T}
    default_nce::T
end

function (m::BinnedMedianNceModel{T})(mz::AbstractFloat, charge::Integer) where {T}
    c = Int(charge)
    if c < 1 || c > 6 || m.offsets[c] == 0x00
        best_c = 0
        best_d = typemax(Int)
        for k in 1:6
            m.offsets[k] == 0x00 && continue
            d = abs(k - c)
            if d < best_d
                best_d = d
                best_c = k
            end
        end
        best_c == 0 && return m.default_nce
        c = best_c
    end
    nb = Int(m.n_bins[c])
    idx = clamp(floor(Int, (T(mz) - m.mz_min[c]) / m.bin_width[c]) + 1, 1, nb)
    return m.medians[Int(m.offsets[c]) + idx - 1]
end

(m::BinnedMedianNceModel)() = m.default_nce

function (m::BinnedMedianNceModel)(x::AbstractVector, charge::AbstractVector)
    return map((xi, ci) -> m(xi, ci), x, charge)
end

function fit_binned_median_nce(
    mz::AbstractVector{T},
    nce::AbstractVector{T},
    charge::AbstractVector,
    default_nce::T;
    min_per_bin::Int = 50
) where {T<:AbstractFloat}
    all_medians = T[]
    offsets = zeros(UInt8, 6)
    n_bins_arr = zeros(UInt8, 6)
    mz_mins = zeros(T, 6)
    bin_widths = zeros(T, 6)

    for c in sort(unique(UInt8.(charge)))
        (c < 0x01 || c > 0x06) && continue
        ci = Int(c)
        mask = UInt8.(charge) .== c
        n_pts = count(mask)
        n_pts < min_per_bin && continue

        c_mz = mz[mask]
        c_nce = nce[mask]
        mz_lo, mz_hi = extrema(c_mz)
        if mz_hi == mz_lo
            mz_hi = mz_lo + one(T)
        end

        for nb in (10, 5, 3, 2, 1)
            bw = (mz_hi - mz_lo) / T(nb)
            ok = true
            meds = Vector{T}(undef, nb)
            for b in 1:nb
                lo = mz_lo + T(b - 1) * bw
                hi = (b == nb) ? mz_hi + one(T) : mz_lo + T(b) * bw
                in_bin = c_nce[(c_mz .>= lo) .& (c_mz .< hi)]
                if length(in_bin) < min_per_bin
                    ok = false
                    break
                end
                meds[b] = median(in_bin)
            end
            if ok
                offsets[ci] = UInt8(length(all_medians) + 1)
                n_bins_arr[ci] = UInt8(nb)
                mz_mins[ci] = mz_lo
                bin_widths[ci] = bw
                append!(all_medians, meds)
                break
            end
        end
    end

    return BinnedMedianNceModel{T}(
        all_medians,
        NTuple{6, UInt8}(offsets),
        NTuple{6, UInt8}(n_bins_arr),
        NTuple{6, T}(mz_mins),
        NTuple{6, T}(bin_widths),
        default_nce
    )
end