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

const CHOLMOD_LOCK = ReentrantLock()

"""
    whitsmddw(x, y, w, λ; d=2)

Whittaker smoother with divided differences and weights. This is a Julia port 
of Paul Eilers’ 2003 MATLAB function. 

# Arguments
- `x::AbstractVector`: Sample positions (must be strictly increasing).
- `y::AbstractVector`: Data values at each position.
- `w::AbstractVector`: Weights for each data point.
- `λ::Real`: Smoothing parameter (larger = smoother).
- `d::Integer=2`: Order of differences (default 2).

# Returns
- `z::Vector`: Smoothed series (same length as `y`).
"""
function whitsmddw(x::AbstractVector{Float32},
                   y::AbstractVector{Float32},
                   w::AbstractVector{Float32},
                   n::Int,
                   λ::Float32 = 0.0002f0;
                   d::Int64=2)

    # create constant‑time views on the active slice
    xs = @view x[1:n]
    ys = @view y[1:n]
    ws = @view w[1:n]

    # Build the dth-order differencing matrix
    D = ddmat(xs, d)  # (n-d) x n

    # Build weights as a sparse diagonal matrix
    W = spdiagm(0 => ws)  # n x n

    # Form the system M*z = w.*y, where M = (W + λ * D' * D)
    M = W + λ * (D' * D)

    # CHOLMODFactorization requires Float64
    b = convert(Vector{Float64}, ws .* ys)

    # Initialize
    prob = LinearProblem(Symmetric(M), b)

    # CHOLMOD is not thread-safe — serialize access to prevent EXCEPTION_ACCESS_VIOLATION.
    # If this lock becomes a bottleneck, replace with pre-allocated dense buffers per thread
    # or a direct pentadiagonal solver (the Whittaker-Henderson system with d=2 is always pentadiagonal).
    sol = lock(CHOLMOD_LOCK) do
        solve(prob, CHOLMODFactorization())
    end

    # Downstream requires Float32
    z = Float32.(sol.u)

    # Don't allow negatives
    @. z = max(z, 0f0)

    return z
end

"""
    ddmat(x, d)

Compute the divided differencing matrix of order `d` for sample positions `x`.

# Arguments
- `x::AbstractVector`: 1D array of sample positions
- `d::Integer`: order of the differences

# Returns
- A sparse matrix `D` of size `(length(x) - d) x length(x)` such that 
  `D * Y` gives the divided differences of order `d` applied to the vector `Y`.

This function is a Julia port of the MATLAB code by Paul Eilers (2003).
"""
function ddmat(x::AbstractVector{Float32}, d::Int64;)
    # Compute divided differencing matrix of order d
    #
    # Input
    #   x:  vector of sampling positions
    #   d:  order of diffferences
    # Output
    #   D:  the matrix; D * Y gives divided differences of order d
    #
    # Paul Eilers, 2003
    
    m = size(x,1);
    if d == 0
        # Sparse identity matrix
        return spdiagm(0 => ones(Float64, m))
    else
        # Differences in the x-array for current order
        dx = x[d+1:end] .- x[1:end-d]

        # Build sparse diagonal matrix with 1/dx along the main diagonal
        V = spdiagm(0 => 1.0f0 ./ dx)  # size: (m-d) x (m-d)

        # Recursively get the (d-1)th differencing matrix
        D_lower = ddmat(x, d - 1)    # size: (m-(d-1)) x m

        # Take row-wise consecutive differences of D_lower
        # `diff(D_lower)` in MATLAB is the same as `D_lower[2:end, :] .- D_lower[1:end-1, :]` in Julia
        D_diff = D_lower[2:end, :] .- D_lower[1:end-1, :]  # size: (m-d) x m

        # Multiply to get the final dth differencing matrix
        return V * D_diff  # size: (m-d) x m
    end
end