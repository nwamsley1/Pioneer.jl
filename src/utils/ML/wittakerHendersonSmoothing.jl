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

"""
    WHWorkspace

Pre-allocated workspace for the pentadiagonal Whittaker-Henderson solver.
Create one per thread, reuse across all calls. Zero allocations per solve.
"""
mutable struct WHWorkspace
    # Pentadiagonal system bands (Float64 for numerical stability)
    d0::Vector{Float64}
    d1::Vector{Float64}
    d2::Vector{Float64}
    # LDL' factors
    ld1::Vector{Float64}
    ld2::Vector{Float64}
    diag::Vector{Float64}
    # Solve temporaries
    rhs::Vector{Float64}
    z_f64::Vector{Float64}
    z::Vector{Float32}
    h::Vector{Float64}
    # WHSmooth! temporaries (eliminate per-call allocations)
    w_tmp::Vector{Float32}
    x_tmp::Vector{Float32}
    n_max::Int
end

function WHWorkspace(n_max::Int)
    WHWorkspace(
        zeros(Float64, n_max),     # d0
        zeros(Float64, n_max - 1), # d1
        zeros(Float64, n_max - 2), # d2
        zeros(Float64, n_max - 1), # ld1
        zeros(Float64, n_max - 2), # ld2
        zeros(Float64, n_max),     # diag
        zeros(Float64, n_max),     # rhs
        zeros(Float64, n_max),     # z_f64
        zeros(Float32, n_max),     # z
        zeros(Float64, n_max - 1), # h
        zeros(Float32, n_max),     # w_tmp
        zeros(Float32, n_max),     # x_tmp
        n_max
    )
end

"""
    build_DtD_bands!(ws, x, n, λ)

Compute the 5 bands of λ*D'D directly from the sample positions x[1:n].
D is the (n-2)xn second divided difference matrix.

Row j of D (j=1,...,n-2) has nonzeros at columns j, j+1, j+2:
  D[j,j]   = 1 / (h_j * (h_j + h_{j+1}))
  D[j,j+1] = -1 / (h_j * h_{j+1})
  D[j,j+2] = 1 / (h_{j+1} * (h_j + h_{j+1}))

D'D = sum_j (row_j)' * (row_j), a sum of rank-1 outer products.
Each contributes to a 3x3 block at (j:j+2, j:j+2).
We accumulate only the upper triangle (symmetric), storing 3 bands.
"""
function build_DtD_bands!(ws::WHWorkspace, x::AbstractVector{Float32}, n::Int, λ::Float32)
    h = ws.h
    d0 = ws.d0
    d1 = ws.d1
    d2 = ws.d2
    λ64 = Float64(λ)

    # Zero the bands
    @inbounds for i in 1:n
        d0[i] = 0.0
    end
    @inbounds for i in 1:n-1
        d1[i] = 0.0
        h[i] = Float64(x[i+1]) - Float64(x[i])
    end
    @inbounds for i in 1:n-2
        d2[i] = 0.0
    end

    # Accumulate D'D from each row of D
    @inbounds for j in 1:n-2
        hj = h[j]
        hj1 = h[j+1]
        s = hj + hj1

        # Row j entries of D (Newton divided difference)
        a = 1.0 / (hj * s)
        b = -1.0 / (hj * hj1)
        c = 1.0 / (hj1 * s)

        # Rank-1 update: D'D += [a; b; c] * [a, b, c]'
        d0[j]   += λ64 * a * a
        d0[j+1] += λ64 * b * b
        d0[j+2] += λ64 * c * c

        d1[j]   += λ64 * a * b
        d1[j+1] += λ64 * b * c

        d2[j]   += λ64 * a * c
    end

    return nothing
end

"""
    whitsmddw!(ws, x, y, w, n, λ)

Zero-allocation Whittaker-Henderson smoother using pentadiagonal LDL' factorization.

Solves (W + λ D'D) z = W y where M = W + λ D'D is symmetric positive definite pentadiagonal.

All storage is pre-allocated in `ws`. Returns a view into `ws.z[1:n]`.
"""
function whitsmddw!(ws::WHWorkspace,
                     x::AbstractVector{Float32},
                     y::AbstractVector{Float32},
                     w::AbstractVector{Float32},
                     n::Int,
                     λ::Float32)

    d0 = ws.d0
    d1 = ws.d1
    d2 = ws.d2
    ld1 = ws.ld1
    ld2 = ws.ld2
    diag = ws.diag
    rhs = ws.rhs
    z64 = ws.z_f64
    z = ws.z

    # Step 1: Build λ*D'D bands (Float64)
    build_DtD_bands!(ws, x, n, λ)

    # Step 2: Add W to main diagonal
    @inbounds for i in 1:n
        d0[i] += Float64(w[i])
    end

    # Step 3: Form RHS = w .* y (Float64)
    @inbounds for i in 1:n
        rhs[i] = Float64(w[i]) * Float64(y[i])
    end

    # Step 4: LDL' factorization of the symmetric pentadiagonal matrix
    #
    # M is symmetric with bands d0 (main), d1 (off-1), d2 (off-2).
    # We factor M = L D L' where L is unit lower triangular with bandwidth 2.
    #
    # L stored as: ld1[i] = L[i+1, i], ld2[i] = L[i+2, i]
    # D stored as: diag[i] = D[i,i]

    # Row 1
    @inbounds diag[1] = d0[1]

    # Row 2
    if n >= 2
        @inbounds begin
            ld1[1] = d1[1] / diag[1]
            diag[2] = d0[2] - ld1[1] * d1[1]
        end
    end

    # Rows 3..n
    @inbounds for i in 3:n
        ld2[i-2] = d2[i-2] / diag[i-2]
        ld1[i-1] = (d1[i-1] - ld2[i-2] * ld1[i-2] * diag[i-2]) / diag[i-1]
        diag[i] = d0[i] - ld1[i-1]^2 * diag[i-1] - ld2[i-2]^2 * diag[i-2]
    end

    # Step 5: Forward solve L*y' = rhs
    @inbounds z64[1] = rhs[1]
    if n >= 2
        @inbounds z64[2] = rhs[2] - ld1[1] * z64[1]
    end
    @inbounds for i in 3:n
        z64[i] = rhs[i] - ld1[i-1] * z64[i-1] - ld2[i-2] * z64[i-2]
    end

    # Step 6: Diagonal solve D*y'' = y'
    @inbounds for i in 1:n
        z64[i] = z64[i] / diag[i]
    end

    # Step 7: Back solve L'*z = y''
    @inbounds for i in n-1:-1:1
        z64[i] -= ld1[i] * z64[i+1]
        if i + 2 <= n
            z64[i] -= ld2[i] * z64[i+2]
        end
    end

    # Step 8: Convert to Float32 and clamp negatives
    @inbounds for i in 1:n
        z[i] = max(Float32(z64[i]), 0f0)
    end

    return @view z[1:n]
end
