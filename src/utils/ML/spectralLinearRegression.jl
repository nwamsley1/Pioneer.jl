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

abstract type RegularizationType end
struct L2Norm <: RegularizationType end
struct NoNorm <: RegularizationType end

function solveOLS!(Hs::SparseArray{Ti, T},
                    r::Vector{T},
                    X₁::Vector{T},
                    colnorm2::Vector{T},
                    max_iter_outer::Int64,
                    relative_convergence_threshold::T) where {Ti<:Integer,T<:AbstractFloat}

    # Precompute column norms (sum of squared entries per column)
    @inbounds for col in 1:Hs.n
        s = zero(T)
        for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
            s += Hs.nzval[i]^2
        end
        colnorm2[col] = s
    end

    max_weight = T(0)
    i = 0
    while i < max_iter_outer
        _diff = T(0)
        weight_floor = i >= 5 ? max_weight * T(1e-7) : T(0)
        max_weight = T(0)

        for col in 1:Hs.n
            cn2 = colnorm2[col]
            iszero(cn2) && continue

            # Compute gradient (X^T r for this column)
            grad = zero(T)
            @inbounds @fastmath for j in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
                grad += Hs.nzval[j] * r[Hs.rowval[j]]
            end

            # Exact coordinate update: new_x = x - grad / ||col||^2, clamped >= 0
            X0 = X₁[col]
            X₁[col] = max(X0 - grad / cn2, zero(T))

            # Update residuals in-place
            δx_raw = X₁[col] - X0
            if !iszero(δx_raw)
                @inbounds @fastmath for j in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
                    r[Hs.rowval[j]] += Hs.nzval[j] * δx_raw
                end
            end

            # Track max weight
            if X₁[col] > max_weight
                max_weight = X₁[col]
            end

            # Convergence tracking
            if X₁[col] > weight_floor
                rel_change = abs(δx_raw) / abs(X₁[col])
                if rel_change > _diff
                    _diff = rel_change
                end
            end
        end

        if _diff < relative_convergence_threshold
            return (true, i)
        end
        i += 1
    end

    return (false, max_iter_outer)
end