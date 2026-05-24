# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# False Transfer Rate (FTR) utilities for match-between-runs.
#
# Ported verbatim from v0.6.6's src/utils/ML/ftrUtilities.jl (commit 081a1f72
# removed it). Self-contained: no Pioneer dependencies.

"""
    get_ftr_threshold(scores, is_bad_transfer, alpha; doSort=true, mask=nothing) -> τ

Return the minimum score threshold τ such that the empirical FTR

```
FTR(τ) = (# of bad transfers with score ≥ τ) / (# of all candidates with score ≥ τ)
```

is ≤ `alpha`. Walks scores descending; the largest τ for which FTR(τ) ≤ α
is returned. Use with the candidate set (e.g., `MBR_transfer_candidate=true`)
either by passing pre-filtered vectors or via the `mask` keyword.
"""
function get_ftr_threshold(scores::AbstractVector{U},
                           is_bad_transfer::AbstractVector{Bool},
                           alpha::Real;
                           doSort::Bool = true,
                           mask::Union{Nothing, AbstractVector{Bool}} = nothing
                          ) where {U <: Real}
    @assert length(scores) == length(is_bad_transfer)

    if mask === nothing
        order = doSort ? sortperm(scores; rev = true, alg = QuickSort) : collect(eachindex(scores))
    else
        selected = findall(mask)
        order = doSort ? selected[sortperm(view(scores, selected); rev = true, alg = QuickSort)] : selected
    end

    if isempty(scores) || isempty(order)
        return U(Inf)
    end

    n_total = 0
    n_bad   = 0
    τ       = U(Inf)
    for idx in order
        n_total += 1
        n_bad   += is_bad_transfer[idx] ? 1 : 0
        if (n_bad / n_total) <= alpha
            τ = scores[idx]
        end
    end
    return τ
end

"""
    get_ftr!(scores, is_target, is_bad_transfer, ftrs; doSort=true, mask=nothing)

Compute the empirical FTR curve. `ftrs[i]` ends up holding the FTR at
threshold `scores[i]`. `is_target` is currently unused; kept for
v0.6.6 compatibility / future per-class extensions.
"""
function get_ftr!(scores::AbstractVector{U},
                  ::AbstractVector{Bool},     # is_target (unused)
                  is_bad_transfer::AbstractVector{Bool},
                  ftrs::AbstractVector{T};
                  doSort::Bool = true,
                  mask::Union{Nothing, AbstractVector{Bool}} = nothing
                 ) where {U <: Real, T <: Real}
    @assert length(scores) == length(is_bad_transfer) == length(ftrs)

    if mask === nothing
        order = doSort ? sortperm(scores; rev = true, alg = QuickSort) : collect(eachindex(scores))
    else
        selected = findall(mask)
        order = doSort ? selected[sortperm(view(scores, selected); rev = true, alg = QuickSort)] : selected
    end

    n_total = 0
    n_bad   = 0
    for idx in order
        n_total += 1
        n_bad   += is_bad_transfer[idx] ? 1 : 0
        ftrs[idx] = n_total > 0 ? T(n_bad / n_total) : T(Inf)
    end
    return ftrs
end
