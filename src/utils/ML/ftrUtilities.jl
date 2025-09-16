"""
False Transfer Rate (FTR) utilities.

This module provides helpers to estimate score thresholds for
match-between-runs workflows using transfer decoys.
"""

"""
    get_ftr_threshold(scores::AbstractVector{U},
                      is_target::AbstractVector{Bool},
                      is_transfer_decoy::AbstractVector{Bool},
                      alpha::Real; doSort::Bool = true,
                      mask::Union{Nothing,AbstractVector{Bool}} = nothing) where {U<:Real}

Return the minimum score threshold `τ` such that

```
FTR(τ) = (# of transfer decoys with score ≥ τ) / (# of targets with score ≥ τ)
```

is less than or equal to `alpha`.  The input vectors must be of equal
length and correspond element-wise to candidate precursor–run pairs.
"""
function get_ftr_threshold(scores::AbstractVector{U},
                           is_target::AbstractVector{Bool},
                           is_bad_transfer::AbstractVector{Bool},
                           alpha::Real; doSort::Bool = true,
                           mask::Union{Nothing,AbstractVector{Bool}} = nothing) where {U<:Real}
    @assert length(scores) == length(is_target) == length(is_bad_transfer)

    if mask === nothing
        order = doSort ? sortperm(scores, rev=true, alg=QuickSort) : collect(eachindex(scores))
    else
        selected = findall(mask)
        order = doSort ? selected[sortperm(view(scores, selected), rev=true, alg=QuickSort)] : selected
    end

    num_transfers = 0
    num_bad_transfers = 0
    τ = maximum(scores)
    best_count = 0
    for idx in order
        num_transfers += 1
        num_bad_transfers += is_bad_transfer[idx] ? 1 : 0
        
        if (num_transfers > 0) && ((num_bad_transfers / num_transfers) <= alpha)
            τ = scores[idx]
            best_count = num_transfers
        end
    end

    println("FTR probability threshold: ", τ, " Num passing candidate transfers: ", best_count, " out of ", num_transfers, "\n\n")

    return τ
end

"""
    get_ftr!(scores::AbstractVector{U},
             is_target::AbstractVector{Bool},
             is_transfer_decoy::AbstractVector{Bool},
             ftrs::AbstractVector{T}; doSort::Bool = true,
             mask::Union{Nothing,AbstractVector{Bool}} = nothing) where {U<:Real,T<:Real}

Compute the empirical FTR curve.  `ftrs[i]` contains the FTR for the
candidate at `scores[i]`.
"""
function get_ftr!(scores::AbstractVector{U},
                  is_target::AbstractVector{Bool},
                  is_transfer_decoy::AbstractVector{Bool},
                  ftrs::AbstractVector{T}; doSort::Bool = true,
                  mask::Union{Nothing,AbstractVector{Bool}} = nothing) where {U<:Real,T<:Real}
    @assert length(scores) == length(is_target) == length(is_transfer_decoy) == length(ftrs)

    if mask === nothing
        order = doSort ? sortperm(scores, rev=true, alg=QuickSort) : collect(eachindex(scores))
    else
        selected = findall(mask)
        order = doSort ? selected[sortperm(view(scores, selected), rev=true, alg=QuickSort)] : selected
    end

    target_cum = 0
    transfer_cum = 0
    for idx in order
        target_cum += is_target[idx] ? 1 : 0
        transfer_cum += is_transfer_decoy[idx] ? 1 : 0
        ftrs[idx] = target_cum > 0 ? (transfer_cum / target_cum) : Inf
    end
    return ftrs
end