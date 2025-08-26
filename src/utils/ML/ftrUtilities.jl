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
                           is_transfer_decoy::AbstractVector{Bool},
                           alpha::Real; doSort::Bool = true,
                           mask::Union{Nothing,AbstractVector{Bool}} = nothing) where {U<:Real}
    @assert length(scores) == length(is_target) == length(is_transfer_decoy)

    if mask === nothing
        order = doSort ? sortperm(scores, rev=true, alg=QuickSort) : collect(eachindex(scores))
    else
        selected = findall(mask)
        order = doSort ? selected[sortperm(view(scores, selected), rev=true, alg=QuickSort)] : selected
    end

    target_cum = 0
    transfer_cum = 0
    τ = maximum(scores)
    best_count = 0
    for idx in order
        target_cum += 1
        transfer_cum += is_transfer_decoy[idx] ? 1 : 0
        
        if target_cum > 0 && (transfer_cum / target_cum) <= alpha
            τ = scores[idx]
            best_count = target_cum
        end
    end

    #println(best_count, " ", τ, " ",  transfer_cum, " ",  target_cum, " ", transfer_cum / target_cum, " ", alpha, "\n\n")

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