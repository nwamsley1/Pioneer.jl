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
                      mask::Union{Nothing,AbstractVector{Bool}} = nothing,
                      ftr_scale_factor::Real = 1.0) where {U<:Real}

Return the minimum score threshold `τ` such that

```
FTR(τ) = (# of transfer decoys with score ≥ τ) / (# of targets with score ≥ τ)
```

is less than or equal to `alpha`.  The input vectors must be of equal
length and correspond element-wise to candidate precursor–run pairs.

# Arguments
- `ftr_scale_factor`: Correction factor for runtime decoy purging.
  When decoys are purged, use `1.0 / decoy_fraction` to correct FTR estimates.
  Default is 1.0 (no correction).
"""
function get_ftr_threshold(scores::AbstractVector{U},
                           is_bad_transfer::AbstractVector{Bool},
                           alpha::Real; doSort::Bool = true,
                           mask::Union{Nothing,AbstractVector{Bool}} = nothing,
                           ftr_scale_factor::Real = 1.0) where {U<:Real}
    @assert length(scores) == length(is_bad_transfer)

    if mask === nothing
        order = doSort ? sortperm(scores, rev=true, alg=QuickSort) : collect(eachindex(scores))
    else
        selected = findall(mask)
        order = doSort ? selected[sortperm(view(scores, selected), rev=true, alg=QuickSort)] : selected
    end

    # Handle empty input case
    if isempty(scores)
        return U(Inf)  # Return threshold that rejects all candidates
    end

    num_transfers = 0
    num_bad_transfers = 0
    τ = maximum(scores)
    best_count = 0
    for idx in order
        num_transfers += 1
        num_bad_transfers += is_bad_transfer[idx] ? 1 : 0

        # Apply FTR scale factor to correct for missing decoys
        if (num_transfers > 0) && (((num_bad_transfers / num_transfers) * ftr_scale_factor) <= alpha)
            τ = scores[idx]
            best_count = num_transfers
        end
    end

    return τ
end

"""
    get_ftr!(scores::AbstractVector{U},
             is_target::AbstractVector{Bool},
             is_transfer_decoy::AbstractVector{Bool},
             ftrs::AbstractVector{T}; doSort::Bool = true,
             mask::Union{Nothing,AbstractVector{Bool}} = nothing,
             ftr_scale_factor::Real = 1.0) where {U<:Real,T<:Real}

Compute the empirical FTR curve.  `ftrs[i]` contains the FTR for the
candidate at `scores[i]`.

# Arguments
- `ftr_scale_factor`: Correction factor for runtime decoy purging.
  When decoys are purged, use `1.0 / decoy_fraction` to correct FTR estimates.
  Default is 1.0 (no correction).
"""
function get_ftr!(scores::AbstractVector{U},
                  is_target::AbstractVector{Bool},
                  is_transfer_decoy::AbstractVector{Bool},
                  ftrs::AbstractVector{T}; doSort::Bool = true,
                  mask::Union{Nothing,AbstractVector{Bool}} = nothing,
                  ftr_scale_factor::Real = 1.0) where {U<:Real,T<:Real}
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
        target_cum += 1  # Count ALL candidates at this threshold level
        transfer_cum += is_transfer_decoy[idx] ? 1 : 0
        # Apply FTR scale factor to correct for missing decoys
        ftrs[idx] = target_cum > 0 ? ((transfer_cum / target_cum) * ftr_scale_factor) : Inf
    end
    return ftrs
end