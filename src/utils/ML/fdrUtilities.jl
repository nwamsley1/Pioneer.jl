"""
FDR (False Discovery Rate) and q-value calculation utilities.

This module provides functions for calculating FDR and q-values using the 
target-decoy approach, with support for library target/decoy ratio correction.
"""

"""
    get_qvalues!(probs::Vector{U}, labels::Vector{Bool}, qvals::Vector{T}; 
                 doSort::Bool=true, fdr_scale_factor::Float32=1.0f0) where {T,U<:AbstractFloat}
    get_qvalues!(PSMs::DataFrame, probs::Vector{Float64}, labels::Vector{Bool})

Calculates q-values (false discovery rate estimates) for PSMs.

# Arguments
- `probs`: Vector of probability scores
- `labels`: Vector of target/decoy labels (true = target, false = decoy)
- `qvals`: Vector to store calculated q-values
- `doSort`: Whether to sort by probability scores (default: true)
- `fdr_scale_factor`: Scale factor to correct for library target/decoy ratio (default: 1.0)

# Process
1. Sorts PSMs by probability score (if doSort=true)
2. Calculates running ratio of decoys to targets
3. Applies FDR scale factor to correct for library imbalance
4. Assigns q-values based on corrected decoy/target ratio
5. Ensures q-values are monotonically non-increasing

Implements target-decoy approach for FDR estimation with library ratio correction.
"""
function get_qvalues!(probs::AbstractVector{U}, labels::AbstractVector{Bool}, qvals::AbstractVector{T}; 
                      doSort::Bool = true, fdr_scale_factor::Float32 = 1.0f0
) where {T,U<:AbstractFloat}

    if doSort
        order = sortperm(probs, rev = true, alg=QuickSort) #Sort class probabilities
    else
        order = eachindex(probs)
    end

    targets = 0
    decoys = 1 # pseudocount to guarantee finite sample control of the FDR
    @inbounds @fastmath for i in order
            targets += labels[i]
            decoys += (1 - labels[i])
            # Apply FDR scale factor to correct for library target/decoy ratio
            qvals[i] = (decoys * fdr_scale_factor) / targets
    end

    fdr = Inf
    @inbounds @fastmath for i in reverse(order)
        if qvals[i] > fdr
            qvals[i] = fdr
        else
            fdr = qvals[i]
        end
    end
end

# DataFrame convenience method
get_qvalues!(PSMs::DataFrame, probs::Vector{Float64}, labels::Vector{Bool}) = get_qvalues!(PSMs, allowmissing(probs), allowmissing(labels))

"""
    get_PEP!(scores::AbstractVector{U}, is_target::AbstractVector{Bool}, fdrs::AbstractVector{T};
                    doSort=true, fdr_scale_factor::Float32=1.0f0) where {T,U<:AbstractFloat}
    Estimate posterior error probability using isotonic regression.
    The function fits a non-decreasing decoy probability curve over the score
    distribution via a weighted Pool Adjacent Violators Algorithm (PAVA). Each
    decoy observation is weighted by `fdr_scale_factor` to correct for library
    target/decoy imbalance. The fitted decoy probabilities are converted to PEP
    values (decoy_prob/(1 - decoy_prob)).
    # Arguments
    - `scores`: Vector of scores (higher = better)
    - `is_target`: Vector of target/decoy labels (true = target, false = decoy)
    - `fdrs`: Vector to store calculated local PEP values
    - `doSort`: Whether to sort by scores before fitting (default: true)
    - `fdr_scale_factor`: Scale factor to correct for library target/decoy ratio
"""
function get_PEP!(scores::AbstractVector{U}, is_target::AbstractVector{Bool}, fdrs::AbstractVector{T};
    doSort::Bool=true, fdr_scale_factor::Float32=1.0f0) where {T,U<:AbstractFloat}

    @assert length(scores) == length(is_target)
    N = length(scores)
    if N == 0
        return
    end

    # sort by score if requested
    order = doSort ? sortperm(scores, rev=true, alg=QuickSort) : collect(eachindex(scores))

    # prepare labels and weights
    labels = Vector{Float64}(undef, N + 1)
    weights = Vector{Float64}(undef, N + 1)
    labels[1] = 0.5              # pseudo observation
    weights[1] = 1.0

    @inbounds for j in 1:N
        idx = order[j]
        # decoy = 1.0, target = 0.0
        labels[j+1] = is_target[idx] ? 0.0 : 1.0
        weights[j+1] = is_target[idx] ? 1.0 : float(fdr_scale_factor)
    end

    fitted = _weighted_pava(labels, weights)
    fitted = fitted[2:end]  # remove pseudo row

    pep = fitted ./ (1 .- fitted)
    pep = clamp.(pep, 0.0, 1.0)

    @inbounds for (j, idx) in enumerate(order)
        fdrs[idx] = T(pep[j])
    end
    return
end

""" Weighted pool adjacent violators algorithm used by `get_PEP!`."""
function _weighted_pava(y::Vector{Float64}, w::Vector{Float64})
    n = length(y)
    # preallocate stacks of maximum size n
    v    = Vector{Float64}(undef, n)
    wt   = Vector{Float64}(undef, n)
    lenv = Vector{Int}(undef,    n)

    m = 0  # current stack height
    for i in 1:n
        # “push” y[i], w[i], 1 onto our stack
        m += 1
        @inbounds begin
            v[m]    = y[i]
            wt[m]   = w[i]
            lenv[m] = 1
        end

        # merge down while out of order
        while m > 1 && v[m-1] > v[m]
            @inbounds begin
                new_w   = wt[m-1] + wt[m]
                new_v   = (v[m-1]*wt[m-1] + v[m]*wt[m]) / new_w
                new_len = lenv[m-1] + lenv[m]
                m -= 1
                v[m]    = new_v
                wt[m]   = new_w
                lenv[m] = new_len
            end
        end
    end

    # now expand back out into result
    result = Vector{Float64}(undef, n)
    idx = 1
    for j in 1:m
        @inbounds for _ in 1:lenv[j]
            result[idx] = v[j]
            idx += 1
        end
    end

    return result
end


# No need to export since this is included directly in the Pioneer module