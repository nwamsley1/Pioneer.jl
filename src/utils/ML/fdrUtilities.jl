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
    get_local_FDR!(scores::AbstractVector{U}, is_target::AbstractVector{Bool}, fdrs::AbstractVector{T}; 
                   window_size::Int=1000, doSort=true, fdr_scale_factor::Float32=1.0f0) where {T,U<:AbstractFloat}

Calculates local FDR (False Discovery Rate) within a sliding window.

# Arguments
- `scores`: Vector of scores (higher = better)
- `is_target`: Vector of target/decoy labels (true = target, false = decoy)
- `fdrs`: Vector to store calculated local FDRs
- `window_size`: Size of sliding window for local FDR calculation (default: 1000)
- `doSort`: Whether to sort by scores (default: true)
- `fdr_scale_factor`: Scale factor to correct for library target/decoy ratio (default: 1.0)

# Process
1. Sorts items by descending score (if doSort=true)
2. Builds prefix sums for efficient window calculations
3. For each position, calculates FDR within next window_size items
4. Applies FDR scale factor to correct for library imbalance

Local FDR provides a more granular estimate of FDR in specific score regions,
useful for filtering when global FDR may be too conservative.
"""
function get_local_FDR!(scores::AbstractVector{U}, is_target::AbstractVector{Bool}, fdrs::AbstractVector{T}; 
    window_size::Int=1000, doSort=true, fdr_scale_factor::Float32=1.0f0) where {T,U<:AbstractFloat}
    @assert length(scores) == length(is_target)
    N = length(scores)
    if N == 0
        return
    end
    # 1) Sort items by descending score
    if doSort
        idxs = sortperm(scores, rev = true, alg = QuickSort) #Sort class probabilities
    else
        idxs = eachindex(scores)
    end
    # We'll define rank i as the i-th item in that sorted order
    # so rank 1 => idxs[1], rank N => idxs[N].

    # 2) Build prefix sums for decoys & targets in sorted order
    #    This lets us quickly count how many decoys/targets are in a range.
    decoy_prefix = zeros(Int, N+1)   # decoy_prefix[i] = # decoys among top i items
    target_prefix = zeros(Int, N+1)  # same for targets

    decoy_prefix[1]  = !is_target[1]
    target_prefix[1] = is_target[1]

    @inbounds @fastmath for rank in 2:N
        i = idxs[rank] 
        decoy_prefix[rank]  = decoy_prefix[rank-1]  + !is_target[i]
        target_prefix[rank] = target_prefix[rank-1] + is_target[i]
    end

    # 3) For each rank i, compute local FDR from this point to the next X entries
    fdrs[idxs[1]] = ((decoy_prefix[window_size] + 1) * fdr_scale_factor) / max(1, target_prefix[window_size])

    for rank in 2:(N-window_size)
        # Count decoys/targets in [L,R] using prefix sums
        decs_in_window   = decoy_prefix[rank + window_size]  - decoy_prefix[rank-1]
        targs_in_window  = target_prefix[rank + window_size] - target_prefix[rank-1]

        # local FDR = (#decoys * scale_factor) / max(1, #targets)
        # adding a pseudocount to guarantee finite sample control of the FDR
        fdrs[idxs[rank]] = ((decs_in_window + 1) * fdr_scale_factor) / max(1, targs_in_window)
    end

    return 
end


# No need to export since this is included directly in the Pioneer module