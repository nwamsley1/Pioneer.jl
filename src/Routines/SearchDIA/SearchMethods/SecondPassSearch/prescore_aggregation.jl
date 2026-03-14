#==========================================================
Prescore Aggregation — Raw Log-Odds Combination

Combines per-file LightGBM prescore probabilities into a
single global score per precursor via log-odds averaging.

Uses raw LightGBM probabilities directly (no PEP calibration).
Simulation showed isotonic regression PEP calibration introduces
systematic anti-conservative FDR bias (mean inflation 1.29x at
1% FDR across 100 seeds) because each sample's own target/decoy
label influences its own PEP estimate. Raw log-odds of CV-scored
LightGBM probabilities is perfectly calibrated (mean inflation
0.999, see scripts/CLAUDE.md for full results).
==========================================================#

"""
    _logodds_combine(probs, top_n, floor) -> Float32

Take top-n values from `probs`, clamp to [floor, 1-ε],
convert to log-odds, average, and convert back to probability.

Optimized to avoid per-call allocations:
- Fast path for single-element vectors (scalar ops, no sort)
- Fast path for n==1 (find max, scalar ops)
- General case uses partialsort! (in-place) + scalar loop
"""
function _logodds_combine(probs::Vector{Float32}, top_n::Int, floor::Float32)::Float32
    isempty(probs) && return 0.0f0
    eps = 1f-6
    n = min(length(probs), top_n)
    if n == 1
        # Single value needed: find max without sorting
        p = length(probs) == 1 ? @inbounds(probs[1]) : maximum(probs)
        c = clamp(p, floor, 1 - eps)
        # logit → average(1 element) → sigmoid is identity
        return c
    end
    # In-place partial sort — safe because probs is a temporary consumed only once
    partialsort!(probs, 1:n; rev=true)
    lo_sum = 0.0f0
    @inbounds for i in 1:n
        c = clamp(probs[i], floor, 1 - eps)
        lo_sum += log(c / (1 - c))
    end
    avg = lo_sum / n
    return 1.0f0 / (1 + exp(-avg))
end
