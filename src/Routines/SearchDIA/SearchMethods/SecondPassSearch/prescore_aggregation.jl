#==========================================================
Prescore Aggregation Strategies

Trait-based system for combining per-file LightGBM prescore
probabilities into a single global score per precursor.
==========================================================#

"""Abstract type for how per-file prescore probabilities are aggregated globally."""
abstract type PrescoreAggregationStrategy end

"""
    RawLogOddsAggregation <: PrescoreAggregationStrategy

Aggregate raw LightGBM probabilities via log-odds averaging of top-√n values.
Fast but assumes cross-file calibration of LightGBM outputs.
"""
struct RawLogOddsAggregation <: PrescoreAggregationStrategy end

"""
    PEPCalibratedAggregation <: PrescoreAggregationStrategy

First convert per-file scores to calibrated PEP, then aggregate
(1 - PEP) via log-odds averaging. Theoretically grounded because
PEP is a calibrated local probability, removing cross-file bias.
"""
struct PEPCalibratedAggregation <: PrescoreAggregationStrategy end

#==========================================================
Per-file score calibration
==========================================================#

"""
    calibrate_file_scores(::RawLogOddsAggregation, probs, targets) -> Vector{Float32}

Identity — return raw LightGBM probabilities unchanged.
"""
function calibrate_file_scores(::RawLogOddsAggregation,
                               probs::Vector{Float32},
                               targets::AbstractVector{Bool})
    return probs
end

"""
    calibrate_file_scores(::PEPCalibratedAggregation, probs, targets) -> Vector{Float32}

Convert per-file LightGBM probabilities to calibrated (1 - PEP).
Falls back to raw probabilities when target/decoy counts are too low
for reliable PEP estimation.
"""
function calibrate_file_scores(::PEPCalibratedAggregation,
                               probs::Vector{Float32},
                               targets::AbstractVector{Bool})
    n_targets = count(targets)
    n_decoys  = count(!, targets)

    # Need enough targets and decoys for reliable PEP estimation
    if n_targets < 50 || n_decoys < 20
        return probs
    end

    peps = Vector{Float32}(undef, length(probs))
    get_PEP!(probs, targets, peps; doSort=true)

    # Return 1 - PEP (calibrated probability of being correct)
    calibrated = Vector{Float32}(undef, length(peps))
    @inbounds for i in eachindex(peps)
        calibrated[i] = 1.0f0 - peps[i]
    end
    return calibrated
end

#==========================================================
Cross-file score combination
==========================================================#

"""
    combine_scores(::RawLogOddsAggregation, probs, top_n) -> Float32

Log-odds average of top-n raw probabilities with floor clamp at 0.1.
"""
function combine_scores(::RawLogOddsAggregation, probs::Vector{Float32}, top_n::Int)::Float32
    return _logodds_combine(probs, top_n, 0.1f0)
end

"""
    combine_scores(::PEPCalibratedAggregation, probs, top_n) -> Float32

Log-odds average of top-n calibrated (1-PEP) values with floor clamp at 0.01.
Lower floor is safe because PEP calibration removes the noise that motivates
the conservative 0.1 floor in the raw approach.
"""
function combine_scores(::PEPCalibratedAggregation, probs::Vector{Float32}, top_n::Int)::Float32
    return _logodds_combine(probs, top_n, 0.01f0)
end

"""
    _logodds_combine(probs, top_n, floor) -> Float32

Shared implementation: take top-n values, clamp to [floor, 1-ε],
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
