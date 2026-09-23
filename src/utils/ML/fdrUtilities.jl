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
FDR (False Discovery Rate) and q-value calculation utilities.

This module provides functions for calculating FDR and q-values using the 
target-decoy approach, with support for library target/decoy ratio correction.
"""

function _score_order(scores)
    sortperm(scores; by=_score_key, rev=true)
end

function _get_qvalues_from_order!(scores, labels, qvalues, order, scale)
    n = length(order)
    targets = decoys = 0
    i = 1
    while i <= n
        score = _score_key(scores[order[i]])
        j = i
        while j <= n && _score_key(scores[order[j]]) == score
            labels[order[j]] ? (targets += 1) : (decoys += 1)
            j += 1
        end
        q = _group_fdr(targets, decoys, scale)
        for k in i:(j-1)
            qvalues[order[k]] = q
        end
        i = j
    end
    minimum_q = Inf
    for k in n:-1:1
        row = order[k]
        minimum_q = min(minimum_q, qvalues[row])
        qvalues[row] = minimum_q
    end
    return nothing
end

function _score_array_calibration(scores, labels; kwargs...)
    build_score_calibration(emit -> _emit_score_arrays(emit, scores, labels); kwargs...)
end

"""
    get_qvalues!(scores, labels, qvalues; doSort=true, fdr_scale_factor=1.0f0,
                 memory_budget_bytes=SCORE_WORKSPACE_BYTES)

Assign each tied score the same q-value: count the whole score group, then take
suffix minima of scaled cumulative decoy/target ratios. `doSort=false` requires
scores ordered from best to worst. Large inputs use bounded external sorting;
the workspace budget excludes caller-owned arrays and fixed I/O overhead.
"""
function get_qvalues!(probs::AbstractVector{U}, labels::AbstractVector{Bool}, qvals::AbstractVector{T};
    doSort::Bool=true, fdr_scale_factor::Float32=1.0f0,
    memory_budget_bytes::Int=SCORE_WORKSPACE_BYTES,
) where {T,U<:AbstractFloat}
    length(probs) == length(labels) == length(qvals) || throw(DimensionMismatch("Score arrays differ in length"))
    _validate_score_options(fdr_scale_factor, memory_budget_bytes)
    isempty(probs) && return nothing
    if length(probs) > memory_budget_bytes ÷ 64
        fit = _score_array_calibration(probs, labels; compute_pep=false,
            fdr_scale_factor, memory_budget_bytes)
        try
            for i in 1:length(probs)
                qvals[i] = fit.qval_spline(probs[i])
            end
        finally
            _close_score_calibration(fit)
        end
    else
        order = doSort ? _score_order(probs) : (1:length(probs))
        _get_qvalues_from_order!(probs, labels, qvals, order, fdr_scale_factor)
    end
    return nothing
end

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
function get_PEP!(scores::AbstractVector{U}, is_target::AbstractVector{Bool}, peps::AbstractVector{T};
    doSort::Bool=true, fdr_scale_factor::Float32=1.0f0,
    memory_budget_bytes::Int=SCORE_WORKSPACE_BYTES,
) where {T,U<:AbstractFloat}
    length(scores) == length(is_target) == length(peps) || throw(DimensionMismatch("Score arrays differ in length"))
    _validate_score_options(fdr_scale_factor, memory_budget_bytes)
    isempty(scores) && return nothing
    if length(scores) > memory_budget_bytes ÷ 64
        fit = _score_array_calibration(scores, is_target; fdr_scale_factor, memory_budget_bytes)
        try
            for i in 1:length(scores)
                peps[i] = fit.pep_interp(scores[i])
            end
        finally
            _close_score_calibration(fit)
        end
    else
        order = doSort ? _score_order(scores) : (1:length(scores))
        _get_PEP_from_order!(scores, is_target, peps, order, fdr_scale_factor; memory_budget_bytes)
    end
    return nothing
end

"""Fit equal-score groups with weighted PAVA, preserving the zero-valued
unit-weight pseudocount before the first group. `order` is descending by score.
"""
function _get_PEP_from_order!(scores, labels, peps, order, scale;
    memory_budget_bytes::Int=SCORE_WORKSPACE_BYTES)
    length(scores) == length(labels) == length(peps) == length(order) ||
        throw(DimensionMismatch("Score arrays differ in length"))
    _validate_score_options(scale, memory_budget_bytes)
    if length(scores) > memory_budget_bytes ÷ 64
        return get_PEP!(scores, labels, peps; fdr_scale_factor=Float32(scale), memory_budget_bytes)
    end
    # Block sizes count observations, including all observations in a tied group.
    blocks = ScorePAVABlock[ScorePAVABlock(0.0, 1.0, 0)]
    i = 1
    while i <= length(order)
        score = _score_key(scores[order[i]])
        targets = decoys = 0
        j = i
        while j <= length(order) && _score_key(scores[order[j]]) == score
            labels[order[j]] ? (targets += 1) : (decoys += 1)
            j += 1
        end
        d = Float64(decoys) * scale
        block = ScorePAVABlock(d, targets + d, j-i)
        while !isempty(blocks) && last(blocks).decoys / last(blocks).weight > block.decoys / block.weight
            top = pop!(blocks)
            block = ScorePAVABlock(top.decoys + block.decoys, top.weight + block.weight, top.groups + block.groups)
        end
        push!(blocks, block)
        i = j
    end
    position = 1
    for block in blocks
        pep = clamp(block.decoys / (block.weight - block.decoys), 0.0, 1.0)
        for _ in 1:block.groups
            peps[order[position]] = pep
            position += 1
        end
    end
    return nothing
end

"""Compute q-values and PEPs using one shared score ordering."""
function get_score_statistics!(scores, labels, qvalues, peps;
    fdr_scale_factor::Float32=1.0f0, memory_budget_bytes::Int=SCORE_WORKSPACE_BYTES)
    length(scores) == length(labels) == length(qvalues) == length(peps) ||
        throw(DimensionMismatch("Score arrays differ in length"))
    _validate_score_options(fdr_scale_factor, memory_budget_bytes)
    isempty(scores) && return nothing
    if length(scores) > memory_budget_bytes ÷ 64
        fit = _score_array_calibration(scores, labels; fdr_scale_factor, memory_budget_bytes)
        try
            for i in 1:length(scores)
                qvalues[i] = fit.qval_spline(scores[i])
                peps[i] = fit.pep_interp(scores[i])
            end
        finally
            _close_score_calibration(fit)
        end
    else
        order = _score_order(scores)
        _get_qvalues_from_order!(scores, labels, qvalues, order, fdr_scale_factor)
        _get_PEP_from_order!(scores, labels, peps, order, fdr_scale_factor; memory_budget_bytes)
    end
    return nothing
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
