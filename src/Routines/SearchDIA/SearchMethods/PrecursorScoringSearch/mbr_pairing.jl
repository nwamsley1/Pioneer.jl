# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# MBR pairing.
#
# Builds experiment-global counterfactual-partner candidates used by MBR
# Batch F. For every unique `precursor_idx` observed across the per-file
# second-pass Arrow tables, stores a short ranked list of OPPOSITE-class
# partners ordered by tier (same cv_fold × prec_mz decile, same fold, global)
# and then nearest predicted iRT.
#
# The result is a zero-padded `Matrix{UInt32}` with rows = candidate rank
# and columns = precursor_idx. MBR streaming chooses the first candidate
# with a donor row from a file other than the receiver file. Nothing here
# mutates the spectral library or any PSM column.
#
# Memory: streams the per-file Arrow tables read-only — never materialises
# a concatenated PSM DataFrame. Only (:precursor_idx, :target, :cv_fold)
# columns are touched; everything else (mz, irt) comes from the library.

using Statistics

const MBR_COUNTERFACTUAL_MAX_PARTNER_CANDIDATES = 8

# Streams the per-file Arrow tables, projects just the three columns we
# need, and returns one entry per unique pid: target flag, library mz +
# irt, and cv_fold (the first cv_fold encountered — it's constant per
# pid because folds are precursor-keyed in the fold-split files).
function _collect_unique_precursors_streaming(
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
)
    prec_mz_full  = getMz(precursors)
    prec_irt_full = getIrt(precursors)

    seen         = Set{UInt32}()
    plist_pids   = UInt32[]
    plist_target = Bool[]
    plist_mz     = Float32[]
    plist_irt    = Float32[]
    plist_fold   = UInt8[]

    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        n = length(tbl.precursor_idx)
        n == 0 && continue
        pid_c    = tbl.precursor_idx
        target_c = tbl.target
        fold_c   = tbl.cv_fold
        @inbounds for i in 1:n
            pid = UInt32(pid_c[i])
            pid in seen && continue
            push!(seen, pid)
            push!(plist_pids, pid)
            push!(plist_target, Bool(target_c[i]))
            push!(plist_mz, Float32(prec_mz_full[pid]))
            push!(plist_irt, Float32(prec_irt_full[pid]))
            push!(plist_fold, UInt8(fold_c[i]))
        end
    end
    return (pids = plist_pids, target = plist_target,
            mz = plist_mz, irt = plist_irt, fold = plist_fold)
end

# 10 mz-quantile deciles. Returns Int8 bin index ∈ 1..10 for each input
# value. Non-finite mz lands in bin 5 (the median bin).
function _compute_mz_deciles(mzs::Vector{Float32})
    finite = filter(isfinite, mzs)
    edges = isempty(finite) ? Float32[0f0, 1f0] :
                Float32.(unique(quantile(finite, 0:0.1:1)))
    if length(edges) < 11
        mn, mx = isempty(finite) ? (0f0, 1f0) : extrema(finite)
        edges = collect(Float32, LinRange(mn, mx, 11))
    end
    bins = Vector{Int}(undef, length(mzs))
    @inbounds for i in eachindex(mzs)
        v = mzs[i]
        bins[i] = isfinite(v) ? clamp(searchsortedlast(edges, v), 1, 10) : 5
    end
    return bins
end

# A pool is a (pids, irts) NamedTuple with irts sorted ascending so we
# can do binary search for iRT-nearest lookups.
const _IrtPool = NamedTuple{(:pids, :irts), Tuple{Vector{UInt32}, Vector{Float32}}}

_empty_pool() = (pids = UInt32[], irts = Float32[])

# Build a sorted pool from a list of (irt, pid) tuples — sort by irt
# (tiebreak by pid for determinism), then split.
function _sort_to_pool(pairs::Vector{Tuple{Float32, UInt32}})::_IrtPool
    sort!(pairs; by = x -> (x[1], x[2]))
    return (pids = UInt32[x[2] for x in pairs],
            irts = Float32[x[1] for x in pairs])
end

# Per-(cv_fold, mz_decile) sorted iRT pools, one per target/decoy class.
function _build_stratum_pools(
    pids::Vector{UInt32}, target::Vector{Bool},
    fold::Vector{UInt8}, bin_mz::Vector{Int}, irt::Vector{Float32},
)
    tmp_t = Dict{Tuple{Int,Int}, Vector{Tuple{Float32, UInt32}}}()
    tmp_d = Dict{Tuple{Int,Int}, Vector{Tuple{Float32, UInt32}}}()
    @inbounds for i in eachindex(pids)
        key = (Int(fold[i]), bin_mz[i])
        dst = target[i] ? tmp_t : tmp_d
        push!(get!(() -> Tuple{Float32, UInt32}[], dst, key), (irt[i], pids[i]))
    end
    pools_t = Dict{Tuple{Int,Int}, _IrtPool}(k => _sort_to_pool(v) for (k, v) in tmp_t)
    pools_d = Dict{Tuple{Int,Int}, _IrtPool}(k => _sort_to_pool(v) for (k, v) in tmp_d)
    return pools_t, pools_d
end

# Collapse all (fold, *) stratum pools for a single fold into one sorted
# fold-wide pool — used as the first fallback when the in-stratum pool
# is empty.
function _build_fold_pool(stratum_pools::Dict{Tuple{Int,Int}, _IrtPool}, fold::Int)::_IrtPool
    pairs = Tuple{Float32, UInt32}[]
    for ((f, _), pool) in stratum_pools
        f == fold || continue
        @inbounds for k in eachindex(pool.pids)
            push!(pairs, (pool.irts[k], pool.pids[k]))
        end
    end
    return _sort_to_pool(pairs)
end

# Experiment-wide fallback pool (ignores fold + mz).
function _build_global_pool(fold_pools::Dict{Int, _IrtPool})::_IrtPool
    pairs = Tuple{Float32, UInt32}[]
    for (_, pool) in fold_pools
        @inbounds for k in eachindex(pool.pids)
            push!(pairs, (pool.irts[k], pool.pids[k]))
        end
    end
    return _sort_to_pool(pairs)
end

# Binary-search a sorted iRT pool, then expand outward to append nearest
# unmatched pids into `candidate_mat[:, col]`.
@inline function _candidate_already_added(
    candidate_mat::Matrix{UInt32},
    col::Int,
    n_added::Int,
    pid::UInt32,
)
    @inbounds for row in 1:n_added
        candidate_mat[row, col] == pid && return true
    end
    return false
end

function _push_nearest_partner_candidates!(
    candidate_mat::Matrix{UInt32},
    col::Int,
    n_added::Int,
    pool::_IrtPool,
    target_irt::Float32,
)
    max_candidates = size(candidate_mat, 1)
    n = length(pool.irts)
    n == 0 && return n_added

    right = searchsortedfirst(pool.irts, target_irt)
    left = right - 1
    @inbounds while n_added < max_candidates && (left >= 1 || right <= n)
        use_left = if right > n
            true
        elseif left < 1
            false
        else
            abs(target_irt - pool.irts[left]) <= abs(pool.irts[right] - target_irt)
        end
        pool_idx = use_left ? left : right
        use_left ? (left -= 1) : (right += 1)
        pid = pool.pids[pool_idx]
        _candidate_already_added(candidate_mat, col, n_added, pid) && continue
        n_added += 1
        candidate_mat[n_added, col] = pid
    end
    return n_added
end

function _assign_partner_candidates_with_fallback!(
    candidate_mat::Matrix{UInt32},
    col::Int,
    my_pid::UInt32,
    my_irt::Float32,
    my_fold::Int,
    my_mz::Int,
    is_target::Bool,
    pools_t::Dict,
    pools_d::Dict,
    fold_pool_t::Dict,
    fold_pool_d::Dict,
    global_pool_t::_IrtPool,
    global_pool_d::_IrtPool,
)
    n_added = 0

    opp_local = is_target ?
        get(pools_d, (my_fold, my_mz), _empty_pool()) :
        get(pools_t, (my_fold, my_mz), _empty_pool())
    n_added = _push_nearest_partner_candidates!(candidate_mat, col, n_added, opp_local, my_irt)

    opp_fold = is_target ? fold_pool_d[my_fold] : fold_pool_t[my_fold]
    n_added = _push_nearest_partner_candidates!(candidate_mat, col, n_added, opp_fold, my_irt)

    opp_global = is_target ? global_pool_d : global_pool_t
    n_added = _push_nearest_partner_candidates!(candidate_mat, col, n_added, opp_global, my_irt)

    n_added == 0 && error(
        "build_counterfactual_partner_candidate_matrix: no opposite-class partner anywhere " *
        "for pid=$my_pid (target=$is_target)."
    )
    return n_added
end

"""
    build_counterfactual_partner_candidate_matrix(file_paths, precursors; max_candidates=8) -> Matrix{UInt32}

MBR Batch F: build a zero-padded precursor-indexed matrix of ranked
counterfactual partner candidates by streaming the per-file Arrow tables
(no DataFrame materialised).

For every unique `precursor_idx` observed across `file_paths`, candidates
are OPPOSITE-class precursors ordered by nearest predicted iRT within the
same (cv_fold × prec_mz decile) stratum, then same-fold-any-mz fallback,
then experiment-global fallback.
"""
function build_counterfactual_partner_candidate_matrix(
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    max_candidates::Int = MBR_COUNTERFACTUAL_MAX_PARTNER_CANDIDATES,
)
    unique = _collect_unique_precursors_streaming(file_paths, precursors)
    n_precs = length(unique.pids)
    max_candidates > 0 || error("max counterfactual partner candidates must be positive")
    n_precs == 0 && return zeros(UInt32, max_candidates, 0)

    bin_mz = _compute_mz_deciles(unique.mz)
    pools_t, pools_d = _build_stratum_pools(
        unique.pids, unique.target, unique.fold, bin_mz, unique.irt,
    )

    # Pre-build per-fold and experiment-global fallback pools. Two folds
    # (cv_fold ∈ {0, 1}) is the convention upstream.
    fold_pool_t = Dict{Int, _IrtPool}(
        0 => _build_fold_pool(pools_t, 0),
        1 => _build_fold_pool(pools_t, 1),
    )
    fold_pool_d = Dict{Int, _IrtPool}(
        0 => _build_fold_pool(pools_d, 0),
        1 => _build_fold_pool(pools_d, 1),
    )
    global_pool_t = _build_global_pool(fold_pool_t)
    global_pool_d = _build_global_pool(fold_pool_d)

    partner_candidates = zeros(UInt32, max_candidates, Int(maximum(unique.pids)))
    @inbounds for i in 1:n_precs
        _assign_partner_candidates_with_fallback!(
            partner_candidates,
            Int(unique.pids[i]),
            unique.pids[i],
            unique.irt[i],
            Int(unique.fold[i]), bin_mz[i], unique.target[i],
            pools_t, pools_d, fold_pool_t, fold_pool_d,
            global_pool_t, global_pool_d,
        )
    end

    return partner_candidates
end
