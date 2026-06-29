# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# MBR pairing.
#
# Builds experiment-global counterfactual-partner pools used by MBR
# Batch F. For every unique `precursor_idx` observed across the per-file
# second-pass Arrow tables, stores iRT pools needed to choose a different
# precursor as the counterfactual partner for each receiver file.
#
# MBR streaming resolves the partner against donor availability: first the
# highest-scoring donor from a different precursor within a local iRT window in
# the same cv_fold × prec_mz decile, then nearest-iRT broader fallbacks. Nothing here
# mutates the spectral library or any PSM column.
#
# Memory: streams the per-file Arrow tables read-only — never materialises
# a concatenated PSM DataFrame. Only (:precursor_idx, :target, :cv_fold)
# columns are touched; everything else (mz, irt) comes from the library.

using Statistics

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
function _compute_mz_decile_edges(mzs::Vector{Float32})
    finite = filter(isfinite, mzs)
    edges = isempty(finite) ? Float32[0f0, 1f0] :
                Float32.(unique(quantile(finite, 0:0.1:1)))
    if length(edges) < 11
        mn, mx = isempty(finite) ? (0f0, 1f0) : extrema(finite)
        edges = collect(Float32, LinRange(mn, mx, 11))
    end
    return edges
end

function _assign_mz_deciles(mzs::Vector{Float32}, edges::Vector{Float32})
    bins = Vector{Int}(undef, length(mzs))
    @inbounds for i in eachindex(mzs)
        v = mzs[i]
        bins[i] = isfinite(v) ? clamp(searchsortedlast(edges, v), 1, 10) : 5
    end
    return bins
end

function _compute_mz_deciles(mzs::Vector{Float32})
    return _assign_mz_deciles(mzs, _compute_mz_decile_edges(mzs))
end

# A pool is a (pids, irts) NamedTuple with irts sorted ascending so we
# can do binary search for local-window and iRT-nearest lookups.
const _IrtPool = NamedTuple{(:pids, :irts), Tuple{Vector{UInt32}, Vector{Float32}}}

_empty_pool() = (pids = UInt32[], irts = Float32[])

struct _CounterfactualPartnerPools
    target_by_pid::Vector{Bool}
    fold_by_pid::Vector{UInt8}
    mz_bin_by_pid::Vector{UInt8}
    irt_by_pid::Vector{Float32}
    pools::Dict{Tuple{Int, Int}, _IrtPool}
    fold_pool::Dict{Int, _IrtPool}
    global_pool::_IrtPool
end

# Build a sorted pool from a list of (irt, pid) tuples — sort by irt
# (tiebreak by pid for determinism), then split.
function _sort_to_pool(pairs::Vector{Tuple{Float32, UInt32}})::_IrtPool
    sort!(pairs; by = x -> (x[1], x[2]))
    return (pids = UInt32[x[2] for x in pairs],
            irts = Float32[x[1] for x in pairs])
end

# Per-(cv_fold, mz_decile) sorted iRT pools.
function _build_stratum_pools(
    pids::Vector{UInt32},
    fold::Vector{UInt8}, bin_mz::Vector{Int}, irt::Vector{Float32},
)
    tmp = Dict{Tuple{Int,Int}, Vector{Tuple{Float32, UInt32}}}()
    @inbounds for i in eachindex(pids)
        key = (Int(fold[i]), bin_mz[i])
        push!(get!(() -> Tuple{Float32, UInt32}[], tmp, key), (irt[i], pids[i]))
    end
    return Dict{Tuple{Int,Int}, _IrtPool}(k => _sort_to_pool(v) for (k, v) in tmp)
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

"""
    build_counterfactual_partner_pools(file_paths, precursors) -> _CounterfactualPartnerPools

MBR Batch F: build the precursor-indexed metadata and sorted iRT pools needed
for deterministic, file-aware counterfactual partner resolution.

For each receiver row, MBR streaming first chooses the highest-scoring donor
from a different precursor inside a local iRT window in the same
(cv_fold × prec_mz decile). If that primary stratum has no eligible cross-file
donor, it falls back through same-bin nearest iRT, same-fold nearest iRT, then
global nearest iRT.
"""
function build_counterfactual_partner_pools(
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    ;
    receiver_file_paths::Vector{String} = file_paths,
)
    pool_unique = _collect_unique_precursors_streaming(file_paths, precursors)
    receiver_unique = receiver_file_paths == file_paths ?
                      pool_unique :
                      _collect_unique_precursors_streaming(receiver_file_paths, precursors)
    n_pool_precs = length(pool_unique.pids)
    n_receiver_precs = length(receiver_unique.pids)
    if n_pool_precs == 0 && n_receiver_precs == 0
        empty_stratum = Dict{Tuple{Int, Int}, _IrtPool}()
        empty_fold = Dict{Int, _IrtPool}(0 => _empty_pool(), 1 => _empty_pool())
        return _CounterfactualPartnerPools(
            Bool[], UInt8[], UInt8[], Float32[],
            empty_stratum,
            empty_fold,
            _empty_pool(),
        )
    end

    mz_edges = _compute_mz_decile_edges(pool_unique.mz)
    pool_bin_mz = _assign_mz_deciles(pool_unique.mz, mz_edges)
    pools = _build_stratum_pools(
        pool_unique.pids, pool_unique.fold, pool_bin_mz, pool_unique.irt,
    )

    # Pre-build per-fold and experiment-global fallback pools. Two folds
    # (cv_fold ∈ {0, 1}) is the convention upstream.
    fold_pool = Dict{Int, _IrtPool}(
        0 => _build_fold_pool(pools, 0),
        1 => _build_fold_pool(pools, 1),
    )
    global_pool = _build_global_pool(fold_pool)

    max_pid = Int(maximum(receiver_unique.pids))
    target_by_pid = falses(max_pid)
    fold_by_pid = zeros(UInt8, max_pid)
    mz_bin_by_pid = zeros(UInt8, max_pid)
    irt_by_pid = zeros(Float32, max_pid)

    receiver_bin_mz = _assign_mz_deciles(receiver_unique.mz, mz_edges)
    @inbounds for i in 1:n_receiver_precs
        pid = Int(receiver_unique.pids[i])
        target_by_pid[pid] = receiver_unique.target[i]
        fold_by_pid[pid] = receiver_unique.fold[i]
        mz_bin_by_pid[pid] = UInt8(receiver_bin_mz[i])
        irt_by_pid[pid] = receiver_unique.irt[i]
    end

    return _CounterfactualPartnerPools(
        target_by_pid,
        fold_by_pid,
        mz_bin_by_pid,
        irt_by_pid,
        pools,
        fold_pool,
        global_pool,
    )
end
