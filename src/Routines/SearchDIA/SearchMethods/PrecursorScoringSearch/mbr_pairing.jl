# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# MBR pairing.
#
# Builds the experiment-global counterfactual-partner map used by MBR
# Batch F. For every unique `precursor_idx` observed across the per-file
# second-pass Arrow tables, assigns exactly one partner of the OPPOSITE
# target/decoy class — the one whose predicted iRT is closest within the
# same (cv_fold × prec_mz decile) stratum.
#
# The result is a `Dict{UInt32, UInt32}` (pid → partner_pid) consumed by
# the MBR streaming feature compute as the lookup for each row's _false
# donor (see mbr_streaming.jl). Nothing here mutates the spectral library
# or any PSM column.
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

# Binary-search the pool for the pid whose irt is closest to `target_irt`.
# Returns UInt32(0) if the pool is empty.
@inline function _nearest_pid(pool::_IrtPool, target_irt::Float32)
    n = length(pool.irts)
    n == 0 && return UInt32(0)
    idx = searchsortedfirst(pool.irts, target_irt)
    if idx > n
        return pool.pids[n]
    elseif idx == 1
        return pool.pids[1]
    else
        d_left  = abs(target_irt - pool.irts[idx - 1])
        d_right = abs(pool.irts[idx] - target_irt)
        return d_left <= d_right ? pool.pids[idx - 1] : pool.pids[idx]
    end
end

# Three-tier nearest-iRT fallback: in-stratum → same fold (any mz) →
# experiment-global. Returns (partner_pid, tier) where tier ∈ (:stratum,
# :fb_fold, :fb_global). Throws if no opposite-class precursor exists
# anywhere — that experiment is degenerate.
function _assign_partner_with_fallback(
    my_pid::UInt32, my_irt::Float32, my_fold::Int, my_mz::Int, is_target::Bool,
    pools_t::Dict, pools_d::Dict,
    fold_pool_t::Dict, fold_pool_d::Dict,
    global_pool_t::_IrtPool, global_pool_d::_IrtPool,
)
    opp_local  = is_target ?
        get(pools_d, (my_fold, my_mz), _empty_pool()) :
        get(pools_t, (my_fold, my_mz), _empty_pool())
    partner = _nearest_pid(opp_local, my_irt)
    partner != UInt32(0) && return (partner, :stratum)

    opp_fold = is_target ? fold_pool_d[my_fold] : fold_pool_t[my_fold]
    partner = _nearest_pid(opp_fold, my_irt)
    partner != UInt32(0) && return (partner, :fb_fold)

    opp_global = is_target ? global_pool_d : global_pool_t
    partner = _nearest_pid(opp_global, my_irt)
    partner == UInt32(0) && error(
        "build_counterfactual_partner_map: no opposite-class partner anywhere " *
        "for pid=$my_pid (target=$is_target)."
    )
    return (partner, :fb_global)
end

"""
    build_counterfactual_partner_map(file_paths, precursors) -> Dict{UInt32, UInt32}

MBR Batch F: build the experiment-global pid → counterfactual_partner_pid
map by streaming the per-file Arrow tables (no DataFrame materialised).

For every unique `precursor_idx` observed across `file_paths`, the partner
is the OPPOSITE-class precursor whose predicted iRT is closest within the
same (cv_fold × prec_mz decile) stratum. Falls back to same-fold-any-mz
and then experiment-global if a stratum is one-sided. Errors if no
opposite-class precursor exists anywhere.

Returns the partner map (consumed by mbr_streaming.jl for each row's
_false donor lookup). Hard requirement: every observed pid gets a
partner (asserted).
"""
function build_counterfactual_partner_map(
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
)
    t0 = time()
    unique = _collect_unique_precursors_streaming(file_paths, precursors)
    n_precs = length(unique.pids)

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

    pid_to_partner = Dict{UInt32, UInt32}()
    sizehint!(pid_to_partner, n_precs)
    n_stratum = 0; n_fb_fold = 0; n_fb_global = 0
    @inbounds for i in 1:n_precs
        partner, tier = _assign_partner_with_fallback(
            unique.pids[i], unique.irt[i],
            Int(unique.fold[i]), bin_mz[i], unique.target[i],
            pools_t, pools_d, fold_pool_t, fold_pool_d,
            global_pool_t, global_pool_d,
        )
        pid_to_partner[unique.pids[i]] = partner
        tier === :stratum   ? (n_stratum   += 1) :
        tier === :fb_fold   ? (n_fb_fold   += 1) : (n_fb_global += 1)
    end

    elapsed = time() - t0
    pct_strat = round(100 * n_stratum   / max(n_precs, 1), digits = 2)
    pct_fbf   = round(100 * n_fb_fold   / max(n_precs, 1), digits = 3)
    pct_fbg   = round(100 * n_fb_global / max(n_precs, 1), digits = 3)

    @user_info "MBR Batch F (iRT-NN) — counterfactual partner map:"
    @user_info "  unique precursors: $n_precs"
    @user_info "  partner = closest-pred-iRT opposite within (cv_fold × mz_decile)"
    @user_info "  in stratum:                     $n_stratum ($pct_strat%)"
    @user_info "  fallback (cv_fold, any mz):     $n_fb_fold ($pct_fbf%)"
    @user_info "  fallback (experiment-global):   $n_fb_global ($pct_fbg%)"
    @user_info "  total partnered: $(length(pid_to_partner)) / $n_precs"
    @user_info "  build_counterfactual_partner_map elapsed: $(round(elapsed, digits = 2))s"

    @assert length(pid_to_partner) == n_precs "Counterfactual partner map missing $(n_precs - length(pid_to_partner)) precursors"

    return pid_to_partner
end
