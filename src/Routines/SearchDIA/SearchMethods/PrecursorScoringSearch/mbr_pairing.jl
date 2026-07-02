# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# MBR pairing.
#
# Builds experiment-global counterfactual-partner pools used by MBR
# Batch F. For every unique `precursor_idx` observed across the per-file
# second-pass Arrow tables, stores shared same-charge, same-length iRT pools
# needed to choose the nearest valid non-self counterfactual partner for each
# receiver file.
#
# MBR streaming resolves the partner against donor availability: first the
# nearest non-self precursor in the same cv_fold × prec_mz decile × charge ×
# length stratum that has a cross-file donor, then same-fold/same-charge/
# same-length, then experiment-wide same-charge/same-length. Nothing here
# mutates the spectral library or any PSM column.
#
# Memory: streams the per-file Arrow tables read-only — never materialises
# a concatenated PSM DataFrame. Only (:precursor_idx, :target, :cv_fold)
# columns are touched; everything else (mz, charge, length, irt) comes from the
# library.

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
    prec_charge_full = getCharge(precursors)
    prec_length_full = getLength(precursors)

    seen         = Set{UInt32}()
    plist_pids   = UInt32[]
    plist_target = Bool[]
    plist_mz     = Float32[]
    plist_charge = UInt8[]
    plist_length = UInt8[]
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
            push!(plist_charge, UInt8(prec_charge_full[pid]))
            push!(plist_length, UInt8(prec_length_full[pid]))
            push!(plist_irt, Float32(prec_irt_full[pid]))
            push!(plist_fold, UInt8(fold_c[i]))
        end
    end
    return (pids = plist_pids, target = plist_target,
            mz = plist_mz, charge = plist_charge, length = plist_length,
            irt = plist_irt, fold = plist_fold)
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

struct _CounterfactualPartnerPools
    target_by_pid::Vector{Bool}
    fold_by_pid::Vector{UInt8}
    mz_bin_by_pid::Vector{UInt8}
    charge_by_pid::Vector{UInt8}
    length_by_pid::Vector{UInt8}
    irt_by_pid::Vector{Float32}
    pools::Dict{Tuple{Int, Int, Int, Int}, _IrtPool}
    fold_charge_length_pool::Dict{Tuple{Int, Int, Int}, _IrtPool}
    charge_length_pool::Dict{Tuple{Int, Int}, _IrtPool}
end

# Build a sorted pool from a list of (irt, pid) tuples — sort by irt
# (tiebreak by pid for determinism), then split.
function _sort_to_pool(pairs::Vector{Tuple{Float32, UInt32}})::_IrtPool
    sort!(pairs; by = x -> (x[1], x[2]))
    return (pids = UInt32[x[2] for x in pairs],
            irts = Float32[x[1] for x in pairs])
end

# Per-(cv_fold, mz_decile, charge, length) sorted iRT pools shared across
# targets and decoys.
function _build_stratum_pools(
    pids::Vector{UInt32}, fold::Vector{UInt8}, bin_mz::Vector{Int},
    charge::Vector{UInt8}, length::Vector{UInt8}, irt::Vector{Float32},
)
    tmp = Dict{Tuple{Int,Int,Int,Int}, Vector{Tuple{Float32, UInt32}}}()
    @inbounds for i in eachindex(pids)
        key = (Int(fold[i]), bin_mz[i], Int(charge[i]), Int(length[i]))
        push!(get!(() -> Tuple{Float32, UInt32}[], tmp, key), (irt[i], pids[i]))
    end
    return Dict{Tuple{Int,Int,Int,Int}, _IrtPool}(k => _sort_to_pool(v) for (k, v) in tmp)
end

# Same-fold, same-charge, same-length fallback pools.
function _build_fold_charge_length_pools(
    pids::Vector{UInt32}, fold::Vector{UInt8}, charge::Vector{UInt8},
    length::Vector{UInt8}, irt::Vector{Float32},
)
    tmp = Dict{Tuple{Int,Int,Int}, Vector{Tuple{Float32, UInt32}}}()
    @inbounds for i in eachindex(pids)
        key = (Int(fold[i]), Int(charge[i]), Int(length[i]))
        push!(get!(() -> Tuple{Float32, UInt32}[], tmp, key), (irt[i], pids[i]))
    end
    return Dict{Tuple{Int,Int,Int}, _IrtPool}(k => _sort_to_pool(v) for (k, v) in tmp)
end

# Experiment-wide same-charge, same-length fallback pools.
function _build_charge_length_pools(
    pids::Vector{UInt32}, charge::Vector{UInt8}, length::Vector{UInt8}, irt::Vector{Float32},
)
    tmp = Dict{Tuple{Int,Int}, Vector{Tuple{Float32, UInt32}}}()
    @inbounds for i in eachindex(pids)
        key = (Int(charge[i]), Int(length[i]))
        push!(get!(() -> Tuple{Float32, UInt32}[], tmp, key), (irt[i], pids[i]))
    end
    return Dict{Tuple{Int,Int}, _IrtPool}(k => _sort_to_pool(v) for (k, v) in tmp)
end

"""
    build_counterfactual_partner_pools(file_paths, precursors) -> _CounterfactualPartnerPools

    MBR Batch F: build the precursor-indexed metadata and sorted shared
    same-charge, same-length iRT pools needed for deterministic, file-aware
    counterfactual partner resolution.

    For each receiver row, MBR streaming chooses the nearest non-self precursor
    by predicted iRT that has a donor in a file other than the receiver file,
    checking same (cv_fold × prec_mz decile × charge × length), then same-fold
    same-charge same-length, then experiment-wide same-charge same-length.
"""
function build_counterfactual_partner_pools(
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
)
    unique = _collect_unique_precursors_streaming(file_paths, precursors)
    n_precs = length(unique.pids)
    if n_precs == 0
        return _CounterfactualPartnerPools(
            Bool[], UInt8[], UInt8[], UInt8[], UInt8[], Float32[],
            Dict{Tuple{Int, Int, Int, Int}, _IrtPool}(),
            Dict{Tuple{Int, Int, Int}, _IrtPool}(),
            Dict{Tuple{Int, Int}, _IrtPool}(),
        )
    end

    bin_mz = _compute_mz_deciles(unique.mz)
    pools = _build_stratum_pools(
        unique.pids, unique.fold, bin_mz, unique.charge, unique.length,
        unique.irt,
    )

    fold_charge_length_pool = _build_fold_charge_length_pools(
        unique.pids, unique.fold, unique.charge, unique.length, unique.irt,
    )
    charge_length_pool = _build_charge_length_pools(
        unique.pids, unique.charge, unique.length, unique.irt,
    )

    max_pid = Int(maximum(unique.pids))
    target_by_pid = falses(max_pid)
    fold_by_pid = zeros(UInt8, max_pid)
    mz_bin_by_pid = zeros(UInt8, max_pid)
    charge_by_pid = zeros(UInt8, max_pid)
    length_by_pid = zeros(UInt8, max_pid)
    irt_by_pid = zeros(Float32, max_pid)

    @inbounds for i in 1:n_precs
        pid = Int(unique.pids[i])
        target_by_pid[pid] = unique.target[i]
        fold_by_pid[pid] = unique.fold[i]
        mz_bin_by_pid[pid] = UInt8(bin_mz[i])
        charge_by_pid[pid] = unique.charge[i]
        length_by_pid[pid] = unique.length[i]
        irt_by_pid[pid] = unique.irt[i]
    end

    return _CounterfactualPartnerPools(
        target_by_pid,
        fold_by_pid,
        mz_bin_by_pid,
        charge_by_pid,
        length_by_pid,
        irt_by_pid,
        pools,
        fold_charge_length_pool,
        charge_length_pool,
    )
end
