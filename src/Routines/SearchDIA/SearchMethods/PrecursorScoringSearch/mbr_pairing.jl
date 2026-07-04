# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# MBR pairing.
#
# Builds experiment-global counterfactual-partner pools used by MBR
# Batch F. For every unique `precursor_idx` observed across the per-file
# second-pass Arrow tables, stores shared same-length/same-charge iRT pools
# needed to choose the nearest valid non-self counterfactual partner for each
# receiver file.
#
# MBR streaming resolves the partner against donor availability by choosing
# the nearest non-self same-length/same-charge precursor by library iRT that
# has a cross-file donor.
# Nothing here
# mutates the spectral library or any PSM column.
#
# Memory: streams the per-file Arrow tables read-only — never materialises
# a concatenated PSM DataFrame. Only :precursor_idx is touched; everything
# else (mz, charge, irt) comes from the library.

# Streams the per-file Arrow tables and returns one entry per unique pid with
# library m/z, charge, and iRT metadata.
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
    plist_mz     = Float32[]
    plist_charge = UInt8[]
    plist_length = UInt8[]
    plist_irt    = Float32[]

    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        n = length(tbl.precursor_idx)
        n == 0 && continue
        pid_c    = tbl.precursor_idx
        @inbounds for i in 1:n
            pid = UInt32(pid_c[i])
            pid in seen && continue
            push!(seen, pid)
            push!(plist_pids, pid)
            push!(plist_mz, Float32(prec_mz_full[pid]))
            push!(plist_charge, UInt8(prec_charge_full[pid]))
            push!(plist_length, UInt8(prec_length_full[pid]))
            push!(plist_irt, Float32(prec_irt_full[pid]))
        end
    end
    return (pids = plist_pids, mz = plist_mz, charge = plist_charge,
            length = plist_length, irt = plist_irt)
end

# A pool is a (pids, mzs) NamedTuple with m/z values sorted ascending so we can
# do binary search for m/z-nearest lookups.
const _MzPool = NamedTuple{(:pids, :mzs), Tuple{Vector{UInt32}, Vector{Float32}}}

# A pool is a (pids, irts) NamedTuple with library iRT values sorted ascending
# so we can do binary search for iRT-nearest lookups.
const _IrtPool = NamedTuple{(:pids, :irts), Tuple{Vector{UInt32}, Vector{Float32}}}

struct _CounterfactualEligibilityByFile
    dataset_passed::BitSet
    run_passed_by_file::Dict{UInt32, BitSet}
end

_empty_pool() = (pids = UInt32[], mzs = Float32[])
_empty_irt_pool() = (pids = UInt32[], irts = Float32[])

struct _CounterfactualPartnerPools
    mz_by_pid::Vector{Float32}
    charge_by_pid::Vector{UInt8}
    length_by_pid::Vector{UInt8}
    irt_by_pid::Vector{Float32}
    charge_pool::Dict{Int, _MzPool}
    charge_length_irt_pool::Dict{Tuple{Int, Int}, _IrtPool}
end

_CounterfactualPartnerPools(
    mz_by_pid::Vector{Float32},
    charge_by_pid::Vector{UInt8},
    irt_by_pid::Vector{Float32},
    charge_pool::Dict{Int, _MzPool},
) = _CounterfactualPartnerPools(
    mz_by_pid,
    charge_by_pid,
    UInt8[],
    irt_by_pid,
    charge_pool,
    Dict{Tuple{Int, Int}, _IrtPool}(),
)

# Build a sorted pool from a list of (mz, pid) tuples — sort by m/z
# (tiebreak by pid for determinism), then split.
function _sort_to_pool(pairs::Vector{Tuple{Float32, UInt32}})::_MzPool
    sort!(pairs; by = x -> (x[1], x[2]))
    return (pids = UInt32[x[2] for x in pairs],
            mzs = Float32[x[1] for x in pairs])
end

function _sort_to_irt_pool(pairs::Vector{Tuple{Float32, UInt32}})::_IrtPool
    sort!(pairs; by = x -> (x[1], x[2]))
    return (pids = UInt32[x[2] for x in pairs],
            irts = Float32[x[1] for x in pairs])
end

# Experiment-wide same-charge fallback pools.
function _build_charge_pools(
    pids::Vector{UInt32}, charge::Vector{UInt8}, mz::Vector{Float32},
)
    tmp = Dict{Int, Vector{Tuple{Float32, UInt32}}}()
    @inbounds for i in eachindex(pids)
        key = Int(charge[i])
        push!(get!(() -> Tuple{Float32, UInt32}[], tmp, key), (mz[i], pids[i]))
    end
    return Dict{Int, _MzPool}(k => _sort_to_pool(v) for (k, v) in tmp)
end

function _build_charge_length_irt_pools(
    pids::Vector{UInt32},
    charge::Vector{UInt8},
    length::Vector{UInt8},
    irt::Vector{Float32},
)
    tmp = Dict{Tuple{Int, Int}, Vector{Tuple{Float32, UInt32}}}()
    @inbounds for i in eachindex(pids)
        key = (Int(charge[i]), Int(length[i]))
        push!(get!(() -> Tuple{Float32, UInt32}[], tmp, key), (irt[i], pids[i]))
    end
    return Dict{Tuple{Int, Int}, _IrtPool}(k => _sort_to_irt_pool(v) for (k, v) in tmp)
end

"""
    build_counterfactual_partner_pools(file_paths, precursors) -> _CounterfactualPartnerPools

    MBR Batch F: build the precursor-indexed metadata and sorted shared
    same-length/same-charge iRT pools needed for deterministic, file-aware counterfactual
    partner resolution.

    For each receiver row, MBR streaming chooses the nearest non-self precursor
    by library iRT that has the same length and charge and has a donor in a
    file other than the receiver file.
"""
function build_counterfactual_partner_pools(
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
)
    unique = _collect_unique_precursors_streaming(file_paths, precursors)
    n_precs = length(unique.pids)
    if n_precs == 0
        return _CounterfactualPartnerPools(
            Float32[], UInt8[], UInt8[], Float32[],
            Dict{Int, _MzPool}(),
            Dict{Tuple{Int, Int}, _IrtPool}(),
        )
    end

    charge_pool = _build_charge_pools(
        unique.pids, unique.charge, unique.mz,
    )
    charge_length_irt_pool = _build_charge_length_irt_pools(
        unique.pids, unique.charge, unique.length, unique.irt,
    )

    max_pid = Int(maximum(unique.pids))
    mz_by_pid = zeros(Float32, max_pid)
    charge_by_pid = zeros(UInt8, max_pid)
    length_by_pid = zeros(UInt8, max_pid)
    irt_by_pid = zeros(Float32, max_pid)

    @inbounds for i in 1:n_precs
        pid = Int(unique.pids[i])
        mz_by_pid[pid] = unique.mz[i]
        charge_by_pid[pid] = unique.charge[i]
        length_by_pid[pid] = unique.length[i]
        irt_by_pid[pid] = unique.irt[i]
    end

    return _CounterfactualPartnerPools(
        mz_by_pid,
        charge_by_pid,
        length_by_pid,
        irt_by_pid,
        charge_pool,
        charge_length_irt_pool,
    )
end

"""
    build_counterfactual_receiver_eligibility(file_paths; q_value_threshold = 0.01f0)
        -> Dict{UInt32, BitSet}

Build the compact membership indexes used to decide whether a precursor can
act as a counterfactual transfer target in a receiver run. A counterfactual
precursor is eligible when it passed the full-dataset global q-value gate but
did not pass the run-level q-value gate in that receiver run.
"""
function build_counterfactual_receiver_eligibility(
    file_paths::Vector{String};
    q_value_threshold::Float32 = 0.01f0,
)::_CounterfactualEligibilityByFile
    dataset_passed = BitSet()
    run_passed_by_file = Dict{UInt32, BitSet}()
    for path in file_paths
        tbl = Arrow.Table(path)
        for col in (:precursor_idx, :ms_file_idx, :qval, :global_qval)
            hasproperty(tbl, col) ||
                error("Post-qvalue MBR counterfactual eligibility requires column $col in $path")
        end
        n = length(tbl.precursor_idx)
        @inbounds for i in 1:n
            pid = Int(tbl.precursor_idx[i])
            file_key = UInt32(tbl.ms_file_idx[i])
            if Float32(tbl.global_qval[i]) <= q_value_threshold
                push!(dataset_passed, pid)
            end
            if Float32(tbl.qval[i]) <= q_value_threshold
                push!(get!(() -> BitSet(), run_passed_by_file, file_key), pid)
            end
        end
    end
    return _CounterfactualEligibilityByFile(dataset_passed, run_passed_by_file)
end
