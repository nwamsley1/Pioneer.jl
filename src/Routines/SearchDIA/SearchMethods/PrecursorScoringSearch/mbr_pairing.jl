# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# MBR pairing.
#
# Builds experiment-global counterfactual-partner pools used by MBR
# Batch F. The shared same-length/same-charge iRT pools contain every finite
# refined-iRT observation observed across the per-file second-pass Arrow tables,
# so a precursor seen in multiple runs can be paired by its minimum refined-iRT
# distance to the receiver row.
#
# MBR streaming resolves the partner against donor availability by choosing
# the nearest non-self same-length/same-charge precursor by refined iRT that
# has a cross-file donor.
# Nothing here mutates the spectral library or any PSM column.
#
# Memory: streams the per-file Arrow tables read-only — never materialises
# a concatenated PSM DataFrame. Pool iRT comes from row-level :irt_pred when
# present and finite, falling back to raw library iRT otherwise; mz, charge,
# and length come from the library.

# Streams the per-file Arrow tables and returns one entry per unique pid with
# target flag, library m/z, charge, length, iRT, and CV-fold metadata.
function _collect_unique_precursors_streaming(
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    ;
    unique_precursors::Bool = true,
    irt_column::Symbol = :irt_pred,
    require_irt_column::Bool = false,
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
    plist_file   = UInt32[]

    for fpath in file_paths
        tbl = Arrow.Table(fpath)
        n = length(tbl.precursor_idx)
        n == 0 && continue
        pid_c    = tbl.precursor_idx
        file_c   = hasproperty(tbl, :ms_file_idx) ? tbl.ms_file_idx : nothing
        target_c = hasproperty(tbl, :target) ? tbl.target : trues(n)
        fold_c   = hasproperty(tbl, :cv_fold) ? tbl.cv_fold : zeros(UInt8, n)
        irt_c = hasproperty(tbl, irt_column) ? getproperty(tbl, irt_column) : nothing
        if irt_c === nothing && require_irt_column
            error("MBR counterfactual pairing requires column $irt_column in $fpath")
        end
        @inbounds for i in 1:n
            pid = UInt32(pid_c[i])
            if unique_precursors
                pid in seen && continue
                push!(seen, pid)
            end
            file_idx = file_c === nothing ? UInt32(0) : UInt32(file_c[i])
            raw_irt = Float32(prec_irt_full[pid])
            refined_irt = raw_irt
            if irt_c !== nothing
                candidate_irt = Float32(irt_c[i])
                if isfinite(candidate_irt)
                    refined_irt = candidate_irt
                elseif require_irt_column
                    error("MBR counterfactual pairing found non-finite $irt_column for precursor $pid in $fpath")
                end
            elseif require_irt_column
                error("MBR counterfactual pairing requires column $irt_column in $fpath")
            end
            push!(plist_pids, pid)
            push!(plist_target, Bool(target_c[i]))
            push!(plist_mz, Float32(prec_mz_full[pid]))
            push!(plist_charge, UInt8(prec_charge_full[pid]))
            push!(plist_length, UInt8(prec_length_full[pid]))
            push!(plist_irt, refined_irt)
            push!(plist_fold, UInt8(fold_c[i]))
            push!(plist_file, file_idx)
        end
    end
    return (pids = plist_pids, target = plist_target, mz = plist_mz,
            charge = plist_charge, length = plist_length, irt = plist_irt,
            fold = plist_fold, file_idx = plist_file)
end

# A pool is a (pids, irts) NamedTuple with refined iRT observations sorted
# ascending so we can do binary search for iRT-nearest lookups. The same pid may
# appear multiple times if it was observed with different refined iRTs in
# different runs.
const _IrtPool = NamedTuple{(:pids, :irts), Tuple{Vector{UInt32}, Vector{Float32}}}

struct _CounterfactualEligibilityByFile
    dataset_passed::BitSet
    run_passed_by_file::Dict{UInt32, BitSet}
end

_empty_irt_pool() = (pids = UInt32[], irts = Float32[])
_empty_pool() = _empty_irt_pool()

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
    file_charge_length_pool::Dict{Tuple{UInt32, Int, Int}, _IrtPool}
end

_CounterfactualPartnerPools(
    target_by_pid,
    fold_by_pid,
    mz_bin_by_pid,
    charge_by_pid,
    length_by_pid,
    irt_by_pid,
    pools,
    fold_charge_length_pool,
    charge_length_pool,
) = _CounterfactualPartnerPools(
    target_by_pid,
    fold_by_pid,
    mz_bin_by_pid,
    charge_by_pid,
    length_by_pid,
    irt_by_pid,
    pools,
    fold_charge_length_pool,
    charge_length_pool,
    Dict{Tuple{UInt32, Int, Int}, _IrtPool}(),
)

_CounterfactualPartnerPools(
    mz_by_pid::Vector{Float32},
    charge_by_pid::Vector{UInt8},
    irt_by_pid::Vector{Float32},
    charge_pool::Dict{Int, <:NamedTuple},
) = _CounterfactualPartnerPools(
    falses(length(mz_by_pid)),
    zeros(UInt8, length(mz_by_pid)),
    ones(UInt8, length(mz_by_pid)),
    charge_by_pid,
    zeros(UInt8, length(mz_by_pid)),
    irt_by_pid,
    Dict{Tuple{Int, Int, Int, Int}, _IrtPool}(),
    Dict{Tuple{Int, Int, Int}, _IrtPool}(),
    Dict{Tuple{Int, Int}, _IrtPool}(),
)

function _compute_mz_decile_edges(mzs::Vector{Float32})
    finite = filter(isfinite, mzs)
    if isempty(finite)
        return Float32.(collect(LinRange(0.0f0, 1.0f0, 11)))
    end
    sorted = sort(finite)
    n = length(sorted)
    edges = Vector{Float32}(undef, 11)
    @inbounds for i in 0:10
        idx = clamp(round(Int, 1 + i * (n - 1) / 10), 1, n)
        edges[i + 1] = sorted[idx]
    end
    if length(unique(edges)) < 11
        mn, mx = extrema(finite)
        edges = Float32.(collect(LinRange(mn, mx, 11)))
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

function _sort_to_pool(pairs::Vector{Tuple{Float32, UInt32}})::_IrtPool
    sort!(pairs; by = x -> (x[1], x[2]))
    return (pids = UInt32[x[2] for x in pairs],
            irts = Float32[x[1] for x in pairs])
end

function _build_stratum_pools(
    pids::Vector{UInt32},
    fold::Vector{UInt8},
    bin_mz::Vector{Int},
    charge::Vector{UInt8},
    length::Vector{UInt8},
    irt::Vector{Float32},
)
    tmp = Dict{Tuple{Int, Int, Int, Int}, Vector{Tuple{Float32, UInt32}}}()
    @inbounds for i in eachindex(pids)
        key = (Int(fold[i]), bin_mz[i], Int(charge[i]), Int(length[i]))
        push!(get!(() -> Tuple{Float32, UInt32}[], tmp, key), (irt[i], pids[i]))
    end
    return Dict{Tuple{Int, Int, Int, Int}, _IrtPool}(k => _sort_to_pool(v) for (k, v) in tmp)
end

function _build_fold_charge_length_pools(
    pids::Vector{UInt32},
    fold::Vector{UInt8},
    charge::Vector{UInt8},
    length::Vector{UInt8},
    irt::Vector{Float32},
)
    tmp = Dict{Tuple{Int, Int, Int}, Vector{Tuple{Float32, UInt32}}}()
    @inbounds for i in eachindex(pids)
        key = (Int(fold[i]), Int(charge[i]), Int(length[i]))
        push!(get!(() -> Tuple{Float32, UInt32}[], tmp, key), (irt[i], pids[i]))
    end
    return Dict{Tuple{Int, Int, Int}, _IrtPool}(k => _sort_to_pool(v) for (k, v) in tmp)
end

function _build_charge_length_pools(
    pids::Vector{UInt32}, charge::Vector{UInt8}, length::Vector{UInt8}, irt::Vector{Float32},
)
    tmp = Dict{Tuple{Int, Int}, Vector{Tuple{Float32, UInt32}}}()
    @inbounds for i in eachindex(pids)
        key = (Int(charge[i]), Int(length[i]))
        push!(get!(() -> Tuple{Float32, UInt32}[], tmp, key), (irt[i], pids[i]))
    end
    return Dict{Tuple{Int, Int}, _IrtPool}(k => _sort_to_pool(v) for (k, v) in tmp)
end

function _build_file_charge_length_pools(
    pids::Vector{UInt32},
    file_idx::Vector{UInt32},
    charge::Vector{UInt8},
    length::Vector{UInt8},
    irt::Vector{Float32},
)
    tmp = Dict{Tuple{UInt32, Int, Int}, Vector{Tuple{Float32, UInt32}}}()
    @inbounds for i in eachindex(pids)
        key = (file_idx[i], Int(charge[i]), Int(length[i]))
        push!(get!(() -> Tuple{Float32, UInt32}[], tmp, key), (irt[i], pids[i]))
    end
    return Dict{Tuple{UInt32, Int, Int}, _IrtPool}(k => _sort_to_pool(v) for (k, v) in tmp)
end

"""
    build_counterfactual_partner_pools(file_paths, precursors) -> _CounterfactualPartnerPools

    MBR Batch F: build the precursor-indexed metadata and sorted shared
    same-length/same-charge iRT pools needed for deterministic, file-aware counterfactual
    partner resolution.

    For each receiver row, MBR streaming chooses the nearest non-self precursor
    by refined iRT that has the same length and charge and has a donor in a file
    other than the receiver file.
"""
function build_counterfactual_partner_pools(
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    ;
    receiver_file_paths::Vector{String} = file_paths,
    irt_column::Symbol = :irt_pred,
    require_irt_column::Bool = false,
)
    pool_unique = _collect_unique_precursors_streaming(
        file_paths,
        precursors;
        unique_precursors = false,
        irt_column = irt_column,
        require_irt_column = require_irt_column,
    )
    receiver_unique = _collect_unique_precursors_streaming(
        receiver_file_paths,
        precursors;
        unique_precursors = true,
        irt_column = irt_column,
        require_irt_column = require_irt_column,
    )
    n_pool_precs = length(pool_unique.pids)
    n_receiver_precs = length(receiver_unique.pids)
    if n_pool_precs == 0 && n_receiver_precs == 0
        return _CounterfactualPartnerPools(
            Bool[], UInt8[], UInt8[], UInt8[], UInt8[], Float32[],
            Dict{Tuple{Int, Int, Int, Int}, _IrtPool}(),
            Dict{Tuple{Int, Int, Int}, _IrtPool}(),
            Dict{Tuple{Int, Int}, _IrtPool}(),
            Dict{Tuple{UInt32, Int, Int}, _IrtPool}(),
        )
    end

    mz_edges = _compute_mz_decile_edges(pool_unique.mz)
    pool_bin_mz = _assign_mz_deciles(pool_unique.mz, mz_edges)
    pools = _build_stratum_pools(
        pool_unique.pids, pool_unique.fold, pool_bin_mz, pool_unique.charge,
        pool_unique.length, pool_unique.irt,
    )
    fold_charge_length_pool = _build_fold_charge_length_pools(
        pool_unique.pids, pool_unique.fold, pool_unique.charge,
        pool_unique.length, pool_unique.irt,
    )
    charge_length_pool = _build_charge_length_pools(
        pool_unique.pids, pool_unique.charge, pool_unique.length, pool_unique.irt,
    )
    file_charge_length_pool = _build_file_charge_length_pools(
        pool_unique.pids,
        pool_unique.file_idx,
        pool_unique.charge,
        pool_unique.length,
        pool_unique.irt,
    )

    max_pid = Int(maximum(receiver_unique.pids))
    target_by_pid = falses(max_pid)
    fold_by_pid = zeros(UInt8, max_pid)
    mz_bin_by_pid = zeros(UInt8, max_pid)
    charge_by_pid = zeros(UInt8, max_pid)
    length_by_pid = zeros(UInt8, max_pid)
    irt_by_pid = zeros(Float32, max_pid)

    receiver_bin_mz = _assign_mz_deciles(receiver_unique.mz, mz_edges)
    @inbounds for i in 1:n_receiver_precs
        pid = Int(receiver_unique.pids[i])
        target_by_pid[pid] = receiver_unique.target[i]
        fold_by_pid[pid] = receiver_unique.fold[i]
        mz_bin_by_pid[pid] = UInt8(receiver_bin_mz[i])
        charge_by_pid[pid] = receiver_unique.charge[i]
        length_by_pid[pid] = receiver_unique.length[i]
        irt_by_pid[pid] = receiver_unique.irt[i]
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
        file_charge_length_pool,
    )
end

function build_post_integration_counterfactual_partner_pools(
    file_paths::Vector{String},
    precursors::LibraryPrecursors;
    receiver_file_paths::Vector{String} = file_paths,
)
    return build_counterfactual_partner_pools(
        file_paths,
        precursors;
        receiver_file_paths = receiver_file_paths,
        irt_column = :irt_pred,
        require_irt_column = true,
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
