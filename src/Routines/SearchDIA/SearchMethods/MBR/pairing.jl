# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

using Statistics

function _sort_mbr_irt_pool(
    pairs::Vector{Tuple{Float32, UInt32}},
)::_MBRIrtPool
    sort!(pairs; by = x -> (x[1], x[2]))
    return (
        pids = UInt32[x[2] for x in pairs],
        irts = Float32[x[1] for x in pairs],
    )
end

"""
    build_mbr_partner_pools(file_paths, precursors)

Build deterministic same-charge/same-length counterfactual pools. The
experiment-wide pool supports sparse-run fallback; the file-specific pool
lets false transfers use the exact same donor run as the real transfer.
"""
function build_mbr_partner_pools(
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
)
    n_precursors = length(precursors)
    charge_by_pid = zeros(UInt8, n_precursors)
    length_by_pid = zeros(UInt8, n_precursors)
    irt_by_pid = Float32.(getIrt(precursors))
    charges = getCharge(precursors)
    sequences = getSequence(precursors)
    @inbounds for pid in 1:n_precursors
        charge_by_pid[pid] = UInt8(charges[pid])
        length_by_pid[pid] = UInt8(min(length(sequences[pid]), typemax(UInt8)))
    end

    global_pairs = Dict{
        Tuple{Int, Int},
        Vector{Tuple{Float32, UInt32}},
    }()
    file_pairs = Dict{
        Tuple{UInt32, Int, Int},
        Vector{Tuple{Float32, UInt32}},
    }()

    for path in file_paths
        tbl = Arrow.Table(path)
        for col in (:precursor_idx, :ms_file_idx)
            hasproperty(tbl, col) ||
                error("MBR counterfactual pools require column $col in $path")
        end
        irt_col = hasproperty(tbl, :irt_pred) ? tbl.irt_pred : nothing
        @inbounds for row in eachindex(tbl.precursor_idx)
            pid = UInt32(tbl.precursor_idx[row])
            pid_int = Int(pid)
            file_idx = UInt32(tbl.ms_file_idx[row])
            irt = irt_col === nothing ?
                irt_by_pid[pid_int] :
                Float32(irt_col[row])
            isfinite(irt) || (irt = irt_by_pid[pid_int])
            key = (
                Int(charge_by_pid[pid_int]),
                Int(length_by_pid[pid_int]),
            )
            push!(
                get!(() -> Tuple{Float32, UInt32}[], global_pairs, key),
                (irt, pid),
            )
            push!(
                get!(
                    () -> Tuple{Float32, UInt32}[],
                    file_pairs,
                    (file_idx, key[1], key[2]),
                ),
                (irt, pid),
            )
        end
    end

    charge_length_pool = Dict{Tuple{Int, Int}, _MBRIrtPool}(
        key => _sort_mbr_irt_pool(values)
        for (key, values) in global_pairs
    )
    file_charge_length_pool =
        Dict{Tuple{UInt32, Int, Int}, _MBRIrtPool}(
            key => _sort_mbr_irt_pool(values)
            for (key, values) in file_pairs
        )
    return _MBRPartnerPools(
        charge_by_pid,
        length_by_pid,
        irt_by_pid,
        charge_length_pool,
        file_charge_length_pool,
    )
end

"""
    build_mbr_counterfactual_eligibility(file_paths; q_value_threshold)

Counterfactual precursor `p` is eligible in receiver run `r` only when `p`
passes the experiment-wide global q-value gate but does not already pass the
run-level q-value gate in `r`. That makes a counterfactual represent a
plausible *false transfer*, rather than an already-established identification.
"""
function build_mbr_counterfactual_eligibility(
    file_paths::Vector{String};
    q_value_threshold::Float32,
)
    global_passed = BitSet()
    run_passed_by_file = Dict{UInt32, BitSet}()
    for path in file_paths
        tbl = Arrow.Table(path)
        for col in (:precursor_idx, :ms_file_idx, :qval, :global_qval)
            hasproperty(tbl, col) ||
                error("MBR eligibility requires column $col in $path")
        end
        @inbounds for row in eachindex(tbl.precursor_idx)
            pid = Int(tbl.precursor_idx[row])
            file_idx = UInt32(tbl.ms_file_idx[row])
            global_qval = Float32(tbl.global_qval[row])
            run_qval = Float32(tbl.qval[row])
            if isfinite(global_qval) && global_qval <= q_value_threshold
                push!(global_passed, pid)
            end
            if isfinite(run_qval) && run_qval <= q_value_threshold
                push!(
                    get!(() -> BitSet(), run_passed_by_file, file_idx),
                    pid,
                )
            end
        end
    end
    return _MBRCounterfactualEligibility(
        global_passed,
        run_passed_by_file,
    )
end

@inline function _mbr_counterfactual_eligible(
    eligibility::_MBRCounterfactualEligibility,
    receiver_file::UInt32,
    pid::UInt32,
)
    pid_int = Int(pid)
    pid_int in eligibility.global_passed || return false
    run_passed = get(eligibility.run_passed_by_file, receiver_file, nothing)
    return run_passed === nothing || !(pid_int in run_passed)
end

@inline function _mbr_nearest_pool_order(
    pool::_MBRIrtPool,
    target_irt::Float32,
)
    n = length(pool.irts)
    n == 0 && return Int[]
    right = searchsortedfirst(pool.irts, target_irt)
    left = right - 1
    order = Int[]
    sizehint!(order, n)
    while left >= 1 || right <= n
        use_left = if right > n
            true
        elseif left < 1
            false
        else
            abs(target_irt - pool.irts[left]) <=
                abs(pool.irts[right] - target_irt)
        end
        push!(order, use_left ? left : right)
        use_left ? (left -= 1) : (right += 1)
    end
    return order
end
