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
    build_partitioned_index_from_lib(spec_lib; partition_width=5.0f0,
        frag_bin_tol_ppm=2.5f0, rt_bin_tol=3.0f0,
        rank_to_score=UInt8[8,4,4,2,2,1,1],
        y_start_index=UInt8(4), b_start_index=UInt8(3))

Build a LocalPartitionedFragmentIndex from scratch using the spectral library's
DetailedFrag data. Each partition gets its own independently-constructed index
with LocalFragment entries (UInt16 local_id + UInt8 score = 4 bytes).

Precursor IDs are remapped to partition-local UInt16 values (1..N, N ≤ 65535).
Partitions that would exceed 65535 unique precursors are automatically split.
"""
function build_partitioned_index_from_lib(
    spec_lib::SpectralLibrary;
    partition_width::Float32 = 5.0f0,
    frag_bin_tol_ppm::Float32 = 2.5f0,
    rt_bin_tol::Float32 = 3.0f0,
    rank_to_score::Vector{UInt8} = UInt8[8, 4, 4, 2, 2, 1, 1],
    y_start_index::UInt8 = UInt8(4),
    b_start_index::UInt8 = UInt8(3),
)
    precursors = getPrecursors(spec_lib)
    frag_lookup = getFragmentLookupTable(spec_lib)
    detailed_frags = getFragments(frag_lookup)
    prec_mzs = getMz(precursors)
    prec_irts = getIrt(precursors)
    n_precursors = length(prec_mzs)
    max_rank = length(rank_to_score)

    # ── Step 1: Assign precursors to initial partitions by prec_mz ───────────
    min_prec_mz = Float32(Inf)
    max_prec_mz = Float32(-Inf)
    for i in 1:n_precursors
        pmz = prec_mzs[i]
        min_prec_mz = min(min_prec_mz, pmz)
        max_prec_mz = max(max_prec_mz, pmz)
    end
    n_initial = max(1, ceil(Int, (max_prec_mz - min_prec_mz) / partition_width))
    @debug_l2 "build_partitioned_index: prec m/z [$(round(min_prec_mz, digits=2)), $(round(max_prec_mz, digits=2))], $(partition_width) Da → $(n_initial) initial partitions"

    # Collect global precursor IDs per initial partition
    initial_partition_pids = [UInt32[] for _ in 1:n_initial]
    for pid in UInt32(1):UInt32(n_precursors)
        pmz = prec_mzs[pid]
        k = clamp(floor(Int, (pmz - min_prec_mz) / partition_width) + 1, 1, n_initial)
        push!(initial_partition_pids[k], pid)
    end

    # ── Step 2: Split partitions exceeding MAX_LOCAL_PRECS ───────────────────
    final_partition_pids = Vector{UInt32}[]
    for pids in initial_partition_pids
        if length(pids) <= MAX_LOCAL_PRECS
            push!(final_partition_pids, pids)
        else
            sort!(pids, by = pid -> prec_mzs[pid])
            for i in 1:MAX_LOCAL_PRECS:length(pids)
                chunk_end = min(i + MAX_LOCAL_PRECS - 1, length(pids))
                push!(final_partition_pids, pids[i:chunk_end])
            end
        end
    end
    n_partitions = length(final_partition_pids)
    if n_partitions != n_initial
        @debug_l2 "build_partitioned_index: split to $(n_partitions) partitions ($(n_partitions - n_initial) extra from UInt16 limit)"
    end

    # ── Step 3: Build SimpleFrags + local ID mapping per partition ────────────
    partition_frags = [SimpleFrag{Float32}[] for _ in 1:n_partitions]
    partition_local_to_global = [UInt32[] for _ in 1:n_partitions]
    global_to_local = Dict{UInt32, UInt16}()

    for k in 1:n_partitions
        pids = final_partition_pids[k]
        empty!(global_to_local)
        local_to_global = zeros(UInt32, length(pids))
        for (i, pid) in enumerate(pids)
            lid = UInt16(i)
            global_to_local[pid] = lid
            local_to_global[i] = pid
        end

        for pid in pids
            pmz = prec_mzs[pid]
            pirt = prec_irts[pid]
            frag_range = getPrecFragRange(frag_lookup, pid)
            lid = global_to_local[pid]

            rank = 0
            for fi in frag_range
                dfrag = detailed_frags[fi]

                # Apply fragment index ion-type filters
                if dfrag.is_y
                    dfrag.ion_position < y_start_index && continue
                elseif dfrag.is_b
                    dfrag.ion_position < b_start_index && continue
                elseif dfrag.is_p
                    continue
                end
                dfrag.is_isotope && continue

                rank += 1
                rank > max_rank && break

                push!(partition_frags[k], SimpleFrag{Float32}(
                    getMz(dfrag),
                    UInt32(lid),
                    pmz,
                    pirt,
                    UInt8(0),
                    rank_to_score[rank],
                ))
            end
        end

        partition_local_to_global[k] = local_to_global
    end

    total_frags = sum(length, partition_frags)
    @debug_l2 "build_partitioned_index: $(total_frags) total fragments across $(n_partitions) partitions"

    # ── Step 4: Build LocalPartition per partition ────────────────────────────
    partitions = Vector{LocalPartition{Float32}}(undef, n_partitions)

    for k in 1:n_partitions
        frags_k = partition_frags[k]
        l2g = partition_local_to_global[k]
        n_local = UInt16(length(l2g))

        if isempty(frags_k)
            partitions[k] = LocalPartition{Float32}(
                SoAFragBins{Float32}(Float32[], Float32[], UInt32[], UInt32[]),
                FragIndexBin{Float32}[],
                LocalFragment[],
                l2g,
                n_local,
                UInt16[],
            )
            continue
        end

        partitions[k] = _build_local_partition(frags_k, l2g, n_local,
                                                frag_bin_tol_ppm, rt_bin_tol)
    end

    # Compute per-partition prec_mz bounds
    partition_bounds = Vector{Tuple{Float32, Float32}}(undef, n_partitions)
    for k in 1:n_partitions
        pids = final_partition_pids[k]
        if isempty(pids)
            partition_bounds[k] = (Float32(Inf), Float32(-Inf))
        else
            pmin = Float32(Inf)
            pmax = Float32(-Inf)
            for pid in pids
                pmz = prec_mzs[pid]
                pmin = min(pmin, pmz)
                pmax = max(pmax, pmz)
            end
            partition_bounds[k] = (pmin, pmax)
        end
    end

    return LocalPartitionedFragmentIndex{Float32}(partitions, partition_bounds, n_partitions)
end

"""
Build a LocalPartition from SimpleFrags whose prec_id field already contains
local UInt16 IDs (stored as UInt32). Produces LocalFragment entries.
"""
function _build_local_partition(
    frag_ions::Vector{SimpleFrag{Float32}},
    local_to_global::Vector{UInt32},
    n_local::UInt16,
    frag_bin_tol_ppm::Float32,
    rt_bin_tol::Float32,
)
    sort!(frag_ions, by = x -> getIRT(x))

    n = length(frag_ions)
    local_fragments = Vector{LocalFragment}(undef, n)
    rt_bins = Vector{FragIndexBin{Float32}}(undef, n)
    soa = SoAFragBins{Float32}(
        Vector{Float32}(undef, n),
        Vector{Float32}(undef, n),
        Vector{UInt32}(undef, n),
        Vector{UInt32}(undef, n),
    )
    rt_bin_idx = 0
    frag_bin_idx = 0

    start_idx = 1
    start_irt = getIRT(frag_ions[1])

    for i in 1:n
        stop_irt = getIRT(frag_ions[i])
        if (stop_irt - start_irt > rt_bin_tol) && (i > start_idx)
            stop_idx = i - 1
            stop_irt_val = getIRT(frag_ions[stop_idx])
            sort!(@view(frag_ions[start_idx:stop_idx]), by = x -> getMZ(x))
            first_fb = frag_bin_idx + 1
            frag_bin_idx = _build_local_frag_bins!(local_fragments, soa,
                frag_bin_idx, frag_ions, start_idx, stop_idx, frag_bin_tol_ppm)
            rt_bin_idx += 1
            rt_bins[rt_bin_idx] = FragIndexBin{Float32}(
                start_irt, stop_irt_val, UInt32(first_fb), UInt32(frag_bin_idx))
            start_idx = i
            start_irt = getIRT(frag_ions[i])
        end
    end

    # Last RT bin
    stop_idx = n
    stop_irt_val = getIRT(frag_ions[stop_idx])
    sort!(@view(frag_ions[start_idx:stop_idx]), by = x -> getMZ(x))
    first_fb = frag_bin_idx + 1
    frag_bin_idx = _build_local_frag_bins!(local_fragments, soa,
        frag_bin_idx, frag_ions, start_idx, stop_idx, frag_bin_tol_ppm)
    rt_bin_idx += 1
    rt_bins[rt_bin_idx] = FragIndexBin{Float32}(
        start_irt, stop_irt_val, UInt32(first_fb), UInt32(frag_bin_idx))

    # Resize SoA arrays to actual count + SIMD padding on highs
    resize!(soa.lows, frag_bin_idx)
    resize!(soa.highs, frag_bin_idx + 7)
    for pad_i in (frag_bin_idx + 1):(frag_bin_idx + 7)
        soa.highs[pad_i] = Float32(Inf)  # safe sentinel for SIMD _vload8
    end
    resize!(soa.first_bins, frag_bin_idx)
    resize!(soa.last_bins, frag_bin_idx)
    rb_final = rt_bins[1:rt_bin_idx]
    skip_hints = _compute_skip_hints(soa, rb_final)

    return LocalPartition{Float32}(
        soa,
        rb_final,
        local_fragments,
        local_to_global,
        n_local,
        skip_hints,
    )
end

"""
Build fragment m/z bins producing LocalFragment entries (UInt16 local IDs).
Writes bin metadata into SoA parallel arrays. Returns updated frag_bin_idx.
"""
function _build_local_frag_bins!(
    local_fragments::Vector{LocalFragment},
    soa::SoAFragBins{Float32},
    frag_bin_idx::Int,
    frag_ions::Vector{SimpleFrag{Float32}},
    start::Int, stop::Int,
    frag_bin_tol_ppm::Float32,
)
    start_idx = start
    start_mz = getMZ(frag_ions[start])

    for i in start:stop
        stop_mz = getMZ(frag_ions[i])
        diff_mz = stop_mz - start_mz
        mean_mz = (stop_mz + start_mz) / 2
        if (diff_mz / (mean_mz / 1.0f6) > frag_bin_tol_ppm) && (i > start_idx)
            bin_stop = i - 1
            bin_stop_mz = getMZ(frag_ions[bin_stop])
            sort!(@view(frag_ions[start_idx:bin_stop]), by = x -> getPrecMZ(x))
            frag_bin_idx += 1
            soa.lows[frag_bin_idx] = start_mz
            soa.highs[frag_bin_idx] = bin_stop_mz
            soa.first_bins[frag_bin_idx] = UInt32(start_idx)
            soa.last_bins[frag_bin_idx] = UInt32(bin_stop)
            for idx in start_idx:bin_stop
                sf = frag_ions[idx]
                local_fragments[idx] = LocalFragment(
                    UInt16(getPrecID(sf)), getScore(sf))
            end
            start_idx = i
            start_mz = getMZ(frag_ions[i])
        end
    end

    # Last frag bin
    stop_mz = getMZ(frag_ions[stop])
    sort!(@view(frag_ions[start_idx:stop]), by = x -> getPrecMZ(x))
    frag_bin_idx += 1
    soa.lows[frag_bin_idx] = start_mz
    soa.highs[frag_bin_idx] = stop_mz
    soa.first_bins[frag_bin_idx] = UInt32(start_idx)
    soa.last_bins[frag_bin_idx] = UInt32(stop)
    for idx in start_idx:stop
        sf = frag_ions[idx]
        local_fragments[idx] = LocalFragment(
            UInt16(getPrecID(sf)), getScore(sf))
    end

    return frag_bin_idx
end

"""
Compute per-frag-bin skip hints for the hinted search.
hints[j] = k where frag_bins.lows[j+k] - frag_bins.lows[j] >= 5.0 Da.
"""
function _compute_skip_hints(
    frag_bins::SoAFragBins{Float32},
    rt_bins::Vector{FragIndexBin{Float32}},
)
    n_fb = length(frag_bins)
    hints = ones(UInt16, n_fb)
    lows = frag_bins.lows

    for rt_bin in rt_bins
        range = getSubBinRange(rt_bin)
        fb_start = Int(first(range))
        fb_end = Int(last(range))
        fb_start > fb_end && continue

        for j in fb_start:fb_end
            target_low = lows[j] + 5.0f0
            max_k = fb_end - j
            max_k <= 0 && continue

            if lows[fb_end] < target_low
                hints[j] = UInt16(max_k)
                continue
            end

            # Binary search for smallest k where lows[j+k] >= target_low
            lo_k = 1
            hi_k = max_k
            result_k = hi_k
            while lo_k <= hi_k
                mid_k = (lo_k + hi_k) >>> 1
                if lows[j + mid_k] >= target_low
                    result_k = mid_k
                    hi_k = mid_k - 1
                else
                    lo_k = mid_k + 1
                end
            end

            hints[j] = UInt16(clamp(result_k, 1, 65535))
        end
    end

    return hints
end
