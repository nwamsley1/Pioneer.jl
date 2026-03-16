#!/usr/bin/env julia
#
# Diagnose score mismatches between NativeFragmentIndex and PartitionedFragmentIndex.
# Reads serialized diagnostic data from benchmark.jl.
#
# Usage:
#   julia --project=. experiments/PartitionedFragIndex/diagnose.jl /path/to/partitioned_diag.jls
#

using Pioneer
using Serialization

include(joinpath(@__DIR__, "partitioned_types.jl"))

if isempty(ARGS)
    error("Usage: julia --project=. experiments/PartitionedFragIndex/diagnose.jl <partitioned_diag.jls>")
end

println("Loading diagnostic data from $(ARGS[1]) ...")
d = deserialize(ARGS[1])

native_frag_bins  = d.native_frag_bins
native_rt_bins    = d.native_rt_bins
native_fragments  = d.native_fragments
pfi               = d.partitioned_index
precursor_mzs     = d.precursor_mzs
mem               = d.mass_err_model
diag_scans        = d.diag_scans

println("Loaded $(length(diag_scans)) mismatched scans")
println("Native index: $(length(native_frag_bins)) frag_bins, $(length(native_rt_bins)) rt_bins, $(length(native_fragments)) fragments")
println("Partitioned index: $(pfi.n_partitions) partitions")
println()

for (scan_idx, info) in sort(collect(diag_scans), by=first)
    # Find precursors with score mismatches
    mismatched_pids = UInt32[]
    for pid in keys(info.baseline_scores)
        if haskey(info.partitioned_scores, pid) && info.baseline_scores[pid] != info.partitioned_scores[pid]
            push!(mismatched_pids, pid)
        end
    end
    isempty(mismatched_pids) && continue

    sort!(mismatched_pids)

    println("="^80)
    println("SCAN $scan_idx — $(length(mismatched_pids)) precursors with score mismatches")
    println("  quad_window=[$(info.prec_min), $(info.prec_max)]  center_mz=$(info.center_mz)  iso_width=$(info.iso_width)")
    println("  irt_window=[$(info.irt_lo), $(info.irt_hi)]  rt_bin_idx=$(info.rt_bin_idx)")
    println("  n_peaks=$(length(info.masses))  iso_bounds=$(info.iso_bounds)")
    println()

    for target_pid in mismatched_pids[1:min(3, length(mismatched_pids))]
        base_score = info.baseline_scores[target_pid]
        part_score = info.partitioned_scores[target_pid]
        target_prec_mz = precursor_mzs[target_pid]
        target_k = get_partition_idx(pfi, target_prec_mz)
        part_lo = pfi.prec_mz_min + (target_k - 1) * pfi.partition_width
        part_hi = part_lo + pfi.partition_width

        println("  --- prec_id=$target_pid  prec_mz=$target_prec_mz  baseline=$base_score  partitioned=$part_score ---")
        println("      partition=$target_k  range=[$part_lo, $part_hi]")

        # Compute corrected masses and tolerances once
        mass_ranges = Tuple{Float32, Float32}[]
        for mass in info.masses
            corrected_mz = Pioneer.getCorrectedMz(mem, mass)
            frag_min, frag_max = Pioneer.getMzBoundsReverse(mem, corrected_mz)
            push!(mass_ranges, (frag_min, frag_max))
        end

        # ── BASELINE: brute-force count of fragment hits per RT bin ──────
        println("\n      BASELINE (brute-force walk):")
        baseline_total = 0
        baseline_hits_detail = Dict{Int, Vector{NamedTuple{(:frag_idx, :frag_mz, :score, :fb_idx, :peak_idx), Tuple{UInt32, Float32, UInt8, Int, Int}}}}()

        local_rt = info.rt_bin_idx
        while local_rt <= length(native_rt_bins) && Pioneer.getLow(native_rt_bins[local_rt]) < info.irt_hi
            rt_bin = native_rt_bins[local_rt]
            sub_range = Pioneer.getSubBinRange(rt_bin)
            hits_this_rt = 0

            for (peak_idx, (frag_min, frag_max)) in enumerate(mass_ranges)
                for fb_idx in sub_range
                    fbin = native_frag_bins[fb_idx]
                    Pioneer.getLow(fbin) > frag_max && break
                    Pioneer.getHigh(fbin) < frag_min && continue

                    frag_range = Pioneer.getSubBinRange(fbin)
                    for fi in frag_range
                        frag = native_fragments[fi]
                        Pioneer.getPrecID(frag) != target_pid && continue
                        pmz = Pioneer.getPrecMZ(frag)
                        (pmz < info.prec_min || pmz > info.prec_max) && continue
                        hits_this_rt += 1
                        if !haskey(baseline_hits_detail, local_rt)
                            baseline_hits_detail[local_rt] = []
                        end
                        push!(baseline_hits_detail[local_rt], (frag_idx=fi, frag_mz=Pioneer.getLow(fbin), score=Pioneer.getScore(frag), fb_idx=fb_idx, peak_idx=peak_idx))
                    end
                end
            end

            if hits_this_rt > 0
                println("        RT bin $local_rt (iRT=[$(Pioneer.getLow(rt_bin)), $(Pioneer.getHigh(rt_bin))]): $hits_this_rt hits")
            end
            baseline_total += hits_this_rt
            local_rt += 1
        end
        println("        TOTAL: $baseline_total")

        # ── PARTITIONED: brute-force count per partition per RT bin ───────
        first_k = get_partition_idx(pfi, info.prec_min)
        last_k  = get_partition_idx(pfi, info.prec_max)
        println("\n      PARTITIONED (searching partitions $first_k:$last_k):")
        partitioned_total = 0

        for k in first_k:last_k
            partition = getPartition(pfi, k)
            p_rt_bins   = Pioneer.getRTBins(partition)
            p_frag_bins = Pioneer.getFragBins(partition)
            p_fragments = Pioneer.getFragments(partition)
            isempty(p_frag_bins) && continue

            p_local_rt = info.rt_bin_idx
            if p_local_rt > length(p_rt_bins)
                p_local_rt = length(p_rt_bins)
            end

            while p_local_rt <= length(p_rt_bins) && Pioneer.getLow(p_rt_bins[p_local_rt]) < info.irt_hi
                p_rt_bin = p_rt_bins[p_local_rt]
                p_sub_range = Pioneer.getSubBinRange(p_rt_bin)

                if first(p_sub_range) <= last(p_sub_range)
                    rt_hits = 0
                    for (peak_idx, (frag_min, frag_max)) in enumerate(mass_ranges)
                        for fb_idx in p_sub_range
                            fbin = p_frag_bins[fb_idx]
                            Pioneer.getLow(fbin) > frag_max && break
                            Pioneer.getHigh(fbin) < frag_min && continue

                            frag_range = Pioneer.getSubBinRange(fbin)
                            for fi in frag_range
                                frag = p_fragments[fi]
                                Pioneer.getPrecID(frag) != target_pid && continue
                                rt_hits += 1
                            end
                        end
                    end

                    if rt_hits > 0
                        println("        Partition $k, RT bin $p_local_rt: $rt_hits hits")
                    end
                    partitioned_total += rt_hits
                end
                p_local_rt += 1
            end
        end
        println("        TOTAL: $partitioned_total")

        # ── Compare totals ───────────────────────────────────────────────
        if baseline_total == partitioned_total
            println("\n      Fragment hit counts MATCH ($baseline_total)")
            println("      Score difference must come from the scoring algorithm (exponential/binary search), not missing fragments.")
            println("      This means the issue is in queryFragmentPartitioned! vs queryFragment! frag-bin lookup logic.")
        else
            diff = baseline_total - partitioned_total
            println("\n      *** FRAGMENT COUNT MISMATCH: baseline=$baseline_total  partitioned=$partitioned_total  diff=$diff ***")

            # Check: how many total fragments exist for this precursor in each index?
            n_native = count(f -> Pioneer.getPrecID(f) == target_pid, native_fragments)
            n_part = 0
            for pk in 1:pfi.n_partitions
                n_part += count(f -> Pioneer.getPrecID(f) == target_pid, Pioneer.getFragments(getPartition(pfi, pk)))
            end
            println("      Total fragments for prec_id=$target_pid: native=$n_native  across_all_partitions=$n_part")

            if n_native != n_part
                println("      *** BUILD BUG: fragments lost or duplicated during partitioning! ***")
            else
                println("      Fragment count OK globally — issue is in which frag bins/RT bins are searched.")
                # Show which RT bins had hits in baseline but not partitioned
                for (rt_idx, hits) in sort(collect(baseline_hits_detail))
                    println("      Baseline RT bin $rt_idx: $(length(hits)) hits")
                    for h in hits[1:min(5, length(hits))]
                        println("        frag_idx=$(h.frag_idx) fb_idx=$(h.fb_idx) peak_idx=$(h.peak_idx) score=$(h.score)")
                    end
                end
            end
        end
        println()
    end
end

println("\nDiagnosis complete.")
