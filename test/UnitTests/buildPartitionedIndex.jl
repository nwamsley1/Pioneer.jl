using Pioneer: SoAFragBins, LocalFragment, LocalPartition, LocalPartitionedFragmentIndex,
    FragIndexBin, MAX_LOCAL_PRECS,
    getFragBins, getRTBins, getFragments, getSkipHints,
    getPartitions, getPartition, getNPartitions, get_partition_range,
    getLow, getHigh, getSubBinRange, getPrecID, getScore,
    _build_local_partition, _compute_skip_hints, _build_local_frag_bins!,
    SimpleFrag

@testset "buildPartitionedIndex" begin

    @testset "_build_local_partition from SimpleFrags" begin
        # Build synthetic SimpleFrags that mimic what build_partitioned_index_from_lib creates
        frag_ions = SimpleFrag{Float32}[]
        n_precs = 10
        for pid in 1:n_precs
            pmz = Float32(500.0 + pid * 5.0)
            pirt = Float32(pid * 1.0)
            for rank in 1:3
                fmz = Float32(100.0 + pid * 10.0 + rank * 0.5)
                push!(frag_ions, SimpleFrag{Float32}(fmz, UInt32(pid), pmz, pirt, UInt8(0), UInt8(4 - rank + 1)))
            end
        end

        l2g = UInt32.(1:n_precs)
        partition = _build_local_partition(frag_ions, l2g, UInt16(n_precs), 2.5f0, 0.0f0, 3.0f0)

        # Structural checks
        @test !isempty(getFragBins(partition))
        @test !isempty(getRTBins(partition))
        @test length(getFragments(partition)) == length(frag_ions)
        @test partition.n_local_precs == UInt16(n_precs)
        @test length(partition.local_to_global) == n_precs

        # SoA highs has SIMD padding
        fb = getFragBins(partition)
        n_fb = length(fb)
        @test length(fb.highs) >= n_fb + 7
        for pad_i in (n_fb + 1):length(fb.highs)
            @test fb.highs[pad_i] == Float32(Inf)
        end

        # Skip hints have correct length and are >= 1
        hints = getSkipHints(partition)
        @test length(hints) == n_fb
        for h in hints
            @test h >= UInt16(1)
        end

        # RT bins have valid sub-bin ranges pointing into frag bins
        for rt_bin in getRTBins(partition)
            r = getSubBinRange(rt_bin)
            @test first(r) >= 1
            @test last(r) <= n_fb
        end

        # Fragment local IDs are valid
        for f in getFragments(partition)
            @test getPrecID(f) >= UInt16(1)
            @test getPrecID(f) <= UInt16(n_precs)
        end
    end

    @testset "_compute_skip_hints correctness" begin
        # Build simple SoA with known lows
        n = 20
        lows = Float32[100.0 + i * 0.5 for i in 0:(n-1)]
        highs = vcat(Float32[100.0 + i * 0.5 + 0.4 for i in 0:(n-1)], fill(Float32(Inf), 7))
        first_bins = UInt32.(1:n)
        last_bins = UInt32.(1:n)
        soa = SoAFragBins{Float32}(lows, highs, first_bins, last_bins)
        rt_bins = [FragIndexBin{Float32}(0.0f0, 100.0f0, UInt32(1), UInt32(n))]

        hints = _compute_skip_hints(soa, rt_bins)
        @test length(hints) == n

        # For uniform bins 0.5 Da wide: +5 Da = 10 bins
        # hint[1] should be ~10 (bins from lows[1]=100.0 to lows[11]=105.0)
        @test hints[1] == UInt16(10)

        # All hints >= 1
        for h in hints
            @test h >= UInt16(1)
        end

        # Hints near end should be <= n - j
        for j in 1:n
            @test hints[j] <= UInt16(n - j) || n - j <= 0
        end
    end

    @testset "Empty partition handling" begin
        soa = SoAFragBins{Float32}(Float32[], Float32[], UInt32[], UInt32[])
        lp = LocalPartition{Float32}(
            soa,
            FragIndexBin{Float32}[],
            LocalFragment[],
            UInt32[],
            UInt16(0),
            UInt16[],
        )
        @test isempty(getFragBins(lp))
        @test isempty(getRTBins(lp))
        @test isempty(getFragments(lp))
        @test isempty(getSkipHints(lp))
    end

    @testset "Balanced partition splitting" begin
        # Simulate the split logic from build_partitioned_index_from_lib
        # with a set of precursor IDs that exceeds MAX_LOCAL_PRECS
        function _split_balanced_test(pids::Vector{UInt32}, mzs::Vector{Float32}, max_precs::Int)
            out = Vector{UInt32}[]
            function _split!(out, pids, mzs, max_precs)
                if length(pids) <= max_precs
                    push!(out, pids)
                else
                    sort!(pids, by = pid -> mzs[pid])
                    mid = length(pids) ÷ 2
                    _split!(out, pids[1:mid], mzs, max_precs)
                    _split!(out, pids[mid+1:end], mzs, max_precs)
                end
            end
            _split!(out, pids, mzs, max_precs)
            return out
        end

        # Test: 100 precursors with max=30 → should split into ~4 balanced partitions
        n = 100
        mzs = Float32[500.0 + i * 0.1 for i in 1:n]
        pids = UInt32.(1:n)
        parts = _split_balanced_test(pids, mzs, 30)

        # All precursors preserved
        @test sort(vcat(parts...)) == UInt32.(1:n)
        # No partition exceeds limit
        for p in parts
            @test length(p) <= 30
        end
        # Balanced: all partitions within 2x of each other
        sizes = length.(parts)
        @test maximum(sizes) <= 2 * minimum(sizes)
        # Partitions are in m/z order (each partition's max < next partition's min)
        for i in 1:length(parts)-1
            max_mz_i = maximum(mzs[pid] for pid in parts[i])
            min_mz_next = minimum(mzs[pid] for pid in parts[i+1])
            @test max_mz_i <= min_mz_next
        end

        # Test: exactly at limit → no split
        parts_exact = _split_balanced_test(UInt32.(1:30), mzs, 30)
        @test length(parts_exact) == 1
        @test length(parts_exact[1]) == 30

        # Test: one over limit → splits into two halves
        parts_one_over = _split_balanced_test(UInt32.(1:31), mzs, 30)
        @test length(parts_one_over) == 2
        @test length(parts_one_over[1]) == 15
        @test length(parts_one_over[2]) == 16

        # Test: large split (200K → should produce ~4 partitions of ~50K)
        n_large = 200_000
        mzs_large = Float32[400.0 + i * 0.001 for i in 1:n_large]
        pids_large = UInt32.(1:n_large)
        parts_large = _split_balanced_test(pids_large, mzs_large, MAX_LOCAL_PRECS)
        @test sort(vcat(parts_large...)) == UInt32.(1:n_large)
        for p in parts_large
            @test length(p) <= MAX_LOCAL_PRECS
        end
        # Should be 4 partitions (200K / 65535 → need 4)
        @test length(parts_large) == 4
        # All roughly balanced (~50K each)
        for p in parts_large
            @test length(p) >= 40_000
        end
    end

    @testset "mDa fragment bins" begin
        # Build fragments with known m/z values spanning 100-200 Da
        frag_ions = SimpleFrag{Float32}[]
        n_precs = 5
        for pid in 1:n_precs
            for fmz in Float32[100.0, 120.0, 150.0, 180.0, 200.0]
                push!(frag_ions, SimpleFrag{Float32}(
                    fmz + pid * 0.0001f0,  # slight offset per precursor
                    UInt32(pid), 500.0f0, Float32(pid), UInt8(0), UInt8(1)))
            end
        end
        l2g = UInt32.(1:n_precs)

        # Build with 2 mDa bins (ppm=0 → mDa mode)
        partition_mda2 = _build_local_partition(frag_ions, l2g, UInt16(n_precs),
            0.0f0, 2.0f0, 100.0f0)
        n_bins_mda2 = length(getFragBins(partition_mda2))

        # Build with 4 mDa bins
        partition_mda4 = _build_local_partition(copy(frag_ions), l2g, UInt16(n_precs),
            0.0f0, 4.0f0, 100.0f0)
        n_bins_mda4 = length(getFragBins(partition_mda4))

        # 4 mDa should have fewer bins than 2 mDa (wider tolerance)
        @test n_bins_mda4 <= n_bins_mda2

        # Build with ppm mode for comparison
        partition_ppm = _build_local_partition(copy(frag_ions), l2g, UInt16(n_precs),
            2.5f0, 0.0f0, 100.0f0)
        n_bins_ppm = length(getFragBins(partition_ppm))

        # All should produce non-empty partitions with all fragments
        @test length(getFragments(partition_mda2)) == length(frag_ions)
        @test length(getFragments(partition_mda4)) == length(frag_ions)
        @test length(getFragments(partition_ppm)) == length(frag_ions)

        # mDa bins should have uniform width (check that the bin widths are consistent)
        fb_mda2 = getFragBins(partition_mda2)
        if length(fb_mda2) > 2
            widths = [fb_mda2.highs[i] - fb_mda2.lows[i] for i in 1:length(fb_mda2)]
            # All bin widths should be <= 2 mDa = 0.002 Da (allowing for rounding)
            for w in widths
                @test w <= 0.003f0  # 3 mDa generous bound
            end
        end

        # ppm bins should have variable width (wider at higher m/z)
        fb_ppm = getFragBins(partition_ppm)
        if length(fb_ppm) > 5
            low_mz_widths = [fb_ppm.highs[i] - fb_ppm.lows[i] for i in 1:min(3, length(fb_ppm))
                             if fb_ppm.lows[i] < 130.0f0]
            high_mz_widths = [fb_ppm.highs[i] - fb_ppm.lows[i] for i in 1:length(fb_ppm)
                              if fb_ppm.lows[i] > 170.0f0]
            if !isempty(low_mz_widths) && !isempty(high_mz_widths)
                # High m/z bins should be wider than low m/z bins in ppm mode
                @test sum(high_mz_widths)/length(high_mz_widths) > sum(low_mz_widths)/length(low_mz_widths)
            end
        end
    end

    @testset "LocalPartitionedFragmentIndex construction" begin
        # Build two small partitions
        frag_ions1 = [SimpleFrag{Float32}(100.0f0, UInt32(1), 500.0f0, 1.0f0, UInt8(0), UInt8(3)),
                      SimpleFrag{Float32}(200.0f0, UInt32(2), 502.0f0, 2.0f0, UInt8(0), UInt8(5))]
        p1 = _build_local_partition(frag_ions1, UInt32[1, 2], UInt16(2), 10.0f0, 0.0f0, 10.0f0)

        frag_ions2 = [SimpleFrag{Float32}(150.0f0, UInt32(1), 510.0f0, 1.5f0, UInt8(0), UInt8(2))]
        p2 = _build_local_partition(frag_ions2, UInt32[3], UInt16(1), 10.0f0, 0.0f0, 10.0f0)

        bounds = [(500.0f0, 502.0f0), (510.0f0, 510.0f0)]
        pfi = LocalPartitionedFragmentIndex{Float32}([p1, p2], bounds, 2)

        @test getNPartitions(pfi) == 2
        @test length(getPartitions(pfi)) == 2

        # Partition range query
        f, l = get_partition_range(pfi, 501.0f0, 505.0f0)
        @test f == 1
        @test l == 1

        f, l = get_partition_range(pfi, 500.0f0, 515.0f0)
        @test f == 1
        @test l == 2
    end
end
