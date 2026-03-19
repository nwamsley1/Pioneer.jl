using Pioneer: SoAFragBins, LocalFragment, LocalPartition, LocalPartitionedFragmentIndex,
    FragIndexBin, MAX_LOCAL_PRECS,
    getFragBins, getRTBins, getFragments, getSkipHints,
    getPartitions, getPartition, getNPartitions, get_partition_range,
    getLow, getHigh, getSubBinRange,
    _build_local_partition, _compute_skip_hints,
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
        partition = _build_local_partition(frag_ions, l2g, UInt16(n_precs), 2.5f0, 3.0f0)

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
            @test Pioneer.getPrecID(f) >= UInt16(1)
            @test Pioneer.getPrecID(f) <= UInt16(n_precs)
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

    @testset "LocalPartitionedFragmentIndex construction" begin
        # Build two small partitions
        frag_ions1 = [SimpleFrag{Float32}(100.0f0, UInt32(1), 500.0f0, 1.0f0, UInt8(0), UInt8(3)),
                      SimpleFrag{Float32}(200.0f0, UInt32(2), 502.0f0, 2.0f0, UInt8(0), UInt8(5))]
        p1 = _build_local_partition(frag_ions1, UInt32[1, 2], UInt16(2), 10.0f0, 10.0f0)

        frag_ions2 = [SimpleFrag{Float32}(150.0f0, UInt32(1), 510.0f0, 1.5f0, UInt8(0), UInt8(2))]
        p2 = _build_local_partition(frag_ions2, UInt32[3], UInt16(1), 10.0f0, 10.0f0)

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
