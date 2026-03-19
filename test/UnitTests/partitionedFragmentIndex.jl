using Pioneer: SoAFragBins, LocalFragment, LocalCounter, LocalPartition,
    LocalPartitionedFragmentIndex, FragIndexBin, MassErrorModel,
    HINT_LINEAR_THRESHOLD, MAX_LOCAL_PRECS,
    getFragBins, getRTBins, getFragments, getSkipHints,
    getPartitions, getPartition, getNPartitions, get_partition_range,
    getPrecID, getScore, getLow, getHigh, getSubBinRange,
    F32x8, _vbroadcast8, _vload8, _vcmpge_mask, _find_first_ge,
    _findFirstFragBin_hybrid, searchFragmentBinUnconditional!,
    queryFragmentHinted!, _score_partition_hinted!, _find_rt_bin_start

# ─── Helpers ─────────────────────────────────────────────────────────────────

"""Extract non-zero scores from a LocalCounter as Dict{UInt16, UInt8}."""
function extract_local_scores(lc::LocalCounter{UInt16, UInt8})
    scores = Dict{UInt16, UInt8}()
    @inbounds for i in 1:(lc.size - 1)
        lid = lc.ids[i]
        sc = lc.counts[lid]
        sc > 0 && (scores[lid] = sc)
    end
    return scores
end

"""Build a SoAFragBins + LocalFragment[] with uniform bins."""
function make_soa_frag_bins(;
        n_bins::Int, start_mz::Float32, bin_width::Float32,
        frags_per_bin::Int, n_precs::Int)
    lows = Vector{Float32}(undef, n_bins)
    highs = Vector{Float32}(undef, n_bins + 7)  # SIMD padding
    first_bins = Vector{UInt32}(undef, n_bins)
    last_bins = Vector{UInt32}(undef, n_bins)
    fragments = LocalFragment[]
    frag_idx = UInt32(0)

    for i in 1:n_bins
        lo = start_mz + (i - 1) * bin_width
        hi = lo + bin_width
        lows[i] = lo
        highs[i] = hi
        first_f = frag_idx + one(UInt32)
        for j in 1:frags_per_bin
            frag_idx += one(UInt32)
            lid = UInt16(((i - 1) * frags_per_bin + j - 1) % n_precs + 1)
            push!(fragments, LocalFragment(lid, UInt8(1)))
        end
        first_bins[i] = first_f
        last_bins[i] = frag_idx
    end

    # Pad highs with Inf sentinels
    for pad_i in (n_bins + 1):(n_bins + 7)
        highs[pad_i] = Float32(Inf)
    end

    soa = SoAFragBins{Float32}(lows, highs, first_bins, last_bins)
    return soa, fragments
end

"""Build skip hints for SoA frag bins."""
function make_soa_hints(soa::SoAFragBins{Float32}, rt_bins::Vector{FragIndexBin{Float32}})
    n = length(soa)
    hints = ones(UInt16, n)
    lows = soa.lows

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

"""Build a synthetic LocalPartition with uniform SoA bins."""
function make_test_partition(;
        n_frag_bins::Int, start_mz::Float32, bin_width::Float32,
        frags_per_bin::Int, n_precs::Int,
        n_rt_bins::Int=1, rt_lo::Float32=0.0f0, rt_hi::Float32=10.0f0)
    soa, fragments = make_soa_frag_bins(;
        n_bins=n_frag_bins, start_mz=start_mz, bin_width=bin_width,
        frags_per_bin=frags_per_bin, n_precs=n_precs)

    rt_bins = FragIndexBin{Float32}[]
    bins_per_rt = max(1, n_frag_bins ÷ n_rt_bins)
    for r in 1:n_rt_bins
        rlo = rt_lo + (r - 1) * (rt_hi - rt_lo) / n_rt_bins
        rhi = rt_lo + r * (rt_hi - rt_lo) / n_rt_bins
        fb_first = UInt32((r - 1) * bins_per_rt + 1)
        fb_last = r == n_rt_bins ? UInt32(n_frag_bins) : UInt32(r * bins_per_rt)
        push!(rt_bins, FragIndexBin{Float32}(rlo, rhi, fb_first, fb_last))
    end

    hints = make_soa_hints(soa, rt_bins)
    l2g = UInt32.(1:n_precs)

    return LocalPartition{Float32}(soa, rt_bins, fragments, l2g, UInt16(n_precs), hints)
end

"""Brute-force reference: score all matching bins by linear scan."""
function brute_force_query!(counter::LocalCounter{UInt16, UInt8},
        soa::SoAFragBins{Float32}, fragments::Vector{LocalFragment},
        frag_mz_min::Float32, frag_mz_max::Float32, max_idx::UInt32)
    for j in UInt32(1):max_idx
        if soa.highs[j] >= frag_mz_min && soa.lows[j] <= frag_mz_max
            for fi in soa.first_bins[j]:soa.last_bins[j]
                frag = fragments[fi]
                Pioneer.inc!(counter, getPrecID(frag), getScore(frag))
            end
        end
    end
end

# ─── Tests ───────────────────────────────────────────────────────────────────

@testset "PartitionedFragmentIndex" begin

    @testset "Type construction" begin
        # LocalFragment
        lf = LocalFragment(UInt16(42), UInt8(7))
        @test getPrecID(lf) == UInt16(42)
        @test getScore(lf) == UInt8(7)

        # SoAFragBins
        soa = SoAFragBins{Float32}([1.0f0, 2.0f0], [1.5f0, 2.5f0], [UInt32(1), UInt32(2)], [UInt32(1), UInt32(3)])
        @test length(soa) == 2
        @test !isempty(soa)
        @test isempty(SoAFragBins{Float32}(Float32[], Float32[], UInt32[], UInt32[]))

        # LocalCounter
        lc = LocalCounter(UInt16, UInt8, 100)
        @test lc.size == 1
        Pioneer.inc!(lc, UInt16(5), UInt8(3))
        @test lc.counts[5] == UInt8(3)
        @test lc.size == 2
        Pioneer.inc!(lc, UInt16(5), UInt8(2))
        @test lc.counts[5] == UInt8(5)
        @test lc.size == 2  # no new encounter
        Pioneer.reset!(lc)
        @test lc.size == 1
        @test lc.counts[5] == UInt8(0)

        # LocalPartition
        soa2, frags = make_soa_frag_bins(n_bins=5, start_mz=100.0f0, bin_width=1.0f0,
                                          frags_per_bin=2, n_precs=3)
        rt = [FragIndexBin{Float32}(0.0f0, 10.0f0, UInt32(1), UInt32(5))]
        hints = ones(UInt16, 5)
        lp = LocalPartition{Float32}(soa2, rt, frags, UInt32[1,2,3], UInt16(3), hints)
        @test length(getFragBins(lp)) == 5
        @test length(getRTBins(lp)) == 1
        @test length(getFragments(lp)) == 10
        @test length(getSkipHints(lp)) == 5
    end

    @testset "SIMD _find_first_ge" begin
        # Basic: find first element >= threshold
        arr = Float32[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                      Inf, Inf, Inf, Inf, Inf, Inf, Inf]  # 7 padding
        @test _find_first_ge(arr, UInt32(1), UInt32(10), 5.0f0) == UInt32(5)
        @test _find_first_ge(arr, UInt32(1), UInt32(10), 1.0f0) == UInt32(1)
        @test _find_first_ge(arr, UInt32(1), UInt32(10), 10.0f0) == UInt32(10)
        @test _find_first_ge(arr, UInt32(1), UInt32(10), 11.0f0) == UInt32(11)  # not found

        # Single element
        @test _find_first_ge(arr, UInt32(3), UInt32(3), 3.0f0) == UInt32(3)
        @test _find_first_ge(arr, UInt32(3), UInt32(3), 4.0f0) == UInt32(4)  # not found

        # Larger array (exercises SIMD path)
        big = vcat(Float32.(1:100), fill(Float32(Inf), 7))
        @test _find_first_ge(big, UInt32(1), UInt32(100), 50.0f0) == UInt32(50)
        @test _find_first_ge(big, UInt32(1), UInt32(100), 99.5f0) == UInt32(100)
    end

    @testset "_findFirstFragBin_hybrid" begin
        arr = vcat(Float32.(1:200), fill(Float32(Inf), 7))

        # Should find same result as _find_first_ge
        for threshold in [UInt32(1), UInt32(8), UInt32(32), UInt32(128)]
            @test _findFirstFragBin_hybrid(arr, UInt32(1), UInt32(200), 100.0f0, threshold) ==
                  _find_first_ge(arr, UInt32(1), UInt32(200), 100.0f0)
            @test _findFirstFragBin_hybrid(arr, UInt32(1), UInt32(200), 1.0f0, threshold) ==
                  _find_first_ge(arr, UInt32(1), UInt32(200), 1.0f0)
            @test _findFirstFragBin_hybrid(arr, UInt32(50), UInt32(150), 120.0f0, threshold) ==
                  _find_first_ge(arr, UInt32(50), UInt32(150), 120.0f0)
        end
    end

    @testset "searchFragmentBinUnconditional!" begin
        frags = [LocalFragment(UInt16(1), UInt8(3)),
                 LocalFragment(UInt16(2), UInt8(5)),
                 LocalFragment(UInt16(1), UInt8(2))]
        lc = LocalCounter(UInt16, UInt8, 5)
        searchFragmentBinUnconditional!(lc, frags, UInt32(1):UInt32(3))
        scores = extract_local_scores(lc)
        @test scores[UInt16(1)] == UInt8(5)  # 3 + 2
        @test scores[UInt16(2)] == UInt8(5)
    end

    @testset "queryFragmentHinted! vs brute force" begin
        soa, frags = make_soa_frag_bins(n_bins=100, start_mz=100.0f0,
            bin_width=0.5f0, frags_per_bin=5, n_precs=20)
        rt = [FragIndexBin{Float32}(0.0f0, 10.0f0, UInt32(1), UInt32(100))]
        hints = make_soa_hints(soa, rt)
        max_idx = UInt32(100)

        peaks = Float32[105.0, 112.0, 125.0, 140.0, 148.0]
        tol = 0.3f0
        prev_mz = 0.0f0
        lb, ub = UInt32(1), UInt32(1)

        for peak in peaks
            frag_mz_min = peak - tol
            frag_mz_max = peak + tol

            # Brute force reference
            lc_ref = LocalCounter(UInt16, UInt8, 21)
            brute_force_query!(lc_ref, soa, frags, frag_mz_min, frag_mz_max, max_idx)
            ref_scores = extract_local_scores(lc_ref)

            # Hinted search
            lc_hnt = LocalCounter(UInt16, UInt8, 21)
            lb, ub = queryFragmentHinted!(lc_hnt, max_idx, lb, ub,
                soa, frags, frag_mz_min, frag_mz_max, hints, prev_mz, HINT_LINEAR_THRESHOLD)
            hnt_scores = extract_local_scores(lc_hnt)

            @test ref_scores == hnt_scores
            prev_mz = frag_mz_min
        end
    end

    @testset "queryFragmentHinted! first peak (prev_mz=0)" begin
        soa, frags = make_soa_frag_bins(n_bins=100, start_mz=100.0f0,
            bin_width=0.5f0, frags_per_bin=5, n_precs=20)
        rt = [FragIndexBin{Float32}(0.0f0, 10.0f0, UInt32(1), UInt32(100))]
        hints = make_soa_hints(soa, rt)
        max_idx = UInt32(100)
        frag_mz_min, frag_mz_max = 125.0f0, 125.5f0

        lc_ref = LocalCounter(UInt16, UInt8, 21)
        brute_force_query!(lc_ref, soa, frags, frag_mz_min, frag_mz_max, max_idx)
        ref = extract_local_scores(lc_ref)

        lc_hnt = LocalCounter(UInt16, UInt8, 21)
        queryFragmentHinted!(lc_hnt, max_idx, UInt32(1), UInt32(1),
            soa, frags, frag_mz_min, frag_mz_max, hints, 0.0f0, HINT_LINEAR_THRESHOLD)
        hnt = extract_local_scores(lc_hnt)

        @test ref == hnt
    end

    @testset "queryFragmentHinted! no match" begin
        soa, frags = make_soa_frag_bins(n_bins=50, start_mz=100.0f0,
            bin_width=1.0f0, frags_per_bin=2, n_precs=5)
        rt = [FragIndexBin{Float32}(0.0f0, 10.0f0, UInt32(1), UInt32(50))]
        hints = make_soa_hints(soa, rt)
        max_idx = UInt32(50)

        lc = LocalCounter(UInt16, UInt8, 6)
        queryFragmentHinted!(lc, max_idx, UInt32(1), UInt32(1),
            soa, frags, 500.0f0, 501.0f0, hints, 90.0f0, HINT_LINEAR_THRESHOLD)
        @test isempty(extract_local_scores(lc))
    end

    @testset "queryFragmentHinted! threshold sweep" begin
        soa, frags = make_soa_frag_bins(n_bins=200, start_mz=100.0f0,
            bin_width=0.5f0, frags_per_bin=3, n_precs=15)
        rt = [FragIndexBin{Float32}(0.0f0, 10.0f0, UInt32(1), UInt32(200))]
        hints = make_soa_hints(soa, rt)
        max_idx = UInt32(200)

        peaks = Float32[110.0, 115.0, 130.0, 160.0, 190.0]
        tol = 0.3f0
        thresholds = [UInt32(1), UInt32(8), UInt32(32), UInt32(128), UInt32(1000000)]

        # Get reference with brute force
        ref_all = Dict{Float32, Dict{UInt16, UInt8}}()
        for peak in peaks
            lc = LocalCounter(UInt16, UInt8, 16)
            brute_force_query!(lc, soa, frags, peak - tol, peak + tol, max_idx)
            ref_all[peak] = extract_local_scores(lc)
        end

        for threshold in thresholds
            lb, ub = UInt32(1), UInt32(1)
            prev_mz = 0.0f0
            for peak in peaks
                lc = LocalCounter(UInt16, UInt8, 16)
                lb, ub = queryFragmentHinted!(lc, max_idx, lb, ub,
                    soa, frags, peak - tol, peak + tol,
                    hints, prev_mz, threshold)
                @test extract_local_scores(lc) == ref_all[peak]
                prev_mz = peak - tol
            end
        end
    end

    @testset "_score_partition_hinted! correctness" begin
        partition = make_test_partition(
            n_frag_bins=500, start_mz=100.0f0, bin_width=0.5f0,
            frags_per_bin=4, n_precs=30,
            n_rt_bins=3, rt_lo=0.0f0, rt_hi=15.0f0)

        masses = Union{Missing, Float32}[
            110.0f0, 125.0f0, 150.0f0, 175.0f0, 200.0f0, 225.0f0, 250.0f0, 300.0f0]
        mem = MassErrorModel(0.0f0, (10.0f0, 10.0f0))

        # Hinted with default threshold
        lc_hnt = LocalCounter(UInt16, UInt8, 31)
        _score_partition_hinted!(lc_hnt, partition, 0.0f0, 10.0f0, masses, mem)
        hnt = extract_local_scores(lc_hnt)
        @test !isempty(hnt)

        # Different thresholds should produce identical results
        for threshold in [UInt32(1), UInt32(32), UInt32(128), UInt32(1000000)]
            lc_t = LocalCounter(UInt16, UInt8, 31)
            _score_partition_hinted!(lc_t, partition, 0.0f0, 10.0f0, masses, mem;
                                      linear_threshold=threshold)
            t_scores = extract_local_scores(lc_t)
            @test t_scores == hnt
        end
    end

    @testset "_score_partition_hinted! PPM tolerance" begin
        partition = make_test_partition(
            n_frag_bins=300, start_mz=200.0f0, bin_width=0.5f0,
            frags_per_bin=3, n_precs=20,
            n_rt_bins=2, rt_lo=0.0f0, rt_hi=10.0f0)

        masses = Union{Missing, Float32}[210.0f0, 230.0f0, 260.0f0, 300.0f0, 340.0f0]
        mem = MassErrorModel(0.0f0, (20.0f0, 20.0f0))

        lc_ref = LocalCounter(UInt16, UInt8, 21)
        _score_partition_hinted!(lc_ref, partition, 0.0f0, 5.0f0, masses, mem)
        ref = extract_local_scores(lc_ref)
        @test !isempty(ref)

        for threshold in [UInt32(1), UInt32(64), UInt32(1000000)]
            lc = LocalCounter(UInt16, UInt8, 21)
            _score_partition_hinted!(lc, partition, 0.0f0, 5.0f0, masses, mem;
                                      linear_threshold=threshold)
            @test extract_local_scores(lc) == ref
        end
    end

    @testset "_find_rt_bin_start" begin
        rt_bins = [
            FragIndexBin{Float32}(0.0f0, 5.0f0, UInt32(1), UInt32(10)),
            FragIndexBin{Float32}(5.0f0, 10.0f0, UInt32(11), UInt32(20)),
            FragIndexBin{Float32}(10.0f0, 15.0f0, UInt32(21), UInt32(30)),
        ]
        @test _find_rt_bin_start(rt_bins, 0.0f0) == 1
        @test _find_rt_bin_start(rt_bins, 5.0f0) == 1  # first bin high=5.0 >= 5.0
        @test _find_rt_bin_start(rt_bins, 7.0f0) == 2
        @test _find_rt_bin_start(rt_bins, 12.0f0) == 3
        @test _find_rt_bin_start(rt_bins, 16.0f0) == 4  # past end
    end

    @testset "get_partition_range" begin
        bounds = [
            (100.0f0, 105.0f0),
            (105.0f0, 110.0f0),
            (110.0f0, 115.0f0),
        ]
        parts = [LocalPartition{Float32}(
            SoAFragBins{Float32}(Float32[], Float32[], UInt32[], UInt32[]),
            FragIndexBin{Float32}[], LocalFragment[], UInt32[], UInt16(0), UInt16[])
            for _ in 1:3]
        pfi = LocalPartitionedFragmentIndex{Float32}(parts, bounds, 3)

        # Query that spans partitions 1-2
        f, l = get_partition_range(pfi, 103.0f0, 108.0f0)
        @test f == 1
        @test l == 2

        # Query that spans all
        f, l = get_partition_range(pfi, 100.0f0, 115.0f0)
        @test f == 1
        @test l == 3

        # Query in single partition
        f, l = get_partition_range(pfi, 111.0f0, 114.0f0)
        @test f == 3
        @test l == 3

        # Query outside range
        f, l = get_partition_range(pfi, 200.0f0, 300.0f0)
        @test f > l  # empty range
    end

    @testset "Stress: many bins, many peaks" begin
        soa, frags = make_soa_frag_bins(n_bins=2000, start_mz=100.0f0,
            bin_width=0.25f0, frags_per_bin=5, n_precs=50)
        rt = [FragIndexBin{Float32}(0.0f0, 10.0f0, UInt32(1), UInt32(2000))]
        hints = make_soa_hints(soa, rt)
        max_idx = UInt32(2000)

        peaks = Float32[100.0 + i * 25.0 for i in 0:19]
        tol = 0.2f0

        lb, ub = UInt32(1), UInt32(1)
        prev_mz = 0.0f0

        for peak in peaks
            frag_mz_min = peak - tol
            frag_mz_max = peak + tol

            lc_ref = LocalCounter(UInt16, UInt8, 51)
            brute_force_query!(lc_ref, soa, frags, frag_mz_min, frag_mz_max, max_idx)
            ref = extract_local_scores(lc_ref)

            lc_hnt = LocalCounter(UInt16, UInt8, 51)
            lb, ub = queryFragmentHinted!(lc_hnt, max_idx, lb, ub,
                soa, frags, frag_mz_min, frag_mz_max, hints, prev_mz, HINT_LINEAR_THRESHOLD)
            hnt = extract_local_scores(lc_hnt)

            @test ref == hnt
            prev_mz = frag_mz_min
        end
    end
end
