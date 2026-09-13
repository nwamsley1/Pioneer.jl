using Pioneer: SoAFragBins, LocalFragment, Counter, LocalPartition,
    LocalPartitionedFragmentIndex, FragIndexBin, MassErrorModel,
    HINT_LINEAR_THRESHOLD, MAX_LOCAL_PRECS,
    getFragBins, getRTBins, getFragments, getSkipHints,
    getPartitions, getPartition, getNPartitions, get_partition_range,
    getPrecID, getScore, getLow, getHigh, getSubBinRange,
    F32x8, _vbroadcast8, _vload8, _vcmpge_mask, _find_first_ge,
    _findFirstFragBin_hybrid, searchFragmentBinUnconditional!,
    queryFragmentHinted!, _score_partition_hinted!, _find_rt_bin_start
using StaticArrays: SVector

# ─── Helpers ─────────────────────────────────────────────────────────────────

"""Extract non-zero scores from a Counter as Dict{UInt16, UInt8}."""
function extract_local_scores(lc::Counter{UInt16, UInt8})
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
function brute_force_query!(counter::Counter{UInt16, UInt8},
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

"""Independent partition oracle: retain RT/peak/bin encounter order and OR masks."""
function brute_force_partition_windows!(counter, partition, irt_low, irt_high,
        masses, intensities, model; intensity_threshold=0.0f0)
    bins = getFragBins(partition)
    fragments = getFragments(partition)
    for rt_bin in getRTBins(partition)
        getHigh(rt_bin) >= irt_low && getLow(rt_bin) < irt_high || continue
        isempty(getSubBinRange(rt_bin)) && continue
        for peak_i in eachindex(masses)
            intensity = intensities[peak_i]::Float32
            intensity < intensity_threshold && continue
            @fastmath _, low, high = Pioneer.getCorrectedMzAndBounds(
                model, masses[peak_i]::Float32, intensity)
            for bin_i in getSubBinRange(rt_bin)
                bins.highs[bin_i] >= low && bins.lows[bin_i] <= high || continue
                for frag_i in bins.first_bins[bin_i]:bins.last_bins[bin_i]
                    fragment = fragments[frag_i]
                    Pioneer.or!(counter, getPrecID(fragment), getScore(fragment))
                end
            end
        end
    end
    return nothing
end

"""Repeat the same m/z bins in three RT bins, with empty bins before/between them."""
function make_window_reuse_partition()
    centers = Float32[100, 110, 120, 130, 140]
    local_ids = ([9, 2, 12, 4, 7], [5, 11, 3, 10, 1], [15, 6, 14, 8, 13])
    rt_bounds = ((1.0f0, 2.0f0), (2.0f0, 3.0f0), (4.0f0, 5.0f0))
    lows, highs = Float32[], Float32[]
    first_bins, last_bins = UInt32[], UInt32[]
    fragments = LocalFragment[]
    rt_bins = [FragIndexBin{Float32}(0.0f0, 1.0f0, UInt32(1), UInt32(0))]
    for r in 1:3
        first_bin = UInt32(length(lows) + 1)
        for (j, mz) in enumerate(centers)
            push!(lows, mz - 0.001f0)
            push!(highs, mz + 0.001f0)
            push!(first_bins, UInt32(length(fragments) + 1))
            mask = UInt8(1) << UInt8(mod(j + r - 2, 8))
            push!(fragments, LocalFragment(UInt16(local_ids[r][j]), mask))
            push!(fragments, LocalFragment(UInt16(16), mask))
            push!(last_bins, UInt32(length(fragments)))
        end
        low_rt, high_rt = rt_bounds[r]
        push!(rt_bins, FragIndexBin{Float32}(low_rt, high_rt, first_bin, UInt32(length(lows))))
        r == 2 && push!(rt_bins, FragIndexBin{Float32}(3.0f0, 4.0f0, UInt32(11), UInt32(10)))
    end
    append!(highs, fill(Inf32, 7))
    bins = SoAFragBins{Float32}(lows, highs, first_bins, last_bins)
    return LocalPartition{Float32}(bins, rt_bins, fragments, UInt32.(101:116),
        UInt16(16), make_soa_hints(bins, rt_bins))
end

function window_test_spline(base::Float32, slope::Float32, low::Float32, high::Float32)
    width = high - low
    change = slope * width
    coefficients = SVector{8, Float32}(base, change, 0, 0, base + change, change, 0, 0)
    return Pioneer.UniformSpline{8, Float32}(coefficients, 3, low, high, width)
end

"""Real intensity/scout models with varying bias, spread, and extrapolated tails."""
function window_reuse_models(; bias_shift=0.0f0)
    mz_bias = window_test_spline(0.002f0 + bias_shift, 0.00002f0, 90.0f0, 150.0f0)
    intensity_bias = window_test_spline(0.0f0, 0.0002f0, 1.0f0, 10.0f0)
    spread = window_test_spline(1.5f0, 0.04f0, 1.0f0, 10.0f0)
    rt_bias = window_test_spline(0.0f0, 0.0f0, 0.0f0, 5.0f0)
    mz_spread = window_test_spline(1.0f0, 0.001f0, 90.0f0, 150.0f0)
    extrap(s) = Pioneer.make_spline_extrap(s, s.first, s.last)
    intensity = Pioneer.IntensityMassErrorModel(
        mz_bias, intensity_bias, spread, rt_bias,
        extrap(mz_bias), extrap(intensity_bias), extrap(spread), extrap(rt_bias),
        2.0f0, 0.03f0, mz_spread, extrap(mz_spread), 1.0f0, 0.0f0, 0.0f0,
        -10.0f0, 0.03f0, 0.0f0, 5.0f0)
    scout = Pioneer.ScoutCalibratedMassErrorModel(mz_bias, extrap(mz_bias),
        intensity_bias, extrap(intensity_bias), true, 0.006f0)
    scout_no_intensity = Pioneer.ScoutCalibratedMassErrorModel(mz_bias, extrap(mz_bias),
        intensity_bias, extrap(intensity_bias), false, 0.006f0)
    return (intensity, scout, scout_no_intensity)
end

struct PartitionWindowCountingModel{M} <: Pioneer.AbstractMassErrorModel
    model::M
    calls::Base.RefValue{Int}
end

function Pioneer.getCorrectedMzAndBounds(model::PartitionWindowCountingModel,
        mz::Float32, intensity::Float32)
    model.calls[] += 1
    return Pioneer.getCorrectedMzAndBounds(model.model, mz, intensity)
end

function partition_counter_state(counter)
    return (ids=copy(counter.ids[1:(counter.size - 1)]), counts=copy(counter.counts))
end

function accepted_partition_windows(masses, intensities, model, threshold)
    lows, highs = Float32[], Float32[]
    for i in eachindex(masses)
        intensity = intensities[i]::Float32
        intensity < threshold && continue
        @fastmath _, low, high = Pioneer.getCorrectedMzAndBounds(model, masses[i]::Float32, intensity)
        push!(lows, low)
        push!(highs, high)
    end
    return lows, highs
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

        # Counter
        lc = Counter(UInt16, UInt8, 100)
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
        lc = Counter(UInt16, UInt8, 5)
        searchFragmentBinUnconditional!(lc, frags, UInt32(1):UInt32(3))
        scores = extract_local_scores(lc)
        # Per-precursor score is bitmask-OR'd across observed fragments
        # (commit f472cd06 replaced additive accumulation with bitmask
        # scoring throughout the index). pid=1 sees ranks 3 & 2 →
        # 3 | 2 = 3; pid=2 sees rank 5 alone → 5.
        @test scores[UInt16(1)] == UInt8(3)  # 3 | 2 = 0b011
        @test scores[UInt16(2)] == UInt8(5)  # 5 alone = 0b101
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
            lc_ref = Counter(UInt16, UInt8, 21)
            brute_force_query!(lc_ref, soa, frags, frag_mz_min, frag_mz_max, max_idx)
            ref_scores = extract_local_scores(lc_ref)

            # Hinted search
            lc_hnt = Counter(UInt16, UInt8, 21)
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

        lc_ref = Counter(UInt16, UInt8, 21)
        brute_force_query!(lc_ref, soa, frags, frag_mz_min, frag_mz_max, max_idx)
        ref = extract_local_scores(lc_ref)

        lc_hnt = Counter(UInt16, UInt8, 21)
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

        lc = Counter(UInt16, UInt8, 6)
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
            lc = Counter(UInt16, UInt8, 16)
            brute_force_query!(lc, soa, frags, peak - tol, peak + tol, max_idx)
            ref_all[peak] = extract_local_scores(lc)
        end

        for threshold in thresholds
            lb, ub = UInt32(1), UInt32(1)
            prev_mz = 0.0f0
            for peak in peaks
                lc = Counter(UInt16, UInt8, 16)
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
        # `_score_partition_hinted!` now takes intensities (added between
        # `masses` and `mass_err_model`); use uniform 1.0 so the default
        # `intensity_threshold = 0.0` admits everything → behavior matches
        # the pre-intensity-arg signature for these tests.
        intensities = Union{Missing, Float32}[ones(Float32, length(masses))...]
        mem = MassErrorModel(0.0f0, (10.0f0, 10.0f0))

        # Hinted with default threshold
        lc_hnt = Counter(UInt16, UInt8, 31)
        _score_partition_hinted!(lc_hnt, partition, 0.0f0, 10.0f0, masses, intensities, mem)
        hnt = extract_local_scores(lc_hnt)
        @test !isempty(hnt)

        # Different thresholds should produce identical results
        for threshold in [UInt32(1), UInt32(32), UInt32(128), UInt32(1000000)]
            lc_t = Counter(UInt16, UInt8, 31)
            _score_partition_hinted!(lc_t, partition, 0.0f0, 10.0f0, masses, intensities, mem;
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
        intensities = Union{Missing, Float32}[ones(Float32, length(masses))...]
        mem = MassErrorModel(0.0f0, (20.0f0, 20.0f0))

        lc_ref = Counter(UInt16, UInt8, 21)
        _score_partition_hinted!(lc_ref, partition, 0.0f0, 5.0f0, masses, intensities, mem)
        ref = extract_local_scores(lc_ref)
        @test !isempty(ref)

        for threshold in [UInt32(1), UInt32(64), UInt32(1000000)]
            lc = Counter(UInt16, UInt8, 21)
            _score_partition_hinted!(lc, partition, 0.0f0, 5.0f0, masses, intensities, mem;
                                      linear_threshold=threshold)
            @test extract_local_scores(lc) == ref
        end
    end

    @testset "Mass windows reused across RT bins" begin
        partition = make_window_reuse_partition()
        masses = Union{Missing, Float32}[100.001f0, 110.006f0, 120.0f0, 130.012f0, 140.0f0]
        intensities = Union{Missing, Float32}[1.0f0, 16.0f0, 8.0f0, 16.0f0, 4096.0f0]
        models = (MassErrorModel(2.0f0, (10.0f0, 20.0f0)), window_reuse_models()...)
        for model in models, threshold in (0.0f0, 16.0f0),
                rt_range in ((0.0f0, 5.0f0), (2.0f0, 4.5f0), (1.0f0, 2.0f0))
            reference = Counter(UInt16, UInt8, 17)
            brute_force_partition_windows!(reference, partition, rt_range..., masses, intensities, model;
                intensity_threshold=threshold)
            # Cover both the compatible no-scratch call and caller-owned scratch.
            for with_scratch in (false, true)
                counter = Counter(UInt16, UInt8, 17)
                low_buf, high_buf = Float32[], Float32[]
                if with_scratch
                    _score_partition_hinted!(counter, partition, rt_range..., masses, intensities, model;
                        intensity_threshold=threshold, mz_low_buf=low_buf, mz_high_buf=high_buf)
                    low, high = accepted_partition_windows(masses, intensities, model, threshold)
                    @test low_buf[1:length(low)] == low
                    @test high_buf[1:length(high)] == high
                else
                    _score_partition_hinted!(counter, partition, rt_range..., masses, intensities, model;
                        intensity_threshold=threshold)
                end
                @test partition_counter_state(counter) == partition_counter_state(reference)
            end
        end

        @testset "Intensity cutoff ties preserve precursor encounter order" begin
            exact_masses = Union{Missing, Float32}[100, 110, 120, 130, 140]
            tied_intensities = Union{Missing, Float32}[5, 10, 10, 2, 10]
            counter = Counter(UInt16, UInt8, 17)
            _score_partition_hinted!(counter, partition, 0.0f0, 5.0f0,
                exact_masses, tied_intensities, MassErrorModel(0.0f0, (0.0f0, 0.0f0));
                intensity_threshold=10.0f0, mz_low_buf=Float32[], mz_high_buf=Float32[])
            @test counter.ids[1:(counter.size - 1)] == UInt16[2, 16, 12, 7, 11, 3, 1, 6, 14, 13]
            @test counter.counts[16] == UInt8(0x7e)

            # A zero-width observed window on either fragment-bin boundary matches.
            for edges in (getFragBins(partition).lows, getFragBins(partition).highs)
                boundary_masses = Union{Missing, Float32}[edges[1:5]...]
                reference = Counter(UInt16, UInt8, 17)
                Pioneer.reset!(counter)
                model = MassErrorModel(0.0f0, (0.0f0, 0.0f0))
                brute_force_partition_windows!(reference, partition, 0.0f0, 5.0f0,
                    boundary_masses, tied_intensities, model; intensity_threshold=10.0f0)
                _score_partition_hinted!(counter, partition, 0.0f0, 5.0f0,
                    boundary_masses, tied_intensities, model; intensity_threshold=10.0f0,
                    mz_low_buf=Float32[], mz_high_buf=Float32[])
                @test partition_counter_state(counter) == partition_counter_state(reference)
                @test counter.size > 1
            end
        end

        @testset "Correct once per accepted peak and replace previous scan/model windows" begin
            low_buf, high_buf = fill(-1.0f0, 12), fill(-2.0f0, 12)
            counter = Counter(UInt16, UInt8, 17)
            changed_models = window_reuse_models(bias_shift=0.015f0)
            scans = ((masses, intensities, models[2], 0.0f0),
                (masses, intensities, changed_models[1], 16.0f0),
                (masses[2:3], intensities[2:3], changed_models[2], 0.0f0),
                (masses, reverse(intensities), models[3], 0.0f0))
            previous_low = Float32[]
            for (scan_masses, scan_intensities, model, threshold) in scans
                counted = PartitionWindowCountingModel(model, Ref(0))
                Pioneer.reset!(counter)
                reference = Counter(UInt16, UInt8, 17)
                brute_force_partition_windows!(reference, partition, 0.0f0, 5.0f0,
                    scan_masses, scan_intensities, model; intensity_threshold=threshold)
                _score_partition_hinted!(counter, partition, 0.0f0, 5.0f0,
                    scan_masses, scan_intensities, counted; intensity_threshold=threshold,
                    mz_low_buf=low_buf, mz_high_buf=high_buf)
                low, high = accepted_partition_windows(scan_masses, scan_intensities, model, threshold)
                @test counted.calls[] == length(low)
                @test low_buf[1:length(low)] == low
                @test high_buf[1:length(high)] == high
                @test partition_counter_state(counter) == partition_counter_state(reference)
                @test low != previous_low
                previous_low = low
            end
        end

        @testset "Lazy windows for empty, disjoint and filtered scans" begin
            empty_partition = LocalPartition{Float32}(
                SoAFragBins{Float32}(Float32[], Float32[], UInt32[], UInt32[]),
                FragIndexBin{Float32}[], LocalFragment[], UInt32[], UInt16(0), UInt16[])
            no_rt_partition = LocalPartition{Float32}(getFragBins(partition),
                FragIndexBin{Float32}[], getFragments(partition), UInt32.(101:116),
                UInt16(16), getSkipHints(partition))
            unreadable = Union{Missing, Float32}[missing]
            empty_scan = Union{Missing, Float32}[]
            cases = ((empty_partition, 0.0f0, 5.0f0, unreadable, unreadable, 0.0f0),
                (no_rt_partition, 0.0f0, 5.0f0, unreadable, unreadable, 0.0f0),
                (partition, 6.0f0, 7.0f0, unreadable, unreadable, 0.0f0),
                (partition, 0.0f0, 0.5f0, unreadable, unreadable, 0.0f0),
                (partition, 0.0f0, 5.0f0, empty_scan, empty_scan, 0.0f0),
                (partition, 0.0f0, 5.0f0, unreadable, Union{Missing, Float32}[1], 2.0f0))
            for (part, rt_low, rt_high, scan_masses, scan_intensities, threshold) in cases
                counter = Counter(UInt16, UInt8, 17)
                low_buf, high_buf = fill(-1.0f0, 8), fill(-2.0f0, 8)
                counted = PartitionWindowCountingModel(models[2], Ref(0))
                _score_partition_hinted!(counter, part, rt_low, rt_high,
                    scan_masses, scan_intensities, counted; intensity_threshold=threshold,
                    mz_low_buf=low_buf, mz_high_buf=high_buf)
                @test counted.calls[] == 0
                @test counter.size == 1
                @test all(iszero, counter.counts)
                @test low_buf == fill(-1.0f0, 8)
                @test high_buf == fill(-2.0f0, 8)
            end
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

            lc_ref = Counter(UInt16, UInt8, 51)
            brute_force_query!(lc_ref, soa, frags, frag_mz_min, frag_mz_max, max_idx)
            ref = extract_local_scores(lc_ref)

            lc_hnt = Counter(UInt16, UInt8, 51)
            lb, ub = queryFragmentHinted!(lc_hnt, max_idx, lb, ub,
                soa, frags, frag_mz_min, frag_mz_max, hints, prev_mz, HINT_LINEAR_THRESHOLD)
            hnt = extract_local_scores(lc_hnt)

            @test ref == hnt
            prev_mz = frag_mz_min
        end
    end
end
