# Synthetic tests for queryFragmentHinted! correctness.
#
# Validates that the hint-informed exponential search produces identical results
# to queryFragmentPartitioned! (the proven exponential+binary search with step_size=4)
# and to a reference step_size=1 search.
#
# Tests cover:
#   - queryFragmentHinted! vs queryFragmentPartitioned! on identical frag bin arrays
#   - Edge cases: first peak (prev_mz=0), zero/negative delta, single bin, boundary bins
#   - Full _score_partition_hinted! vs _score_partition! on synthetic partitions
#   - Various linear_threshold values
#
# Run:  julia --project=. experiments/PartitionedFragIndex/test_hinted_search.jl

using Pioneer

const SCRIPT_DIR = @__DIR__
include(joinpath(SCRIPT_DIR, "partitioned_types.jl"))
include(joinpath(SCRIPT_DIR, "build_partitioned_index.jl"))
include(joinpath(SCRIPT_DIR, "search_partitioned_index.jl"))

# ─── Helpers ─────────────────────────────────────────────────────────────────

"""Build a uniform frag bin array: N bins of width `bin_width` Da starting at `start_mz`.
Each bin has `frags_per_bin` fragments with local_ids cycling 1:n_precs and score=1."""
function make_uniform_frag_bins(;
        n_bins::Int, start_mz::Float32, bin_width::Float32,
        frags_per_bin::Int, n_precs::Int)
    frag_bins = Vector{Pioneer.FragIndexBin{Float32}}(undef, n_bins)
    fragments = LocalFragment[]
    frag_idx = UInt32(0)
    for i in 1:n_bins
        lo = start_mz + (i - 1) * bin_width
        hi = lo + bin_width
        first_f = frag_idx + one(UInt32)
        for j in 1:frags_per_bin
            frag_idx += one(UInt32)
            lid = UInt16(((i - 1) * frags_per_bin + j - 1) % n_precs + 1)
            push!(fragments, LocalFragment(lid, UInt8(1)))
        end
        frag_bins[i] = Pioneer.FragIndexBin{Float32}(lo, hi, first_f, frag_idx)
    end
    return frag_bins, fragments
end

"""Build skip hints for a uniform bin array (constant density)."""
function make_uniform_hints(frag_bins::Vector{Pioneer.FragIndexBin{Float32}})
    n = length(frag_bins)
    hints = ones(UInt8, n)
    for j in 1:n
        look_ahead = min(50, n - j)
        if look_ahead > 0
            mz_span = Pioneer.getLow(frag_bins[j + look_ahead]) - Pioneer.getLow(frag_bins[j])
            if mz_span > 0.0f0
                bpd = Float64(look_ahead) / Float64(mz_span)
                hints[j] = UInt8(clamp(round(Int, bpd), 1, 255))
            end
        end
    end
    return hints
end

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

"""Build a synthetic LocalPartition with uniform bins."""
function make_test_partition(;
        n_frag_bins::Int, start_mz::Float32, bin_width::Float32,
        frags_per_bin::Int, n_precs::Int,
        n_rt_bins::Int=1, rt_lo::Float32=0.0f0, rt_hi::Float32=10.0f0)
    frag_bins, fragments = make_uniform_frag_bins(;
        n_bins=n_frag_bins, start_mz=start_mz, bin_width=bin_width,
        frags_per_bin=frags_per_bin, n_precs=n_precs)

    # Build RT bins that each span equal portions of the frag bins
    rt_bins = Pioneer.FragIndexBin{Float32}[]
    bins_per_rt = max(1, n_frag_bins ÷ n_rt_bins)
    for r in 1:n_rt_bins
        rlo = rt_lo + (r - 1) * (rt_hi - rt_lo) / n_rt_bins
        rhi = rt_lo + r * (rt_hi - rt_lo) / n_rt_bins
        fb_first = UInt32((r - 1) * bins_per_rt + 1)
        fb_last = r == n_rt_bins ? UInt32(n_frag_bins) : UInt32(r * bins_per_rt)
        push!(rt_bins, Pioneer.FragIndexBin{Float32}(rlo, rhi, fb_first, fb_last))
    end

    hints = make_uniform_hints(frag_bins)
    l2g = UInt32.(1:n_precs)  # identity mapping

    return LocalPartition{Float32}(
        frag_bins, rt_bins, fragments, l2g, UInt16(n_precs), hints)
end

# ─── Test: queryFragmentHinted! vs queryFragmentPartitioned! ─────────────────

"""
Test 1: Single query — first peak (prev_mz=0).
Hinted search should degenerate to step_size=1 and produce identical results.
"""
function test_first_peak()
    frag_bins, fragments = make_uniform_frag_bins(;
        n_bins=100, start_mz=100.0f0, bin_width=0.5f0,
        frags_per_bin=5, n_precs=20)
    hints = make_uniform_hints(frag_bins)
    max_idx = UInt32(length(frag_bins))

    # Query near the middle
    frag_mz_min, frag_mz_max = 125.0f0, 125.5f0

    # Reference: queryFragmentPartitioned!
    lc_ref = LocalCounter(UInt16, UInt8, 21)
    lb_ref, ub_ref = queryFragmentPartitioned!(lc_ref, max_idx,
        UInt32(1), UInt32(1), frag_bins, fragments, frag_mz_min, frag_mz_max)
    ref_scores = extract_local_scores(lc_ref)

    # Hinted: prev_mz=0 → est_step=1
    lc_hnt = LocalCounter(UInt16, UInt8, 21)
    lb_hnt, ub_hnt = queryFragmentHinted!(lc_hnt, max_idx,
        UInt32(1), UInt32(1), frag_bins, fragments, frag_mz_min, frag_mz_max,
        hints, 0.0f0, 16)
    hnt_scores = extract_local_scores(lc_hnt)

    @assert ref_scores == hnt_scores "First peak: scores differ\n  ref=$ref_scores\n  hnt=$hnt_scores"
end

"""
Test 2: Consecutive peaks — different hint values produce identical results.
Since lb is preserved (not advanced by exponential search), the step_size/hint
only affects how quickly we find the upper bound, not which bins are scored.
"""
function test_consecutive_peaks()
    frag_bins, fragments = make_uniform_frag_bins(;
        n_bins=200, start_mz=100.0f0, bin_width=0.5f0,
        frags_per_bin=3, n_precs=15)
    max_idx = UInt32(length(frag_bins))

    peaks = Float32[110.0, 112.0, 115.0, 120.0, 130.0, 150.0, 180.0, 195.0]
    tol = 0.3f0

    # Build two different hint vectors — one accurate, one all-ones (worst case)
    hints_good = make_uniform_hints(frag_bins)
    hints_bad = ones(UInt8, length(frag_bins))  # no useful hint info

    lb_good, ub_good = UInt32(1), UInt32(1)
    lb_bad, ub_bad = UInt32(1), UInt32(1)
    prev_mz = 0.0f0

    for peak in peaks
        frag_mz_min = peak - tol
        frag_mz_max = peak + tol

        lc_good = LocalCounter(UInt16, UInt8, 16)
        lb_good, ub_good = queryFragmentHinted!(lc_good, max_idx,
            lb_good, ub_good, frag_bins, fragments, frag_mz_min, frag_mz_max,
            hints_good, prev_mz, 16)
        scores_good = extract_local_scores(lc_good)

        lc_bad = LocalCounter(UInt16, UInt8, 16)
        lb_bad, ub_bad = queryFragmentHinted!(lc_bad, max_idx,
            lb_bad, ub_bad, frag_bins, fragments, frag_mz_min, frag_mz_max,
            hints_bad, prev_mz, 16)
        scores_bad = extract_local_scores(lc_bad)

        @assert scores_good == scores_bad (
            "Consecutive peaks: good vs bad hints differ at peak=$peak\n" *
            "  good=$scores_good\n  bad=$scores_bad")

        prev_mz = frag_mz_min
    end
end

"""
Test 3: Single frag bin — both should find and score exactly one bin.
"""
function test_single_bin()
    frag_bins = [Pioneer.FragIndexBin{Float32}(500.0f0, 500.5f0, UInt32(1), UInt32(3))]
    fragments = [LocalFragment(UInt16(1), UInt8(2)),
                 LocalFragment(UInt16(2), UInt8(3)),
                 LocalFragment(UInt16(1), UInt8(1))]
    hints = UInt8[1]
    max_idx = UInt32(1)

    # Query that hits the bin
    frag_mz_min, frag_mz_max = 500.0f0, 500.5f0

    lc_ref = LocalCounter(UInt16, UInt8, 3)
    queryFragmentPartitioned!(lc_ref, max_idx,
        UInt32(1), UInt32(1), frag_bins, fragments, frag_mz_min, frag_mz_max)
    ref = extract_local_scores(lc_ref)

    lc_hnt = LocalCounter(UInt16, UInt8, 3)
    queryFragmentHinted!(lc_hnt, max_idx,
        UInt32(1), UInt32(1), frag_bins, fragments, frag_mz_min, frag_mz_max,
        hints, 0.0f0, 16)
    hnt = extract_local_scores(lc_hnt)

    @assert ref == hnt "Single bin: scores differ\n  ref=$ref\n  hnt=$hnt"
    @assert ref == Dict{UInt16,UInt8}(1 => 3, 2 => 3) "Single bin: unexpected scores $ref"
end

"""
Test 4: Query that misses all bins — both should return empty scores.
"""
function test_no_match()
    frag_bins, fragments = make_uniform_frag_bins(;
        n_bins=50, start_mz=100.0f0, bin_width=1.0f0,
        frags_per_bin=2, n_precs=5)
    hints = make_uniform_hints(frag_bins)
    max_idx = UInt32(length(frag_bins))

    # Query way above all bins
    frag_mz_min, frag_mz_max = 500.0f0, 501.0f0

    lc_ref = LocalCounter(UInt16, UInt8, 6)
    queryFragmentPartitioned!(lc_ref, max_idx,
        UInt32(1), UInt32(1), frag_bins, fragments, frag_mz_min, frag_mz_max)
    ref = extract_local_scores(lc_ref)

    lc_hnt = LocalCounter(UInt16, UInt8, 6)
    queryFragmentHinted!(lc_hnt, max_idx,
        UInt32(1), UInt32(1), frag_bins, fragments, frag_mz_min, frag_mz_max,
        hints, 90.0f0, 16)
    hnt = extract_local_scores(lc_hnt)

    @assert isempty(ref) "No match ref: expected empty, got $ref"
    @assert isempty(hnt) "No match hnt: expected empty, got $hnt"
end

"""
Test 5: Negative delta_mz (prev_mz > frag_mz_min) — should use step_size=1.
"""
function test_negative_delta()
    frag_bins, fragments = make_uniform_frag_bins(;
        n_bins=50, start_mz=100.0f0, bin_width=1.0f0,
        frags_per_bin=3, n_precs=10)
    hints = make_uniform_hints(frag_bins)
    max_idx = UInt32(length(frag_bins))

    frag_mz_min, frag_mz_max = 110.0f0, 110.5f0
    prev_mz = 120.0f0  # higher than query — negative delta

    lc_ref = LocalCounter(UInt16, UInt8, 11)
    queryFragmentPartitioned!(lc_ref, max_idx,
        UInt32(1), UInt32(1), frag_bins, fragments, frag_mz_min, frag_mz_max)
    ref = extract_local_scores(lc_ref)

    lc_hnt = LocalCounter(UInt16, UInt8, 11)
    queryFragmentHinted!(lc_hnt, max_idx,
        UInt32(1), UInt32(1), frag_bins, fragments, frag_mz_min, frag_mz_max,
        hints, prev_mz, 16)
    hnt = extract_local_scores(lc_hnt)

    @assert ref == hnt "Negative delta: scores differ\n  ref=$ref\n  hnt=$hnt"
end

"""
Test 6: Large gap — hint should give large step, binary search path after exponential.
"""
function test_large_gap()
    frag_bins, fragments = make_uniform_frag_bins(;
        n_bins=1000, start_mz=100.0f0, bin_width=0.5f0,
        frags_per_bin=2, n_precs=10)
    hints = make_uniform_hints(frag_bins)
    max_idx = UInt32(length(frag_bins))

    # First peak near start
    frag_mz_min1, frag_mz_max1 = 102.0f0, 102.3f0
    lc1_ref = LocalCounter(UInt16, UInt8, 11)
    lb_ref, ub_ref = queryFragmentPartitioned!(lc1_ref, max_idx,
        UInt32(1), UInt32(1), frag_bins, fragments, frag_mz_min1, frag_mz_max1)

    lc1_hnt = LocalCounter(UInt16, UInt8, 11)
    lb_hnt, ub_hnt = queryFragmentHinted!(lc1_hnt, max_idx,
        UInt32(1), UInt32(1), frag_bins, fragments, frag_mz_min1, frag_mz_max1,
        hints, 0.0f0, 16)

    # Second peak far away (~400 bins)
    frag_mz_min2, frag_mz_max2 = 300.0f0, 300.3f0
    lc2_ref = LocalCounter(UInt16, UInt8, 11)
    lb_ref, ub_ref = queryFragmentPartitioned!(lc2_ref, max_idx,
        lb_ref, ub_ref, frag_bins, fragments, frag_mz_min2, frag_mz_max2)
    ref = extract_local_scores(lc2_ref)

    lc2_hnt = LocalCounter(UInt16, UInt8, 11)
    lb_hnt, ub_hnt = queryFragmentHinted!(lc2_hnt, max_idx,
        lb_hnt, ub_hnt, frag_bins, fragments, frag_mz_min2, frag_mz_max2,
        hints, frag_mz_min1, 16)
    hnt = extract_local_scores(lc2_hnt)

    @assert ref == hnt "Large gap: scores differ\n  ref=$ref\n  hnt=$hnt"
end

"""
Test 7: All linear_threshold values produce identical results.
Since lb is preserved, threshold only affects whether we use linear scan
or binary search — both must find the same first matching bin.
"""
function test_linear_threshold_extremes()
    frag_bins, fragments = make_uniform_frag_bins(;
        n_bins=200, start_mz=100.0f0, bin_width=0.5f0,
        frags_per_bin=3, n_precs=15)
    hints = make_uniform_hints(frag_bins)
    max_idx = UInt32(length(frag_bins))

    peaks = Float32[110.0, 115.0, 130.0, 160.0, 190.0]
    tol = 0.3f0

    thresholds = [0, 1, 2, 4, 8, 16, 32, 64, 128, 1000000]

    # Collect results for threshold=0 as reference
    ref_all = Dict{Float32, Dict{UInt16, UInt8}}()
    lb_ref, ub_ref = UInt32(1), UInt32(1)
    prev_mz = 0.0f0
    for peak in peaks
        frag_mz_min = peak - tol
        frag_mz_max = peak + tol
        lc = LocalCounter(UInt16, UInt8, 16)
        lb_ref, ub_ref = queryFragmentHinted!(lc, max_idx,
            lb_ref, ub_ref, frag_bins, fragments, frag_mz_min, frag_mz_max,
            hints, prev_mz, 0)
        ref_all[peak] = extract_local_scores(lc)
        prev_mz = frag_mz_min
    end

    # Every other threshold must match
    for threshold in thresholds[2:end]
        lb, ub = UInt32(1), UInt32(1)
        prev_mz = 0.0f0
        for peak in peaks
            frag_mz_min = peak - tol
            frag_mz_max = peak + tol
            lc = LocalCounter(UInt16, UInt8, 16)
            lb, ub = queryFragmentHinted!(lc, max_idx,
                lb, ub, frag_bins, fragments, frag_mz_min, frag_mz_max,
                hints, prev_mz, threshold)
            hnt = extract_local_scores(lc)

            @assert ref_all[peak] == hnt (
                "Threshold=$threshold vs 0: scores differ at peak=$peak\n" *
                "  ref=$(ref_all[peak])\n  hnt=$hnt")
            prev_mz = frag_mz_min
        end
    end
end

"""
Test 8: Full _score_partition_hinted! vs _score_partition! on a synthetic partition.
Exercises the complete RT bin + frag bin + scoring pipeline.
"""
function test_score_partition_full()
    partition = make_test_partition(;
        n_frag_bins=500, start_mz=100.0f0, bin_width=0.5f0,
        frags_per_bin=4, n_precs=30,
        n_rt_bins=3, rt_lo=0.0f0, rt_hi=15.0f0)

    masses = Union{Missing, Float32}[
        110.0f0, 125.0f0, 150.0f0, 175.0f0, 200.0f0, 225.0f0, 250.0f0, 300.0f0]
    mem = Pioneer.MassErrorModel(0.0f0, (10.0f0, 10.0f0))

    # Reference
    lc_ref = LocalCounter(UInt16, UInt8, 31)
    _score_partition!(lc_ref, partition, 0.0f0, 10.0f0, masses, mem)
    ref = extract_local_scores(lc_ref)
    reset!(lc_ref)

    # Hinted (default threshold)
    lc_hnt = LocalCounter(UInt16, UInt8, 31)
    _score_partition_hinted!(lc_hnt, partition, 0.0f0, 10.0f0, masses, mem)
    hnt = extract_local_scores(lc_hnt)
    reset!(lc_hnt)

    @assert !isempty(ref) "Full partition: reference found no scores"
    @assert ref == hnt "Full partition: scores differ\n  ref IDs: $(sort(collect(keys(ref))))\n  hnt IDs: $(sort(collect(keys(hnt))))"

    # Also test with different thresholds
    for threshold in [0, 4, 32, 128, 1000000]
        lc_t = LocalCounter(UInt16, UInt8, 31)
        _score_partition_hinted!(lc_t, partition, 0.0f0, 10.0f0, masses, mem;
                                  linear_threshold=threshold)
        t_scores = extract_local_scores(lc_t)
        reset!(lc_t)
        @assert t_scores == ref (
            "Full partition threshold=$threshold: scores differ from reference")
    end
end

"""
Test 9: Non-uniform bin widths — bins get wider at higher m/z.
Hints should reflect varying density correctly.
"""
function test_nonuniform_bins()
    # Build bins with increasing width: 0.1, 0.2, 0.3, ... Da
    frag_bins = Pioneer.FragIndexBin{Float32}[]
    fragments = LocalFragment[]
    mz = 100.0f0
    frag_idx = UInt32(0)
    n_precs = 10
    for i in 1:50
        width = Float32(0.1 * i)
        lo = mz
        hi = mz + width
        first_f = frag_idx + one(UInt32)
        for j in 1:3
            frag_idx += one(UInt32)
            push!(fragments, LocalFragment(UInt16((i + j - 1) % n_precs + 1), UInt8(1)))
        end
        push!(frag_bins, Pioneer.FragIndexBin{Float32}(lo, hi, first_f, frag_idx))
        mz = hi
    end

    hints = make_uniform_hints(frag_bins)
    max_idx = UInt32(length(frag_bins))

    # Query a few peaks across the range
    peaks = Float32[105.0, 120.0, 150.0, 200.0]
    tol = 1.0f0

    lb_ref, ub_ref = UInt32(1), UInt32(1)
    lb_hnt, ub_hnt = UInt32(1), UInt32(1)
    prev_mz = 0.0f0

    for peak in peaks
        frag_mz_min = peak - tol
        frag_mz_max = peak + tol

        lc_ref = LocalCounter(UInt16, UInt8, n_precs + 1)
        lb_ref, ub_ref = queryFragmentPartitioned!(lc_ref, max_idx,
            lb_ref, ub_ref, frag_bins, fragments, frag_mz_min, frag_mz_max)
        ref = extract_local_scores(lc_ref)

        lc_hnt = LocalCounter(UInt16, UInt8, n_precs + 1)
        lb_hnt, ub_hnt = queryFragmentHinted!(lc_hnt, max_idx,
            lb_hnt, ub_hnt, frag_bins, fragments, frag_mz_min, frag_mz_max,
            hints, prev_mz, 16)
        hnt = extract_local_scores(lc_hnt)

        @assert ref == hnt (
            "Non-uniform peak=$peak: scores differ\n  ref=$ref\n  hnt=$hnt")
        prev_mz = frag_mz_min
    end
end

"""
Test 10: Stress — 2000 bins, 20 sorted peaks spanning the full range.
Compare queryFragmentPartitioned! (step_size=4) against hinted.
"""
function test_stress_many_peaks()
    frag_bins, fragments = make_uniform_frag_bins(;
        n_bins=2000, start_mz=100.0f0, bin_width=0.25f0,
        frags_per_bin=5, n_precs=50)
    hints = make_uniform_hints(frag_bins)
    max_idx = UInt32(length(frag_bins))

    # 20 sorted peaks spanning the range
    peaks = Float32[100.0 + i * 25.0 for i in 0:19]
    tol = 0.2f0

    lb_ref, ub_ref = UInt32(1), UInt32(1)
    lb_hnt, ub_hnt = UInt32(1), UInt32(1)
    prev_mz = 0.0f0

    for peak in peaks
        frag_mz_min = peak - tol
        frag_mz_max = peak + tol

        # Reference: queryFragmentPartitioned! (step_size=4, the default)
        lc_ref = LocalCounter(UInt16, UInt8, 51)
        lb_ref, ub_ref = queryFragmentPartitioned!(lc_ref, max_idx,
            lb_ref, ub_ref, frag_bins, fragments, frag_mz_min, frag_mz_max)
        ref = extract_local_scores(lc_ref)

        # Hinted
        lc_hnt = LocalCounter(UInt16, UInt8, 51)
        lb_hnt, ub_hnt = queryFragmentHinted!(lc_hnt, max_idx,
            lb_hnt, ub_hnt, frag_bins, fragments, frag_mz_min, frag_mz_max,
            hints, prev_mz, 16)
        hnt = extract_local_scores(lc_hnt)

        @assert ref == hnt (
            "Stress peak=$peak: scores differ\n  ref=$ref\n  hnt=$hnt")
        prev_mz = frag_mz_min
    end
end

"""
Test 11: Boundary — queries at first, last, and middle bin edges.
Each query is standalone from (1,1) with prev_mz=0 (first-peak semantics).
"""
function test_boundary_bins()
    frag_bins, fragments = make_uniform_frag_bins(;
        n_bins=10, start_mz=100.0f0, bin_width=1.0f0,
        frags_per_bin=2, n_precs=5)
    hints = make_uniform_hints(frag_bins)
    max_idx = UInt32(length(frag_bins))

    # Single-bin queries from fresh (1,1) bounds
    for (fmz_min, fmz_max, label) in [
        (100.0f0, 101.0f0, "first bin"),
        (109.0f0, 110.0f0, "last bin"),
        (104.5f0, 105.5f0, "middle bin"),
    ]
        lc_ref = LocalCounter(UInt16, UInt8, 6)
        queryFragmentPartitioned!(lc_ref, max_idx,
            UInt32(1), UInt32(1), frag_bins, fragments, fmz_min, fmz_max)
        ref = extract_local_scores(lc_ref)

        lc_hnt = LocalCounter(UInt16, UInt8, 6)
        queryFragmentHinted!(lc_hnt, max_idx,
            UInt32(1), UInt32(1), frag_bins, fragments, fmz_min, fmz_max,
            hints, 0.0f0, 16)
        hnt = extract_local_scores(lc_hnt)

        @assert ref == hnt "Boundary ($label): scores differ\n  ref=$ref\n  hnt=$hnt"
    end

    # Sequential queries with carry-forward (the real usage pattern)
    sorted_queries = [
        (100.0f0, 101.0f0),
        (103.0f0, 104.0f0),
        (106.0f0, 107.0f0),
        (109.0f0, 110.0f0),
    ]
    lb_ref, ub_ref = UInt32(1), UInt32(1)
    lb_hnt, ub_hnt = UInt32(1), UInt32(1)
    prev_mz = 0.0f0

    for (fmz_min, fmz_max) in sorted_queries
        lc_ref = LocalCounter(UInt16, UInt8, 6)
        lb_ref, ub_ref = queryFragmentPartitioned!(lc_ref, max_idx,
            lb_ref, ub_ref, frag_bins, fragments, fmz_min, fmz_max)
        ref = extract_local_scores(lc_ref)

        lc_hnt = LocalCounter(UInt16, UInt8, 6)
        lb_hnt, ub_hnt = queryFragmentHinted!(lc_hnt, max_idx,
            lb_hnt, ub_hnt, frag_bins, fragments, fmz_min, fmz_max,
            hints, prev_mz, 16)
        hnt = extract_local_scores(lc_hnt)

        @assert ref == hnt (
            "Boundary sequential ($fmz_min-$fmz_max): scores differ\n  ref=$ref\n  hnt=$hnt")
        prev_mz = fmz_min
    end
end

"""
Test 12: Full partition with real-ish MassErrorModel — ppm-based tolerances.
"""
function test_partition_with_ppm_tolerance()
    partition = make_test_partition(;
        n_frag_bins=300, start_mz=200.0f0, bin_width=0.5f0,
        frags_per_bin=3, n_precs=20,
        n_rt_bins=2, rt_lo=0.0f0, rt_hi=10.0f0)

    # PPM-based mass error (20 ppm)
    masses = Union{Missing, Float32}[
        210.0f0, 230.0f0, 260.0f0, 300.0f0, 340.0f0]
    mem = Pioneer.MassErrorModel(0.0f0, (20.0f0, 20.0f0))

    lc_ref = LocalCounter(UInt16, UInt8, 21)
    _score_partition!(lc_ref, partition, 0.0f0, 5.0f0, masses, mem)
    ref = extract_local_scores(lc_ref)
    reset!(lc_ref)

    for threshold in [0, 16, 64, 1000000]
        lc_hnt = LocalCounter(UInt16, UInt8, 21)
        _score_partition_hinted!(lc_hnt, partition, 0.0f0, 5.0f0, masses, mem;
                                  linear_threshold=threshold)
        hnt = extract_local_scores(lc_hnt)
        reset!(lc_hnt)
        @assert ref == hnt (
            "PPM partition threshold=$threshold: scores differ\n" *
            "  ref: $(length(ref)) IDs\n  hnt: $(length(hnt)) IDs\n" *
            "  diff: $(setdiff(keys(ref), keys(hnt)))")
    end
end

# ─── Run all ─────────────────────────────────────────────────────────────────

function run_all_tests()
    println("Running queryFragmentHinted! correctness tests...\n")
    tests = [
        test_first_peak,
        test_consecutive_peaks,
        test_single_bin,
        test_no_match,
        test_negative_delta,
        test_large_gap,
        test_linear_threshold_extremes,
        test_score_partition_full,
        test_nonuniform_bins,
        test_stress_many_peaks,
        test_boundary_bins,
        test_partition_with_ppm_tolerance,
    ]
    n_pass = 0
    n_fail = 0
    for t in tests
        try
            t()
            println("  PASS: $(nameof(t))")
            n_pass += 1
        catch e
            println("  FAIL: $(nameof(t))")
            for (exc, bt) in Base.catch_stack()
                showerror(stdout, exc, bt)
                println()
            end
            n_fail += 1
        end
    end
    println("\n$n_pass passed, $n_fail failed out of $(length(tests)) tests.")
    n_fail > 0 && exit(1)
end

run_all_tests()
