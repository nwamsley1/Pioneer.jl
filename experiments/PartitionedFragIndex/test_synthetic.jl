# Synthetic data tests for the partitioned fragment index.
#
# Validates correctness of:
#   - searchFragmentBinUnconditional!
#   - queryFragmentPartitioned!
#   - searchScanPartitioned!
#   - build_partitioned_index
#   - Post-filter for false positives outside the quad window
#
# Run:  julia --project=. experiments/PartitionedFragIndex/test_synthetic.jl

using Pioneer

const SCRIPT_DIR = @__DIR__
include(joinpath(SCRIPT_DIR, "partitioned_types.jl"))
include(joinpath(SCRIPT_DIR, "build_partitioned_index.jl"))
include(joinpath(SCRIPT_DIR, "search_partitioned_index.jl"))

# ─── Helpers ─────────────────────────────────────────────────────────────────

"""
    make_test_index(rt_specs)

Build a `NativeFragmentIndex{Float32}` from a compact specification.

`rt_specs` is a vector of `(rt_lo, rt_hi, frag_bin_specs)` where each
`frag_bin_specs` is `[(fb_lo, fb_hi, [(prec_id, prec_mz, score), ...]), ...]`.

Fragments within each frag bin are sorted by `prec_mz` (required by
`searchFragmentBin!`'s binary search).
"""
function make_test_index(rt_specs::Vector)
    frags     = Pioneer.IndexFragment{Float32}[]
    frag_bins = Pioneer.FragIndexBin{Float32}[]
    rt_bins   = Pioneer.FragIndexBin{Float32}[]

    for (rt_lo, rt_hi, fb_specs) in rt_specs
        first_fb = UInt32(length(frag_bins) + 1)
        for (fb_lo, fb_hi, frag_specs) in fb_specs
            first_f = UInt32(length(frags) + 1)
            sorted  = sort(collect(frag_specs), by = x -> x[2])
            for (pid, pmz, sc) in sorted
                push!(frags, Pioneer.IndexFragment{Float32}(
                    UInt32(pid), Float32(pmz), UInt8(sc), UInt8(1)))
            end
            last_f = UInt32(length(frags))
            push!(frag_bins, Pioneer.FragIndexBin{Float32}(
                Float32(fb_lo), Float32(fb_hi), first_f, last_f))
        end
        last_fb = UInt32(length(frag_bins))
        push!(rt_bins, Pioneer.FragIndexBin{Float32}(
            Float32(rt_lo), Float32(rt_hi), first_fb, last_fb))
    end

    Pioneer.NativeFragmentIndex{Float32}(frag_bins, rt_bins, frags)
end

"""Extract non-zero scores from a Counter as Dict{UInt32, UInt8}."""
function extract_scores(counter::Pioneer.Counter{UInt32, UInt8})
    scores = Dict{UInt32, UInt8}()
    for i in 1:(Pioneer.getSize(counter) - 1)
        id = Pioneer.getID(counter, i)
        id == 0 && continue
        sc = Pioneer.getCount(counter, id)
        sc > 0 && (scores[id] = sc)
    end
    scores
end

"""Run baseline `searchScan!` and return scores dict."""
function run_baseline(nfi, masses_f32, rt_bin_idx, irt_high, quad_func, iso_bounds, max_prec_id)
    counter     = Pioneer.Counter(UInt32, UInt8, max_prec_id)
    masses      = Union{Missing, Float32}[m for m in masses_f32]
    intensities = Union{Missing, Float32}[1000.0f0 for _ in masses_f32]
    mem         = Pioneer.MassErrorModel(0.0f0, (10.0f0, 10.0f0))

    Pioneer.searchScan!(counter,
        Pioneer.getRTBins(nfi), Pioneer.getFragBins(nfi), Pioneer.getFragments(nfi),
        masses, intensities, Int64(rt_bin_idx), Float32(irt_high),
        mem, quad_func, iso_bounds)

    extract_scores(counter)
end

"""
Run partitioned `searchScanPartitioned!` and return scores dict.
When `precursor_mzs` is provided, applies post-filter (mirrors the logic in
`searchFragmentIndexPartitioned`).
"""
function run_partitioned(pfi, masses_f32, irt_low, irt_high, quad_func, iso_bounds, max_prec_id;
                         precursor_mzs=nothing)
    counter = Pioneer.Counter(UInt32, UInt8, max_prec_id)
    masses  = Union{Missing, Float32}[m for m in masses_f32]
    mem     = Pioneer.MassErrorModel(0.0f0, (10.0f0, 10.0f0))

    searchScanPartitioned!(counter, pfi, Float32(irt_low), Float32(irt_high),
                           masses, mem, quad_func, iso_bounds)

    if precursor_mzs !== nothing
        prec_min = Float32(Pioneer.getPrecMinBound(quad_func) -
                           Pioneer.NEUTRON * first(iso_bounds) / 2)
        prec_max = Float32(Pioneer.getPrecMaxBound(quad_func) +
                           Pioneer.NEUTRON * last(iso_bounds) / 2)
        @inbounds for idx in 1:(Pioneer.getSize(counter) - 1)
            pid = Pioneer.getID(counter, idx)
            pid == 0 && continue
            pmz = precursor_mzs[pid]
            if pmz < prec_min || pmz > prec_max
                counter.counts[pid] = zero(UInt8)
            end
        end
    end

    extract_scores(counter)
end

# ─── Tests ───────────────────────────────────────────────────────────────────

"""
Test 1: Single RT bin, single frag bin, all precursors in window.
Expected: counter scores = {1→3, 2→5, 3→2}.
"""
function test_single_bin()
    nfi = make_test_index([
        (0.0, 10.0, [
            (499.99, 500.01, [(1, 505.0, 3), (2, 505.0, 5), (3, 505.0, 2)])
        ])
    ])
    masses = Float32[500.0]
    quad   = Pioneer.SquareQuadFunction(500.0f0, 510.0f0, 505.0f0)
    iso    = (UInt8(0), UInt8(0))
    expected = Dict{UInt32,UInt8}(1=>3, 2=>5, 3=>2)

    baseline = run_baseline(nfi, masses, 1, 10.0f0, quad, iso, 10)
    @assert baseline == expected "Baseline wrong: $baseline (expected $expected)"

    pfi  = build_partitioned_index(nfi; partition_width=10.0f0)
    part = run_partitioned(pfi, masses, 0.0f0, 10.0f0, quad, iso, 10)
    @assert part == expected "Partitioned wrong: $part (expected $expected)"
end

"""
Test 2: Multiple frag bins — peaks must hit the correct bins.
prec_id=1 gets score contributions from 3 different fragment m/z bins.
Expected total for prec_id=1 = 1+2+4 = 7.
"""
function test_multiple_frag_bins()
    nfi = make_test_index([
        (0.0, 10.0, [
            ( 99.99, 100.01, [(1, 505.0, 1)]),
            (199.99, 200.01, [(1, 505.0, 2)]),
            (299.99, 300.01, [(1, 505.0, 4)])
        ])
    ])
    masses   = Float32[100.0, 200.0, 300.0]
    quad     = Pioneer.SquareQuadFunction(500.0f0, 510.0f0, 505.0f0)
    iso      = (UInt8(0), UInt8(0))
    expected = Dict{UInt32,UInt8}(1 => 7)

    baseline = run_baseline(nfi, masses, 1, 10.0f0, quad, iso, 10)
    @assert baseline == expected "Baseline wrong: $baseline (expected $expected)"

    pfi  = build_partitioned_index(nfi; partition_width=20.0f0)
    part = run_partitioned(pfi, masses, 0.0f0, 10.0f0, quad, iso, 10)
    @assert part == expected "Partitioned wrong: $part (expected $expected)"
end

"""
Test 3: Precursors at a partition boundary.
partition_width=10, anchors at prec_mz=500 and 520 force 2 partitions.
prec_id=1 (prec_mz=509.9) lands in partition 1; prec_id=2 (510.1) in partition 2.
Quad window [505, 515] spans both — both should be scored after post-filter.
"""
function test_partition_boundary()
    nfi = make_test_index([
        (0.0, 10.0, [
            (499.99, 500.01, [
                (1,   509.9, 3),
                (2,   510.1, 4),
                (99,  500.0, 1),   # anchor low
                (100, 520.0, 1),   # anchor high
            ])
        ])
    ])
    masses = Float32[500.0]
    quad   = Pioneer.SquareQuadFunction(505.0f0, 515.0f0, 510.0f0)
    iso    = (UInt8(0), UInt8(0))

    # Baseline: binary search on prec_mz filters to [505, 515]
    baseline = run_baseline(nfi, masses, 1, 10.0f0, quad, iso, 200)
    @assert baseline[UInt32(1)]  == 3 "Baseline prec_id=1 wrong"
    @assert baseline[UInt32(2)]  == 4 "Baseline prec_id=2 wrong"
    @assert !haskey(baseline, UInt32(99))  "Baseline should exclude prec_id=99"
    @assert !haskey(baseline, UInt32(100)) "Baseline should exclude prec_id=100"

    pfi = build_partitioned_index(nfi; partition_width=10.0f0)
    @assert getNPartitions(pfi) == 2 "Expected 2 partitions, got $(getNPartitions(pfi))"

    precursor_mzs = zeros(Float32, 200)
    precursor_mzs[1]   = 509.9f0
    precursor_mzs[2]   = 510.1f0
    precursor_mzs[99]  = 500.0f0
    precursor_mzs[100] = 520.0f0

    part = run_partitioned(pfi, masses, 0.0f0, 10.0f0, quad, iso, 200;
                           precursor_mzs=precursor_mzs)
    @assert part[UInt32(1)] == 3 "Partitioned prec_id=1 wrong: $(get(part, UInt32(1), 0))"
    @assert part[UInt32(2)] == 4 "Partitioned prec_id=2 wrong: $(get(part, UInt32(2), 0))"
    @assert !haskey(part, UInt32(99))  "Post-filter should remove prec_id=99"
    @assert !haskey(part, UInt32(100)) "Post-filter should remove prec_id=100"
end

"""
Test 4: Post-filter removes false positives.
Wide partition (30 Da) keeps prec_id=1 (prec_mz=505) and prec_id=2 (prec_mz=522)
in the same partition. Quad window [500, 510] should only pass prec_id=1.
"""
function test_post_filter()
    nfi = make_test_index([
        (0.0, 10.0, [
            (499.99, 500.01, [(1, 505.0, 3), (2, 522.0, 5)])
        ])
    ])
    masses = Float32[500.0]
    quad   = Pioneer.SquareQuadFunction(500.0f0, 510.0f0, 505.0f0)
    iso    = (UInt8(0), UInt8(0))

    # Baseline: prec_id=2 outside [500,510] → not scored
    baseline = run_baseline(nfi, masses, 1, 10.0f0, quad, iso, 200)
    @assert baseline == Dict{UInt32,UInt8}(1 => 3) "Baseline wrong: $baseline"

    pfi = build_partitioned_index(nfi; partition_width=30.0f0)
    @assert getNPartitions(pfi) == 1 "Expected 1 partition"

    # Without post-filter: unconditional scoring includes prec_id=2
    raw = run_partitioned(pfi, masses, 0.0f0, 10.0f0, quad, iso, 200)
    @assert haskey(raw, UInt32(2)) "Without filter, prec_id=2 should be scored"

    # With post-filter: prec_id=2 removed
    precursor_mzs = zeros(Float32, 200)
    precursor_mzs[1] = 505.0f0
    precursor_mzs[2] = 522.0f0

    filtered = run_partitioned(pfi, masses, 0.0f0, 10.0f0, quad, iso, 200;
                               precursor_mzs=precursor_mzs)
    @assert filtered == Dict{UInt32,UInt8}(1 => 3) "Post-filter failed: $filtered"
end

"""
Test 5: Multiple RT bins — verify RT windowing works.
3 RT bins with scores 1, 2, 3. With irt_high=9.99, only bins 1 and 2 are searched
(total=3). With irt_high=15.0, all three (total=6).
"""
function test_multiple_rt_bins()
    nfi = make_test_index([
        ( 0.0,  5.0, [(499.99, 500.01, [(1, 505.0, 1)])]),
        ( 5.0, 10.0, [(499.99, 500.01, [(1, 505.0, 2)])]),
        (10.0, 15.0, [(499.99, 500.01, [(1, 505.0, 3)])]),
    ])
    masses = Float32[500.0]
    quad   = Pioneer.SquareQuadFunction(500.0f0, 510.0f0, 505.0f0)
    iso    = (UInt8(0), UInt8(0))

    pfi = build_partitioned_index(nfi; partition_width=20.0f0)

    # irt_high=9.99 → bins 1 & 2 → total 3
    bl_2 = run_baseline(nfi, masses, 1, 9.99f0, quad, iso, 10)
    @assert bl_2 == Dict{UInt32,UInt8}(1 => 3) "Baseline (2 bins) wrong: $bl_2"
    pt_2 = run_partitioned(pfi, masses, 0.0f0, 9.99f0, quad, iso, 10)
    @assert pt_2 == bl_2 "Partitioned (2 bins) wrong: $pt_2 vs $bl_2"

    # irt_high=15.0 → all 3 bins → total 6
    bl_3 = run_baseline(nfi, masses, 1, 15.0f0, quad, iso, 10)
    @assert bl_3 == Dict{UInt32,UInt8}(1 => 6) "Baseline (3 bins) wrong: $bl_3"
    pt_3 = run_partitioned(pfi, masses, 0.0f0, 15.0f0, quad, iso, 10)
    @assert pt_3 == bl_3 "Partitioned (3 bins) wrong: $pt_3 vs $bl_3"
end

"""
Test 6: `build_partitioned_index` structural correctness.
Index with 11 precursors spanning prec_mz [500, 600], partition_width=25 → 4 partitions.
Checks: total fragment count, partition assignment, RT bin count, frag bin bounds.
"""
function test_build_structural()
    frag_specs = [(Int(pmz - 499), pmz, 2) for pmz in 500.0:10.0:600.0]  # 11 entries

    nfi = make_test_index([
        (0.0, 5.0, [
            ( 99.99, 100.01, frag_specs),
            (199.99, 200.01, frag_specs),
        ]),
        (5.0, 10.0, [
            ( 99.99, 100.01, frag_specs),
        ]),
    ])

    pfi = build_partitioned_index(nfi; partition_width=25.0f0)

    # Total fragment count preserved
    orig_count = length(Pioneer.getFragments(nfi))
    part_count = sum(length(Pioneer.getFragments(p)) for p in getPartitions(pfi))
    @assert part_count == orig_count "Fragment count: $part_count vs $orig_count"

    # Each fragment's prec_mz falls within its partition's range
    for (k, part) in enumerate(getPartitions(pfi))
        lo = pfi.prec_mz_min + (k - 1) * pfi.partition_width
        hi = pfi.prec_mz_min + k * pfi.partition_width
        for frag in Pioneer.getFragments(part)
            pmz = Pioneer.getPrecMZ(frag)
            if k < pfi.n_partitions
                @assert lo <= pmz < hi "P$k: prec_mz=$pmz outside [$lo, $hi)"
            else
                @assert lo <= pmz <= hi "P$k (last): prec_mz=$pmz outside [$lo, $hi]"
            end
        end
    end

    # Every partition has the same number of RT bins as the original
    orig_n_rt = length(Pioneer.getRTBins(nfi))
    for (k, part) in enumerate(getPartitions(pfi))
        @assert length(Pioneer.getRTBins(part)) == orig_n_rt (
            "P$k: $(length(Pioneer.getRTBins(part))) RT bins vs $orig_n_rt")
    end

    # Frag bin lb/ub values come from the original
    orig_lbs = Set(Pioneer.getLow(fb) for fb in Pioneer.getFragBins(nfi))
    orig_ubs = Set(Pioneer.getHigh(fb) for fb in Pioneer.getFragBins(nfi))
    for part in getPartitions(pfi)
        for fb in Pioneer.getFragBins(part)
            @assert Pioneer.getLow(fb) in orig_lbs "Unknown frag bin lb: $(Pioneer.getLow(fb))"
            @assert Pioneer.getHigh(fb) in orig_ubs "Unknown frag bin ub: $(Pioneer.getHigh(fb))"
        end
    end

    # Empty RT bins have empty sub-bin ranges (first > last)
    for (k, part) in enumerate(getPartitions(pfi))
        for rt_bin in Pioneer.getRTBins(part)
            r = Pioneer.getSubBinRange(rt_bin)
            if isempty(r)
                @assert first(r) > last(r) "P$k: empty range but first <= last"
            end
        end
    end
end

"""
Test 7: Empty partitions / empty frag bins.
RT bin 1 has fragments only for prec_mz=505 (partition 1); RT bin 2 only for
prec_mz=555 and 606 (partitions 2, 3). Partitions with empty RT bins must not
crash and must produce correct scores.
"""
function test_empty_partitions()
    nfi = make_test_index([
        (0.0, 5.0, [
            (499.99, 500.01, [(1, 505.0, 3)])
        ]),
        (5.0, 10.0, [
            (499.99, 500.01, [(2, 555.0, 4), (3, 606.0, 1)])
        ]),
    ])

    pfi = build_partitioned_index(nfi; partition_width=50.0f0)
    @assert getNPartitions(pfi) >= 2 "Expected ≥ 2 partitions, got $(getNPartitions(pfi))"

    masses = Float32[500.0]
    iso    = (UInt8(0), UInt8(0))

    precursor_mzs = zeros(Float32, 200)
    precursor_mzs[1] = 505.0f0
    precursor_mzs[2] = 555.0f0
    precursor_mzs[3] = 606.0f0

    # Narrow quad [500, 510] + irt_high covering only RT bin 1
    quad_narrow = Pioneer.SquareQuadFunction(500.0f0, 510.0f0, 505.0f0)
    bl = run_baseline(nfi, masses, 1, 4.99f0, quad_narrow, iso, 200)
    @assert bl == Dict{UInt32,UInt8}(1 => 3) "Narrow baseline wrong: $bl"
    pt = run_partitioned(pfi, masses, 0.0f0, 4.99f0, quad_narrow, iso, 200;
                         precursor_mzs=precursor_mzs)
    @assert pt == bl "Narrow partitioned ($pt) != baseline ($bl)"

    # Narrow quad but irt_high=6.0 — walks into partition 1's empty RT bin 2
    pt2 = run_partitioned(pfi, masses, 0.0f0, 6.0f0, quad_narrow, iso, 200;
                          precursor_mzs=precursor_mzs)
    @assert pt2 == Dict{UInt32,UInt8}(1 => 3) "Empty RT bin caused wrong scores: $pt2"

    # Wide quad [500, 610] + both RT bins
    quad_wide = Pioneer.SquareQuadFunction(500.0f0, 610.0f0, 555.0f0)
    bl_w = run_baseline(nfi, masses, 1, 10.0f0, quad_wide, iso, 200)
    pt_w = run_partitioned(pfi, masses, 0.0f0, 10.0f0, quad_wide, iso, 200;
                           precursor_mzs=precursor_mzs)
    @assert pt_w == bl_w "Wide search: partitioned ($pt_w) != baseline ($bl_w)"
end

"""
Test 8: Stress test with many frag bins per RT bin.
200 frag bins, 20 precursors per bin spanning prec_mz [500, 690], 2 RT bins,
8 query masses. Exercises exponential search carry-over between masses and
verifies partitioning doesn't break frag-bin lookup.
"""
function test_many_frag_bins_stress()
    n_fb = 200
    n_precs = 20

    fb_specs = []
    for i in 1:n_fb
        fb_lo = 100.0 + (i - 1) * 0.5
        fb_hi = fb_lo + 0.5
        sc = (i % 3) + 1
        frags = [(j, 500.0 + (j - 1) * 10.0, sc) for j in 1:n_precs]
        push!(fb_specs, (fb_lo, fb_hi, frags))
    end

    nfi = make_test_index([
        (0.0, 5.0, fb_specs),
        (5.0, 10.0, fb_specs),
    ])

    masses = Float32[102.25, 115.75, 128.25, 141.75, 155.25, 168.75, 182.25, 195.75]
    quad = Pioneer.SquareQuadFunction(490.0f0, 710.0f0, 600.0f0)
    iso  = (UInt8(0), UInt8(0))

    baseline = run_baseline(nfi, masses, 1, 10.0f0, quad, iso, n_precs + 5)
    @assert !isempty(baseline) "Stress: baseline found no scores"

    pfi = build_partitioned_index(nfi; partition_width=10.0f0)

    precursor_mzs = zeros(Float32, n_precs + 5)
    for j in 1:n_precs
        precursor_mzs[j] = Float32(500.0 + (j - 1) * 10.0)
    end

    part = run_partitioned(pfi, masses, 0.0f0, 10.0f0, quad, iso, n_precs + 5;
                           precursor_mzs=precursor_mzs)
    @assert part == baseline (
        "Stress: partitioned != baseline\n" *
        "  Mismatched keys: $(setdiff(keys(baseline), keys(part)) ∪ setdiff(keys(part), keys(baseline)))\n" *
        "  Sample diffs: $([(k, baseline[k], get(part, k, 0x00)) for k in first(collect(keys(baseline)), 5)])")
end

"""
Test 9: Sparse partitions — after partitioning, each partition's frag bins have gaps.
Odd frag bins → prec_id=1 (prec_mz=505), even → prec_id=2 (prec_mz=530).
partition_width=10 puts them in different partitions. Verifies the exponential search
handles sparse frag bin distributions correctly.
"""
function test_sparse_partitions()
    fb_specs = []
    for i in 1:10
        fb_lo = 100.0 + (i - 1)
        fb_hi = fb_lo + 1.0
        if isodd(i)
            frags = [(1, 505.0, 2)]
        else
            frags = [(2, 530.0, 3)]
        end
        push!(fb_specs, (fb_lo, fb_hi, frags))
    end

    nfi = make_test_index([(0.0, 10.0, fb_specs)])

    # 10 masses, one per frag bin
    masses = Float32[Float32(100.0 + i - 0.5) for i in 1:10]
    quad = Pioneer.SquareQuadFunction(500.0f0, 535.0f0, 517.5f0)
    iso  = (UInt8(0), UInt8(0))

    baseline = run_baseline(nfi, masses, 1, 10.0f0, quad, iso, 10)
    # prec_id=1: 5 odd bins × score 2 = 10
    # prec_id=2: 5 even bins × score 3 = 15
    @assert baseline == Dict{UInt32,UInt8}(1 => 10, 2 => 15) "Sparse baseline wrong: $baseline"

    pfi = build_partitioned_index(nfi; partition_width=10.0f0)

    precursor_mzs = zeros(Float32, 10)
    precursor_mzs[1] = 505.0f0
    precursor_mzs[2] = 530.0f0

    part = run_partitioned(pfi, masses, 0.0f0, 10.0f0, quad, iso, 10;
                           precursor_mzs=precursor_mzs)
    @assert part == baseline "Sparse: partitioned ($part) != baseline ($baseline)"
end

"""
Test 10: Narrow quad window with isotope error bounds.
Isotope error bounds expand the effective precursor window by ±NEUTRON*bounds/2.
Verifies both paths handle this consistently.
"""
function test_isotope_err_bounds()
    nfi = make_test_index([
        (0.0, 10.0, [
            (499.99, 500.01, [
                (1, 504.0, 3),   # just below quad min=505 but within iso expansion
                (2, 507.0, 5),   # solidly in window
                (3, 511.5, 2),   # just above quad max=510 but within iso expansion
            ])
        ])
    ])

    masses = Float32[500.0]
    quad = Pioneer.SquareQuadFunction(505.0f0, 510.0f0, 507.5f0)
    # iso_bounds=(2,2) → prec_min = 505 - NEUTRON*1 ≈ 504.0, prec_max = 510 + NEUTRON*1 ≈ 511.0
    iso  = (UInt8(2), UInt8(2))

    baseline = run_baseline(nfi, masses, 1, 10.0f0, quad, iso, 10)
    # prec_id=1 (504.0): prec_min ≈ 504.0 → borderline, depends on NEUTRON exact value
    # prec_id=2 (507.0): definitely in
    # prec_id=3 (511.5): prec_max ≈ 511.003 → 511.5 > 511.003 → might be out
    @assert haskey(baseline, UInt32(2)) "Isotope: prec_id=2 must be scored"

    pfi = build_partitioned_index(nfi; partition_width=20.0f0)

    precursor_mzs = zeros(Float32, 10)
    precursor_mzs[1] = 504.0f0
    precursor_mzs[2] = 507.0f0
    precursor_mzs[3] = 511.5f0

    part = run_partitioned(pfi, masses, 0.0f0, 10.0f0, quad, iso, 10;
                           precursor_mzs=precursor_mzs)
    @assert part == baseline "Isotope: partitioned ($part) != baseline ($baseline)"
end

# ─── Run all ─────────────────────────────────────────────────────────────────

function run_all_tests()
    println("Running synthetic partitioned fragment index tests…\n")
    tests = [
        test_single_bin,
        test_multiple_frag_bins,
        test_partition_boundary,
        test_post_filter,
        test_multiple_rt_bins,
        test_build_structural,
        test_empty_partitions,
        test_many_frag_bins_stress,
        test_sparse_partitions,
        test_isotope_err_bounds,
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
                showerror(stdout, exc, bt; backtrace_type=:abbreviated)
                println()
            end
            n_fail += 1
        end
    end
    println("\n$n_pass passed, $n_fail failed out of $(length(tests)) tests.")
    n_fail > 0 && exit(1)
end

run_all_tests()
