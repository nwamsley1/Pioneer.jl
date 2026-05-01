# Unit tests for Counter (and related bitmask filters: CountFilter, LUTFilter,
# PatternAccumulator).
#
# Counter is the array-backed sparse counter used in fragment index search
# (searchFragmentIndexPartitionMajorHinted) and in the fused-path scratch
# buffers (FragIndexScratch). inc!/or! support both additive scoring and
# bitmask-OR scoring; reset! is O(n_active).
#
# These tests are extracted from the former test_arraydict_counter.jl —
# the ArrayDict portion was deleted alongside the classic search path
# (replaced by AbstractPrecursorMap). Counter remains alive and load-bearing.
#
# Run: julia --project=. test/UnitTests/test_counter.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
end

# ═══════════════════════════════════════════════════════════════════════════════
# Counter
# ═══════════════════════════════════════════════════════════════════════════════
#
# Branch map (Counter.jl):
#   constructor           — size=1                                   [test 1]
#   inc! first encounter  — counts[id]==0, so size incremented       [test 2a]
#   inc! repeat encounter — counts[id]>0, size NOT incremented       [test 2b]
#   reset!(c)             — full reset, loop over all ids            [test 3]
#   reset!(c) empty       — size==1, loop body never executes        [test 4]

@testset "Counter" begin

    # ── Test 1: Constructor ──────────────────────────────────────────────
    @testset "constructor" begin
        c = Pioneer.Counter(UInt32, UInt8, 10)
        @test c.size == 1   # 1-indexed: next insertion slot
        @test length(c.ids) == 10
        @test length(c.counts) == 10
        @test all(c.ids .== 0)
        @test all(c.counts .== 0)
    end

    # ── Test 1b: Constructor with UInt16 keys (partitioned index path) ───
    @testset "constructor UInt16 keys" begin
        c = Pioneer.Counter(UInt16, UInt8, 100)
        @test c.size == 1
        @test length(c.ids) == 100
        Pioneer.inc!(c, UInt16(50), UInt8(3))
        @test c.counts[50] == UInt8(3)
        @test c.size == 2
    end

    # ── Test 2: inc! — first encounter vs repeat ─────────────────────────
    @testset "inc! first and repeat encounters" begin
        c = Pioneer.Counter(UInt32, UInt8, 10)

        # First encounter with id=3: count goes 0→2, size increments 1→2
        Pioneer.inc!(c, UInt32(3), UInt8(2))
        @test c.counts[3] == UInt8(2)
        @test c.size == 2         # incremented (was first encounter)
        @test c.ids[1] == UInt32(3)

        # Repeat encounter with id=3: count goes 2→5, size stays 2
        Pioneer.inc!(c, UInt32(3), UInt8(3))
        @test c.counts[3] == UInt8(5)
        @test c.size == 2         # NOT incremented (repeat)

        # First encounter with id=7: count goes 0→1, size increments 2→3
        Pioneer.inc!(c, UInt32(7), UInt8(1))
        @test c.counts[7] == UInt8(1)
        @test c.size == 3
        @test c.ids[2] == UInt32(7)
    end

    # ── Test 3: reset!(c) — full reset ───────────────────────────────────
    @testset "reset! full" begin
        c = Pioneer.Counter(UInt32, UInt8, 10)
        Pioneer.inc!(c, UInt32(2), UInt8(5))
        Pioneer.inc!(c, UInt32(6), UInt8(3))
        Pioneer.inc!(c, UInt32(9), UInt8(7))
        @test c.size == 4  # 3 unique ids + 1

        Pioneer.reset!(c)
        @test c.size == 1
        @test c.counts[2] == 0
        @test c.counts[6] == 0
        @test c.counts[9] == 0
        @test c.ids[1] == 0
        @test c.ids[2] == 0
        @test c.ids[3] == 0
    end

    # ── Test 4: reset!(c) on fresh counter — empty loop ──────────────────
    @testset "reset! on fresh" begin
        c = Pioneer.Counter(UInt32, UInt8, 10)
        @test c.size == 1
        Pioneer.reset!(c)  # loop range is 1:0, body never executes
        @test c.size == 1
    end

    # ── Test 5: Full round-trip — inc, reset, reuse ──────────────────────
    @testset "full round-trip" begin
        c = Pioneer.Counter(UInt32, UInt8, 20)

        # Simulate fragment index search: multiple fragments per precursor
        for _ in 1:4
            Pioneer.inc!(c, UInt32(10), UInt8(1))  # 4 fragments
        end
        for _ in 1:2
            Pioneer.inc!(c, UInt32(15), UInt8(1))  # 2 fragments
        end
        for _ in 1:6
            Pioneer.inc!(c, UInt32(8), UInt8(1))   # 6 fragments
        end
        @test c.size == 4  # 3 unique + 1
        @test c.counts[10] == UInt8(4)
        @test c.counts[15] == UInt8(2)
        @test c.counts[8] == UInt8(6)

        # Reset and reuse for next scan
        Pioneer.reset!(c)
        @test c.size == 1
        @test c.counts[10] == 0
        @test c.counts[15] == 0
        @test c.counts[8] == 0

        # Can insert again
        Pioneer.inc!(c, UInt32(5), UInt8(3))
        @test c.counts[5] == UInt8(3)
        @test c.size == 2
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Counter — or! (bitmask fragment index scoring)
# ═══════════════════════════════════════════════════════════════════════════════
#
# or! is the bitmask variant of inc!. Instead of adding scores, it OR's
# bitmasks so each bit tracks whether a specific fragment rank was matched.
#
# Branch map:
#   or! first encounter   — counts[id]==0, size incremented          [test 1]
#   or! repeat encounter  — counts[id]>0, size NOT incremented       [test 2]
#   or! accumulates via OR, not addition                             [test 3]
#   or! idempotent for same bit                                      [test 4]
#   or! multiple IDs don't interfere                                 [test 5]
#   or! + reset! round-trip                                          [test 6]
#   bitmask scoring: rank → bit position mapping                     [test 7]
#   bitmask scoring: count_ones threshold check                      [test 8]
#   or! vs inc! produce different results                            [test 9]
#   or! performance: no slower than inc! on same workload            [test 10]

@testset "Counter or! (bitmask)" begin

    # ── Test 1: or! first encounter — size increments ────────────────────
    @testset "or! first encounter" begin
        c = Pioneer.Counter(UInt16, UInt8, 100)
        Pioneer.or!(c, UInt16(5), UInt8(0x01))
        @test c.counts[5] == UInt8(0x01)
        @test c.size == 2  # incremented from 1→2
        @test c.ids[1] == UInt16(5)
    end

    # ── Test 2: or! repeat encounter — size does NOT increment ───────────
    @testset "or! repeat encounter" begin
        c = Pioneer.Counter(UInt16, UInt8, 100)
        Pioneer.or!(c, UInt16(5), UInt8(0x01))
        @test c.size == 2
        Pioneer.or!(c, UInt16(5), UInt8(0x02))
        @test c.size == 2  # NOT incremented (repeat)
        @test c.counts[5] == UInt8(0x03)  # 0x01 | 0x02
    end

    # ── Test 3: or! accumulates via OR, not addition ─────────────────────
    @testset "or! uses bitwise OR" begin
        c = Pioneer.Counter(UInt16, UInt8, 100)
        # OR'ing the same mask twice should be idempotent
        Pioneer.or!(c, UInt16(1), UInt8(0x04))
        Pioneer.or!(c, UInt16(1), UInt8(0x04))
        @test c.counts[1] == UInt8(0x04)  # NOT 0x08

        # Whereas inc! would have added: 0x04 + 0x04 = 0x08
        c2 = Pioneer.Counter(UInt16, UInt8, 100)
        Pioneer.inc!(c2, UInt16(1), UInt8(0x04))
        Pioneer.inc!(c2, UInt16(1), UInt8(0x04))
        @test c2.counts[1] == UInt8(0x08)  # addition
    end

    # ── Test 4: or! builds up bitmask from individual rank bits ──────────
    @testset "or! builds bitmask from ranks" begin
        c = Pioneer.Counter(UInt16, UInt8, 100)
        # Simulate fragment matching: each fragment rank sets one bit
        # Rank 1 → bit 0 (0x01), rank 2 → bit 1 (0x02), etc.
        Pioneer.or!(c, UInt16(1), UInt8(0x01))  # rank 1
        Pioneer.or!(c, UInt16(1), UInt8(0x04))  # rank 3
        Pioneer.or!(c, UInt16(1), UInt8(0x10))  # rank 5
        @test c.counts[1] == UInt8(0x15)  # 0b00010101 = ranks 1,3,5
        @test count_ones(c.counts[1]) == 3
        @test c.size == 2  # still just 1 unique ID
    end

    # ── Test 5: multiple IDs don't interfere ─────────────────────────────
    @testset "or! multiple IDs independent" begin
        c = Pioneer.Counter(UInt16, UInt8, 100)
        Pioneer.or!(c, UInt16(1), UInt8(0x0F))  # ID 1: bits 0-3
        Pioneer.or!(c, UInt16(2), UInt8(0xF0))  # ID 2: bits 4-7
        Pioneer.or!(c, UInt16(3), UInt8(0xFF))  # ID 3: all bits
        @test c.counts[1] == UInt8(0x0F)
        @test c.counts[2] == UInt8(0xF0)
        @test c.counts[3] == UInt8(0xFF)
        @test c.size == 4  # 3 unique IDs + 1
    end

    # ── Test 6: or! + reset! round-trip ──────────────────────────────────
    @testset "or! reset round-trip" begin
        c = Pioneer.Counter(UInt16, UInt8, 100)
        Pioneer.or!(c, UInt16(10), UInt8(0xFF))
        Pioneer.or!(c, UInt16(20), UInt8(0x55))
        @test c.size == 3
        @test c.counts[10] == 0xFF
        @test c.counts[20] == 0x55

        Pioneer.reset!(c)
        @test c.size == 1
        @test c.counts[10] == 0x00
        @test c.counts[20] == 0x00

        # Reuse after reset
        Pioneer.or!(c, UInt16(10), UInt8(0x01))
        @test c.counts[10] == 0x01
        @test c.size == 2
    end

    # ── Test 7: bitmask scoring — rank to bit position mapping ───────────
    @testset "rank to bitmask mapping" begin
        # Verify the UInt8(1) << UInt8(rank-1) mapping used in build_poin_lib.jl
        for rank in 1:8
            mask = UInt8(1) << UInt8(rank - 1)
            @test count_ones(mask) == 1
            @test mask == UInt8(2^(rank - 1))
            # Verify the correct bit is set
            @test (mask >> (rank - 1)) & UInt8(1) == UInt8(1)
            # Verify no other bits are set
            for other_rank in 1:8
                if other_rank != rank
                    @test (mask >> (other_rank - 1)) & UInt8(1) == UInt8(0)
                end
            end
        end
    end

    # ── Test 8: count_ones threshold check ───────────────────────────────
    @testset "count_ones threshold" begin
        # Simulate the threshold check from search.jl:
        # count_ones(score) < min_score && continue

        # 3of8 (min_score=3): need at least 3 bits set
        @test count_ones(UInt8(0b00000111)) == 3  # passes
        @test count_ones(UInt8(0b00000011)) == 2  # fails
        @test count_ones(UInt8(0b10100001)) == 3  # passes (non-contiguous)
        @test count_ones(UInt8(0b11111111)) == 8  # passes

        # 7of8 for parameter tuning (min_score=7)
        @test count_ones(UInt8(0b11111110)) == 7  # passes
        @test count_ones(UInt8(0b11111100)) == 6  # fails
        @test count_ones(UInt8(0b11111111)) == 8  # passes

        # Edge cases
        @test count_ones(UInt8(0)) == 0
        @test count_ones(UInt8(0xFF)) == 8
    end

    # ── Test 9: or! vs inc! produce different results ────────────────────
    @testset "or! vs inc! semantics differ" begin
        c_or = Pioneer.Counter(UInt16, UInt8, 100)
        c_inc = Pioneer.Counter(UInt16, UInt8, 100)

        # Same sequence of operations
        masks = UInt8[0x01, 0x02, 0x01, 0x04, 0x02]
        for m in masks
            Pioneer.or!(c_or, UInt16(1), m)
            Pioneer.inc!(c_inc, UInt16(1), m)
        end

        # or! gives the union of all bits: 0x01|0x02|0x04 = 0x07
        @test c_or.counts[1] == UInt8(0x07)
        @test count_ones(c_or.counts[1]) == 3

        # inc! gives the sum: 1+2+1+4+2 = 10
        @test c_inc.counts[1] == UInt8(10)

        # Both track the same unique ID count
        @test c_or.size == c_inc.size
    end

    # ── Test 10: Full bitmask search simulation ──────────────────────────
    @testset "full bitmask search simulation" begin
        # Simulate one scan of fragment index search with bitmask scoring.
        # 3 precursors, each with some fragments matching.
        c = Pioneer.Counter(UInt16, UInt8, 100)

        # Precursor 1: fragments at ranks 1,2,3,5 matched → 0b00010111
        Pioneer.or!(c, UInt16(1), UInt8(1) << UInt8(0))  # rank 1
        Pioneer.or!(c, UInt16(1), UInt8(1) << UInt8(1))  # rank 2
        Pioneer.or!(c, UInt16(1), UInt8(1) << UInt8(2))  # rank 3
        Pioneer.or!(c, UInt16(1), UInt8(1) << UInt8(4))  # rank 5

        # Precursor 2: fragments at ranks 1,4 matched → 0b00001001
        Pioneer.or!(c, UInt16(2), UInt8(1) << UInt8(0))  # rank 1
        Pioneer.or!(c, UInt16(2), UInt8(1) << UInt8(3))  # rank 4

        # Precursor 3: all 8 ranks matched → 0xFF
        for rank in 1:8
            Pioneer.or!(c, UInt16(3), UInt8(1) << UInt8(rank - 1))
        end

        @test c.size == 4  # 3 unique + 1

        # Verify bitmask patterns
        @test c.counts[1] == UInt8(0b00010111)
        @test c.counts[2] == UInt8(0b00001001)
        @test c.counts[3] == UInt8(0xFF)

        # Verify count_ones matches expected fragment counts
        @test count_ones(c.counts[1]) == 4
        @test count_ones(c.counts[2]) == 2
        @test count_ones(c.counts[3]) == 8

        # Apply threshold: min_score=3 (at least 3 fragments matched)
        min_score = UInt8(3)
        passing = UInt16[]
        for i in 1:(c.size - 1)
            lid = c.ids[i]
            score = c.counts[lid]
            count_ones(score) >= min_score && push!(passing, lid)
        end
        @test Set(passing) == Set([UInt16(1), UInt16(3)])  # prec 2 fails (only 2 bits)
        @test UInt16(2) ∉ passing

        # Reset and verify clean state
        Pioneer.reset!(c)
        @test c.counts[1] == 0x00
        @test c.counts[2] == 0x00
        @test c.counts[3] == 0x00
    end

    # ── Test 11: filter_counter! with CountFilter ─────────────────────────
    @testset "filter_counter! CountFilter" begin
        c = Pioneer.Counter(UInt16, UInt8, 100)
        # 3 precursors with different bitmask scores
        Pioneer.or!(c, UInt16(1), UInt8(0x07))  # 3 bits → passes min=3
        Pioneer.or!(c, UInt16(2), UInt8(0x03))  # 2 bits → fails min=3
        Pioneer.or!(c, UInt16(3), UInt8(0xFF))  # 8 bits → passes min=3
        @test c.size == 4

        f = Pioneer.CountFilter(UInt8(3))
        n_pass = Pioneer.filter_counter!(c, f)
        @test n_pass == 2
        # Passing IDs compacted to front
        passing_ids = Set([c.ids[1], c.ids[2]])
        @test passing_ids == Set([UInt16(1), UInt16(3)])
        # Counts still intact for reading
        @test c.counts[1] == UInt8(0x07)
        @test c.counts[3] == UInt8(0xFF)
    end

    # ── Test 12: filter_counter! with LUTFilter ─────────────────────────
    @testset "filter_counter! LUTFilter" begin
        c = Pioneer.Counter(UInt16, UInt8, 100)
        Pioneer.or!(c, UInt16(1), UInt8(0x07))  # pattern 7
        Pioneer.or!(c, UInt16(2), UInt8(0x03))  # pattern 3
        Pioneer.or!(c, UInt16(3), UInt8(0xFF))  # pattern 255

        # Build a LUT that passes only patterns with 3+ bits
        table = Bool[count_ones(UInt8(i)) >= 3 for i in 0:255]
        f = Pioneer.LUTFilter(table)
        n_pass = Pioneer.filter_counter!(c, f)
        @test n_pass == 2  # patterns 7 (3 bits) and 255 (8 bits) pass; 3 (2 bits) fails
        passing_ids = Set([c.ids[1], c.ids[2]])
        @test passing_ids == Set([UInt16(1), UInt16(3)])
    end

    # ── Test 13: filter_counter! returns 0 when nothing passes ──────────
    @testset "filter_counter! nothing passes" begin
        c = Pioneer.Counter(UInt16, UInt8, 100)
        Pioneer.or!(c, UInt16(1), UInt8(0x01))  # 1 bit
        Pioneer.or!(c, UInt16(2), UInt8(0x02))  # 1 bit

        f = Pioneer.CountFilter(UInt8(5))  # need 5 bits
        @test Pioneer.filter_counter!(c, f) == 0
    end

    # ── Test 14: filter_counter! + reset round-trip ──────────────────────
    @testset "filter_counter! + reset round-trip" begin
        c = Pioneer.Counter(UInt16, UInt8, 100)
        Pioneer.or!(c, UInt16(5), UInt8(0x1F))  # 5 bits
        Pioneer.or!(c, UInt16(10), UInt8(0x01)) # 1 bit

        f = Pioneer.CountFilter(UInt8(3))
        n_pass = Pioneer.filter_counter!(c, f)
        @test n_pass == 1
        @test c.ids[1] == UInt16(5)
        @test c.counts[5] == UInt8(0x1F)  # still readable

        # Reset clears everything
        Pioneer.reset!(c)
        @test c.size == 1
        @test c.counts[5] == 0x00
    end

    # ── Test 15: passes_filter dispatch ──────────────────────────────────
    @testset "passes_filter dispatch" begin
        cf = Pioneer.CountFilter(UInt8(3))
        @test Pioneer.passes_filter(cf, UInt8(0x07)) == true   # 3 bits
        @test Pioneer.passes_filter(cf, UInt8(0x03)) == false  # 2 bits
        @test Pioneer.passes_filter(cf, UInt8(0xFF)) == true   # 8 bits

        table = fill(false, 256)
        table[0xFF + 1] = true  # only all-8 passes
        lf = Pioneer.LUTFilter(table)
        @test Pioneer.passes_filter(lf, UInt8(0xFF)) == true
        @test Pioneer.passes_filter(lf, UInt8(0x07)) == false
        @test Pioneer.passes_filter(lf, UInt8(0x00)) == false
    end

    # ── Test 16: regression — reset after filter must clear ALL entries ──
    @testset "reset clears stale entries after filter" begin
        c = Pioneer.Counter(UInt16, UInt8, 20)

        # Scan 1: IDs 5, 3, 8 — only ID 3 passes (8 bits)
        Pioneer.or!(c, UInt16(5), UInt8(0x01))   # 1 bit
        Pioneer.or!(c, UInt16(3), UInt8(0xFF))   # 8 bits
        Pioneer.or!(c, UInt16(8), UInt8(0x03))   # 2 bits
        @test c.size == 4

        # Simulate what the hot loop does: check filter inline, then reset
        f = Pioneer.CountFilter(UInt8(3))
        for i in 1:(c.size - 1)
            lid = c.ids[i]
            Pioneer.passes_filter(f, c.counts[lid]) || continue
            # emit passing ID (just verify it's correct)
            @test lid == UInt16(3)
        end

        # Reset — must clear ALL entries including non-passing ones
        Pioneer.reset!(c)

        # Scan 2: only ID 7 gets scored
        Pioneer.or!(c, UInt16(7), UInt8(0x0F))  # 4 bits
        @test c.size == 2  # just ID 7

        # CRITICAL: stale counts from scan 1 must be zero
        @test c.counts[5] == 0x00
        @test c.counts[3] == 0x00
        @test c.counts[8] == 0x00
        @test c.counts[7] == UInt8(0x0F)  # current scan's entry
    end

    # ── Test 17: LUTFilter 256-bit bitset ───────────────────────────────
    @testset "LUTFilter bitset construction" begin
        # From Vector{Bool}
        v = Bool[count_ones(UInt8(i)) >= 3 for i in 0:255]
        lf = Pioneer.LUTFilter(v)
        @test lf.bits isa NTuple{4, UInt64}
        @test Pioneer.passes_filter(lf, UInt8(0x07)) == true   # 3 bits
        @test Pioneer.passes_filter(lf, UInt8(0x03)) == false  # 2 bits
        @test Pioneer.passes_filter(lf, UInt8(0xFF)) == true   # 8 bits
        @test Pioneer.passes_filter(lf, UInt8(0x00)) == false  # 0 bits

        # Exhaustive check: every pattern matches the original Bool vector
        for i in 0:255
            @test Pioneer.passes_filter(lf, UInt8(i)) == v[i + 1]
        end

        # Wrong length should error
        @test_throws ErrorException Pioneer.LUTFilter(Bool[true, false])
    end

    # ── Test 18: PatternAccumulator ─────────────────────────────────────
    @testset "PatternAccumulator" begin
        is_decoy = Bool[false, true, false, true, false]  # 5 precursors
        acc = Pioneer.PatternAccumulator(2, is_decoy; min_score=UInt8(2))
        @test length(acc.thread_tc) == 2
        @test length(acc.thread_dc) == 2
        @test acc.min_score == UInt8(2)

        # Accumulate: target prec 1 with score 0x07, decoy prec 2 with score 0x03
        Pioneer.accumulate!(acc, 1, UInt32(1), UInt8(0x07))
        Pioneer.accumulate!(acc, 1, UInt32(2), UInt8(0x03))
        Pioneer.accumulate!(acc, 2, UInt32(3), UInt8(0xFF))
        Pioneer.accumulate!(acc, 2, UInt32(4), UInt8(0x0F))

        # Thread 1: prec 1 (target, score 7+1=8) → tc[8], prec 2 (decoy, score 3+1=4) → dc[4]
        @test acc.thread_tc[1][Int(0x07) + 1] == 1
        @test acc.thread_dc[1][Int(0x03) + 1] == 1
        # Thread 2: prec 3 (target, score 0xFF+1=256) → tc[256], prec 4 (decoy, score 0x0F+1=16) → dc[16]
        @test acc.thread_tc[2][Int(0xFF) + 1] == 1
        @test acc.thread_dc[2][Int(0x0F) + 1] == 1

        # Merge
        tc, dc = Pioneer.merge_accumulator(acc)
        @test length(tc) == 256
        @test tc[Int(0x07) + 1] == 1
        @test tc[Int(0xFF) + 1] == 1
        @test dc[Int(0x03) + 1] == 1
        @test dc[Int(0x0F) + 1] == 1
        @test sum(tc) == 2  # 2 target hits
        @test sum(dc) == 2  # 2 decoy hits
    end

    # ── Test 19: or! performance parity with inc! ────────────────────────
    @testset "or! performance parity" begin
        # n_ids=1000, n_ops=100_000 means each id sees 100 increments. With
        # `inc!`'s `score ∈ 1:8`, a UInt8 counter would wrap to 0 after
        # ~32 iterations of the same id (e.g. 8 × 32 = 256 ≡ 0 mod 256).
        # `inc!` interprets `counts[id] == 0` as "first encounter" and
        # ghost-increments `c.size`, eventually pushing it past the array
        # bound. Production never hits this because real searches don't
        # reuse a Counter across that many additions per id, but the test
        # does — fail surfaces only under `Pkg.test`'s
        # `--check-bounds=yes` (CI), where `@inbounds` in inc! is ignored.
        # Using UInt16 counts widens the wrap horizon to 65 535,
        # comfortably above 100 × 8.
        n_ids = 1000
        n_ops = 100_000
        c_or = Pioneer.Counter(UInt16, UInt16, n_ids + 1)
        c_inc = Pioneer.Counter(UInt16, UInt16, n_ids + 1)

        # Warmup
        for _ in 1:1000
            Pioneer.or!(c_or, UInt16(1), UInt16(0x01))
            Pioneer.inc!(c_inc, UInt16(1), UInt16(0x01))
        end
        Pioneer.reset!(c_or); Pioneer.reset!(c_inc)

        # Benchmark or!
        t_or = @elapsed begin
            for i in 1:n_ops
                id = UInt16(mod(i, n_ids) + 1)
                mask = UInt16(1) << UInt16(mod(i, 8))
                Pioneer.or!(c_or, id, mask)
            end
        end

        Pioneer.reset!(c_or)

        # Benchmark inc!
        t_inc = @elapsed begin
            for i in 1:n_ops
                id = UInt16(mod(i, n_ids) + 1)
                score = UInt16(mod(i, 8) + 1)
                Pioneer.inc!(c_inc, id, score)
            end
        end

        # or! should be within 2x of inc! (same algorithm, just |= vs +=)
        @test t_or < t_inc * 2.0
        # Print for visibility (not a hard assertion)
        @info "Counter performance: or!=$(round(t_or*1e6/n_ops, digits=2))ns/op, inc!=$(round(t_inc*1e6/n_ops, digits=2))ns/op, ratio=$(round(t_or/t_inc, digits=2))"
    end
end
