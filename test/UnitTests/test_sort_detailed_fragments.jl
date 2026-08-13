# Unit tests for sort_detailed_fragments_by_mz! (BuildSpecLib helper).
#
# Target: src/Routines/BuildSpecLib/build/build_poin_lib.jl
# This helper is called between partitioned-index construction and the
# detailed_fragments.jls write so the on-disk file satisfies run_fused!'s
# m/z-sorted-within-precursor pre-condition (verify_mz_sorted in fusedScan.jl).
#
# Run: julia --project=. test/UnitTests/test_sort_detailed_fragments.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
end

const sort_detailed_fragments_by_mz! = Pioneer.sort_detailed_fragments_by_mz!
const CompactFrag = Pioneer.CompactFrag

# Helper: build a CompactFrag with given (prec_id, mz, rank).
function _frag(pid::Integer, mz::Real, rank::Integer = 0)
    CompactFrag(
        UInt32(pid), Float32(mz), Float16(1000),
        true, false, false, false,                  # is_y, is_b, is_p, is_isotope
        UInt8(1), UInt8(1), UInt8(2), UInt8(rank), UInt8(0),
    )
end

@testset "sort_detailed_fragments_by_mz!" begin

    @testset "single precursor sorted in place" begin
        # Rank-ordered (intensity-descending) input → m/z-ascending output.
        frags = [_frag(1, 300.0, 0), _frag(1, 100.0, 1), _frag(1, 200.0, 2)]
        prec_ranges = UInt64[1, 4]
        n_mut = sort_detailed_fragments_by_mz!(frags, prec_ranges)
        @test n_mut == 1                          # this precursor was re-sorted
        @test [Float32(f.mz) for f in frags] == Float32[100.0, 200.0, 300.0]
    end

    @testset "already-sorted precursor not mutated" begin
        frags = [_frag(1, 100.0, 0), _frag(1, 200.0, 1), _frag(1, 300.0, 2)]
        prec_ranges = UInt64[1, 4]
        n_mut = sort_detailed_fragments_by_mz!(frags, prec_ranges)
        @test n_mut == 0                          # was already m/z-sorted
        @test [Float32(f.mz) for f in frags] == Float32[100.0, 200.0, 300.0]
    end

    @testset "multiple precursors sorted independently" begin
        # prec 1: [400, 100, 250]  →  [100, 250, 400]
        # prec 2: [600, 500]       →  [500, 600]
        # prec 3: [700]            →  unchanged
        frags = [_frag(1, 400.0), _frag(1, 100.0), _frag(1, 250.0),
                 _frag(2, 600.0), _frag(2, 500.0),
                 _frag(3, 700.0)]
        prec_ranges = UInt64[1, 4, 6, 7]
        n_mut = sort_detailed_fragments_by_mz!(frags, prec_ranges)
        @test n_mut == 2                          # prec 1 + prec 2 mutated; prec 3 single

        @test [Float32(f.mz) for f in frags[1:3]] == Float32[100.0, 250.0, 400.0]
        @test [Float32(f.mz) for f in frags[4:5]] == Float32[500.0, 600.0]
        @test Float32(frags[6].mz) === 700.0f0
    end

    @testset "empty precursor range (zero fragments) skipped" begin
        # Precursor 1 has fragment at index 1; precursor 2 has empty range.
        frags = [_frag(1, 100.0)]
        prec_ranges = UInt64[1, 2, 2]              # prec 2: [2, 2) empty
        n_mut = sort_detailed_fragments_by_mz!(frags, prec_ranges)
        @test n_mut == 0
    end

    @testset "preserves rank, ion type, charge, sulfur after sort" begin
        # A fragment's packed metadata (rank, charge, etc.) must NOT change
        # when sorted — only its position in the array does.
        frags = [
            CompactFrag(UInt32(1), 300.0f0, Float16(500), true, false, false, false,
                        UInt8(2), UInt8(3), UInt8(2), UInt8(7), UInt8(1)),
            CompactFrag(UInt32(1), 100.0f0, Float16(800), false, true, false, false,
                        UInt8(1), UInt8(5), UInt8(2), UInt8(2), UInt8(0)),
        ]
        prec_ranges = UInt64[1, 3]
        sort_detailed_fragments_by_mz!(frags, prec_ranges)

        # After sort: m/z-100 fragment first, m/z-300 second.
        @test Float32(frags[1].mz) === 100.0f0
        @test Pioneer.getRank(frags[1]) == UInt8(2)
        @test Pioneer.isB(frags[1]) == true
        @test Pioneer.isY(frags[1]) == false

        @test Float32(frags[2].mz) === 300.0f0
        @test Pioneer.getRank(frags[2]) == UInt8(7)
        @test Pioneer.isY(frags[2]) == true
    end

    @testset "stable: equal-mz fragments keep relative order" begin
        # Two fragments with identical m/z (different ranks). Sort should
        # not reorder them (InsertionSort is stable).
        frags = [_frag(1, 200.0, 1), _frag(1, 200.0, 2), _frag(1, 100.0, 3)]
        prec_ranges = UInt64[1, 4]
        sort_detailed_fragments_by_mz!(frags, prec_ranges)
        # m/z-100 came last in input but first by m/z order
        @test Float32(frags[1].mz) === 100.0f0
        @test Pioneer.getRank(frags[1]) == UInt8(3)
        # The two m/z-200 fragments retain their original (rank 1, then rank 2) order
        @test Float32(frags[2].mz) === 200.0f0
        @test Float32(frags[3].mz) === 200.0f0
        @test Pioneer.getRank(frags[2]) == UInt8(1)
        @test Pioneer.getRank(frags[3]) == UInt8(2)
    end

end
