# Unit tests for filter_low_quality_precursors!
#
# Target: Pioneer.filter_low_quality_precursors! (SecondPassSearch/utils.jl:206-244)
# Two filters: (1) dynamic range — zero weights below threshold,
# (2) fragment count — require ≥ min_frag_count monoisotopic OR M+1 matched fragments.
#
# Run: julia --project=. test/UnitTests/test_filter_precursors.jl

using Test
using Pioneer

# ---------------------------------------------------------------------------
# Aliases
# ---------------------------------------------------------------------------
const SA                = Pioneer.SparseArray
const sortSparse!       = Pioneer.sortSparse!
const filter_low_quality_precursors! = Pioneer.filter_low_quality_precursors!

# ---------------------------------------------------------------------------
# Helper: build SparseArray from (col, row, nzval, x, matched, isotope) tuples
# ---------------------------------------------------------------------------
function build_test_sparse(entries)
    n = length(entries)
    sa = SA(UInt32(max(n + 10, 50)))
    sa.n_vals = n
    for (i, (col, row, nz, xv, m, iso)) in enumerate(entries)
        sa.colval[i] = UInt16(col)
        sa.rowval[i] = UInt32(row)
        sa.nzval[i]  = Float32(nz)
        sa.x[i]      = Float32(xv)
        sa.matched[i] = m
        sa.isotope[i] = UInt8(iso)
    end
    sortSparse!(sa)
    for i in 1:sa.n_vals
        sa.m = max(sa.m, Int64(sa.rowval[i]))
    end
    return sa
end

# ═══════════════════════════════════════════════════════════════════════════
# Tests
# ═══════════════════════════════════════════════════════════════════════════

@testset "filter_low_quality_precursors!" begin

    # -------------------------------------------------------------------
    @testset "2.1 Dynamic range filter" begin
        # weights = [100, 1, 0.05], dynamic_range = 0.01
        # threshold = 100 * 0.01 = 1.0
        # Prec 3 (0.05 < 1.0) → zeroed
        # Give all precursors ≥3 mono fragments so only range matters
        entries = [
            # Prec 1: 3 mono matched
            (1, 1, 1.0, 10.0, true, 0x00),
            (1, 2, 1.0, 10.0, true, 0x00),
            (1, 3, 1.0, 10.0, true, 0x00),
            # Prec 2: 3 mono matched
            (2, 4, 1.0, 10.0, true, 0x00),
            (2, 5, 1.0, 10.0, true, 0x00),
            (2, 6, 1.0, 10.0, true, 0x00),
            # Prec 3: 3 mono matched
            (3, 7, 1.0, 10.0, true, 0x00),
            (3, 8, 1.0, 10.0, true, 0x00),
            (3, 9, 1.0, 10.0, true, 0x00),
        ]
        H = build_test_sparse(entries)
        weights = Float32[100.0, 1.0, 0.05]

        filter_low_quality_precursors!(weights, H, 3; dynamic_range=0.01f0)

        @test weights[1] == 100.0f0
        @test weights[2] == 1.0f0
        @test weights[3] == 0.0f0
    end

    # -------------------------------------------------------------------
    @testset "2.2 Fragment count filter (mono)" begin
        # Prec 1: 3 mono matched → passes (min_frag_count=3)
        # Prec 2: 2 mono + 1 M+2 → fails (mono=2<3, m1=0<3)
        entries = [
            # Prec 1: 3 monoisotopic
            (1, 1, 1.0, 10.0, true, 0x00),
            (1, 2, 1.0, 10.0, true, 0x00),
            (1, 3, 1.0, 10.0, true, 0x00),
            # Prec 2: 2 mono + 1 M+2
            (2, 4, 1.0, 10.0, true, 0x00),
            (2, 5, 1.0, 10.0, true, 0x00),
            (2, 6, 1.0, 10.0, true, 0x02),  # M+2 isotope
        ]
        H = build_test_sparse(entries)
        weights = Float32[50.0, 50.0]  # both well above dynamic range

        filter_low_quality_precursors!(weights, H, 3; dynamic_range=Float32(1e-3))

        @test weights[1] == 50.0f0   # passes
        @test weights[2] == 0.0f0    # zeroed: mono=2<3, m1=0<3
    end

    # -------------------------------------------------------------------
    @testset "2.3 M+1 isotopes save a precursor" begin
        # 1 mono + 4 M+1 matched, min_frag_count=3 → m1≥3 → passes
        entries = [
            (1, 1, 1.0, 10.0, true, 0x00),  # mono
            (1, 2, 1.0, 10.0, true, 0x01),  # M+1
            (1, 3, 1.0, 10.0, true, 0x01),  # M+1
            (1, 4, 1.0, 10.0, true, 0x01),  # M+1
            (1, 5, 1.0, 10.0, true, 0x01),  # M+1
        ]
        H = build_test_sparse(entries)
        weights = Float32[10.0]

        filter_low_quality_precursors!(weights, H, 3; dynamic_range=Float32(1e-3))

        @test weights[1] == 10.0f0   # passes: m1=4≥3
    end

    # -------------------------------------------------------------------
    @testset "2.4 All zero weights — no crash" begin
        entries = [
            (1, 1, 1.0, 10.0, true, 0x00),
            (1, 2, 1.0, 10.0, true, 0x00),
            (1, 3, 1.0, 10.0, true, 0x00),
            (2, 4, 1.0, 10.0, true, 0x00),
            (2, 5, 1.0, 10.0, true, 0x00),
            (2, 6, 1.0, 10.0, true, 0x00),
        ]
        H = build_test_sparse(entries)
        weights = Float32[0.0, 0.0]

        # max_weight=0, threshold=0 → 0 < 0 is false → no zeroing from range filter
        # But fragment counts pass, so weights remain 0
        filter_low_quality_precursors!(weights, H, 3; dynamic_range=Float32(1e-3))

        @test weights == Float32[0.0, 0.0]
    end

    # -------------------------------------------------------------------
    @testset "2.5 Single precursor passes both filters" begin
        entries = [
            (1, 1, 1.0, 10.0, true, 0x00),
            (1, 2, 1.0, 10.0, true, 0x00),
            (1, 3, 1.0, 10.0, true, 0x00),
        ]
        H = build_test_sparse(entries)
        weights = Float32[50.0]

        filter_low_quality_precursors!(weights, H, 3; dynamic_range=Float32(1e-3))

        @test weights[1] == 50.0f0
    end

    # -------------------------------------------------------------------
    @testset "2.6 Boundary — weight exactly at threshold NOT zeroed" begin
        # weights = [100, 1], dynamic_range=0.01, threshold=1.0
        # weight=1 is NOT < 1.0 (strict <) → NOT zeroed
        entries = [
            (1, 1, 1.0, 10.0, true, 0x00),
            (1, 2, 1.0, 10.0, true, 0x00),
            (1, 3, 1.0, 10.0, true, 0x00),
            (2, 4, 1.0, 10.0, true, 0x00),
            (2, 5, 1.0, 10.0, true, 0x00),
            (2, 6, 1.0, 10.0, true, 0x00),
        ]
        H = build_test_sparse(entries)
        weights = Float32[100.0, 1.0]

        filter_low_quality_precursors!(weights, H, 3; dynamic_range=0.01f0)

        @test weights[1] == 100.0f0
        @test weights[2] == 1.0f0   # exactly at threshold → not zeroed
    end

    # -------------------------------------------------------------------
    @testset "2.7 Unmatched entries → zero fragment counts → zeroed" begin
        # 1 precursor with 5 entries, all matched=false
        # mono_count=0, m1_count=0 → both < min_frag_count → zeroed
        entries = [
            (1, 1, 1.0, 0.0, false, 0x00),
            (1, 2, 1.0, 0.0, false, 0x00),
            (1, 3, 1.0, 0.0, false, 0x01),
            (1, 4, 1.0, 0.0, false, 0x01),
            (1, 5, 1.0, 0.0, false, 0x00),
        ]
        H = build_test_sparse(entries)
        weights = Float32[50.0]

        filter_low_quality_precursors!(weights, H, 3; dynamic_range=Float32(1e-3))

        @test weights[1] == 0.0f0
    end

end  # filter_low_quality_precursors! testset

println("\nAll filter_low_quality_precursors! tests passed!")
