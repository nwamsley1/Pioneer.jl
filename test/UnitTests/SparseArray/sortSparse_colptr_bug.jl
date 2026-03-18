# sortSparse! colptr off-by-one bug
#
# Bug: sortSparse! set colptr[max_col+1] = n_vals instead of n_vals+1.
# Downstream consumers (initResiduals!, multiply!) iterate columns as:
#   start = colptr[col], stop = colptr[col+1] - 1
# So the last column's range was colptr[last]:n_vals-1, DROPPING the last entry.
#
# This test constructs minimal SparseArrays, calls sortSparse!, and verifies
# that every entry is reachable through the colptr ranges.
#
# Run: julia --project=. test/UnitTests/SparseArray/sortSparse_colptr_bug.jl

using Test
using Pioneer

const SA = Pioneer.SparseArray
const sortSparse! = Pioneer.sortSparse!

"""
Helper: iterate all columns using colptr (the way initResiduals!/multiply! do)
and collect every (col, row, nzval) triple that is reachable.
"""
function collect_via_colptr(sa)
    entries = Tuple{UInt16, eltype(sa.rowval), eltype(sa.nzval)}[]
    for col in 1:sa.n
        start = sa.colptr[col]
        stop  = sa.colptr[col+1] - 1   # ← this is how initResiduals!/multiply! read it
        for idx in start:stop
            push!(entries, (sa.colval[idx], sa.rowval[idx], sa.nzval[idx]))
        end
    end
    return entries
end

"""
Helper: fill a SparseArray with explicit (col, row, nzval) data.
"""
function fill_sa!(sa, data::Vector{Tuple{UInt16, UInt32, Float32}})
    sa.n_vals = length(data)
    for (i, (c, r, v)) in enumerate(data)
        sa.colval[i] = c
        sa.rowval[i] = r
        sa.nzval[i] = v
        sa.x[i] = Float32(r)  # use row as observed intensity for easy verification
        sa.matched[i] = true
        sa.isotope[i] = UInt8(0)
    end
end

@testset "sortSparse! colptr off-by-one" begin

    # =========================================================================
    @testset "Test 1: Single column — last entry must not be dropped" begin
        # All 3 entries belong to column 1. After sort, colptr should give
        # access to all 3 entries.
        sa = SA(UInt32(20))
        fill_sa!(sa, [
            (UInt16(1), UInt32(1), 1.0f0),
            (UInt16(1), UInt32(2), 2.0f0),
            (UInt16(1), UInt32(3), 3.0f0),  # ← this is the LAST entry
        ])
        sortSparse!(sa)

        @test sa.n == 1          # 1 column
        @test sa.n_vals == 3     # 3 entries total

        entries = collect_via_colptr(sa)
        @test length(entries) == 3  # All 3 must be reachable

        # Verify the actual nzval values are all present
        nzvals_found = sort([e[3] for e in entries])
        @test nzvals_found == [1.0f0, 2.0f0, 3.0f0]
    end

    # =========================================================================
    @testset "Test 2: Two columns — last entry of last column" begin
        # Column 1 has 2 entries, column 2 has 1 entry.
        # After sort by column, column 2's single entry is at position 3 (the last).
        sa = SA(UInt32(20))
        fill_sa!(sa, [
            (UInt16(2), UInt32(1), 10.0f0),  # col 2 — will be sorted to position 3
            (UInt16(1), UInt32(1), 1.0f0),   # col 1
            (UInt16(1), UInt32(2), 2.0f0),   # col 1
        ])
        sortSparse!(sa)

        @test sa.n == 2
        entries = collect_via_colptr(sa)
        @test length(entries) == 3  # All 3 must be reachable

        # Column 1 should have nzvals {1, 2}, column 2 should have {10}
        col1 = [e[3] for e in entries if e[1] == UInt16(1)]
        col2 = [e[3] for e in entries if e[1] == UInt16(2)]
        @test sort(col1) == [1.0f0, 2.0f0]
        @test col2 == [10.0f0]
    end

    # =========================================================================
    @testset "Test 3: Three columns — each column's entries are complete" begin
        # Simulates a realistic case: 3 precursors, 7 entries total.
        # Entries are given in scrambled column order.
        sa = SA(UInt32(20))
        fill_sa!(sa, [
            (UInt16(3), UInt32(5), 30.0f0),
            (UInt16(1), UInt32(1), 1.0f0),
            (UInt16(2), UInt32(2), 20.0f0),
            (UInt16(1), UInt32(3), 2.0f0),
            (UInt16(3), UInt32(4), 31.0f0),
            (UInt16(2), UInt32(1), 21.0f0),
            (UInt16(3), UInt32(2), 32.0f0),  # ← last entry, belongs to last column
        ])
        sortSparse!(sa)

        @test sa.n == 3
        entries = collect_via_colptr(sa)
        @test length(entries) == 7  # All 7 must be reachable

        col1 = sort([e[3] for e in entries if e[1] == UInt16(1)])
        col2 = sort([e[3] for e in entries if e[1] == UInt16(2)])
        col3 = sort([e[3] for e in entries if e[1] == UInt16(3)])
        @test col1 == [1.0f0, 2.0f0]
        @test col2 == [20.0f0, 21.0f0]
        @test col3 == [30.0f0, 31.0f0, 32.0f0]
    end

    # =========================================================================
    @testset "Test 4: initResiduals! uses all entries (not n_vals-1)" begin
        # The actual downstream impact: initResiduals! computes r = Hw - x.
        # If the last entry is dropped, the residual at that row is wrong.
        #
        # Setup: 2 precursors, 4 entries.
        #   col 1: rows 1,2 with nzval 1.0 each
        #   col 2: rows 3,4 with nzval 1.0 each
        #   x (observed): [10, 20, 30, 40] at rows 1-4
        #   weights: w = [2.0, 3.0]
        #
        # Expected residuals r = Hw - x:
        #   r[1] = 2.0*1.0 - 10 = -8.0
        #   r[2] = 2.0*1.0 - 20 = -18.0
        #   r[3] = 3.0*1.0 - 30 = -27.0
        #   r[4] = 3.0*1.0 - 40 = -37.0  ← this is the last entry of the last column
        #
        # With the bug, r[4] = 0 - 40 = -40.0 (the nzval contribution is missed)

        sa = SA(UInt32(20))
        sa.n_vals = 4
        # Scramble: put col 2's entries first so they end up last after sort
        sa.colval[1] = UInt16(2); sa.rowval[1] = UInt32(3); sa.nzval[1] = 1.0f0; sa.x[1] = 30.0f0
        sa.colval[2] = UInt16(2); sa.rowval[2] = UInt32(4); sa.nzval[2] = 1.0f0; sa.x[2] = 40.0f0
        sa.colval[3] = UInt16(1); sa.rowval[3] = UInt32(1); sa.nzval[3] = 1.0f0; sa.x[3] = 10.0f0
        sa.colval[4] = UInt16(1); sa.rowval[4] = UInt32(2); sa.nzval[4] = 1.0f0; sa.x[4] = 20.0f0
        for i in 1:4; sa.matched[i] = true; sa.isotope[i] = UInt8(0); end

        sortSparse!(sa)

        w = Float32[2.0, 3.0]
        r = zeros(Float32, 4)
        Pioneer.initResiduals!(r, sa, w)

        # Expected: r = Hw - x
        @test r[1] ≈ -8.0f0   atol=1e-5   # col 1, row 1: 2*1 - 10
        @test r[2] ≈ -18.0f0  atol=1e-5   # col 1, row 2: 2*1 - 20
        @test r[3] ≈ -27.0f0  atol=1e-5   # col 2, row 3: 3*1 - 30
        @test r[4] ≈ -37.0f0  atol=1e-5   # col 2, row 4: 3*1 - 40  ← FAILS with bug (-40 instead)
    end

    # =========================================================================
    @testset "Test 5: colptr sentinel value is n_vals+1" begin
        # Direct check: colptr[n_cols+1] should equal n_vals+1 so that
        # colptr[last_col]:(colptr[last_col+1]-1) includes the last entry.
        sa = SA(UInt32(20))
        fill_sa!(sa, [
            (UInt16(1), UInt32(1), 1.0f0),
            (UInt16(2), UInt32(1), 2.0f0),
            (UInt16(3), UInt32(1), 3.0f0),
        ])
        sortSparse!(sa)

        @test sa.n == 3
        @test sa.colptr[sa.n + 1] == sa.n_vals + 1

        # Each column has exactly 1 entry
        for col in 1:sa.n
            col_size = sa.colptr[col+1] - sa.colptr[col]
            @test col_size == 1
        end
    end

end
