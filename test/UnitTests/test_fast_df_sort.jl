# Unit tests for fast_df_sort! (src/utils/sortUtils.jl).
#
# Tests the fast DataFrame sorting utility that uses zip+sortperm instead of
# DataFrames' DFPerm wrapper. Covers single/multi-column sorting, reverse
# ordering with the negation tricks (signed, unsigned, float, bool), and
# the Vector{Symbol} convenience overload.
#
# Run standalone: julia --project=. test/UnitTests/test_fast_df_sort.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
    using DataFrames
end

@testset "fast_df_sort!" begin

# ═══════════════════════════════════════════════════════════════════════════════
# Single column sorting
# ═══════════════════════════════════════════════════════════════════════════════

@testset "single column ascending" begin
    df = DataFrame(x = [3, 1, 2], y = ["c", "a", "b"])
    Pioneer.fast_df_sort!(df, (:x,))
    @test df.x == [1, 2, 3]
    @test df.y == ["a", "b", "c"]
end

@testset "single column descending Int" begin
    df = DataFrame(x = [1, 3, 2], y = ["a", "c", "b"])
    Pioneer.fast_df_sort!(df, (:x,), rev=(true,))
    @test df.x == [3, 2, 1]
    @test df.y == ["c", "b", "a"]
end

@testset "single column descending Float32" begin
    df = DataFrame(x = Float32[1.5, 3.5, 2.5])
    Pioneer.fast_df_sort!(df, (:x,), rev=(true,))
    @test df.x == Float32[3.5, 2.5, 1.5]
end

@testset "single column descending UInt32" begin
    df = DataFrame(x = UInt32[10, 30, 20])
    Pioneer.fast_df_sort!(df, (:x,), rev=(true,))
    @test df.x == UInt32[30, 20, 10]
end

@testset "single column descending Bool" begin
    df = DataFrame(x = [false, true, false, true], y = [1, 2, 3, 4])
    Pioneer.fast_df_sort!(df, (:x,), rev=(true,))
    # true sorts first when reversed
    @test df.x[1] == true
    @test df.x[2] == true
    @test df.x[3] == false
    @test df.x[4] == false
end

# ═══════════════════════════════════════════════════════════════════════════════
# Multi-column sorting
# ═══════════════════════════════════════════════════════════════════════════════

@testset "two columns ascending" begin
    df = DataFrame(a = [2, 1, 1, 2], b = [2, 2, 1, 1])
    Pioneer.fast_df_sort!(df, (:a, :b))
    @test df.a == [1, 1, 2, 2]
    @test df.b == [1, 2, 1, 2]
end

@testset "two columns mixed rev" begin
    df = DataFrame(a = [1, 1, 2, 2], b = Float32[10.0, 20.0, 10.0, 20.0])
    Pioneer.fast_df_sort!(df, (:a, :b), rev=(false, true))
    @test df.a == [1, 1, 2, 2]
    @test df.b == Float32[20.0, 10.0, 20.0, 10.0]
end

@testset "three columns" begin
    df = DataFrame(
        a = UInt16[1, 1, 1, 2],
        b = [true, false, true, true],
        c = Float32[3.0, 1.0, 2.0, 4.0]
    )
    Pioneer.fast_df_sort!(df, (:a, :b, :c), rev=(false, true, false))
    # a ascending, b descending (true first), c ascending
    @test df.a == UInt16[1, 1, 1, 2]
    @test df.b == [true, true, false, true]
    @test df.c == Float32[2.0, 3.0, 1.0, 4.0]
end

# ═══════════════════════════════════════════════════════════════════════════════
# String column reverse — should fall back to sort!
# ═══════════════════════════════════════════════════════════════════════════════

@testset "string column reverse falls back" begin
    df = DataFrame(s = ["banana", "apple", "cherry"], x = [2, 1, 3])
    Pioneer.fast_df_sort!(df, (:s,), rev=(true,))
    @test df.s == ["cherry", "banana", "apple"]
    @test df.x == [3, 2, 1]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Vector{Symbol} convenience overload
# ═══════════════════════════════════════════════════════════════════════════════

@testset "Vector{Symbol} overload" begin
    df = DataFrame(a = [2, 1, 1, 2], b = [2, 2, 1, 1])
    Pioneer.fast_df_sort!(df, [:a, :b], rev=[false, true])
    @test df.a == [1, 1, 2, 2]
    @test df.b == [2, 1, 2, 1]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Edge cases
# ═══════════════════════════════════════════════════════════════════════════════

@testset "empty DataFrame" begin
    df = DataFrame(x = Int[])
    result = Pioneer.fast_df_sort!(df, (:x,))
    @test nrow(result) == 0
end

@testset "single row" begin
    df = DataFrame(x = [42], y = ["only"])
    Pioneer.fast_df_sort!(df, (:x,))
    @test df.x == [42]
    @test df.y == ["only"]
end

@testset "already sorted" begin
    df = DataFrame(x = [1, 2, 3])
    Pioneer.fast_df_sort!(df, (:x,))
    @test df.x == [1, 2, 3]
end

@testset "all same values" begin
    df = DataFrame(x = [5, 5, 5], y = [3, 1, 2])
    Pioneer.fast_df_sort!(df, (:x, :y))
    @test df.y == [1, 2, 3]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Negation helpers
# ═══════════════════════════════════════════════════════════════════════════════

@testset "negation correctness" begin
    # Unsigned: bitwise complement reverses ordering
    @test Pioneer._rev_val(UInt32(0)) > Pioneer._rev_val(UInt32(1))
    @test Pioneer._rev_val(UInt32(0)) == typemax(UInt32)

    # Signed: negation reverses ordering
    @test Pioneer._rev_val(Int64(5)) < Pioneer._rev_val(Int64(3))

    # Float: negation reverses ordering
    @test Pioneer._rev_val(1.0f0) < Pioneer._rev_val(0.5f0)

    # Bool: flip reverses ordering
    @test Pioneer._rev_val(true) == false
    @test Pioneer._rev_val(false) == true
end

@testset "_supports_rev" begin
    @test Pioneer._supports_rev(Int64) == true
    @test Pioneer._supports_rev(Float32) == true
    @test Pioneer._supports_rev(UInt32) == true
    @test Pioneer._supports_rev(Bool) == true
    @test Pioneer._supports_rev(String) == false
end

# ═══════════════════════════════════════════════════════════════════════════════
# Consistency with Base.sort!
# ═══════════════════════════════════════════════════════════════════════════════

@testset "matches Base.sort! result" begin
    n = 100
    df1 = DataFrame(
        a = rand(UInt16, n),
        b = rand(Float32, n),
        c = rand(Bool, n)
    )
    df2 = copy(df1)

    Pioneer.fast_df_sort!(df1, (:a, :b), rev=(true, false))
    sort!(df2, [:a, :b], rev=[true, false])

    @test df1.a == df2.a
    @test df1.b == df2.b
    @test df1.c == df2.c
end

end # top-level testset
