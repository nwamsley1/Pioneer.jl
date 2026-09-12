# Run standalone with: julia --project=. test/UnitTests/test_isotope_interval_reuse.jl
module IsotopeIntervalReuseTests

using Test
using Pioneer
using StaticArrays: SVector

# Above its lower clamp, each synthetic spline has the known linear value
# P(S, i, mass) = (S + 1) * (i + 1) * mass / 128.
# Sulfur rows use different grids so fragment/complement positions must differ.
function synthetic_rows(count=10)
    rows = Vector{Pioneer.CubicSpline{40,Float32}}[]
    for sulfur in 0:5
        row = Pioneer.CubicSpline{40,Float32}[]
        first = 128f0 * (sulfur + 1)
        for isotope_idx in 0:count - 1
            scale = Float32((sulfur + 1) * (isotope_idx + 1))
            coefficients = Float32[]
            for bin in 0:9
                append!(coefficients, (scale * (sulfur + 1 + bin), scale / 128f0, 0f0, 0f0))
            end
            push!(row, Pioneer.CubicSpline(SVector{40,Float32}(coefficients),
                first, first + 9f0 * 128f0, 128f0, 1f0 / 128f0))
        end
        push!(rows, row)
    end
    return rows
end

function check_abundance(model, transmission, frag, prec, requested, expected)
    output = fill(-77f0, length(expected))
    @test Pioneer.getFragAbundance!(output, transmission, model, frag, prec, requested) === nothing
    @test isequal(output, expected)
end

@testset "isotope spline interval sharing" begin
    xml = joinpath(@__DIR__, "..", "..", "assets", "IsotopeSplines_10kDa_10isotopes.xml")
    loaded = Pioneer.parseIsoXML(xml)
    rows = synthetic_rows()
    model = Pioneer.IsotopeSplineModel(rows)

    @testset "M0–9 models require a common grid in each sulfur row" begin
        @test length(loaded.splines) == 6
        @test all(row -> length(row) == 10, loaded.splines)
        @test model(2, 3, 1024f0) === 96f0
        @test Pioneer.IsotopeSplineModel{Float32}(rows)(2, 3, 1024f0) === 96f0
        for field in (:first, :last, :bin_width, :inv_bin_width)
            incompatible = deepcopy(rows)
            s = incompatible[2][10]
            grid = (first=s.first, last=s.last, bin_width=s.bin_width, inv_bin_width=s.inv_bin_width)
            changed = merge(grid, NamedTuple{(field,)}((getfield(grid, field) + 1f0,)))
            incompatible[2][10] = Pioneer.CubicSpline(s.coeffs,
                changed.first, changed.last, changed.bin_width, changed.inv_bin_width)
            @test_throws ArgumentError Pioneer.IsotopeSplineModel(incompatible)
        end
        for count in (11, 21)
            @test_throws ArgumentError Pioneer.IsotopeSplineModel(synthetic_rows(count))
        end
        @test_throws BoundsError loaded(0, 10, 800f0)
    end

    @testset "known cubic values, lower clamp, and knot neighbors" begin
        cubic = Pioneer.CubicSpline(
            SVector{40,Float32}(repeat(Float32[1, 2, 3, 4], 10)),
            128f0, 1280f0, 128f0, 1f0 / 128f0)
        @test cubic(127f0) === 1f0
        @test cubic(128f0) === 1f0
        @test cubic(129f0) === 10f0
        @test cubic(130f0) === 49f0

        # Deliberately different constants across intervals make a wrong knot
        # selection visible even at the exact boundary.
        coefficients = Float32[]
        for bin in 0:9
            append!(coefficients, (Float32(10bin), 1f0, 0f0, 0f0))
        end
        piecewise = Pioneer.CubicSpline(SVector{40,Float32}(coefficients),
            128f0, 1280f0, 128f0, 1f0 / 128f0)
        for (mass, expected) in (
                (prevfloat(128f0), 0f0), (128f0, 0f0),
                (nextfloat(128f0), nextfloat(128f0) - 128f0),
                (prevfloat(256f0), prevfloat(256f0) - 128f0),
                (256f0, 10f0), (nextfloat(256f0), 10f0 + (nextfloat(256f0) - 256f0)),
                (320f0, 74f0), (1280f0, 90f0))
            @test piecewise(mass) === expected
        end
    end

    @testset "all supported lengths retain higher precursor contributions" begin
        sulfur_pairs = ((0, 0), (0, 3), (2, 0), (2, 2), (5, 5), (7, 6))
        for count in 1:10
            sf, sc = sulfur_pairs[mod1(count, length(sulfur_pairs))]
            frag = Pioneer.isotope(1024f0, Int64(sf), Int64(0))
            prec = Pioneer.isotope(2176f0, Int64(sf + sc), Int64(0))
            scale = 72 * (min(sf, 5) + 1) * (min(sc, 5) + 1)
            # Unit transmission gives the triangular sum 1+...+(count-f).
            expected = Float32[scale * (f + 1) * (count - f) * (count - f + 1) / 2
                for f in 0:count - 1]
            check_abundance(model, ones(Float32, count), frag, prec, 0:count - 1, expected)
            high_only = zeros(Float32, count)
            high_only[end] = 1f0
            requested = 0:min(1, count - 1)
            expected = Float32[scale * (f + 1) * (count - f) for f in requested]
            check_abundance(model, high_only, frag, prec, requested, expected)
        end
    end

    @testset "requested ranges preserve unwritten values" begin
        frag = Pioneer.isotope(128f0, Int64(0), Int64(0))
        prec = Pioneer.isotope(384f0, Int64(0), Int64(0))
        transmission = Float32[1, 0.5, 0, 0.25]
        # Analytic values are [6, 5, 3, 2] for this transmission and model.
        for (requested, expected) in ((0:3, Float32[6, 5, 3, 2, -77]),
                (1:2, Float32[-77, 5, 3, -77, -77]),
                (-2:1, Float32[6, 5, -77, -77, -77]),
                (3:6, Float32[-77, -77, -77, 2, -77]))
            check_abundance(model, transmission, frag, prec, requested, expected)
        end
        @test transmission == Float32[1, 0.5, 0, 0.25]
        output = Float32[-0.0, -77, -77, -77]
        Pioneer.getFragAbundance!(output, zeros(Float32, 4), model, frag, prec, 1:3)
        @test reinterpret(UInt32, output) == UInt32[0x80000000, 0, 0, 0]
        empty_model = Pioneer.IsotopeSplineModel(Vector{Pioneer.CubicSpline{40,Float32}}[])
        for (count, requested) in ((0, 0:-1), (5, 8:9), (11, 1:0), (11, 12:14))
            check_abundance(empty_model, ones(Float32, count), frag, prec, requested, Float32[-77])
        end
    end

    @testset "unsupported or missing isotopes fail before output writes" begin
        frag = Pioneer.isotope(1024f0, Int64(0), Int64(0))
        prec = Pioneer.isotope(2176f0, Int64(2), Int64(0))
        short = Pioneer.IsotopeSplineModel(synthetic_rows(5))
        check_abundance(short, ones(Float32, 5), frag, prec, 0:4, Float32[3240, 4320, 3888, 2592, 1080])
        for (tested_model, count, requested, error_type) in (
                (model, 11, 0:1, ArgumentError), (model, 21, 0:20, ArgumentError),
                (short, 6, 0:0, BoundsError), (short, 6, 5:5, BoundsError))
            output = fill(-77f0, count)
            @test_throws error_type Pioneer.getFragAbundance!(
                output, ones(Float32, count), tested_model, frag, prec, requested)
            @test all(==(-77f0), output)
        end
        for output_length in (5, 11)
            output = fill(-77f0, output_length)
            @test_throws BoundsError Pioneer.getFragAbundance!(output, model, frag, prec, (0, 10))
            @test all(==(-77f0), output)
        end
        @test_throws BoundsError Pioneer._precursor_fraction_transmitted(
            ones(Float32, 11), model, (1, 11), 800f0, UInt8(3), UInt8(2))
    end
end

end # module
