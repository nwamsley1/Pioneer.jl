# Unit tests for normalizeQuant (src/utils/normalizeQuant.jl).
#
# Tests RT-dependent quantification normalization used by MaxLFQSearch.
# Creates synthetic Arrow files with known intensity offsets and verifies
# the normalization pipeline corrects them.
#
# Run standalone: julia --project=. test/UnitTests/test_normalizeQuant.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
    using DataFrames
    using Arrow
    using Statistics
    using Dictionaries
end

# Helper: create a synthetic PSM Arrow file with RT-dependent abundance.
# Files meant to get a spline need enough PSMs to fill
# `_min_bins_for_spline(spline_n_knots)` bins of `QUANT_MIN_BIN_OCCUPANCY` each:
# at spline_n_knots=5 that is 8 bins × 100 = 800 PSMs, hence n=1000 below.
function _write_quant_arrow(path::String, n::Int, offset::Float64;
                            rt_range::Tuple{Float64,Float64}=(0.0, 100.0))
    rts = Float32.(collect(LinRange(first(rt_range), last(rt_range), n)))
    # Base abundance: Gaussian peak centered at RT=50
    base = Float64.(exp.(-(rts .- 50.0f0) .^ 2 ./ 500.0f0) .* 1000.0f0 .+ 100.0f0)
    # Apply multiplicative offset (in log2 space)
    abundance = Float32.(base .* 2.0^offset)
    df = DataFrame(
        irt_obs = rts,
        abundance = abundance,
    )
    Arrow.write(path, df)
end

@testset "normalizeQuant" begin

# ═══════════════════════════════════════════════════════════════════════════════
# getQuantSplines — fit per-file splines
# ═══════════════════════════════════════════════════════════════════════════════

@testset "getQuantSplines" begin
    mktempdir() do dir
        # Two files with different intensity levels
        path1 = joinpath(dir, "file1.arrow")
        path2 = joinpath(dir, "file2.arrow")
        _write_quant_arrow(path1, 1000, 0.0)
        _write_quant_arrow(path2, 1000, 2.0)  # 4× higher

        splines, rt_range = Pioneer.getQuantSplines(
            [path1, path2], :abundance; N=50, spline_n_knots=5)

        @testset "returns spline per file" begin
            @test length(splines) == 2
            @test haskey(splines, path1)
            @test haskey(splines, path2)
        end

        @testset "rt_range spans data" begin
            @test rt_range[1] ≈ 0.0f0 atol=1.0
            @test rt_range[2] ≈ 100.0f0 atol=1.0
        end

        @testset "splines reflect offset" begin
            # At RT=50, file2 should be ~2 log2 units higher than file1
            v1 = splines[path1](50.0f0)
            v2 = splines[path2](50.0f0)
            @test v2 - v1 ≈ 2.0 atol=0.5
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# getQuantCorrections — compute offset corrections
# ═══════════════════════════════════════════════════════════════════════════════

@testset "getQuantCorrections" begin
    mktempdir() do dir
        path1 = joinpath(dir, "file1.arrow")
        path2 = joinpath(dir, "file2.arrow")
        _write_quant_arrow(path1, 1000, 0.0)
        _write_quant_arrow(path2, 1000, 2.0)

        splines, rt_range = Pioneer.getQuantSplines(
            [path1, path2], :abundance; N=50, spline_n_knots=5)

        corrections = Pioneer.getQuantCorrections(splines, rt_range; N=50)

        @testset "returns correction per file" begin
            @test length(corrections) == 2
            @test haskey(corrections, path1)
            @test haskey(corrections, path2)
        end

        @testset "corrections are opposite sign" begin
            # File1 is below median → negative correction
            # File2 is above median → positive correction
            c1 = corrections[path1](50.0)
            c2 = corrections[path2](50.0)
            @test c1 < 0.0 || c2 > 0.0
            # They should roughly cancel: c1 + c2 ≈ 0
            @test abs(c1 + c2) < 0.5
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Full pipeline: normalizeQuant end-to-end
# ═══════════════════════════════════════════════════════════════════════════════

@testset "normalizeQuant end-to-end" begin
    mktempdir() do dir
        # Three files with different systematic offsets
        paths = String[]
        offsets = [0.0, 1.5, -1.0]
        for (i, offset) in enumerate(offsets)
            path = joinpath(dir, "file$i.arrow")
            _write_quant_arrow(path, 1000, offset)
            push!(paths, path)
        end

        Pioneer.normalizeQuant(dir, :abundance; N=50, spline_n_knots=5)

        @testset "normalized column added" begin
            for path in paths
                df = DataFrame(Arrow.Table(path))
                @test hasproperty(df, :abundance_normalized)
            end
        end

        @testset "normalized abundances are more similar across files" begin
            # Read normalized abundances at comparable RT positions
            medians_raw = Float64[]
            medians_norm = Float64[]
            for path in paths
                df = DataFrame(Arrow.Table(path))
                mid = div(nrow(df), 2)
                rng = max(1, mid-10):min(nrow(df), mid+10)
                push!(medians_raw, median(df.abundance[rng]))
                push!(medians_norm, median(df.abundance_normalized[rng]))
            end
            # CV (coefficient of variation) should decrease after normalization
            cv_raw = std(medians_raw) / mean(medians_raw)
            cv_norm = std(medians_norm) / mean(medians_norm)
            @test cv_norm < cv_raw
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Edge case: single file — normalization should be near-identity
# ═══════════════════════════════════════════════════════════════════════════════

@testset "single file" begin
    mktempdir() do dir
        path = joinpath(dir, "only.arrow")
        _write_quant_arrow(path, 1000, 0.0)

        Pioneer.normalizeQuant(dir, :abundance; N=50, spline_n_knots=5)

        df = DataFrame(Arrow.Table(path))
        @test hasproperty(df, :abundance_normalized)
        # With one file, correction should be ~zero → normalized ≈ raw
        ratio = median(df.abundance_normalized ./ df.abundance)
        @test ratio ≈ 1.0 atol=0.1
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Edge case: too few PSMs to fit any spline — fall back to identity
# ═══════════════════════════════════════════════════════════════════════════════

@testset "no file supports a spline" begin
    mktempdir() do dir
        # 200 PSMs fills only 2 bins at >= 100 each, short of the 8 a
        # 7-coefficient spline needs, so every file is skipped.
        paths = String[]
        for (i, offset) in enumerate([0.0, 2.0])
            path = joinpath(dir, "file$i.arrow")
            _write_quant_arrow(path, 200, offset)
            push!(paths, path)
        end

        splines, rt_range = Pioneer.getQuantSplines(
            paths, :abundance; N=50, spline_n_knots=5)
        @test isempty(splines)

        # No cross-file median exists; corrections must be empty, not an error.
        corrections = Pioneer.getQuantCorrections(splines, rt_range; N=50)
        @test isempty(corrections)

        Pioneer.normalizeQuant(dir, :abundance; N=50, spline_n_knots=5)

        for path in paths
            df = DataFrame(Arrow.Table(path))
            @test hasproperty(df, :abundance_normalized)
            # Correction factor 1.0: abundances pass through untouched.
            @test all(df.abundance_normalized .≈ Float64.(df.abundance))
        end
    end
end

end # top-level testset
