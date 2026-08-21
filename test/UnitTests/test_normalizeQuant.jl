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
        precursor_idx = UInt32.(1:n),
        irt_obs = rts,
        abundance = abundance,
    )
    Arrow.write(path, df)
end

# Create the left-censored loading example that motivated the matched-precursor
# normalizer. Run 2 is shifted down in log2 space, and values that fall below
# the limit of detection are absent rather than represented as zeros.
function _write_censored_quant_pair(
    dir::String;
    n::Int = 10_000,
    loading_offset::Float32 = 3.0f0,
    log2_lod::Float32 = 4.0f0,
)
    precursor_idx = UInt32.(1:n)
    rts = Float32.(collect(LinRange(0.0, 100.0, n)))
    base_log2 = Float32[4.0f0 + 6.0f0 * ((i - 1) % 101) / 100 for i in 1:n]

    path1 = joinpath(dir, "complete.arrow")
    Arrow.write(path1, DataFrame(
        precursor_idx = precursor_idx,
        irt_obs = rts,
        abundance = 2.0f0 .^ base_log2,
    ))

    run2_log2 = base_log2 .- loading_offset
    observed = run2_log2 .>= log2_lod
    path2 = joinpath(dir, "censored.arrow")
    Arrow.write(path2, DataFrame(
        precursor_idx = precursor_idx[observed],
        irt_obs = rts[observed],
        abundance = 2.0f0 .^ run2_log2[observed],
    ))
    return path1, path2, observed
end

function _write_rt_bias_quant_pair(dir::String; n::Int = 5_000)
    precursor_idx = UInt32.(1:n)
    rts = Float32.(collect(LinRange(0.0, 100.0, n)))
    base_log2 = Float32[8.0f0 + 0.5f0 * sin(Float32(i) / 17.0f0) for i in 1:n]
    rt_bias = @. -2.0f0 + 0.04f0 * rts

    path1 = joinpath(dir, "rt_reference.arrow")
    path2 = joinpath(dir, "rt_biased.arrow")
    Arrow.write(path1, DataFrame(
        precursor_idx = precursor_idx,
        irt_obs = rts,
        abundance = 2.0f0 .^ base_log2,
    ))
    Arrow.write(path2, DataFrame(
        precursor_idx = precursor_idx,
        irt_obs = rts,
        abundance = 2.0f0 .^ (base_log2 .+ rt_bias),
    ))
    return path1, path2
end

# Four runs with two anchor populations interleaved throughout the gradient.
# Stable anchors occur in every run and have residuals ±1 in runs 1 and 2.
# Less-complete anchors occur only in runs 1 and 2 and have residuals ±4. They
# are slightly more numerous, so an unweighted median chooses ±4, while their
# 1/2 completeness weight lets the stable anchors determine the bin summary.
function _write_completeness_weighted_runs(dir::String; n::Int = 2_000)
    precursor_idx = UInt32.(1:n)
    rts = Float32.(collect(LinRange(0.0, 100.0, n)))
    complete = BitVector(mod(i - 1, 20) < 9 for i in 1:n)
    observed = (trues(n), trues(n), complete, complete)

    run_log2 = (
        ifelse.(complete, 9.0f0, 12.0f0),
        ifelse.(complete, 7.0f0, 4.0f0),
        fill(8.0f0, n),
        fill(8.0f0, n),
    )
    paths = String[]
    for run in eachindex(run_log2)
        keep = observed[run]
        path = joinpath(dir, "run$(run).arrow")
        Arrow.write(path, DataFrame(
            precursor_idx = precursor_idx[keep],
            irt_obs = rts[keep],
            abundance = 2.0f0 .^ run_log2[run][keep],
        ))
        push!(paths, path)
    end
    return paths, complete
end

@testset "normalizeQuant" begin

@testset "weighted median" begin
    @test Pioneer._weighted_median([1.0, 2.0, 10.0], ones(3)) == 2.0
    @test Pioneer._weighted_median([1.0, 2.0, 10.0, 20.0], ones(4)) == 6.0
    @test Pioneer._weighted_median([1.0, 2.0, 10.0, 20.0], fill(2.0f0 / 3.0f0, 4)) == 6.0
    @test Pioneer._weighted_median([1.0, 2.0, 10.0], [1.0, 1.0, 3.0]) == 10.0
    @test Pioneer._weighted_median([10.0, 1.0, 2.0], [3.0, 1.0, 1.0]) == 10.0
end

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
# Completeness weights — favor anchors seen consistently across runs
# ═══════════════════════════════════════════════════════════════════════════════

@testset "completeness-weighted matched residuals" begin
    mktempdir() do dir
        paths, complete = _write_completeness_weighted_runs(dir)
        _, weights = Pioneer.getPrecursorQuantReference(paths, :abundance)

        complete_pid = UInt32(findfirst(complete))
        partial_pid = UInt32(findfirst(.!complete))
        @test weights[complete_pid] == 1.0f0
        @test weights[partial_pid] == 0.5f0

        splines, _ = Pioneer.getQuantSplines(
            paths,
            :abundance;
            N=8,
            spline_n_knots=5,
        )
        @test length(splines) == 4
        for rt in (10.0, 50.0, 90.0)
            @test splines[paths[1]](rt) ≈ 1.0 atol=0.05
            @test splines[paths[2]](rt) ≈ -1.0 atol=0.05
            @test splines[paths[3]](rt) ≈ 0.0 atol=0.05
            @test splines[paths[4]](rt) ≈ 0.0 atol=0.05
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
# Differential missingness — normalize matched precursors, not marginals
# ═══════════════════════════════════════════════════════════════════════════════

@testset "left-censored run uses matched precursor ratios" begin
    mktempdir() do dir
        path1, path2, observed = _write_censored_quant_pair(dir)

        before1 = copy(DataFrame(Arrow.Table(path1)).abundance)
        before2 = copy(DataFrame(Arrow.Table(path2)).abundance)
        raw_matched_shift = median(
            log2.(before1[observed]) .- log2.(before2)
        )
        @test raw_matched_shift ≈ 3.0 atol=0.02

        Pioneer.normalizeQuant(
            [path1, path2],
            :abundance;
            N=100,
            spline_n_knots=7,
        )

        after1 = DataFrame(Arrow.Table(path1))
        after2 = DataFrame(Arrow.Table(path2))
        matched_shift = median(
            log2.(after1.abundance_normalized[observed]) .-
            log2.(after2.abundance_normalized)
        )

        # The shared precursors align even though the observed marginal
        # distributions do not: run 2 is still missing its low-abundance tail.
        @test abs(matched_shift) < 0.05
        @test abs(
            median(log2.(after1.abundance_normalized)) -
            median(log2.(after2.abundance_normalized))
        ) > 1.0
        @test (
            maximum(log2.(after2.abundance_normalized)) -
            minimum(log2.(after2.abundance_normalized))
        ) < (
            maximum(log2.(after1.abundance_normalized)) -
            minimum(log2.(after1.abundance_normalized))
        )
    end
end

@testset "matched residual spline corrects RT-dependent bias" begin
    mktempdir() do dir
        path1, path2 = _write_rt_bias_quant_pair(dir)
        Pioneer.normalizeQuant(
            [path1, path2],
            :abundance;
            N=50,
            spline_n_knots=7,
        )

        run1 = DataFrame(Arrow.Table(path1))
        run2 = DataFrame(Arrow.Table(path2))
        normalized_difference = log2.(run2.abundance_normalized) .-
                                log2.(run1.abundance_normalized)

        # Check local agreement throughout the gradient, not just the global
        # median. The injected bias ranges from -2 to +2 log2 units.
        for bin in Pioneer._occupancy_bins(length(normalized_difference), 500, 10)
            @test abs(median(@view(normalized_difference[bin]))) < 0.05
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

@testset "no matched precursor overlap" begin
    mktempdir() do dir
        path1 = joinpath(dir, "run1.arrow")
        path2 = joinpath(dir, "run2.arrow")
        _write_quant_arrow(path1, 1000, 0.0)
        _write_quant_arrow(path2, 1000, 2.0)

        # Copy before replacing the backing Arrow file; Arrow columns may be
        # memory-mapped on macOS.
        run2 = DataFrame(Arrow.Table(path2); copycols=true)
        run2.precursor_idx .+= UInt32(10_000)
        Pioneer.writeArrow(path2, run2)

        before1 = copy(DataFrame(Arrow.Table(path1)).abundance)
        before2 = copy(DataFrame(Arrow.Table(path2)).abundance)
        Pioneer.normalizeQuant([path1, path2], :abundance; N=50, spline_n_knots=5)

        after1 = DataFrame(Arrow.Table(path1))
        after2 = DataFrame(Arrow.Table(path2))
        @test after1.abundance_normalized ≈ before1
        @test after2.abundance_normalized ≈ before2
    end
end

end # top-level testset
