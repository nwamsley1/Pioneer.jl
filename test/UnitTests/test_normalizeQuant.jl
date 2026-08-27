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

function _write_pairwise_chain_runs(
    dir::String;
    n_per_edge::Int = 1_000,
    offsets::Vector{Float32} = Float32[0.0, 1.0, -2.0, 2.0, -1.0],
)
    n_runs = length(offsets)
    paths = String[]
    observed_ids_by_run = Dict{UInt32, Vector{UInt32}}()
    for run in 1:n_runs
        precursor_idx = UInt32[]
        irt_obs = Float32[]
        base_log2 = Float32[]
        first_edge = max(1, run - 1)
        last_edge = min(n_runs - 1, run)
        for edge in first_edge:last_edge
            append!(
                precursor_idx,
                UInt32.(((edge - 1) * n_per_edge + 1):(edge * n_per_edge)),
            )
            append!(irt_obs, Float32.(LinRange(0.0, 100.0, n_per_edge)))
            append!(
                base_log2,
                Float32[8.0f0 + 0.25f0 * sin(Float32(i) / 17.0f0)
                        for i in 1:n_per_edge],
            )
        end
        # Give endpoint runs the same total direct-ID count as interior runs.
        # This fixture isolates tree propagation; count imbalance is exercised
        # separately below.
        if run == 1 || run == n_runs
            filler_start = (n_runs - 1 + run) * n_per_edge + 1
            append!(
                precursor_idx,
                UInt32.(filler_start:(filler_start + n_per_edge - 1)),
            )
            append!(irt_obs, Float32.(LinRange(0.0, 100.0, n_per_edge)))
            append!(
                base_log2,
                Float32[8.0f0 + 0.25f0 * sin(Float32(i) / 17.0f0)
                        for i in 1:n_per_edge],
            )
        end
        path = joinpath(dir, "chain_run$(run).arrow")
        Arrow.write(path, DataFrame(
            precursor_idx = precursor_idx,
            irt_obs = irt_obs,
            abundance = 2.0f0 .^ (base_log2 .+ offsets[run]),
            target = trues(length(precursor_idx)),
            mbr_recovered = falses(length(precursor_idx)),
        ))
        push!(paths, path)
        observed_ids_by_run[UInt32(run)] = precursor_idx
    end
    return paths, Pioneer.build_run_similarity(observed_ids_by_run), offsets
end

function _write_mbr_anchor_pair(
    dir::String;
    n_observed::Int = 1_000,
    n_recovered::Int = 2_000,
)
    precursor_idx = UInt32.(1:(n_observed + n_recovered))
    irt_obs = Float32.(LinRange(0.0, 100.0, n_observed + n_recovered))
    recovered = BitVector(vcat(falses(n_observed), trues(n_recovered)))
    run1_log2 = vcat(fill(8.0f0, n_observed), fill(14.0f0, n_recovered))
    run2_log2 = vcat(fill(10.0f0, n_observed), fill(4.0f0, n_recovered))
    paths = String[]
    for (run, log2_quant) in enumerate((run1_log2, run2_log2))
        path = joinpath(dir, "mbr_run$(run).arrow")
        Arrow.write(path, DataFrame(
            precursor_idx = precursor_idx,
            irt_obs = irt_obs,
            abundance = 2.0f0 .^ log2_quant,
            target = trues(length(precursor_idx)),
            mbr_recovered = recovered,
        ))
        push!(paths, path)
    end
    observed_ids_by_run = Dict(
        UInt32(1) => copy(precursor_idx),
        UInt32(2) => copy(precursor_idx),
    )
    return paths, Pioneer.build_run_similarity(observed_ids_by_run), n_observed
end

function _write_abundance_extreme_biased_quant_pair(dir::String; n::Int = 4_000)
    precursor_idx = UInt32.(1:n)
    irt_obs = Float32.(LinRange(0.0, 100.0, n))
    abundance_decile = Int[mod(i - 1, 10) for i in 1:n]
    biased = BitVector(decile < 3 || decile >= 7 for decile in abundance_decile)
    pair_mean_log2 = Float32[7.0f0 + decile for decile in abundance_decile]

    # The true loading shift is three log2 units. Ratios at both abundance
    # extremes are compressed to one unit, so the ordinary median and a
    # monotonically increasing abundance weight both select the biased value.
    # Middle-peaked abundance weights select the reliable central anchors.
    observed_ratio = Float32[biased[i] ? 1.0f0 : 3.0f0 for i in 1:n]
    reference_log2 = pair_mean_log2 .+ observed_ratio ./ 2.0f0
    observed_log2 = pair_mean_log2 .- observed_ratio ./ 2.0f0

    paths = String[]
    for (run, log2_quant) in enumerate((reference_log2, observed_log2))
        path = joinpath(dir, "detection_biased_run$(run).arrow")
        Arrow.write(path, DataFrame(
            precursor_idx = precursor_idx,
            irt_obs = irt_obs,
            abundance = 2.0f0 .^ log2_quant,
            target = trues(n),
            mbr_recovered = falses(n),
        ))
        push!(paths, path)
    end

    observed_ids_by_run = Dict(
        UInt32(1) => copy(precursor_idx),
        UInt32(2) => copy(precursor_idx),
    )
    return paths, Pioneer.build_run_similarity(observed_ids_by_run), biased
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
# Pairwise maximum spanning tree — local exact matches without global anchors
# ═══════════════════════════════════════════════════════════════════════════════

@testset "pairwise quant anchors" begin
    frame = DataFrame(
        precursor_idx = UInt32[1, 2, 3, 4, 2, 5, 6],
        irt_obs = Float32.(1:7),
        abundance = Float32[1, 4, 4, 16, 2, 1_024, 2_048],
        target = Bool[true, true, true, true, true, false, true],
        mbr_recovered = Bool[false, false, false, false, false, false, true],
    )

    anchors = Pioneer._file_precursor_quant_anchors(
        frame,
        :abundance,
        :precursor_idx,
    )

    @test length(anchors) == 4
    @test anchors[UInt32(1)].log2_quant ≈ 0.0f0
    @test anchors[UInt32(2)].log2_quant ≈ 2.0f0
    @test anchors[UInt32(3)].log2_quant ≈ 2.0f0
    @test anchors[UInt32(4)].log2_quant ≈ 4.0f0
end

@testset "pairwise rank weights require middle abundance in both runs" begin
    left_log2 = Float64[10, 20, 30, 40, 100, 60, 70, 80, 90]
    right_log2 = Float64[10, 20, 30, 40, 0, 60, 70, 80, 90]

    pair_weights = Pioneer._pairwise_middle_rank_weights(
        left_log2,
        right_log2,
    )

    # The fifth anchor has a perfectly middle pairwise mean, but it is the
    # highest-abundance anchor in the left run and the lowest in the right.
    # Separate run ranks therefore downweight it instead of assigning the
    # maximum weight that pair-mean ranking would give it.
    pair_means = (left_log2 .+ right_log2) ./ 2.0
    pair_mean_weight = Pioneer._middle_rank_weight(
        pair_means[5],
        sort(pair_means),
    )
    @test pair_mean_weight ≈ 1.0
    @test pair_weights[5] ≈ 17 / 81
    @test pair_weights[5] < pair_weights[4]
    @test pair_weights[4] > 0.9

    # The estimator is symmetric in edge orientation.
    @test pair_weights ≈ Pioneer._pairwise_middle_rank_weights(
        right_log2,
        left_log2,
    )
    @test isempty(Pioneer._pairwise_middle_rank_weights(Float64[], Float64[]))
    @test_throws DimensionMismatch Pioneer._pairwise_middle_rank_weights(
        left_log2,
        right_log2[1:end-1],
    )
end

@testset "containment-shrunk precursor-count adjustment" begin
    @test Pioneer._pairwise_containment_count_adjustment(100, 100, 80) == 0.0
    @test Pioneer._pairwise_containment_count_adjustment(100, 25, 20) ≈ 0.8
    @test Pioneer._pairwise_containment_count_adjustment(25, 100, 20) ≈ -0.8
    @test Pioneer._pairwise_containment_count_adjustment(100, 25, 0) == 0.0
    @test_throws ArgumentError Pioneer._pairwise_containment_count_adjustment(
        10,
        5,
        6,
    )

    left_anchors = Dict{UInt32, Pioneer.QuantRunAnchor}()
    right_anchors = Dict{UInt32, Pioneer.QuantRunAnchor}()
    for idx in 1:1_000
        rt = Float32(100.0 * (idx - 1) / 999)
        left_anchors[UInt32(idx)] = Pioneer.QuantRunAnchor(10.0f0, rt)
        if idx <= 800
            right_anchors[UInt32(idx)] = Pioneer.QuantRunAnchor(8.0f0, rt)
        end
    end

    edge = Pioneer._fit_pairwise_quant_spline(
        1,
        2,
        1.0f0,
        left_anchors,
        right_anchors;
        N = 20,
        spline_n_knots = 5,
        min_bin_occupancy = 100,
    )
    expected_adjustment = Pioneer._pairwise_containment_count_adjustment(
        length(left_anchors),
        length(right_anchors),
        800,
    )
    @test edge !== nothing
    @test edge.nanchors == 800
    @test edge.spline(50.0) ≈ 2.0 + expected_adjustment atol=0.01
end

@testset "pairwise maximum spanning tree" begin
    mktempdir() do dir
        paths, atlas, offsets = _write_pairwise_chain_runs(dir)
        run_ids = UInt32.(eachindex(paths))

        # Every precursor occurs in exactly two of five runs, below the current
        # experiment-wide 50% anchor requirement. The original matched method
        # therefore cannot fit any run, while adjacent pairwise edges can.
        global_splines, _ = Pioneer.getQuantSplines(
            paths,
            :abundance;
            N=20,
            spline_n_knots=5,
        )
        @test isempty(global_splines)

        tree = Pioneer.getPairwiseQuantTree(
            paths,
            :abundance,
            run_ids,
            atlas;
            N=20,
            spline_n_knots=5,
        )
        @test length(tree) == length(paths) - 1
        @test Set(
            (min(edge.left_position, edge.right_position),
             max(edge.left_position, edge.right_position))
            for edge in tree
        ) == Set((run, run + 1) for run in 1:(length(paths) - 1))
        @test all(edge -> edge.nanchors == 1_000, tree)

        corrections, components = Pioneer.getPairwiseQuantCorrections(
            paths,
            tree;
            N=20,
        )
        @test components == [collect(eachindex(paths))]
        for edge in tree
            left = edge.left_position
            right = edge.right_position
            @test (
                corrections[paths[left]](50.0) -
                corrections[paths[right]](50.0)
            ) ≈ offsets[left] - offsets[right] atol=0.05
        end

        Pioneer.normalizeQuant(
            paths,
            :abundance;
            N=20,
            spline_n_knots=5,
            run_ids=run_ids,
            run_similarity_atlas=atlas,
        )
        normalized_medians = Float64[]
        for path in paths
            frame = DataFrame(Arrow.Table(path))
            push!(normalized_medians, median(log2.(frame.abundance_normalized)))
        end
        @test maximum(normalized_medians) - minimum(normalized_medians) < 0.05
    end
end

@testset "middle-abundance weighting limits ratio compression" begin
    mktempdir() do dir
        paths, atlas, biased = _write_abundance_extreme_biased_quant_pair(dir)
        run_ids = UInt32[1, 2]

        reference = DataFrame(Arrow.Table(paths[1]))
        observed = DataFrame(Arrow.Table(paths[2]))
        unweighted_ratios = log2.(reference.abundance) .-
                            log2.(observed.abundance)

        # Most matched precursors are biased, so the unweighted median
        # underestimates the known three-unit loading difference.
        @test median(unweighted_ratios) ≈ 1.0 atol=1e-5

        tree = Pioneer.getPairwiseQuantTree(
            paths,
            :abundance,
            run_ids,
            atlas;
            N=20,
            spline_n_knots=5,
        )

        @test length(tree) == 1
        @test only(tree).nanchors == length(biased)
        @test only(tree).spline(50.0) ≈ 3.0 atol=0.05

        Pioneer.normalizeQuant(
            paths,
            :abundance;
            N=20,
            spline_n_knots=5,
            run_ids=run_ids,
            run_similarity_atlas=atlas,
        )
        normalized_reference = DataFrame(Arrow.Table(paths[1]))
        normalized_observed = DataFrame(Arrow.Table(paths[2]))
        reliable_ratios = log2.(
            normalized_reference.abundance_normalized[.!biased],
        ) .- log2.(
            normalized_observed.abundance_normalized[.!biased],
        )
        @test median(reliable_ratios) ≈ 0.0 atol=0.05
    end
end

@testset "MBR-recovered rows are excluded from pairwise anchors" begin
    mktempdir() do dir
        paths, atlas, n_observed = _write_mbr_anchor_pair(dir)
        run_ids = UInt32[1, 2]
        consensus = Pioneer.getPrecursorQuantConsensus(paths, :abundance)
        @test length(consensus) == n_observed

        tree = Pioneer.getPairwiseQuantTree(
            paths,
            :abundance,
            run_ids,
            atlas;
            N=20,
            spline_n_knots=5,
        )

        @test length(tree) == 1
        @test only(tree).nanchors == n_observed
        @test only(tree).spline(50.0) ≈ -2.0 atol=0.05

        Pioneer.normalizeQuant(
            paths,
            :abundance;
            N=20,
            spline_n_knots=5,
            run_ids=run_ids,
            run_similarity_atlas=atlas,
        )
        run1 = DataFrame(Arrow.Table(paths[1]))
        run2 = DataFrame(Arrow.Table(paths[2]))
        observed_difference = median(
            log2.(run1.abundance_normalized[1:n_observed]) .-
            log2.(run2.abundance_normalized[1:n_observed])
        )
        @test abs(observed_difference) < 0.05
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
