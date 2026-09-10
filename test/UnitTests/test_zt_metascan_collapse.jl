using Test
using Arrow
using DataFrames

using Pioneer
using Pioneer: BasicMassSpecData, getMsOrder, getCenterMz, collapse_to_metascans,
               zt_transmission_template, ZTGeometry, _zt_profile_features,
               ZT_PROFILE_FEATURES, ZT_SHAPE_FEATURES

const CSTEP = 1.0f0
const CSTART = 400.0f0
const CBINS = 20

"""Single-cycle ZT arrow: MS1 head + `CBINS` contiguous MS2 bins of width `CSTEP`."""
function write_collapse_arrow(path::String; n_cycles::Int = 1)
    per = CBINS + 1
    n = n_cycles * per
    ms_order  = Vector{UInt8}(undef, n)
    center_mz = Vector{Union{Missing,Float32}}(undef, n)
    isol_w    = Vector{Union{Missing,Float32}}(undef, n)
    cycle_idx = Vector{UInt32}(undef, n)
    for c in 1:n_cycles, b in 0:CBINS
        i = (c - 1) * per + b + 1
        cycle_idx[i] = UInt32(c)
        if b == 0
            ms_order[i] = UInt8(1); center_mz[i] = missing; isol_w[i] = missing
        else
            ms_order[i] = UInt8(2)
            center_mz[i] = CSTART + Float32(b - 1) * CSTEP
            isol_w[i] = CSTEP
        end
    end
    mz  = [Vector{Union{Missing,Float32}}(Float32[100f0]) for _ in 1:n]
    int = [Vector{Union{Missing,Float32}}(Float32[10f0])  for _ in 1:n]
    Arrow.write(path, DataFrame(
        mz_array = mz, intensity_array = int,
        scanHeader = fill("s", n), scanNumber = Int32.(1:n),
        basePeakMz = fill(100f0, n), basePeakIntensity = fill(1000f0, n),
        retentionTime = Float32.(range(0, length = n, step = 0.01f0)),
        lowMz = fill(100f0, n), highMz = fill(1000f0, n), TIC = fill(5f3, n),
        centerMz = center_mz, isolationWidthMz = isol_w,
        msOrder = ms_order, packetType = Int32.(fill(0, n)), cycle_idx = cycle_idx))
    return BasicMassSpecData(path)
end

struct TestPrecursors; mz::Vector{Float32}; end
Pioneer.getMz(p::TestPrecursors) = p.mz

"""Build a PSM table from explicit (precursor, scan, weight) triples."""
function psm_table(rows::Vector{<:Tuple}; with_frags::Bool = false)
    df = DataFrame(
        precursor_idx = UInt32[r[1] for r in rows],
        scan_idx      = UInt32[r[2] for r in rows],
        weight        = Float32[r[3] for r in rows],
        tag           = collect(1:length(rows)),
    )
    if with_frags
        for b in 1:8
            df[!, Symbol("frag$(b)_int")] = Float32[r[3] / Float32(b) for r in rows]
        end
    end
    return df
end

"""Brute-force center test: |prec_mz - centerMz| <= isolationWidth/2."""
brute_centers(df, spectra, pmz) = [i for i in 1:nrow(df) if
    abs(pmz[df.precursor_idx[i]] - Float32(getCenterMz(spectra, Int(df.scan_idx[i])))) <= CSTEP / 2]

@testset "ZT collapse_to_metascans" begin
    spectra = write_collapse_arrow(joinpath(mktempdir(), "c.arrow"))
    k = 3
    L = 2k + 1
    geom = ZTGeometry(CSTEP, CSTEP, Int32(CBINS), Int32(k), 3.0f0)
    tri, tnorm = zt_transmission_template(geom, k)

    # bin b (1-based) -> scan index, MS1 at the head
    bscan(b) = b + 1
    # precursor 1 sits exactly on bin 10's center; precursor 2 half a bin above bin 10
    pmz = Float32[CSTART + 9 * CSTEP, CSTART + 9 * CSTEP + 0.45f0 * CSTEP]
    precs = TestPrecursors(pmz)

    @testset "profile matching the transmission template -> zt_tri_cosine == 1" begin
        rows = Tuple{Int,Int,Float32}[]
        for d in -k:k
            push!(rows, (1, bscan(10 + d), tri[d + k + 1]))
        end
        df = psm_table(rows)
        meta = collapse_to_metascans(df, spectra, precs, geom)
        @test nrow(meta) == 1
        @test isapprox(meta.zt_tri_cosine[1], 1.0f0; atol = 1e-5)
        # entropy of that profile, computed independently
        w = Float32[tri[j] for j in 1:L]
        W = sum(w); expected_ent = -sum(p -> p > 0 ? p * log(p) : 0f0, w ./ W)
        @test isapprox(meta.zt_entropy[1], expected_ent; atol = 1e-5)
    end

    @testset "flat profile -> entropy == log(2k+1)" begin
        rows = [(1, bscan(10 + d), 1.0f0) for d in -k:k]
        df = psm_table(rows)
        meta = collapse_to_metascans(df, spectra, precs, geom)
        @test nrow(meta) == 1
        @test isapprox(meta.zt_entropy[1], Float32(log(L)); atol = 1e-5)
        @test meta.zt_tri_cosine[1] < 1.0f0        # flat is a poor triangle match
    end

    @testset "center selection matches brute force" begin
        rows = Tuple{Int,Int,Float32}[]
        for b in 4:16, p in 1:2
            push!(rows, (p, bscan(b), Float32(b)))
        end
        df = psm_table(rows)
        want = brute_centers(df, spectra, pmz)
        meta = collapse_to_metascans(df, spectra, precs, geom)
        @test nrow(meta) == length(want)
        @test sort(meta.tag) == sort(df.tag[want])
    end

    @testset "k <= 0 returns the input unchanged" begin
        df = psm_table([(1, bscan(10), 1.0f0)])
        out = collapse_to_metascans(df, spectra, precs, ZTGeometry(CSTEP, CSTEP, Int32(CBINS), Int32(0), 3.0f0))
        @test out === df
    end

    @testset "empty input" begin
        df = psm_table(Tuple{Int,Int,Float32}[])
        @test nrow(collapse_to_metascans(df, spectra, precs, geom)) == 0
    end

    @testset "all feature columns present and finite" begin
        rows = Tuple{Int,Int,Float32}[]
        for d in -k:k, p in 1:2
            push!(rows, (p, bscan(10 + d), tri[d + k + 1] * Float32(p)))
        end
        df = psm_table(rows; with_frags = true)
        meta = collapse_to_metascans(df, spectra, precs, geom)
        @test nrow(meta) == 2
        for f in vcat(ZT_PROFILE_FEATURES, ZT_SHAPE_FEATURES)
            @test hasproperty(meta, f)
            @test all(isfinite, Float64.(meta[!, f]))
        end
        # fragments here are exact scalings of the weight profile, so they track it perfectly
        @test all(meta.n_correlated_fragments_shape .== UInt8(8))
        @test all(meta.frag_apex_dispersion_shape .== 0f0)
    end

    @testset "passenger columns survive the collapse" begin
        rows = [(1, bscan(10 + d), tri[d + k + 1]) for d in -k:k]
        df = psm_table(rows)
        meta = collapse_to_metascans(df, spectra, precs, geom)
        @test hasproperty(meta, :tag)
        @test hasproperty(meta, :weight)
        @test nrow(meta) == 1
    end
end
