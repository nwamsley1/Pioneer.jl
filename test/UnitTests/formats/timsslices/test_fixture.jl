# The truncated E. coli diaPASEF fixture (Zenodo; made by test/fixtures/tools/make_d_fixture.jl): raw frame decode,
# .d -> .tdfs conversion, reading the .tdfs back, and convertBruker end to end.

using SHA

@testset "E. coli .d fixture" begin
    # raw frames: integer decode, identical on every platform
    g = TS.open_tdf(FIXTURE_D)
    buf = TS.FrameBuffer()
    lines = collect(Iterators.drop(eachline(joinpath(@__DIR__, "fixtures", "ecoli_fixture_frames_checksums.csv")), 1))
    @test length(lines) == TS.n_frames(g) == 279
    for line in lines
        row, fid, mt, ns, np, stof, sint, sscan, mtof = parse.(Int, split(line, ','))
        TS.read_frame!(buf, g, row)
        @test (g.frames.id[row], g.frames.msms_type[row], buf.n_scans, buf.n_peaks) == (fid, mt, ns, np)
        @test sum(Int, view(buf.tof, 1:np); init = 0) == stof && sum(Int, view(buf.intensity, 1:np); init = 0) == sint
        @test sum(Int(s) * (buf.scan_start[s + 2] - buf.scan_start[s + 1]) for s in 0:ns-1) == sscan
        @test (np == 0 ? 0 : Int(maximum(view(buf.tof, 1:np)))) == mtof
    end
    # the fixture keeps one window per MS2 frame (537.5 or 562.5 m/z)
    @test all(length(TS.windows(g, j)) == 1 for j in 1:TS.n_frames(g) if TS.is_dia(g, j))

    # conversion: byte-identical to the reference where floating point agrees; otherwise the same slices and
    # centroid totals within a tight tolerance (vectorised smoothing may round differently on another CPU)
    out = mktempdir()
    c = TS.convert_run(FIXTURE_D, out; name = "fixture", log = devnull)
    expected = Dict(split(l, ',')[1] => split(l, ',')[2]
                    for l in Iterators.drop(eachline(joinpath(@__DIR__, "fixtures", "ecoli_fixture_tdfs_sha256.csv")), 1))
    identical = all(bytes2hex(open(sha256, joinpath(c.tdfs, fn))) == h for (fn, h) in expected)
    identical || @warn "fixture .tdfs not byte-identical to the reference on this platform; checking totals instead"
    t = TS.open_tdfs(c.tdfs)
    @test TS.n_frames(t) == 279 && TS.n_slices(t) == 12_414
    d = Pioneer.TdfsMassSpecData(c.tdfs)
    pbuf = Pioneer.PeakDecodeBuffer()
    npk = 0; sint = 0.0; sorted = true; bp_ok = true
    for i in 1:length(d)
        mz, it = Pioneer.getPeaks!(pbuf, d, i)
        npk += length(mz); sint += sum(Float64, it; init = 0.0); sorted &= issorted(mz)
        bp_ok &= Pioneer.getBasePeakIntensity(d, i) === (isempty(it) ? 0f0 : maximum(it))   # format 3: the slice's largest intensity
    end
    @test length(d) == 12_414 && sorted && bp_ok
    @test isapprox(npk, 10_951_008; rtol = identical ? 0 : 1e-4)
    @test isapprox(sint, 1.522929679e9; rtol = identical ? 1e-9 : 1e-4)

    # a format-2 .tdfs (no base_peak_intensity column) still reads: base peak missing, expand recomputes it
    v2 = joinpath(out, "v2.tdfs"); cp(c.tdfs, v2)
    sl = TS.Arrow.Table(joinpath(c.tdfs, "slices.arrow"))
    rm(joinpath(v2, "slices.arrow"))
    TS.Arrow.write(joinpath(v2, "slices.arrow"), NamedTuple(k => getproperty(sl, k) for k in propertynames(sl) if k != :base_peak_intensity))
    m = read(joinpath(v2, "meta.json"), String)
    write(joinpath(v2, "meta.json"), replace(m, r"\"format_version\": *3" => "\"format_version\": 2"))
    @test TS.open_tdfs(v2).meta["format_version"] == 2
    @test ismissing(Pioneer.getBasePeakIntensity(Pioneer.TdfsMassSpecData(v2), 1))
    a3 = TS.Arrow.Table(TS.expand(c.tdfs, joinpath(out, "v3.arrow"); log = devnull))
    a2 = TS.Arrow.Table(TS.expand(v2, joinpath(out, "v2.arrow"); log = devnull))
    @test collect(a2.basePeakIntensity) == collect(a3.basePeakIntensity) == Float32[Pioneer.getBasePeakIntensity(d, i) for i in 1:length(d)]

    # convertBruker (the pioneer convert-bruker entry point) on the fixture
    paths = convertBruker(FIXTURE_D; output_dir = joinpath(out, "cb"))
    @test length(paths) == 1 && isdir(only(paths)) && TS.n_slices(TS.open_tdfs(only(paths))) == 12_414
end
