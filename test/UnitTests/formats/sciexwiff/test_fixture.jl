# The truncated SCIEX SWATH fixture (Zenodo; made by test/fixtures/tools/make_wiff_fixture.jl): opening the patched
# .wiff + truncated .wiff.scan, .wiff -> .scxs conversion, reading the .scxs back, and convertSciex end to end.

using SHA

@testset "SCIEX .wiff fixture" begin
    r = SciexWiff.WiffRun(FIXTURE_WIFF)
    @test length(r) == 159_558 && length(r.windows) == 173                  # all slots kept, most empty
    kept = [k for k in 1:length(r) if !SciexWiff.isempty_scan(r, k)]
    @test length(kept) == 8_916
    @test all(30 <= SciexWiff.retention_time_min(r, k) <= 50 for k in kept)
    buf = SciexWiff.ScanBuffer()
    SciexWiff.read_scan!(buf, r, kept[1]); @test buf.n > 0 && issorted(view(buf.bin, 1:buf.n))

    # conversion: byte-identical to the reference, or the same scans and totals within a tight tolerance on a
    # platform whose floating-point centroiding differs
    out = mktempdir()
    SciexWiff.convert_run(FIXTURE_WIFF, out; params = SciexWiff.ConvertParams(format = :scxs), name = "fixture", log = devnull)
    scxs = joinpath(out, "fixture.scxs")
    expected = Dict(split(l, ',')[1] => split(l, ',')[2]
                    for l in Iterators.drop(eachline(joinpath(@__DIR__, "fixtures", "fixture_scxs_sha256.csv")), 1))
    identical = all(bytes2hex(open(sha256, joinpath(scxs, fn))) == h for (fn, h) in expected)
    identical || @warn "fixture .scxs not byte-identical to the reference on this platform; checking totals instead"
    d = Pioneer.loadMassSpecData(scxs)
    @test d isa Pioneer.ScxsMassSpecData && length(d) == 8_916
    pbuf = Pioneer.PeakDecodeBuffer(); npk = 0; sint = 0.0; sorted = true
    for i in 1:length(d)
        mz, it = Pioneer.getPeaks!(pbuf, d, i)
        npk += length(mz); sint += sum(Float64, it; init = 0.0); sorted &= issorted(mz)
    end
    @test sorted
    @test isapprox(npk, 5_054_136; rtol = identical ? 0 : 1e-4)
    @test isapprox(sint, 5.356238101387197e8; rtol = identical ? 1e-9 : 1e-4)

    # convertSciex (the pioneer convert-sciex entry point) on the fixture
    paths = convertSciex(FIXTURE_WIFF; output_dir = joinpath(out, "cs"))
    @test length(paths) == 1 && length(Pioneer.loadMassSpecData(only(paths))) == 8_916
end
