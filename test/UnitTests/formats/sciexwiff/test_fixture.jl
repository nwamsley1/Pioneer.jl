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
    f = SciexWiff.open_scxs(scxs)
    # not declared ZT: plain name, acquisition_type swath; declaring it adds the Q1 bin grid
    @test f.meta["acquisition_type"] == "swath" && f.meta["acquisition_type_source"] == "user"
    zm = SciexWiff.acquisition_metadata(r; zt_scan = true)
    @test zm["acquisition_type"] == "zt_scan_dia" && parse(Int, zm["q1_bins_per_cycle"]) == 173
    d = Pioneer.loadMassSpecData(scxs)
    @test d isa Pioneer.ScxsMassSpecData && length(d) == 8_916
    pbuf = Pioneer.PeakDecodeBuffer(); npk = 0; sint = 0.0; sorted = true; bp_int_ok = true; bp_mz_ok = true
    for i in 1:length(d)
        mz, it = Pioneer.getPeaks!(pbuf, d, i)
        npk += length(mz); sint += sum(Float64, it; init = 0.0); sorted &= issorted(mz)
        # base peak: the scan's most intense centroid, as read back (centroid units, like the .tdfs base peak)
        if isempty(it)
            bp_int_ok &= Pioneer.getBasePeakIntensity(d, i) === 0f0; bp_mz_ok &= isnan(Pioneer.getBasePeakMz(d, i))
        else
            k = argmax(it)
            bp_int_ok &= Pioneer.getBasePeakIntensity(d, i) === it[k]
            bp_mz_ok &= isapprox(Pioneer.getBasePeakMz(d, i), mz[k]; rtol = 1e-6)
        end
    end
    @test sorted && bp_int_ok && bp_mz_ok
    @test isapprox(npk, 5_054_136; rtol = identical ? 0 : 1e-4)
    @test isapprox(sint, 5.356238101387197e8; rtol = identical ? 1e-9 : 1e-4)

    # convertSciex (the pioneer convert-sciex entry point) on the fixture
    paths = redirect_stdout(devnull) do
        convertSciex(FIXTURE_WIFF; output_dir = joinpath(out, "cs"), zt_scan = false)
    end
    @test length(paths) == 1 && length(Pioneer.loadMassSpecData(only(paths))) == 8_916
    @test basename(only(paths)) == "BenchSample_B_nswath4_25ng.scxs"
end
