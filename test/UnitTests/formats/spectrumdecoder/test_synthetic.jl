# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Reading, checking and writing synthetic acquisition files (write_synthetic_acquisition).

using Arrow

synthetic_dir = mktempdir()

@testset "synthetic centroid acquisition" begin
    p = write_synthetic_acquisition(joinpath(synthetic_dir, "centroid.raw"); cycles = 5, windows = 10)
    f = open_acquisition(p)
    S = scans(f)
    @test f.version == 66
    @test length(S) == 55
    @test count(s -> s.ms_order == 0x01, S) == 5 && count(s -> s.ms_order == 0x02, S) == 50
    @test issorted([s.retention_time for s in S])
    @test all(s -> s.analyser == 0x04 && s.packet_type == 0x00000014, S)
    ms2 = filter(s -> s.ms_order == 0x02, S)
    @test sort(unique([s.center_mz for s in ms2])) == [395.0 + 10w for w in 1:10]
    @test all(s -> s.isolation_width == 10.0 && s.collision_energy == 27.0, ms2)
    @test all(s -> isnan(s.center_mz), filter(s -> s.ms_order == 0x01, S))

    mz, it = peaks(f, 1)
    @test (mz, it) == SpectrumDecoder.synthetic_ms1_peaks(1)
    mz, it = peaks(f, 2)
    @test (mz, it) == SpectrumDecoder.synthetic_ms2_peaks(1)
    @test S[1].tic == sum(Float64, peaks(f, 1)[2])

    R = check_consistency(f)
    @test all(SpectrumDecoder.ok, R)

    out = write_pioneer_arrow(f, joinpath(synthetic_dir, "centroid.arrow"))
    t = Arrow.Table(out)
    @test length(t[:msOrder]) == 55
    @test collect(t[:mz_array][2]) == peaks(f, 2)[1]
    @test collect(t[:intensity_array][2]) == peaks(f, 2)[2]
    @test t[:scanNumber] == Int32.(1:55)
    @test ismissing(t[:centerMz][1]) && t[:centerMz][2] == 405f0
    # a cycle starts at the first MS2 window, so each MS1 scan ends the previous cycle
    @test t[:cycle_idx] == Int32[max(1, cld(i - 1, 11)) for i in 1:55]

    # threads and record-batch size do not change the output
    one = write_pioneer_arrow(f, joinpath(synthetic_dir, "one.arrow"); threads = 1)
    small = write_pioneer_arrow(f, joinpath(synthetic_dir, "small.arrow"); batch_peaks = 50, threads = 4)
    for q in (one, small)
        u = Arrow.Table(q)
        @test [collect(x) for x in u[:mz_array]] == [collect(x) for x in t[:mz_array]]
        @test [collect(x) for x in u[:intensity_array]] == [collect(x) for x in t[:intensity_array]]
        @test u[:retentionTime] == t[:retentionTime]
    end
end

@testset "synthetic profile-only ion-trap MS1" begin
    p = write_synthetic_acquisition(joinpath(synthetic_dir, "profile.raw"); cycles = 4, windows = 10, ms1_profile = true)
    f = open_acquisition(p)
    S = scans(f)
    @test all(s -> s.packet_type == 0x00000013 && s.analyser == 0x00, filter(s -> s.ms_order == 0x01, S))
    @test isempty(peaks(f, 1)[1])
    pmz, pit = profile(f, 1)
    @test length(pmz) == 9199 && pmz[1] == 350.125 && S[1].tic == sum(Float64, pit)
    @test all(SpectrumDecoder.ok, check_consistency(f))

    t = Arrow.Table(write_pioneer_arrow(f, joinpath(synthetic_dir, "profile.arrow")))
    # every isotope envelope (three peaks, partly merged at FWHM 0.5) has a centroid within its span, and every
    # centroid lies within an envelope
    mz1 = collect(t[:mz_array][1])
    tmz, _ = SpectrumDecoder.synthetic_ms1_peaks(1)
    spans = [(tmz[3*e+1] - 0.1, tmz[3*e+3] + 0.1) for e in 0:7]
    @test issorted(mz1)
    @test all(sp -> any(m -> sp[1] <= m <= sp[2], mz1), spans)
    @test all(m -> any(sp -> sp[1] <= m <= sp[2], spans), mz1)
    # without centroiding, profile-only scans are written with no peaks
    u = Arrow.Table(write_pioneer_arrow(f, joinpath(synthetic_dir, "profile_raw.arrow"); centroid_profile = false))
    @test isempty(u[:mz_array][1]) && collect(u[:mz_array][2]) == peaks(f, 2)[1]
end

@testset "files that are refused" begin
    p = write_synthetic_acquisition(joinpath(synthetic_dir, "ok.raw"); cycles = 3, windows = 4)
    b = read(p)
    # truncated: the device block at the end is missing
    write(joinpath(synthetic_dir, "truncated.raw"), b[1:end-8000])
    @test_throws ErrorException open_acquisition(joinpath(synthetic_dir, "truncated.raw"))
    # unknown format version
    v = copy(b); v[0x24+1:0x24+4] .= reinterpret(UInt8, [htol(UInt32(65))])
    write(joinpath(synthetic_dir, "v65.raw"), v)
    @test_throws ErrorException open_acquisition(joinpath(synthetic_dir, "v65.raw"))
    # wrong leading marker
    m = copy(b); m[1] = 0x00
    write(joinpath(synthetic_dir, "marker.raw"), m)
    @test_throws ErrorException open_acquisition(joinpath(synthetic_dir, "marker.raw"))
end
