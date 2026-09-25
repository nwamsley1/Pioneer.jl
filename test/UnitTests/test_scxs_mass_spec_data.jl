# ScxsMassSpecData: SCIEX scans from a `.scxs` (SciexWiff.jl). The synthetic part builds a small `.scxs` with
# TimsSlices' block encoder (the blocks are byte-identical) and checks the getters and the per-scan calibration.
# The real-data part runs when PIONEER_TEST_SCXS points at a directory holding `<name>.scxs` and the `<name>.arrow`
# SciexWiff wrote from the same conversion: every scan's peaks and metadata must match.
# Run: julia --project=. test/UnitTests/test_scxs_mass_spec_data.jl

using Test, Arrow
using TimsSlices
using Pioneer
using Pioneer: ScxsMassSpecData, PeakDecodeBuffer, getPeaks!, loadMassSpecData, is_ms_data_path, is_scxs_path,
               getMzArray, getIntensityArray, getRetentionTime, getLowMz, getHighMz, getTIC, getCenterMz,
               getIsolationWidthMz, getMsOrder, getCycleIdx, getCollisionEnergyEv, getPeakCount, getPeakCounts,
               getRetentionTimes, getTICs, getCenterMzs, getMsOrders, getCycleIdxs, getImScans

"Write a .scxs with the given per-scan (bins, intensities); calibration and metadata are fixed test values."
function write_test_scxs(dir, scans; bin_scale = 4, int_scale = 10.0, cal_a = 4.9e-4, cal_b = -13.7)
    mkpath(dir)
    codec = TimsSlices.BlockCodec(); offs = Int64[]; sizes = Int32[]; off = 0
    open(joinpath(dir, "blocks.bin"), "w") do io
        for (bins, ints) in scans
            blk = TimsSlices.SliceBlock(1, Int32[1, length(bins) + 1], bins, ints)
            nb = TimsSlices.encode_slice!(codec, blk, 1, 3)
            write(io, view(codec.zbuf, 1:nb)); push!(offs, off); push!(sizes, nb); off += nb
        end
    end
    n = length(scans)
    ms1 = [isodd(i) for i in 1:n]
    Arrow.write(joinpath(dir, "scans.arrow"), (
        record = Int32.(1:n), cycle = Int32.((0:n-1) .÷ 2 .+ 1), experiment = Int16.(ifelse.(ms1, 1, 2)),
        ms_order = UInt8.(ifelse.(ms1, 1, 2)), retention_time = Float32.(0.1 .* (1:n)),
        low_mz = Float32.(ifelse.(ms1, 400, 100)), high_mz = Float32.(ifelse.(ms1, 1250, 1500)),
        tic = Float32.(100 .* (1:n)), center_mz = Float32.(ifelse.(ms1, NaN32, 401.45f0)),
        isolation_width = Float32.(ifelse.(ms1, NaN32, 2.9f0)), cal_a = fill(cal_a, n), cal_b = fill(cal_b, n),
        n_peaks = Int32[length(s[1]) for s in scans], block_offset = offs, block_size = sizes))
    write(joinpath(dir, "meta.json"), """{"format_version": 1, "bin_scale": $bin_scale, "int_scale": $int_scale}""")
    dir
end

@testset "ScxsMassSpecData" begin
    @testset "synthetic" begin
        dir = joinpath(mktempdir(), "run.scxs")
        scans = [(UInt32[4_000_000, 4_000_040, 5_000_000], UInt32[10, 25, 7]),   # MS1
                 (UInt32[], UInt32[]),                                          # empty MS2
                 (UInt32[3_000_000], UInt32[123_456]),                          # MS1
                 (UInt32[2_500_000, 2_600_000], UInt32[1, 2])]                  # MS2
        write_test_scxs(dir, scans)
        @test is_scxs_path(dir) && is_ms_data_path(dir) && !is_scxs_path(joinpath(dir, "meta.json"))
        d = loadMassSpecData(dir)
        @test d isa ScxsMassSpecData && length(d) == 4
        @test getImScans(d) === nothing                    # no ion mobility: the Arrow code paths
        @test getPeakCounts(d) == Int32[3, 0, 1, 2]
        @test getMsOrders(d) == UInt8[1, 2, 1, 2] && getCycleIdxs(d) == UInt32[1, 1, 2, 2]
        @test ismissing(getCenterMz(d, 1)) && getCenterMz(d, 2) == 401.45f0 && getIsolationWidthMz(d, 4) == 2.9f0
        @test getLowMz(d, 1) == 400f0 && getHighMz(d, 2) == 1500f0 && getTIC(d, 3) == 300f0
        @test getRetentionTime(d, 4) == 0.4f0 && getCollisionEnergyEv(d, 2) == 0f0
        buf = PeakDecodeBuffer()
        for (i, (bins, ints)) in enumerate(scans)
            mz, it = getPeaks!(buf, d, i)
            @test length(mz) == getPeakCount(d, i) == length(bins)
            @test mz == Float32[(4.9e-4 * (b / 4 / 40 + 13.7))^2 for b in bins]
            @test it == Float32.(ints ./ 10)
        end
        # the buffer caches (file, scan); another file's scan of the same index must decode again
        other = loadMassSpecData(write_test_scxs(joinpath(mktempdir(), "b.scxs"), [(UInt32[8_000_000], UInt32[5])]))
        getPeaks!(buf, d, 1)
        @test length(first(getPeaks!(buf, other, 1))) == 1
        @test_throws BoundsError getPeaks!(buf, d, 5)
        @test_throws ErrorException getMzArray(d, 1)
    end

    dir = get(ENV, "PIONEER_TEST_SCXS", "")
    if !isempty(dir)
        @testset "real data: .scxs equals its Arrow" begin
            scxs = only(filter(p -> endswith(p, ".scxs"), readdir(dir; join = true)))
            s = loadMassSpecData(scxs); a = loadMassSpecData(splitext(scxs)[1] * ".arrow")
            @test length(s) == length(a)
            @test getPeakCounts(s) == Pioneer.getPeakCounts(a)
            for f in (getRetentionTimes, Pioneer.getLowMzs, Pioneer.getHighMzs, getTICs, getMsOrders)
                @test collect(f(s)) == collect(f(a))
            end
            @test all(isequal.(collect(getCenterMzs(s)), collect(getCenterMzs(a))))
            buf = PeakDecodeBuffer()
            @test all(i -> (p = getPeaks!(buf, s, i); p[1] == getMzArray(a, i) && p[2] == getIntensityArray(a, i)), 1:length(s))
        end
    end
end
