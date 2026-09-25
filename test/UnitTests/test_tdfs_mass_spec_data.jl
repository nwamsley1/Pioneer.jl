# TdfsMassSpecData against BasicMassSpecData on the same slices: TimsSlices writes a .tdfs and its expanded Arrow
# from the same quantised slices, so every getter must agree scan for scan, in any access order, from any task.
# Peaks are decoded into caller-owned PeakDecodeBuffers; the task tests check that a task's peaks survive yields
# and thread migration while other tasks decode.
#
# Needs the HeLa timsTOF bundle (PRIDE PXD027359) at TIMSSLICES_TEST_DATA (default ~/BrukerTims/pride); skipped
# when it is absent. Run: julia --project=. test/UnitTests/test_tdfs_mass_spec_data.jl

using Test, Random
using TimsSlices
using Pioneer
using Pioneer: BasicMassSpecData, TdfsMassSpecData, PeakDecodeBuffer, getPeaks!, loadMassSpecData, is_ms_data_path,
               getMzArray, getIntensityArray, getRetentionTime, getLowMz, getHighMz, getTIC, getCenterMz,
               getIsolationWidthMz, getMsOrder, getCycleIdx, getCollisionEnergyEv, getPeakCount, getPeakCounts,
               getRetentionTimes, getTICs, getCenterMzs, getIsolationWidthMzs, getMsOrders, getCycleIdxs,
               getImScans, getImSlope, getFrameIds, getCollisionEnergyEvs, getMzArrays, ArrowTableReference, getMSData

const TDFS_TEST_DATA = get(ENV, "TIMSSLICES_TEST_DATA", expanduser("~/BrukerTims/pride"))
const TDFS_HELA = joinpath(TDFS_TEST_DATA, "20210510_TIMS03_EVO03_PaSk_SA_HeLa_50ng_5_6min_DIA_high_speed_S1-B2_1_25186.d")

if !isdir(TDFS_HELA)
    @warn "HeLa .d bundle not found at $TDFS_HELA; skipping TdfsMassSpecData tests"
else
@testset "TdfsMassSpecData" begin
    out_dir = mktempdir()
    # 120 frames (13 MS1 + 107 MS2), integer bins, default MS2 peak cap
    params = TimsSlices.ConvertParams(format = :both, frames = collect(1:120))
    paths = TimsSlices.convert(TDFS_HELA, out_dir; params = params, name = "hela120", log = devnull)
    a = BasicMassSpecData(paths.arrow)
    t = TdfsMassSpecData(paths.tdfs)
    n = length(a)
    @test length(t) == n > 1000
    @test is_ms_data_path(paths.tdfs) && is_ms_data_path(paths.arrow) && !is_ms_data_path(out_dir)
    @test loadMassSpecData(paths.tdfs) isa TdfsMassSpecData && loadMassSpecData(paths.arrow) isa Pioneer.NonIonMobilityData
    # the instrument's scan-to-1/K0 slope comes from the .tdfs calibration; Arrow data carries none
    @test getImSlope(t) == Float32(abs(t.file.meta["im_slope_1overK0_per_scan"])) > 0
    @test getImSlope(a) === nothing

    same_peaks(buf, i) = begin
        mz, it = getPeaks!(buf, t, i)
        isequal(mz, getMzArray(a, i)) && isequal(it, getIntensityArray(a, i))
    end
    buf = PeakDecodeBuffer()
    same_scan(i) = begin
        getRetentionTime(t, i) == getRetentionTime(a, i) && getLowMz(t, i) == getLowMz(a, i) &&
        getHighMz(t, i) == getHighMz(a, i) && getTIC(t, i) == getTIC(a, i) &&
        isequal(getCenterMz(t, i), getCenterMz(a, i)) && isequal(getIsolationWidthMz(t, i), getIsolationWidthMz(a, i)) &&
        getMsOrder(t, i) == getMsOrder(a, i) && getCycleIdx(t, i) == getCycleIdx(a, i) &&
        getCollisionEnergyEv(t, i) == getCollisionEnergyEv(a, i) && getPeakCount(t, i) == length(getMzArray(a, i)) &&
        same_peaks(buf, i)
    end

    @testset "every scan, in order" begin
        @test all(same_scan, 1:n)
        @test any(i -> getMsOrder(a, i) == 1, 1:n) && any(i -> getMsOrder(a, i) == 2, 1:n)
    end

    @testset "plural getters and timsTOF columns" begin
        @test getRetentionTimes(t) == collect(getRetentionTimes(a))
        @test getTICs(t) == collect(getTICs(a))
        @test isequal(getCenterMzs(t), collect(getCenterMzs(a)))
        @test isequal(getIsolationWidthMzs(t), collect(getIsolationWidthMzs(a)))
        @test getMsOrders(t) == collect(getMsOrders(a))
        @test getCycleIdxs(t) == UInt32.(collect(getCycleIdxs(a)))
        @test getImScans(t) == collect(getImScans(a))
        @test getFrameIds(t) == collect(getFrameIds(a))
        @test getCollisionEnergyEvs(t) == collect(getCollisionEnergyEvs(a))
        @test getPeakCounts(t) == getPeakCounts(a)
        @test_throws ErrorException getMzArrays(t)
        # no implicit shared buffer: per-scan peaks only through getPeaks!
        @test_throws ErrorException getMzArray(t, 1)
        @test_throws ErrorException getIntensityArray(t, 1)
        # Arrow data ignores the buffer and returns its own arrays
        @test getPeaks!(PeakDecodeBuffer(), a, 5) == (getMzArray(a, 5), getIntensityArray(a, 5))
    end

    @testset "random order, and the view contract" begin
        rng = MersenneTwister(7)
        order = shuffle(rng, 1:n)
        @test all(same_scan, order[1:min(end, 3000)])
        # a view is valid until its buffer decodes a different scan: after that it shows the other scan
        i, j = findfirst(i -> getPeakCount(t, i) > 0, 1:n), findlast(i -> getPeakCount(t, i) > 0, 1:n)
        b = PeakDecodeBuffer()
        vi, _ = getPeaks!(b, t, i)
        @test vi[1] == getMzArray(a, i)[1]
        vj, _ = getPeaks!(b, t, j)
        @test vj[1] == getMzArray(a, j)[1]
        @test vi[1] == vj[1]              # vi now aliases the buffer holding scan j
        # a second buffer is independent of the first
        b2 = PeakDecodeBuffer()
        vi2, _ = getPeaks!(b2, t, i)
        @test vi2[1] == getMzArray(a, i)[1] && vj[1] == getMzArray(a, j)[1]
    end

    @testset "buffer cache is keyed on (file, scan), not scan alone" begin
        # A second file whose scan i holds different peaks: a buffer that just decoded scan i of `t` must decode again.
        p2 = TimsSlices.convert(TDFS_HELA, joinpath(out_dir, "second"); params = TimsSlices.ConvertParams(
            format = :both, frames = collect(121:240)), name = "hela121", log = devnull)
        t2 = TdfsMassSpecData(p2.tdfs); a2 = BasicMassSpecData(p2.arrow)
        i = findfirst(k -> getPeakCount(t, k) > 0 && getPeakCount(t2, k) > 0 &&
                           getMzArray(a, k) != getMzArray(a2, k), 1:min(n, length(t2)))
        @test i !== nothing
        b = PeakDecodeBuffer()
        m1, _ = getPeaks!(b, t, i);  @test isequal(m1, getMzArray(a, i))
        m2, _ = getPeaks!(b, t2, i); @test isequal(m2, getMzArray(a2, i))
    end

    @testset "task stress: each task owns a buffer; peaks survive yields" begin
        # More tasks than threads, and a yield between fetching and checking, so tasks interleave on threads and
        # migrate between them while holding views. With one buffer per task nothing else can overwrite them.
        nt = Threads.nthreads()
        ok = Threads.Atomic{Int}(0); bad = Threads.Atomic{Int}(0)
        tasks = [Threads.@spawn begin
            rng = MersenneTwister(k); b = PeakDecodeBuffer()
            for _ in 1:2000
                i = rand(rng, 1:n)
                mz, it = getPeaks!(b, t, i)
                yield()
                good = isequal(mz, getMzArray(a, i)) && isequal(it, getIntensityArray(a, i))
                Threads.atomic_add!(good ? ok : bad, 1)
            end
        end for k in 1:4nt]
        foreach(wait, tasks)
        @test bad[] == 0 && ok[] == 2000 * 4nt
    end

    @testset "reference and getMSData" begin
        ms_dir = joinpath(out_dir, "ms"); mkpath(ms_dir)
        mv(paths.tdfs, joinpath(ms_dir, "hela120.tdfs"))
        ref = ArrowTableReference(ms_dir)
        @test length(ref) == 1 && endswith(ref.file_paths[1], "hela120.tdfs")
        @test getMSData(ref, 1) isa TdfsMassSpecData
        ref2 = ArrowTableReference([joinpath(ms_dir, "hela120.tdfs"), paths.arrow, joinpath(out_dir, "nothing.txt")])
        @test length(ref2) == 2 && ref2.file_id_to_name == ["hela120", "hela120"]
    end
    rm(out_dir; recursive = true)
end
end
