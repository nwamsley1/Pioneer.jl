# TdfsMassSpecData against BasicMassSpecData on the same slices: TimsSlices writes a .tdfs and its expanded Arrow
# from the same quantised slices, so every getter must agree scan for scan, in any access order, from any thread.
#
# Needs the HeLa timsTOF bundle (PRIDE PXD027359) at TIMSSLICES_TEST_DATA (default ~/BrukerTims/pride); skipped
# when it is absent. Run: julia --project=. test/UnitTests/test_tdfs_mass_spec_data.jl

using Test, Random
using TimsSlices
using Pioneer
using Pioneer: BasicMassSpecData, TdfsMassSpecData, loadMassSpecData, is_ms_data_path,
               getMzArray, getIntensityArray, getRetentionTime, getLowMz, getHighMz, getTIC, getCenterMz,
               getIsolationWidthMz, getMsOrder, getCycleIdx, getCollisionEnergyEv, getPeakCount, getPeakCounts,
               getRetentionTimes, getTICs, getCenterMzs, getIsolationWidthMzs, getMsOrders, getCycleIdxs,
               getImScans, getFrameIds, getCollisionEnergyEvs, getMzArrays, ArrowTableReference, getMSData

const TDFS_TEST_DATA = get(ENV, "TIMSSLICES_TEST_DATA", expanduser("~/BrukerTims/pride"))
const TDFS_HELA = joinpath(TDFS_TEST_DATA, "20210510_TIMS03_EVO03_PaSk_SA_HeLa_50ng_5_6min_DIA_high_speed_S1-B2_1_25186.d")

if !isdir(TDFS_HELA)
    @warn "HeLa .d bundle not found at $TDFS_HELA; skipping TdfsMassSpecData tests"
else
@testset "TdfsMassSpecData" begin
    out_dir = mktempdir()
    # 120 frames (13 MS1 + 107 MS2), integer bins, per-level cull so both branches of the metadata are exercised
    params = TimsSlices.ConvertParams(format = :both, frames = collect(1:120), cull_q = 0.01, ms1_cull_q = 0.0)
    paths = TimsSlices.convert(TDFS_HELA, out_dir; params = params, name = "hela120", log = devnull)
    a = BasicMassSpecData(paths.arrow)
    t = TdfsMassSpecData(paths.tdfs)
    n = length(a)
    @test length(t) == n > 1000
    @test is_ms_data_path(paths.tdfs) && is_ms_data_path(paths.arrow) && !is_ms_data_path(out_dir)
    @test loadMassSpecData(paths.tdfs) isa TdfsMassSpecData && loadMassSpecData(paths.arrow) isa Pioneer.NonIonMobilityData

    same_scan(i) = begin
        getRetentionTime(t, i) == getRetentionTime(a, i) && getLowMz(t, i) == getLowMz(a, i) &&
        getHighMz(t, i) == getHighMz(a, i) && getTIC(t, i) == getTIC(a, i) &&
        isequal(getCenterMz(t, i), getCenterMz(a, i)) && isequal(getIsolationWidthMz(t, i), getIsolationWidthMz(a, i)) &&
        getMsOrder(t, i) == getMsOrder(a, i) && getCycleIdx(t, i) == getCycleIdx(a, i) &&
        getCollisionEnergyEv(t, i) == getCollisionEnergyEv(a, i) && getPeakCount(t, i) == length(getMzArray(a, i)) &&
        isequal(getMzArray(t, i), getMzArray(a, i)) && isequal(getIntensityArray(t, i), getIntensityArray(a, i))
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
    end

    @testset "random order, and the view contract" begin
        rng = MersenneTwister(7)
        order = shuffle(rng, 1:n)
        @test all(same_scan, order[1:min(end, 3000)])
        # a view is valid until this thread fetches a different scan: after that it shows the other scan
        i, j = findfirst(i -> getPeakCount(t, i) > 0, 1:n), findlast(i -> getPeakCount(t, i) > 0, 1:n)
        vi = getMzArray(t, i)
        @test vi[1] == getMzArray(a, i)[1]
        vj = getMzArray(t, j)
        @test vj[1] == getMzArray(a, j)[1]
        @test vi[1] == vj[1]              # vi now aliases the scratch holding scan j
        # fetching the same scan twice does not decode twice: mz then intensity share the decode
        vm = getMzArray(t, i); vint = getIntensityArray(t, i)
        @test isequal(vm, getMzArray(a, i)) && isequal(vint, getIntensityArray(a, i))
    end

    @testset "thread stress: each task compares its views immediately" begin
        nt = Threads.nthreads()
        ok = Threads.Atomic{Int}(0); bad = Threads.Atomic{Int}(0)
        tasks = [Threads.@spawn begin
            rng = MersenneTwister(k)
            for _ in 1:2000
                i = rand(rng, 1:n)
                good = isequal(getMzArray(t, i), getMzArray(a, i)) && isequal(getIntensityArray(t, i), getIntensityArray(a, i))
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
