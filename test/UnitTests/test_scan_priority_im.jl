# Tests for the ion-mobility (packet) scan priority order used by ParameterTuningSearch:
# get_ms2_scan_priority_order_im stratifies each RT bin over (isolation window, IM bin) cells,
# and getImScans detects the `imScan` column.

using Test
using Arrow
using DataFrames
using Pioneer: BasicMassSpecData, getImScans, get_ms2_scan_priority_order,
               get_ms2_scan_priority_order_im, SCAN_PRIORITY_N_RT_BINS,
               getMsOrder, getCenterMzs, getRetentionTimes

# Minimal packet-style Arrow file: `n_frames` frames, each with one MS1 packet row followed by
# MS2 packets for every (window, scan) pair. All rows of a frame share the frame's RT.
function write_packet_arrow(path::String; n_frames::Int, windows::Vector{Float32}, scans::AbstractVector{<:Integer},
                            with_im::Bool = true, tic_of = (f, w, s) -> 1.0f0)
    rows_per_frame = 1 + length(windows) * length(scans)
    n = n_frames * rows_per_frame
    mz = Vector{Vector{Union{Missing,Float32}}}(undef, n); it = similar(mz)
    rt = Vector{Float32}(undef, n); ms_order = Vector{UInt8}(undef, n); tic = Vector{Float32}(undef, n)
    cmz = Vector{Union{Missing,Float32}}(undef, n); width = Vector{Union{Missing,Float32}}(undef, n)
    im = Vector{UInt16}(undef, n); frame = Vector{Int32}(undef, n)
    r = 0
    for f in 1:n_frames
        r += 1
        mz[r] = Union{Missing,Float32}[400.0f0]; it[r] = Union{Missing,Float32}[1.0f0]
        rt[r] = Float32(f); ms_order[r] = 0x01; tic[r] = 1.0f0; cmz[r] = missing; width[r] = missing
        im[r] = 0x0000; frame[r] = Int32(f)
        for (wi, w) in enumerate(windows), s in scans
            r += 1
            mz[r] = Union{Missing,Float32}[500.0f0]; it[r] = Union{Missing,Float32}[1.0f0]
            rt[r] = Float32(f); ms_order[r] = 0x02; tic[r] = tic_of(f, wi, s); cmz[r] = w; width[r] = 25.0f0
            im[r] = UInt16(s); frame[r] = Int32(f)
        end
    end
    df = DataFrame(mz_array = mz, intensity_array = it, scanHeader = fill("", n), scanNumber = Int32.(1:n),
                   packetType = zeros(Int32, n), retentionTime = rt, lowMz = fill(100.0f0, n), highMz = fill(1700.0f0, n),
                   TIC = tic, centerMz = cmz, isolationWidthMz = width, collisionEnergyField = width,
                   collisionEnergyEvField = zeros(Float32, n), msOrder = ms_order, cycle_idx = Int32.(frame))
    if with_im
        df[!, :frameId] = frame; df[!, :imScan] = im
    end
    Arrow.write(path, df)
    return path
end

@testset "getImScans detects the packet column" begin
    mktempdir() do d
        p = write_packet_arrow(joinpath(d, "im.arrow"); n_frames = 2, windows = Float32[412.5, 612.5], scans = 0:9)
        ms = BasicMassSpecData(p)
        @test getImScans(ms) !== nothing
        @test length(getImScans(ms)) == length(ms)
        q = write_packet_arrow(joinpath(d, "noim.arrow"); n_frames = 2, windows = Float32[412.5], scans = 0:9, with_im = false)
        @test getImScans(BasicMassSpecData(q)) === nothing
    end
end

@testset "get_ms2_scan_priority_order_im — permutation of MS2 rows, cell rotation" begin
    mktempdir() do d
        windows = Float32[412.5, 612.5, 812.5]
        scans = 0:79                       # 8 IM bins of 10 scans each
        n_im = 8
        # TIC rises with scan within a cell so the highest-TIC packet of each cell is its last scan;
        # window 3 is made much denser to check it cannot dominate the early draws.
        tic_of = (f, w, s) -> Float32(1000 * (w == 3 ? 100 : 1) + s)
        p = write_packet_arrow(joinpath(d, "im.arrow"); n_frames = 1, windows = windows, scans = scans, tic_of = tic_of)
        ms = BasicMassSpecData(p)
        order = get_ms2_scan_priority_order_im(ms, n_im)
        ms2 = [i for i in 1:length(ms) if getMsOrder(ms, i) == 0x02]
        @test sort(order) == ms2                       # every MS2 row exactly once, no MS1 rows
        @test eltype(order) == Int32

        ims = getImScans(ms); cmz = getCenterMzs(ms)
        cellof(i) = (cmz[i], div(Int(ims[i]), 10))     # (window, IM bin)
        n_cells = length(windows) * n_im
        first_cells = [cellof(i) for i in order[1:n_cells]]
        @test length(unique(first_cells)) == n_cells   # one draw per cell before any cell repeats
        # each first draw is the highest-TIC packet of its cell (scan 9 of the bin)
        @test all(Int(ims[i]) % 10 == 9 for i in order[1:n_cells])
        # the dense window supplies exactly a third of the first rotation, not all of it
        @test count(cmz[i] == 812.5f0 for i in order[1:n_cells]) == n_im
        # second rotation again covers every cell
        @test length(unique(cellof(i) for i in order[n_cells+1:2n_cells])) == n_cells
    end
end

@testset "get_ms2_scan_priority_order_im — round-robin over RT bins" begin
    mktempdir() do d
        n_frames = SCAN_PRIORITY_N_RT_BINS          # one frame per RT bin
        p = write_packet_arrow(joinpath(d, "im.arrow"); n_frames = n_frames, windows = Float32[412.5], scans = 0:15)
        ms = BasicMassSpecData(p)
        order = get_ms2_scan_priority_order_im(ms, 4)
        rts = getRetentionTimes(ms)
        # the first n_frames draws come from n_frames distinct RT values
        @test length(unique(rts[i] for i in order[1:n_frames])) == n_frames
        # the plain order also covers the same MS2 rows (same contract, different ordering)
        @test sort(get_ms2_scan_priority_order(ms)) == sort(order)
    end
end
