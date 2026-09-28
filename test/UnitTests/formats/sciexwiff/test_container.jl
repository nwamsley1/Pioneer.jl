# Block codec and .scxs container round trips.

@testset "block codec" begin
    c = S.BlockCodec()
    bins = UInt32[5, 6, 100, 70_000, 4_000_000_000]; ints = UInt32[1, 0, 65_536, 7, typemax(UInt32)]
    nb = S.encode_block!(c, bins, ints, 5, 3)
    z = c.zbuf[1:nb]
    b2 = UInt32[]; i2 = UInt32[]
    S.decode_block!(b2, i2, S.BlockCodec(), z, 5)
    @test b2[1:5] == bins && i2[1:5] == ints
    @test S.encode_block!(c, bins, ints, 0, 3) == 0
    @test_throws ErrorException S.decode_block!(b2, i2, S.BlockCodec(), z, 4)   # wrong peak count
end

@testset "scxs container" begin
    dir = joinpath(mktempdir(), "t.scxs")
    w = S.ScxsWriter(dir, Dict{String, Any}("bin_scale" => 4, "int_scale" => 10.0))
    c = S.BlockCodec()
    rows = [(record = Int32(r), cycle = Int32(1), experiment = Int16(r), ms_order = UInt8(r == 1 ? 1 : 2),
             retention_time = 1f0, low_mz = 100f0, high_mz = 1500f0, tic = 10f0,
             center_mz = r == 1 ? NaN32 : 500f0, isolation_width = r == 1 ? NaN32 : 3f0,
             cal_a = 4.9e-4, cal_b = -13.7, n_peaks = Int32(r == 2 ? 0 : 3)) for r in 1:3]
    peaks = [(UInt32[10, 20, 30], UInt32[1, 2, 3]), (UInt32[], UInt32[]), (UInt32[7, 8, 9], UInt32[4, 5, 6])]
    for (row, (b, i)) in zip(rows, peaks)
        nb = S.encode_block!(c, b, i, length(b), 3)
        S.write_scan!(w, row, c.zbuf[1:nb])
    end
    close(w)
    f = open_scxs(dir)
    @test S.n_scans(f) == 3 && f.meta["n_peaks"] == 6
    b = UInt32[]; i = UInt32[]
    for s in 1:3
        np = S.read_block!(b, i, c, f, s)
        @test np == length(peaks[s][1]) && b[1:np] == peaks[s][1] && i[1:np] == peaks[s][2]
    end
    @test S.stored_mz(f, 1, UInt32(4_000_000)) ≈ S.bin_to_mz(1_000_000.0, 4.9e-4, -13.7)
end

if !isempty(DATA_DIR)
    @testset "real data: convert, both outputs agree" begin
        wiff = first(filter(endswith(".wiff"), readdir(DATA_DIR; join = true)))
        out = mktempdir()
        r = S.convert_run(wiff, out; log = devnull)
        f = open_scxs(r.scxs)
        t = S.Arrow.Table(r.arrow)
        @test S.n_scans(f) == length(t.msOrder)
        b = UInt32[]; i = UInt32[]; c = S.BlockCodec(); bad = 0
        for s in 1:S.n_scans(f)
            np = S.read_block!(b, i, c, f, s)
            bad += !(np == length(t.mz_array[s]) &&
                     all(k -> Float32(S.stored_mz(f, s, b[k])) == t.mz_array[s][k], 1:np))
        end
        @test bad == 0
    end
end
