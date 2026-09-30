# Block decoding and centroiding on synthetic data, plus real-data checks when SCIEXWIFF_DATA is set.

const S = SciexWiff

"Encode peaks (step index m ≥ 1 strictly increasing, intensity) as a .wiff.scan block with metadata."
function synthetic_block(steps, ints; start = 1_000_000, a = 4.9e-4, b = -13.7)
    tok = UInt8[]
    prev = 0
    for (m, v) in zip(steps, ints)
        d = m - prev; prev = m
        if d == 1
        elseif d <= 124; push!(tok, UInt8(0x7f + d))
        elseif d <= 256; push!(tok, 0xfc, UInt8(d - 1))
        elseif d <= 65536; push!(tok, 0xfd, UInt8((d - 1) & 0xff), UInt8((d - 1) >> 8))
        else; push!(tok, 0xfe); append!(tok, reinterpret(UInt8, [UInt32(d - 1)]))
        end
        if v <= 0x7b; push!(tok, UInt8(v))
        elseif v <= 0xff; push!(tok, 0x7c, UInt8(v))
        elseif v <= 0xffff; push!(tok, 0x7d, UInt8(v & 0xff), UInt8(v >> 8))
        else; push!(tok, 0x7e); append!(tok, reinterpret(UInt8, [UInt32(v)]))
        end
    end
    cal = vcat(UInt8[0x0a, 0x12, 0x09], reinterpret(UInt8, [a]), UInt8[0x11], reinterpret(UInt8, [b]))
    seg = UInt8[0x12, 0x06, 0x08]
    v = start                                       # varint start bin
    while v >= 0x80; push!(seg, UInt8(v & 0x7f) | 0x80); v >>= 7; end
    push!(seg, UInt8(v)); push!(seg, 0x10, 0x08)
    seg[2] = UInt8(length(seg) - 2)
    meta = vcat(cal, seg)
    block = vcat(UInt8[0xff, 0xff, 0xff, 0xff], reinterpret(UInt8, [UInt32(start)]), UInt8[0x00], tok)
    pad = 4 - length(block) % 4
    append!(block, fill(0xff, pad))
    file = vcat(zeros(UInt8, S.SCAN_FILE_HEADER), UInt8[length(meta)], meta, block)
    ms = S.SCAN_FILE_HEADER; ff = ms + 1 + length(meta)
    file, Int64(ms), Int64(ff), length(block)
end

@testset "block decode" begin
    steps = [1, 2, 3, 10, 200, 300, 70_000, 70_001, 200_000]
    ints = [5, 123, 124, 255, 256, 65_535, 65_536, 1, 40_000_000]
    file, ms, ff, sz = synthetic_block(steps, ints)
    sb = ScanBuffer()
    S.decode_block!(sb, file, ms, ff, sz)
    @test sb.n == length(steps)
    @test sb.bin[1:sb.n] == UInt32.(1_000_000 .+ 8 .* steps)
    @test sb.intensity[1:sb.n] == UInt32.(ints)
    @test sb.cal_a == 4.9e-4 && sb.cal_b == -13.7
    @test S.bin_to_mz(S.mz_to_bin(512.25, sb.cal_a, sb.cal_b), sb.cal_a, sb.cal_b) ≈ 512.25

    # the last data byte may be 0xff (intensity 255 = `7c ff`), directly before the 0xff padding
    for last in (255, 0xffff, 0x01ff, 0xffffffff)
        f2, ms2, ff2, sz2 = synthetic_block([1, 5, 9], [7, 300, last])
        S.decode_block!(sb, f2, ms2, ff2, sz2)
        @test sb.n == 3 && sb.intensity[1:3] == UInt32[7, 300, last]
    end
    # 0xff inside the data (not followed only by padding) is still an error
    f3, ms3, ff3, sz3 = synthetic_block([1, 2], [3, 4])
    f3[ff3 + 10] = 0xff
    @test_throws S.ScanFormatError S.decode_block!(sb, f3, ms3, ff3, sz3)

    bad = copy(file); bad[ff+1] = 0x00
    @test_throws S.ScanFormatError S.decode_block!(sb, bad, ms, ff, sz)
    @test_throws S.ScanFormatError S.decode_block!(sb, file, ms + 1, ff, sz)     # metadata misaligned
    @test_throws S.ScanFormatError S.decode_block!(sb, file[1:end-8], ms, ff, sz)  # truncated
end

@testset "centroid" begin
    # one Gaussian peak (σ = 1 step) centred between steps, on a sparse single-ion background
    μ = 500.3
    steps = collect(495:506); ints = [round(Int, 5000exp(-0.5 * (s - μ)^2)) for s in steps]
    keep = ints .> 0
    steps = vcat([100, 300], steps[keep], [900]); ints = vcat([40, 42], ints[keep], [41])
    file, ms, ff, sz = synthetic_block(steps, ints)
    sb = ScanBuffer(); S.decode_block!(sb, file, ms, ff, sz)
    cb = CentroidBuffer(); centroid!(cb, sb)
    @test cb.n == 1                                   # single-ion events dropped (min_bins = 2)
    true_mz = S.bin_to_mz(1_000_000 + 8μ, sb.cal_a, sb.cal_b)
    @test abs(cb.mz[1] - true_mz) / true_mz * 1e6 < 0.5
    area = 100 * sum(ints[k] * (S.bin_to_mz(Float64(sb.bin[k]) + 8, sb.cal_a, sb.cal_b) -
                                S.bin_to_mz(Float64(sb.bin[k]), sb.cal_a, sb.cal_b)) for k in 3:sb.n-1 if abs(steps[k] - μ) <= 3.5)
    @test cb.intensity[1] ≈ area rtol = 1e-3

    centroid!(cb, sb, CentroidParams(min_bins = 1))
    @test cb.n == 4
    @test issorted(cb.mz[1:cb.n])

    # two resolved peaks 8 steps apart stay separate
    steps2 = collect(0:20) .+ 1000
    ints2 = [round(Int, 3000exp(-0.5 * (s - 1006)^2) + 2000exp(-0.5 * (s - 1014)^2)) for s in steps2]
    k2 = ints2 .> 0
    file, ms, ff, sz = synthetic_block(steps2[k2], ints2[k2])
    S.decode_block!(sb, file, ms, ff, sz); centroid!(cb, sb)
    @test cb.n == 2
end

if !isempty(DATA_DIR)
    @testset "real data: decode agrees with Idx" begin
        for f in filter(endswith(".wiff"), readdir(DATA_DIR; join = true))
            run = WiffRun(f)
            sb = ScanBuffer(); bpi_ok = bin_ok = n = 0
            for r in 1:length(run)
                read_scan!(sb, run, r); sb.n == 0 && continue
                n += 1
                k = argmax(view(sb.intensity, 1:sb.n))
                bpi_ok += sb.intensity[k] == run.index.base_peak_intensity[r]
                bin_ok += sb.bin[k] == run.index.base_peak_bin[r]
            end
            @test bpi_ok == n
            @test bin_ok >= 0.999n
        end
    end
end

@testset "Idx offsets past 4 GiB" begin
    # three blocks; the third's u32 offset wraps (the .wiff.scan is over 4 GiB)
    rec(off, size) = vcat(reinterpret(UInt8, [UInt32(off), UInt32(size)]), zeros(UInt8, 46))
    # block 2 ends at 44 + 0xfffff000 + 0x2000 = 2^32 + 0x102c; block 3's raw offset 0x1000 is 2^32 + 0x1000 + 44
    idx = vcat(zeros(UInt8, 32), rec(100, 1000), rec(0xfffff000, 0x2000), rec(0x00001000, 64))
    ix = S.ScanIndex(idx)
    @test ix.block_offset == Int64[144, 44 + 0xfffff000, (Int64(1) << 32) + 0x1000 + 44]
    @test ix.meta_start[3] == 44 + 0xfffff000 + 0x2000
end
