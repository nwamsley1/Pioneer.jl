# Synthetic CFB containers for testing the reader without vendor data.

const C = SciexWiff.CFB

putu16!(b, o, v) = (b[o+1:o+2] = reinterpret(UInt8, [htol(UInt16(v))]); b)
putu32!(b, o, v) = (b[o+1:o+4] = reinterpret(UInt8, [htol(UInt32(v))]); b)

function dir_entry(name, type; left = C.NOSTREAM, right = C.NOSTREAM, child = C.NOSTREAM,
        start = C.ENDOFCHAIN, size = 0)
    e = zeros(UInt8, 128)
    u = transcode(UInt16, name)
    for (k, c) in enumerate(u)
        putu16!(e, 2(k - 1), c)
    end
    putu16!(e, 64, 2 * (length(u) + 1))
    e[67] = type
    e[68] = 0x01
    putu32!(e, 68, left); putu32!(e, 72, right); putu32!(e, 76, child)
    putu32!(e, 116, start); putu32!(e, 120, size)
    e
end

"""
Version-3 CFB (512-byte sectors) holding storage `S` with a 5000-byte stream `big`
(regular sectors) and a 100-byte stream `small` (mini stream).
Sectors: 0 FAT, 1 directory, 2 mini FAT, 3 mini stream, 4..13 `big`.
"""
function synthetic_cfb(big, small)
    hdr = zeros(UInt8, 512)
    hdr[1:8] .= collect(C.MAGIC)
    putu16!(hdr, 0x18, 0x3e); putu16!(hdr, 0x1a, 3); putu16!(hdr, 0x1c, 0xfffe)
    putu16!(hdr, 0x1e, 9); putu16!(hdr, 0x20, 6)
    putu32!(hdr, 0x2c, 1)             # one FAT sector
    putu32!(hdr, 0x30, 1)             # directory starts at sector 1
    putu32!(hdr, 0x38, 4096)
    putu32!(hdr, 0x3c, 2); putu32!(hdr, 0x40, 1)
    putu32!(hdr, 0x44, C.ENDOFCHAIN); putu32!(hdr, 0x48, 0)
    putu32!(hdr, 0x4c, 0)
    for k in 1:108
        putu32!(hdr, 0x4c + 4k, C.FREESECT)
    end

    fat = fill(C.FREESECT, 128)
    fat[1] = C.FATSECT
    fat[2] = fat[3] = fat[4] = C.ENDOFCHAIN
    for s in 4:12
        fat[s+1] = s + 1
    end
    fat[14] = C.ENDOFCHAIN
    fatsec = collect(reinterpret(UInt8, htol.(fat)))

    dirsec = vcat(
        dir_entry("Root Entry", C.TYPE_ROOT; child = 1, start = 3, size = 128),
        dir_entry("S", C.TYPE_STORAGE; child = 2),
        dir_entry("big", C.TYPE_STREAM; right = 3, start = 4, size = length(big)),
        dir_entry("small", C.TYPE_STREAM; start = 0, size = length(small)))

    minifat = fill(C.FREESECT, 128)
    minifat[1] = 1
    minifat[2] = C.ENDOFCHAIN
    minifatsec = collect(reinterpret(UInt8, htol.(minifat)))

    ministream = zeros(UInt8, 512)
    ministream[1:length(small)] .= small

    bigsecs = zeros(UInt8, 10 * 512)
    bigsecs[1:length(big)] .= big

    vcat(hdr, fatsec, dirsec, minifatsec, ministream, bigsecs)
end

@testset "CFB" begin
    big = UInt8[(k * 7) % 251 for k in 1:5000]
    small = UInt8[(k * 13) % 256 for k in 1:100]
    cf = C.CompoundFile(synthetic_cfb(big, small))

    @test C.stream_paths(cf) == ["S/big", "S/small"]
    @test C.has_stream(cf, "S/big") && C.has_stream(cf, "/S/small") && C.has_stream(cf, "S\\small")
    @test !C.has_stream(cf, "S")
    @test C.read_stream(cf, "S/big") == big
    @test C.read_stream(cf, "S/small") == small
    @test C.stream_size(cf, "S/big") == 5000
    @test_throws KeyError C.read_stream(cf, "S/missing")

    bad = synthetic_cfb(big, small); bad[1] = 0x00
    @test_throws C.CFBError C.CompoundFile(bad)
    @test_throws C.CFBError C.CompoundFile(zeros(UInt8, 100))

    # a FAT cycle in `big` must raise, not hang
    cyc = synthetic_cfb(big, small)
    putu32!(cyc, 512 + 4 * 13, 4)
    @test_throws C.CFBError C.read_stream(C.CompoundFile(cyc), "S/big")

    # truncated file: big stream runs past EOF
    @test_throws C.CFBError C.read_stream(C.CompoundFile(synthetic_cfb(big, small)[1:4096]), "S/big")
end

if !isempty(DATA_DIR)
    @testset "CFB on real .wiff" begin
        for f in filter(endswith(".wiff"), readdir(DATA_DIR; join = true))
            cf = C.CompoundFile(f)
            @test C.has_stream(cf, "SampleSubtree/Sample1/Idx")
            @test (C.stream_size(cf, "SampleSubtree/Sample1/Idx") - 32) % 54 == 0
        end
    end
end
