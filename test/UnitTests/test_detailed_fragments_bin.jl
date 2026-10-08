# detailed_fragments.bin: raw, memory-mapped fragment table (src/structs/SpectralLibrary/detailed_fragments_bin.jl).

using Test
using Pioneer

@testset "detailed_fragments.bin" begin
    dir = mktempdir()
    spline(i, N) = Pioneer.SplineCompactFrag(UInt32(i), Float32(100 + i), ntuple(k -> Float32(i * k), N),
                                             isodd(i), iseven(i), false, false, 0x01, UInt8(i % 60), 0x02,
                                             UInt8(i % 31), 0x00)
    compact(i) = Pioneer.CompactFrag(UInt32(i), Float32(100 + i), Float16(0.5), true, false, false, false,
                                     0x01, UInt8(i % 60), 0x02, UInt8(i % 31), 0x01)
    ranges = UInt64[1, 4, 4, 11]                      # 3 precursors, the second without fragments
    for (name, frags) in (("spline4", [spline(i, 4) for i in 1:10]),
                          ("spline3", [spline(i, 3) for i in 1:10]),
                          ("compact", [compact(i) for i in 1:10]))   # CompactFrag has padding bytes
        a = joinpath(dir, "$name-a.bin"); b = joinpath(dir, "$name-b.bin")
        Pioneer.write_detailed_frags(a, frags, ranges)
        Pioneer.write_detailed_frags(b, copy(frags), copy(ranges))
        f, r = Pioneer.mmap_detailed_frags(a)
        @test f isa Vector{eltype(frags)}
        @test f == frags
        @test r == ranges
        @test read(a) == read(b)                       # bytes depend only on the contents
        @test_throws ReadOnlyMemoryError (f[1] = f[2])
    end
    # streamed in batches == written at once
    frags = [spline(i, 4) for i in 1:10]
    w = Pioneer.DetailedFragsWriter{eltype(frags)}(joinpath(dir, "streamed.bin"))
    Pioneer.append_frags!(w, frags[1:3]); Pioneer.append_frags!(w, view(frags, 4:10))
    Pioneer.finish_detailed_frags!(w, ranges)
    Pioneer.write_detailed_frags(joinpath(dir, "whole.bin"), frags, ranges)
    @test read(joinpath(dir, "streamed.bin")) == read(joinpath(dir, "whole.bin"))
    # a library directory: .bin preferred, legacy .jls pair otherwise
    lib = mkpath(joinpath(dir, "lib.poin"))
    Pioneer.serialize_to_jls(joinpath(lib, "detailed_fragments.jls"), frags)
    Pioneer.serialize_to_jls(joinpath(lib, "precursor_to_fragment_indices.jls"), ranges)
    f, r = Pioneer.load_detailed_frags_and_ranges(lib)
    @test f == frags && r == ranges
    Pioneer.write_detailed_frags(joinpath(lib, "detailed_fragments.bin"), frags, ranges)
    f, r = Pioneer.load_detailed_frags_and_ranges(lib)
    @test f == frags && r == ranges
    @test_throws ReadOnlyMemoryError (f[1] = f[2])     # came from the .bin
    # not a detailed_fragments.bin
    write(joinpath(dir, "junk.bin"), zeros(UInt8, 128))
    @test_throws ErrorException Pioneer.mmap_detailed_frags(joinpath(dir, "junk.bin"))
end
