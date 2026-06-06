using Test
using Pioneer

@testset "MBR annotation-matched fragment Hellinger" begin
    empty_keys = ntuple(_ -> UInt32(0), 8)
    empty_ints = ntuple(_ -> 0.0f0, 8)

    @testset "matches by annotation rather than rank position" begin
        obs_keys = (UInt32(10), UInt32(20), empty_keys[3:8]...)
        ref_keys = (UInt32(20), UInt32(10), empty_keys[3:8]...)
        obs_ints = (4.0f0, 1.0f0, empty_ints[3:8]...)
        ref_ints = (1.0f0, 4.0f0, empty_ints[3:8]...)

        @test Pioneer._mbr_annotation_hellinger(obs_keys, obs_ints, ref_keys, ref_ints) ≈ 0.0f0
    end

    @testset "uses the union of annotations with zero-filled missing labels" begin
        obs_keys = (UInt32(10), empty_keys[2:8]...)
        ref_keys = (UInt32(20), empty_keys[2:8]...)
        obs_ints = (10.0f0, empty_ints[2:8]...)
        ref_ints = (10.0f0, empty_ints[2:8]...)

        @test Pioneer._mbr_annotation_hellinger(obs_keys, obs_ints, ref_keys, ref_ints) ≈ 1.0f0
    end

    @testset "extracts compact annotation keys by top-8 fragment rank" begin
        y5 = Pioneer.CompactFrag(
            UInt32(1), 500.0f0, Float16(1.0),
            true, false, false, false,
            UInt8(1), UInt8(5), UInt8(2), UInt8(1), UInt8(0),
        )
        b7 = Pioneer.CompactFrag(
            UInt32(1), 600.0f0, Float16(0.8),
            false, true, false, false,
            UInt8(2), UInt8(7), UInt8(2), UInt8(2), UInt8(0),
        )
        isotope = Pioneer.CompactFrag(
            UInt32(1), 601.0f0, Float16(0.5),
            false, true, false, true,
            UInt8(2), UInt8(7), UInt8(2), UInt8(3), UInt8(0),
        )
        lookup = Pioneer.StandardFragmentLookup{Float32}(
            Pioneer.CompactFrag{Float32}[y5, b7, isotope],
            UInt64[1, 4],
        )

        keys = Pioneer._mbr_top8_fragment_keys(lookup, UInt32(1))

        @test keys[1] == Pioneer._mbr_fragment_annotation_key(y5)
        @test keys[2] == Pioneer._mbr_fragment_annotation_key(b7)
        @test keys[3] == UInt32(0)
    end
end
