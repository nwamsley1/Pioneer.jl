using Test
using Pioneer

@testset "MS2 cycle indexing" begin
    @testset "scored PSM schema keeps cycle separate from file index" begin
        main_fields = fieldnames(Pioneer.MainSearchScoredPSM{Float32, Float32})
        tuning_fields = fieldnames(Pioneer.TuningScoredPSM{Float32, Float32})

        @test :ms_file_idx in main_fields
        @test :cycle_idx in main_fields
        @test :scan_idx in main_fields
        @test findfirst(==(:ms_file_idx), main_fields) <
              findfirst(==(:cycle_idx), main_fields) <
              findfirst(==(:scan_idx), main_fields)

        @test :ms_file_idx in tuning_fields
        @test :cycle_idx in tuning_fields
        @test :scan_idx in tuning_fields
        @test findfirst(==(:ms_file_idx), tuning_fields) <
              findfirst(==(:cycle_idx), tuning_fields) <
              findfirst(==(:scan_idx), tuning_fields)
    end

    @testset "cycle advances when the MS2 method returns to its first window" begin
        anchor = (Int32(0), Int32(0))
        has_anchor = false
        cycle_idx = 0

        function advance(ms_order, key)
            cycle_idx, anchor, has_anchor = Pioneer._advance_ms2_cycle_idx(
                cycle_idx,
                anchor,
                has_anchor,
                Int32(ms_order),
                key,
            )
            return cycle_idx, anchor, has_anchor
        end

        @test advance(1, (Int32(0), Int32(0))) == (0, (Int32(0), Int32(0)), false)
        @test advance(2, (Int32(500000), Int32(24000))) ==
              (1, (Int32(500000), Int32(24000)), true)
        @test advance(1, (Int32(0), Int32(0))) ==
              (1, (Int32(500000), Int32(24000)), true)
        @test advance(2, (Int32(525000), Int32(24000))) ==
              (1, (Int32(500000), Int32(24000)), true)
        @test advance(2, (Int32(550000), Int32(24000))) ==
              (1, (Int32(500000), Int32(24000)), true)
        @test advance(2, (Int32(500000), Int32(24000))) ==
              (2, (Int32(500000), Int32(24000)), true)
        @test advance(2, (Int32(525000), Int32(24000))) ==
              (2, (Int32(500000), Int32(24000)), true)
    end
end
