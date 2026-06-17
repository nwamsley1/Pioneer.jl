using Test
using Pioneer

@testset "top-6 fragment traces sum contributing isotope intensities" begin
    frag = Pioneer.CompactFrag(
        UInt32(1), 200.0f0, Float16(1000),
        true, false, false, false,
        UInt8(1), UInt8(4), UInt8(2), UInt8(1), UInt8(0))
    unscored = [Pioneer.MainUnscoredPSM{Float32}()]

    Pioneer.apply_main_scoring!(
        unscored, 1, frag, UInt8(0), 100.0f0, 0.0f0, 3, UInt32(1), 1.0f0, true)
    Pioneer.apply_main_scoring!(
        unscored, 1, frag, UInt8(1), 50.0f0, 0.0f0, 3, UInt32(1), 0.3f0, true)
    Pioneer.apply_main_scoring!(
        unscored, 1, frag, UInt8(2), 80.0f0, 0.0f0, 3, UInt32(1), 0.1f0, false)

    psm = unscored[1]
    @test psm.frag1_int === 150.0f0
    @test psm.isotope_count == UInt8(2)
    @test psm.y_count_iso == UInt8(2)
    @test psm.y_count == UInt8(1)
    @test psm.pred_int_sum_m0 === 1.0f0
end

