using Test
using Pioneer

@testset "top-8 fragment traces sum contributing isotope intensities" begin
    frag = Pioneer.CompactFrag(
        UInt32(1), 200.0f0, Float16(1000),
        true, false, false, false,
        UInt8(1), UInt8(4), UInt8(2), UInt8(1), UInt8(0))
    unscored = [Pioneer.MainUnscoredPSM{Float32}()]

    Pioneer.apply_main_scoring!(
        unscored, 1, frag, UInt8(0), 100.0f0, 0.0f0, 3, UInt32(1), 1.0f0, 10.0f0, true)
    Pioneer.apply_main_scoring!(
        unscored, 1, frag, UInt8(1), 50.0f0, 0.0f0, 3, UInt32(1), 0.3f0, 10.0f0, true)
    Pioneer.apply_main_scoring!(
        unscored, 1, frag, UInt8(2), 80.0f0, 0.0f0, 3, UInt32(1), 0.1f0, 10.0f0, false)

    psm = unscored[1]
    @test psm.frag1_int === 150.0f0
    @test psm.isotope_count == UInt8(2)
    @test psm.y_count_iso == UInt8(2)
    @test psm.y_count == UInt8(1)
    @test psm.pred_int_sum_m0 === 1.0f0
end

@testset "expected predicted missing fraction uses M0 mass-error risk" begin
    frag = Pioneer.CompactFrag(
        UInt32(1), 250.0f0, Float16(1000),
        true, false, false, false,
        UInt8(1), UInt8(4), UInt8(2), UInt8(1), UInt8(0))
    unscored = [Pioneer.MainUnscoredPSM{Float32}()]

    Pioneer.apply_main_scoring!(
        unscored, 1, frag, UInt8(0), 100.0f0, 5.0f0, 3, UInt32(1), 4.0f0, 10.0f0, true)
    Pioneer.apply_main_scoring!(
        unscored, 1, frag, UInt8(0), 100.0f0, 20.0f0, 3, UInt32(1), 6.0f0, 10.0f0, true)
    Pioneer.apply_main_scoring!(
        unscored, 1, frag, UInt8(1), 100.0f0, 10.0f0, 3, UInt32(1), 100.0f0, 10.0f0, true)

    scratch = Pioneer.FusedScratch(2)
    Pioneer.push_match!(scratch, 1, 4.0f0, 100.0f0, 0, 1, 0.25f0)
    Pioneer.push_match!(scratch, 2, 6.0f0, 100.0f0, 0, 1, 1.0f0)
    Pioneer.record_mass_missing_hellinger!(
        Pioneer.FusedStandard(Pioneer.FullPrecCapture(), UInt8(7)), unscored, 1, scratch)

    psm = unscored[1]
    @test psm.pred_int_sum_m0 === 10.0f0
    @test psm.expected_missing_pred_int_sum_m0 ≈ 7.0f0
    @test psm.expected_missing_pred_int_sum_m0 / psm.pred_int_sum_m0 ≈ 0.7f0

    scored = Vector{Pioneer.MainSearchScoredPSM{Float32, Float16}}(undef, 1)
    spectral_scores = [
        Pioneer.SpectralScoresMainSearch{Float16, Float32}(
            zero(Float16), zero(Float16), zero(Float16), zero(Float16), zero(Float16),
            zeros(Float32, 16)...,
        )
    ]
    id_to_col = Pioneer.DensePrecMap{UInt16}(8)
    id_to_col[UInt32(1)] = UInt16(1)
    result = Pioneer.Score!(
        scored, unscored, spectral_scores, Float32[1], id_to_col,
        1, 1, 1.0, 0, 1, 1000.0f0, 42)

    @test result.last_val == 1
    @test scored[1].expected_predicted_missing_fraction ≈ 0.7f0

    pred = Float32[4, 6]
    obs = Float32[100, 100]
    risk = Float32[0.25, 1]
    pred_sum = sum(pred)
    obs_sum = sum(obs)
    bc_sum = sum(sqrt.(pred .* obs))
    base_hsq = 1f0 - bc_sum / sqrt(pred_sum * obs_sum)
    expected_hsq = base_hsq
    for i in eachindex(pred)
        obs_without = obs_sum - obs[i]
        bc_without = bc_sum - sqrt(pred[i] * obs[i])
        hsq_without = obs_without > 0f0 ? 1f0 - bc_without / sqrt(pred_sum * obs_without) : 1f0
        expected_hsq += risk[i] * max(hsq_without - base_hsq, 0f0)
    end
    expected_hsq = min(max(expected_hsq, 0f0), 1f0)
    expected_score = -log2(max(expected_hsq, 1f-10))
    @test scored[1].mass_missing_hellinger ≈ expected_score
end

@testset "rank 8 fragment trace intensity is captured" begin
    frag = Pioneer.CompactFrag(
        UInt32(1), 250.0f0, Float16(1000),
        true, false, false, false,
        UInt8(1), UInt8(8), UInt8(2), UInt8(8), UInt8(0))
    unscored = [Pioneer.MainUnscoredPSM{Float32}()]

    Pioneer.apply_main_scoring!(
        unscored, 1, frag, UInt8(0), 42.0f0, 0.0f0, 3, UInt32(1), 1.0f0, 10.0f0, true)

    @test unscored[1].frag8_int === 42.0f0
end
