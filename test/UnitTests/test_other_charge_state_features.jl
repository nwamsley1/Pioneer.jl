using Test
using Pioneer

@testset "other charge state features stay disabled" begin
    for feature in (
        :other_charge_present,
        :other_charge_irt_apex_delta,
        :other_charge_log2_ms2_weight_ratio,
        :other_charge_log2_ms1_m0_ratio,
        :other_charge_ms1_pair_present,
        :other_charge_ratio_abs_delta,
    )
        @test !(feature in Pioneer.PRESCORE_FEATURES)
        @test !(feature in Pioneer.ADVANCED_FEATURE_SET)
    end
end
