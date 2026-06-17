using Test
using Pioneer
using DataFrames

function _test_fragment_rank_weights(n_frags::Integer)
    weights = [1.0f0 / sqrt(Float32(r)) for r in 1:n_frags]
    scale = Float32(n_frags) / sum(weights)
    return Float32.(weights .* scale)
end

@testset "fragment chromatogram positive correlation summaries" begin
    psms = DataFrame(
        precursor_idx = UInt32[10, 10, 10],
        frag1_int = Float32[1, 2, 3],
        frag2_int = Float32[2, 4, 6],
        frag3_int = Float32[3, 2, 1],
        frag4_int = Float32[0, 0, 0],
        frag5_int = Float32[1, 1, 1],
        frag6_int = Float32[0, 2, 4],
        weight = Float32[1, 2, 3],
        irt_obs = Float32[0, 1, 2],
        ms1_m0_intensity = Float32[1, 2, 3],
    )

    Pioneer._add_fragment_chromatogram_features!(psms)

    rank_weights = _test_fragment_rank_weights(6)
    expected_strength = rank_weights[1] + rank_weights[2] + rank_weights[6]
    expected_effective_n = expected_strength^2 /
        (rank_weights[1]^2 + rank_weights[2]^2 + rank_weights[6]^2)

    @test psms.n_correlated_fragments == UInt8[3, 3, 3]
    @test psms.frag_corr_strength ≈ fill(expected_strength, 3)
    @test psms.frag_corr_effective_n ≈ fill(expected_effective_n, 3)
end
