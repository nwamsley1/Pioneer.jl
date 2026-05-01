# Unit tests for FDR utilities (src/utils/ML/fdrUtilities.jl).
#
# Tests the core q-value and PEP computation functions that underpin
# all FDR control in the Pioneer scoring pipeline.
#
# Run standalone: julia --project=. test/UnitTests/test_fdr_utilities.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
end

@testset "FDR Utilities" begin

# ═══════════════════════════════════════════════════════════════════════════════
# get_qvalues! — target-decoy FDR estimation
# ═══════════════════════════════════════════════════════════════════════════════

@testset "get_qvalues!" begin

    @testset "perfect separation" begin
        # All targets have higher scores than all decoys
        probs  = Float32[0.9, 0.8, 0.7, 0.3, 0.2, 0.1]
        labels = Bool[true, true, true, false, false, false]
        qvals  = zeros(Float32, 6)
        Pioneer.get_qvalues!(probs, labels, qvals)

        # First 3 are targets with 0 decoys seen -> q=0
        @test qvals[1] == 0.0f0
        @test qvals[2] == 0.0f0
        @test qvals[3] == 0.0f0
        # After first decoy: 1 decoy / 3 targets
        @test qvals[4] ≈ 1.0f0 / 3.0f0
    end

    @testset "monotonic non-increasing" begin
        # Interleaved targets and decoys
        probs  = Float32[0.9, 0.8, 0.7, 0.6, 0.5, 0.4]
        labels = Bool[true, false, true, false, true, false]
        qvals  = zeros(Float32, 6)
        Pioneer.get_qvalues!(probs, labels, qvals)

        # Q-values must be monotonically non-increasing when sorted by score
        order = sortperm(probs, rev=true)
        sorted_qvals = qvals[order]
        for i in 2:length(sorted_qvals)
            @test sorted_qvals[i] >= sorted_qvals[i-1] || sorted_qvals[i] ≈ sorted_qvals[i-1]
        end
    end

    @testset "all targets" begin
        probs  = Float32[0.9, 0.8, 0.7]
        labels = Bool[true, true, true]
        qvals  = zeros(Float32, 3)
        Pioneer.get_qvalues!(probs, labels, qvals)
        @test all(q -> q == 0.0f0, qvals)
    end

    @testset "all decoys" begin
        probs  = Float32[0.9, 0.8, 0.7]
        labels = Bool[false, false, false]
        qvals  = zeros(Float32, 3)
        Pioneer.get_qvalues!(probs, labels, qvals)
        # With 0 targets, decoys/targets → Inf
        @test all(q -> isinf(q), qvals)
    end

    @testset "unsorted input with doSort=true" begin
        # Scores deliberately out of order
        probs  = Float32[0.3, 0.9, 0.1, 0.7]
        labels = Bool[false, true, false, true]
        qvals  = zeros(Float32, 4)
        Pioneer.get_qvalues!(probs, labels, qvals)

        # Highest score (0.9) is target -> q=0
        @test qvals[2] == 0.0f0
        # Lowest score (0.1) is decoy -> should have highest q-value
        @test qvals[3] >= qvals[2]
    end

    @testset "fdr_scale_factor" begin
        # With scale factor 2.0, decoy contribution doubles
        probs  = Float32[0.9, 0.8, 0.7, 0.6]
        labels = Bool[true, true, false, false]
        qvals1 = zeros(Float32, 4)
        qvals2 = zeros(Float32, 4)
        Pioneer.get_qvalues!(probs, labels, qvals1; fdr_scale_factor=1.0f0)
        Pioneer.get_qvalues!(probs, labels, qvals2; fdr_scale_factor=2.0f0)

        # With higher scale factor, q-values should be larger (more conservative)
        # At the first decoy position (score=0.7): 1*factor/2 targets
        @test qvals2[3] > qvals1[3]
        @test qvals2[3] ≈ 2.0f0 * qvals1[3]
    end

    @testset "Float64 types" begin
        probs  = Float64[0.9, 0.5, 0.1]
        labels = Bool[true, false, false]
        qvals  = zeros(Float64, 3)
        Pioneer.get_qvalues!(probs, labels, qvals)
        @test qvals[1] == 0.0
        @test qvals[2] > 0.0
    end

    @testset "single element" begin
        probs  = Float32[0.5]
        labels = Bool[true]
        qvals  = zeros(Float32, 1)
        Pioneer.get_qvalues!(probs, labels, qvals)
        @test qvals[1] == 0.0f0
    end

    @testset "doSort=false preserves order" begin
        # Pre-sorted descending
        probs  = Float32[0.9, 0.8, 0.7]
        labels = Bool[true, false, true]
        qvals1 = zeros(Float32, 3)
        qvals2 = zeros(Float32, 3)
        Pioneer.get_qvalues!(probs, labels, qvals1; doSort=true)
        Pioneer.get_qvalues!(probs, labels, qvals2; doSort=false)
        # Pre-sorted input should give same result either way
        @test qvals1 ≈ qvals2
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# get_PEP! — posterior error probability via isotonic regression
# ═══════════════════════════════════════════════════════════════════════════════

@testset "get_PEP!" begin

    @testset "perfect separation" begin
        scores    = Float32[0.9, 0.8, 0.7, 0.3, 0.2, 0.1]
        is_target = Bool[true, true, true, false, false, false]
        pep       = zeros(Float32, 6)
        Pioneer.get_PEP!(scores, is_target, pep)

        # Targets (high scores) should have low PEP
        @test pep[1] < 0.5f0
        @test pep[2] < 0.5f0
        # Decoys (low scores) should have high PEP
        @test pep[4] > 0.0f0
        @test pep[5] > 0.0f0
    end

    @testset "PEP in [0, 1]" begin
        scores    = Float32[0.9, 0.8, 0.3, 0.7, 0.2, 0.5]
        is_target = Bool[true, false, true, false, true, false]
        pep       = zeros(Float32, 6)
        Pioneer.get_PEP!(scores, is_target, pep)

        @test all(p -> 0.0f0 <= p <= 1.0f0, pep)
    end

    @testset "monotonic: higher score → lower PEP" begin
        n = 20
        scores    = Float32.(range(0.0, 1.0, length=n))
        # Low scores are decoys, high scores are targets
        is_target = Bool[i > n÷2 for i in 1:n]
        pep       = zeros(Float32, n)
        Pioneer.get_PEP!(scores, is_target, pep)

        # PEP should generally decrease with increasing score
        # (isotonic regression ensures this via PAVA)
        order = sortperm(scores, rev=true)
        sorted_pep = pep[order]
        for i in 2:length(sorted_pep)
            @test sorted_pep[i] >= sorted_pep[i-1] - 1f-6  # allow tiny float noise
        end
    end

    @testset "all targets" begin
        scores    = Float32[0.9, 0.8, 0.7]
        is_target = Bool[true, true, true]
        pep       = zeros(Float32, 3)
        Pioneer.get_PEP!(scores, is_target, pep)
        # With no decoys, PEP should be 0
        @test all(p -> p == 0.0f0, pep)
    end

    @testset "empty input" begin
        scores    = Float32[]
        is_target = Bool[]
        pep       = Float32[]
        # Should not error
        Pioneer.get_PEP!(scores, is_target, pep)
        @test length(pep) == 0
    end

    @testset "fdr_scale_factor" begin
        scores    = Float32[0.9, 0.8, 0.7, 0.3, 0.2, 0.1]
        is_target = Bool[true, true, true, false, false, false]
        pep1      = zeros(Float32, 6)
        pep2      = zeros(Float32, 6)
        Pioneer.get_PEP!(scores, is_target, pep1; fdr_scale_factor=1.0f0)
        Pioneer.get_PEP!(scores, is_target, pep2; fdr_scale_factor=2.0f0)

        # Higher scale factor → decoys contribute more → higher PEP for targets
        @test sum(pep2) >= sum(pep1)
    end

    @testset "unsorted input" begin
        scores    = Float32[0.2, 0.9, 0.1, 0.8]
        is_target = Bool[false, true, false, true]
        pep       = zeros(Float32, 4)
        Pioneer.get_PEP!(scores, is_target, pep)
        # Highest-scoring PSMs (indices 2, 4) should have lower PEP
        @test pep[2] <= pep[1]
        @test pep[4] <= pep[3]
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# _weighted_pava — isotonic regression helper
# ═══════════════════════════════════════════════════════════════════════════════

@testset "_weighted_pava" begin

    @testset "already monotonic" begin
        y = Float64[0.1, 0.3, 0.5, 0.7]
        w = Float64[1.0, 1.0, 1.0, 1.0]
        result = Pioneer._weighted_pava(y, w)
        @test result ≈ y
    end

    @testset "constant input" begin
        y = Float64[0.5, 0.5, 0.5]
        w = Float64[1.0, 1.0, 1.0]
        result = Pioneer._weighted_pava(y, w)
        @test result ≈ y
    end

    @testset "single violation" begin
        # [0.6, 0.4] violates monotonicity → merged to weighted average
        y = Float64[0.6, 0.4]
        w = Float64[1.0, 1.0]
        result = Pioneer._weighted_pava(y, w)
        @test result[1] ≈ 0.5
        @test result[2] ≈ 0.5
    end

    @testset "weighted merge" begin
        # [0.8, 0.2] with weights [3, 1] → weighted avg = (0.8*3 + 0.2*1)/4 = 0.65
        y = Float64[0.8, 0.2]
        w = Float64[3.0, 1.0]
        result = Pioneer._weighted_pava(y, w)
        @test result[1] ≈ 0.65
        @test result[2] ≈ 0.65
    end

    @testset "chain of violations" begin
        y = Float64[0.5, 0.3, 0.1, 0.9]
        w = Float64[1.0, 1.0, 1.0, 1.0]
        result = Pioneer._weighted_pava(y, w)
        # Must be non-decreasing
        for i in 2:length(result)
            @test result[i] >= result[i-1] - 1e-10
        end
        # Weighted average preserved
        @test sum(result .* w) ≈ sum(y .* w)
    end

    @testset "single element" begin
        y = Float64[0.42]
        w = Float64[1.0]
        result = Pioneer._weighted_pava(y, w)
        @test result ≈ [0.42]
    end
end

end # top-level testset
