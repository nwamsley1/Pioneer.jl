# Unit tests for ScoreFragmentMatches! + ModifyFeatures!
#
# Target: Pioneer.ScoreFragmentMatches! (UnscoredPSMs.jl:75-87)
#         Pioneer.ModifyFeatures!       (UnscoredPSMs.jl:141-217)
# Accumulates per-match features into ComplexUnscoredPSM objects.
#
# Run: julia --project=. test/UnitTests/test_score_fragment_matches.jl

using Test
using Pioneer

# ---------------------------------------------------------------------------
# Aliases
# ---------------------------------------------------------------------------
const FragmentMatch          = Pioneer.FragmentMatch
const ComplexUnscoredPSM     = Pioneer.ComplexUnscoredPSM
const ArrayDict              = Pioneer.ArrayDict
const MassErrorModel         = Pioneer.MassErrorModel
const ScoreFragmentMatches!  = Pioneer.ScoreFragmentMatches!

# ---------------------------------------------------------------------------
# Helper: build a FragmentMatch with specified scoring-relevant fields
# ---------------------------------------------------------------------------
function make_scored_match(;
    prec_id::UInt32,
    intensity::Float32 = 100.0f0,
    theoretical_mz::Float32 = 500.0f0,
    match_mz::Float32 = 500.0f0,
    ion_type::UInt8 = UInt8(2),      # default: y-ion
    is_isotope::Bool = false,
    predicted_rank::UInt8 = UInt8(1),
    frag_index::UInt8 = UInt8(1),
)
    FragmentMatch{Float32}(
        1.0f0,            # predicted_intensity (not used by ModifyFeatures!)
        intensity,
        theoretical_mz,
        match_mz,
        UInt32(1),        # peak_ind (not used by ModifyFeatures!)
        prec_id,
        frag_index,
        UInt8(1),         # frag_charge
        ion_type,
        is_isotope,
        predicted_rank,
    )
end

"""Dummy MassErrorModel — not actually used by ModifyFeatures!."""
const DUMMY_ERR = MassErrorModel(0.0f0, (10.0f0, 10.0f0))

"""Set up ArrayDict mapping prec_ids to 1-based column indices."""
function setup_id_to_col(prec_ids::Vector{UInt32})
    max_id = maximum(prec_ids)
    id_to_col = ArrayDict(UInt32, UInt16, Int(max_id + 1))
    for (col, pid) in enumerate(prec_ids)
        Pioneer.update!(id_to_col, pid, UInt16(col))
    end
    return id_to_col
end

# ═══════════════════════════════════════════════════════════════════════════
# Tests
# ═══════════════════════════════════════════════════════════════════════════

@testset "ScoreFragmentMatches!" begin

    # -------------------------------------------------------------------
    @testset "3.1 Single y-ion" begin
        matches = [
            make_scored_match(
                prec_id = UInt32(1),
                ion_type = UInt8(2),          # y-ion
                predicted_rank = UInt8(3),
                frag_index = UInt8(4),
                theoretical_mz = 500.0f0,
                match_mz = 500.005f0,         # 10 ppm error
                intensity = 1000.0f0,
            ),
        ]
        id_to_col = setup_id_to_col([UInt32(1)])
        results = [ComplexUnscoredPSM{Float32}()]
        m_rank = 5

        ScoreFragmentMatches!(results, id_to_col, matches, 1, DUMMY_ERR, m_rank)

        psm = results[1]
        @test psm.y_count == 1
        @test psm.b_count == 0
        @test psm.isotope_count == 0
        @test psm.best_rank == UInt8(3)
        @test psm.topn == 1             # rank 3 ≤ m_rank 5
        @test psm.longest_y == UInt8(4) # frag_index
        @test psm.precursor_idx == UInt32(1)

        # ppm_err = (500.0 - 500.005) / (500.0 / 1e6) = -10.0
        @test psm.error ≈ 10.0f0 atol=0.5
        @test psm.y_int ≈ 1000.0f0
    end

    # -------------------------------------------------------------------
    @testset "3.2 Three match types (y, b, isotope)" begin
        matches = [
            # y-ion: rank 2, frag_index 5, intensity 100
            make_scored_match(
                prec_id = UInt32(1), ion_type = UInt8(2),
                predicted_rank = UInt8(2), frag_index = UInt8(5),
                intensity = 100.0f0,
            ),
            # b-ion: rank 4, frag_index 3, intensity 200
            make_scored_match(
                prec_id = UInt32(1), ion_type = UInt8(1),
                predicted_rank = UInt8(4), frag_index = UInt8(3),
                intensity = 200.0f0,
            ),
            # isotope: rank 1
            make_scored_match(
                prec_id = UInt32(1), is_isotope = true,
                predicted_rank = UInt8(1), intensity = 50.0f0,
            ),
        ]
        id_to_col = setup_id_to_col([UInt32(1)])
        results = [ComplexUnscoredPSM{Float32}()]
        m_rank = 3

        ScoreFragmentMatches!(results, id_to_col, matches, 3, DUMMY_ERR, m_rank)

        psm = results[1]
        @test psm.y_count == 1
        @test psm.b_count == 1
        @test psm.isotope_count == 1
        @test psm.best_rank == UInt8(2)      # best non-isotope: min(2, 4)
        @test psm.best_rank_iso == UInt8(1)  # isotope rank
        @test psm.topn == 1                  # only y (rank 2 ≤ 3); b rank 4 > 3
        @test psm.topn_iso == 1              # isotope rank 1 ≤ 3
        @test psm.longest_y == UInt8(5)
        @test psm.longest_b == UInt8(3)
        @test psm.y_int ≈ 100.0f0
        @test psm.b_int ≈ 200.0f0
    end

    # -------------------------------------------------------------------
    @testset "3.3 Error accumulation" begin
        # 3 y-ions with known ppm errors
        matches = [
            # ppm = (500 - 500.005) / (500/1e6) = -10.0 → |err| = 10.0
            make_scored_match(
                prec_id = UInt32(1), theoretical_mz = 500.0f0,
                match_mz = 500.005f0,
            ),
            # ppm = (1000 - 999.99) / (1000/1e6) = 10.0 → |err| = 10.0
            make_scored_match(
                prec_id = UInt32(1), theoretical_mz = 1000.0f0,
                match_mz = 999.99f0,
            ),
            # ppm = (800 - 800.004) / (800/1e6) = -5.0 → |err| = 5.0
            make_scored_match(
                prec_id = UInt32(1), theoretical_mz = 800.0f0,
                match_mz = 800.004f0,
            ),
        ]
        id_to_col = setup_id_to_col([UInt32(1)])
        results = [ComplexUnscoredPSM{Float32}()]

        ScoreFragmentMatches!(results, id_to_col, matches, 3, DUMMY_ERR, 10)

        # error = |−10| + |10| + |−5| = 25
        @test results[1].error ≈ 25.0f0 atol=1.0
    end

    # -------------------------------------------------------------------
    @testset "3.4 Intensity accumulation" begin
        matches = [
            make_scored_match(
                prec_id = UInt32(1), ion_type = UInt8(2),
                intensity = 100.0f0,
            ),
            make_scored_match(
                prec_id = UInt32(1), ion_type = UInt8(2),
                intensity = 250.0f0,
            ),
        ]
        id_to_col = setup_id_to_col([UInt32(1)])
        results = [ComplexUnscoredPSM{Float32}()]

        ScoreFragmentMatches!(results, id_to_col, matches, 2, DUMMY_ERR, 10)

        @test results[1].y_int ≈ 350.0f0
        @test results[1].y_count == 2
    end

    # -------------------------------------------------------------------
    @testset "3.5 Multiple precursors — independent accumulation" begin
        matches = [
            make_scored_match(
                prec_id = UInt32(1), ion_type = UInt8(2),
                intensity = 100.0f0, predicted_rank = UInt8(2),
            ),
            make_scored_match(
                prec_id = UInt32(3), ion_type = UInt8(1),
                intensity = 500.0f0, predicted_rank = UInt8(5),
            ),
        ]
        id_to_col = setup_id_to_col([UInt32(1), UInt32(3)])
        results = [ComplexUnscoredPSM{Float32}(), ComplexUnscoredPSM{Float32}()]

        ScoreFragmentMatches!(results, id_to_col, matches, 2, DUMMY_ERR, 10)

        # Prec 1 → col 1
        @test results[1].y_count == 1
        @test results[1].b_count == 0
        @test results[1].y_int ≈ 100.0f0
        @test results[1].best_rank == UInt8(2)

        # Prec 3 → col 2
        @test results[2].b_count == 1
        @test results[2].y_count == 0
        @test results[2].b_int ≈ 500.0f0
        @test results[2].best_rank == UInt8(5)
    end

    # -------------------------------------------------------------------
    @testset "3.6 Default initial state" begin
        psm = ComplexUnscoredPSM{Float32}()

        @test psm.best_rank == UInt8(255)
        @test psm.best_rank_iso == UInt8(255)
        @test psm.topn == 0
        @test psm.topn_iso == 0
        @test psm.y_count == 0
        @test psm.b_count == 0
        @test psm.isotope_count == 0
        @test psm.error == 0.0f0
        @test psm.y_int == 0.0f0
        @test psm.b_int == 0.0f0
        @test psm.precursor_idx == UInt32(0)
    end

    # -------------------------------------------------------------------
    @testset "3.7 best_rank monotonicity" begin
        # Feed rank 5 then rank 2 → best_rank should be 2
        matches_a = [
            make_scored_match(
                prec_id = UInt32(1), predicted_rank = UInt8(5),
            ),
            make_scored_match(
                prec_id = UInt32(1), predicted_rank = UInt8(2),
            ),
        ]
        id_to_col = setup_id_to_col([UInt32(1)])
        results_a = [ComplexUnscoredPSM{Float32}()]
        ScoreFragmentMatches!(results_a, id_to_col, matches_a, 2, DUMMY_ERR, 10)
        @test results_a[1].best_rank == UInt8(2)

        # Feed rank 2 then rank 5 → best_rank should still be 2
        matches_b = [
            make_scored_match(
                prec_id = UInt32(1), predicted_rank = UInt8(2),
            ),
            make_scored_match(
                prec_id = UInt32(1), predicted_rank = UInt8(5),
            ),
        ]
        results_b = [ComplexUnscoredPSM{Float32}()]
        ScoreFragmentMatches!(results_b, id_to_col, matches_b, 2, DUMMY_ERR, 10)
        @test results_b[1].best_rank == UInt8(2)
    end

end  # ScoreFragmentMatches! testset

println("\nAll ScoreFragmentMatches! tests passed!")
