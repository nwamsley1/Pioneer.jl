# Unit tests for matchPeaks!, buildDesignMatrix!, and sortSparse!
#
# These test the core spectral matching pipeline with hand-verifiable inputs.
# Each matchPeaks! test has a diagram showing expected behavior.
#
# Run: julia --project=. test/UnitTests/SparseArray/test_matching_pipeline.jl

using Test
using Pioneer

# ---------------------------------------------------------------------------
# Aliases for internal types
# ---------------------------------------------------------------------------
const SA            = Pioneer.SparseArray
const sortSparse!   = Pioneer.sortSparse!
const matchPeaks!   = Pioneer.matchPeaks!
const buildDesignMatrix! = Pioneer.buildDesignMatrix!
const ArrayDict     = Pioneer.ArrayDict

const FragmentMatch = Pioneer.FragmentMatch
const DetailedFrag  = Pioneer.DetailedFrag
const UnmatchedIon  = Pioneer.UnmatchedIon
const MassErrorModel = Pioneer.MassErrorModel

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

"""
Construct a DetailedFrag{Float32} with only the fields that matter for matching.
All other fields get sensible defaults.
"""
function make_frag(; prec_id::UInt32, mz::Float32, intensity::Float16 = Float16(1.0),
                     is_y::Bool = true, is_isotope::Bool = false)
    DetailedFrag{Float32}(
        prec_id,
        mz,
        intensity,
        UInt16(0),      # ion_type
        is_y,           # is_y
        false,          # is_b
        false,          # is_p
        is_isotope,     # is_isotope
        UInt8(1),       # frag_charge
        UInt8(1),       # ion_position
        UInt8(2),       # prec_charge
        UInt8(1),       # rank
        UInt8(0),       # sulfur_count
    )
end

"""
Construct a FragmentMatch{Float32} with only the fields buildDesignMatrix! reads.
"""
function make_match(; predicted_intensity::Float32, intensity::Float32,
                      peak_ind::UInt32, prec_id::UInt32, is_isotope::Bool = false)
    FragmentMatch{Float32}(
        predicted_intensity,
        intensity,
        Float32(0),     # theoretical_mz (unused by buildDesignMatrix!)
        Float32(0),     # match_mz
        peak_ind,
        prec_id,
        UInt8(0),       # frag_index
        UInt8(1),       # frag_charge
        UInt8(0),       # ion_type
        is_isotope,
        UInt8(0),       # predicted_rank
    )
end

"""
Construct an UnmatchedIon.
"""
function make_miss(; prec_id::UInt32, predicted_intensity::Float32,
                     is_isotope::Bool = false)
    UnmatchedIon(prec_id, predicted_intensity, is_isotope)
end

"""Wide-tolerance MassErrorModel: ±100 ppm, no offset."""
const WIDE_TOL = MassErrorModel(0.0f0, (100.0f0, 100.0f0))

"""Convert SparseArray to dense Matrix{Float32} using nzval."""
function to_dense(H)
    M = zeros(Float32, H.m, H.n)
    for i in 1:H.n_vals
        r = H.rowval[i]
        c = H.colval[i]
        M[r, c] += H.nzval[i]
    end
    return M
end

"""Convert SparseArray x values to dense Matrix{Float32}."""
function to_dense_x(H)
    M = zeros(Float32, H.m, H.n)
    for i in 1:H.n_vals
        r = H.rowval[i]
        c = H.colval[i]
        M[r, c] = H.x[i]
    end
    return M
end

"""
Iterate all columns using colptr (the way initResiduals!/multiply! do)
and collect every (col, row, nzval) triple that is reachable.
"""
function collect_via_colptr(sa)
    entries = Tuple{UInt16, eltype(sa.rowval), eltype(sa.nzval)}[]
    for col in 1:sa.n
        start = sa.colptr[col]
        stop  = sa.colptr[col+1] - 1
        for idx in start:stop
            push!(entries, (sa.colval[idx], sa.rowval[idx], sa.nzval[idx]))
        end
    end
    return entries
end

# ═══════════════════════════════════════════════════════════════════════════
# matchPeaks! tests
# ═══════════════════════════════════════════════════════════════════════════

@testset "matchPeaks!" begin

    # -------------------------------------------------------------------
    @testset "Test 1: Perfect 1-to-1 matching" begin
        #
        #  ions (sorted by mz):  200 ──── 400 ──── 600
        #  peaks (sorted by mz): 200 ──── 400 ──── 600
        #                         ↕         ↕        ↕
        #                       match     match    match
        #
        ions = [
            make_frag(prec_id = UInt32(1), mz = 200.0f0),
            make_frag(prec_id = UInt32(2), mz = 400.0f0),
            make_frag(prec_id = UInt32(3), mz = 600.0f0),
        ]
        masses      = Vector{Union{Missing,Float32}}([200.0f0, 400.0f0, 600.0f0])
        intensities = Vector{Union{Missing,Float32}}([1f5,     2f5,     3f5])

        matches   = [FragmentMatch{Float32}() for _ in 1:100]
        unmatched = [UnmatchedIon() for _ in 1:100]

        nmatched, nunmatched = matchPeaks!(
            matches, unmatched, ions, length(ions),
            masses, intensities, WIDE_TOL, typemax(Float32))

        @test nmatched   == 3
        @test nunmatched == 0

        # Verify each match has the correct peak_ind and prec_id
        @test Pioneer.getPeakInd(matches[1]) == UInt32(1)
        @test Pioneer.getPrecID(matches[1])  == UInt32(1)
        @test Pioneer.getPeakInd(matches[2]) == UInt32(2)
        @test Pioneer.getPrecID(matches[2])  == UInt32(2)
        @test Pioneer.getPeakInd(matches[3]) == UInt32(3)
        @test Pioneer.getPrecID(matches[3])  == UInt32(3)

        # Verify observed intensities are carried through
        @test Pioneer.getIntensity(matches[1]) == 1f5
        @test Pioneer.getIntensity(matches[2]) == 2f5
        @test Pioneer.getIntensity(matches[3]) == 3f5
    end

    # -------------------------------------------------------------------
    @testset "Test 2: All ions miss (gap between ions and peaks)" begin
        #
        #  ions:  200 ── 300
        #  peaks:                500 ── 600
        #                  (no overlap at ±100 ppm)
        #
        ions = [
            make_frag(prec_id = UInt32(1), mz = 200.0f0),
            make_frag(prec_id = UInt32(2), mz = 300.0f0),
        ]
        masses      = Vector{Union{Missing,Float32}}([500.0f0, 600.0f0])
        intensities = Vector{Union{Missing,Float32}}([1f5,   1f5])

        matches   = [FragmentMatch{Float32}() for _ in 1:100]
        unmatched = [UnmatchedIon() for _ in 1:100]

        nmatched, nunmatched = matchPeaks!(
            matches, unmatched, ions, length(ions),
            masses, intensities, WIDE_TOL, typemax(Float32))

        @test nmatched   == 0
        @test nunmatched == 2

        # Verify misses have correct prec_ids
        @test Pioneer.getPrecID(unmatched[1]) == UInt32(1)
        @test Pioneer.getPrecID(unmatched[2]) == UInt32(2)
    end

    # -------------------------------------------------------------------
    @testset "Test 3: Multiple ions match same peak" begin
        #
        #  ions:  500.000 (prec 1) ── 500.001 (prec 2)
        #  peaks: 500.000
        #            ↕                    ↕
        #         match (peak 1)       match (peak 1)
        #
        # Both ions are within ±100 ppm of peak 500.0.
        # 500.001 - 500.0 = 0.001, which is 2 ppm → well within 100 ppm.
        #
        ions = [
            make_frag(prec_id = UInt32(1), mz = 500.0f0),
            make_frag(prec_id = UInt32(2), mz = 500.001f0),
        ]
        masses      = Vector{Union{Missing,Float32}}([500.0f0])
        intensities = Vector{Union{Missing,Float32}}([5f5])

        matches   = [FragmentMatch{Float32}() for _ in 1:100]
        unmatched = [UnmatchedIon() for _ in 1:100]

        nmatched, nunmatched = matchPeaks!(
            matches, unmatched, ions, length(ions),
            masses, intensities, WIDE_TOL, typemax(Float32))

        @test nmatched   == 2
        @test nunmatched == 0

        # Both match the same peak (index 1)
        @test Pioneer.getPeakInd(matches[1]) == UInt32(1)
        @test Pioneer.getPeakInd(matches[2]) == UInt32(1)
        @test Pioneer.getPrecID(matches[1])  == UInt32(1)
        @test Pioneer.getPrecID(matches[2])  == UInt32(2)
    end

    # -------------------------------------------------------------------
    @testset "Test 4: Nearest peak selection" begin
        #
        #  ion:   500.0
        #  peaks: 499.95 ───── 500.01
        #           ↑              ↑
        #         0.05 away     0.01 away  ← closer!
        #
        # Both peaks within ±100 ppm. Ion should match the CLOSER peak (500.01).
        #
        ions = [
            make_frag(prec_id = UInt32(1), mz = 500.0f0),
        ]
        masses      = Vector{Union{Missing,Float32}}([499.95f0, 500.01f0])
        intensities = Vector{Union{Missing,Float32}}([1f5,    2f5])

        matches   = [FragmentMatch{Float32}() for _ in 1:100]
        unmatched = [UnmatchedIon() for _ in 1:100]

        nmatched, nunmatched = matchPeaks!(
            matches, unmatched, ions, length(ions),
            masses, intensities, WIDE_TOL, typemax(Float32))

        @test nmatched   == 1
        @test nunmatched == 0

        # Should match peak 2 (mz=500.01, closer to 500.0 than 499.95)
        @test Pioneer.getPeakInd(matches[1]) == UInt32(2)
        @test Pioneer.getIntensity(matches[1]) == 2f5
    end

    # -------------------------------------------------------------------
    @testset "Test 5a: Empty ions" begin
        masses      = Vector{Union{Missing,Float32}}([100.0f0, 200.0f0])
        intensities = Vector{Union{Missing,Float32}}([1f5,   1f5])
        ions = DetailedFrag{Float32}[]

        matches   = [FragmentMatch{Float32}() for _ in 1:100]
        unmatched = [UnmatchedIon() for _ in 1:100]

        nmatched, nunmatched = matchPeaks!(
            matches, unmatched, ions, 0,
            masses, intensities, WIDE_TOL, typemax(Float32))

        @test nmatched   == 0
        @test nunmatched == 0
    end

    @testset "Test 5b: Empty peaks" begin
        ions = [make_frag(prec_id = UInt32(1), mz = 500.0f0)]
        masses      = Vector{Union{Missing,Float32}}()
        intensities = Vector{Union{Missing,Float32}}()

        matches   = [FragmentMatch{Float32}() for _ in 1:100]
        unmatched = [UnmatchedIon() for _ in 1:100]

        nmatched, nunmatched = matchPeaks!(
            matches, unmatched, ions, length(ions),
            masses, intensities, WIDE_TOL, typemax(Float32))

        @test nmatched   == 0
        @test nunmatched == 0
    end

    # -------------------------------------------------------------------
    @testset "Test 6: Ions beyond high_mass are skipped" begin
        #
        #  ion:   2000.0
        #  peaks: 100.0
        #  high_mass = 1500.0
        #
        # Ion is above high_mass → skipped entirely (not a miss, not a match).
        #
        ions = [
            make_frag(prec_id = UInt32(1), mz = 2000.0f0),
        ]
        masses      = Vector{Union{Missing,Float32}}([100.0f0])
        intensities = Vector{Union{Missing,Float32}}([1f5])

        matches   = [FragmentMatch{Float32}() for _ in 1:100]
        unmatched = [UnmatchedIon() for _ in 1:100]

        nmatched, nunmatched = matchPeaks!(
            matches, unmatched, ions, length(ions),
            masses, intensities, WIDE_TOL, 1500.0f0)

        @test nmatched   == 0
        @test nunmatched == 0   # NOT written to misses — skipped
    end

    # -------------------------------------------------------------------
    @testset "Test 7: Interleaved matches and misses" begin
        #
        #  ions (prec_id):  200(1) ── 350(2) ── 500(3) ── 700(4)
        #  peaks:           200    ──────────── 500
        #                    ↕         miss       ↕         miss
        #                  match                match
        #
        ions = [
            make_frag(prec_id = UInt32(1), mz = 200.0f0),
            make_frag(prec_id = UInt32(2), mz = 350.0f0),
            make_frag(prec_id = UInt32(3), mz = 500.0f0),
            make_frag(prec_id = UInt32(4), mz = 700.0f0),
        ]
        masses      = Vector{Union{Missing,Float32}}([200.0f0, 500.0f0])
        intensities = Vector{Union{Missing,Float32}}([1f5,   2f5])

        matches   = [FragmentMatch{Float32}() for _ in 1:100]
        unmatched = [UnmatchedIon() for _ in 1:100]

        nmatched, nunmatched = matchPeaks!(
            matches, unmatched, ions, length(ions),
            masses, intensities, WIDE_TOL, typemax(Float32))

        @test nmatched   == 2
        @test nunmatched == 2

        # Matched ions
        @test Pioneer.getPrecID(matches[1]) == UInt32(1)
        @test Pioneer.getPrecID(matches[2]) == UInt32(3)
        @test Pioneer.getPeakInd(matches[1]) == UInt32(1)   # peak at 200
        @test Pioneer.getPeakInd(matches[2]) == UInt32(2)   # peak at 500

        # Missed ions
        @test Pioneer.getPrecID(unmatched[1]) == UInt32(2)   # ion at 350
        @test Pioneer.getPrecID(unmatched[2]) == UInt32(4)   # ion at 700
    end

    # -------------------------------------------------------------------
    @testset "Test 8: Mass error model offset enables match" begin
        #
        # Ion at mz 1000.0, peak at mz 1000.01 (10 ppm away).
        #
        # Without offset (±5 ppm tolerance):
        #   corrected peak = 1000.01, bounds = [999.995, 1000.005]
        #   1000.01 > 1000.005 → MISS
        #
        # With +10 ppm offset (±5 ppm tolerance):
        #   corrected peak = 1000.01 - 10*(1000.01/1e6) ≈ 999.99999
        #   bounds = [999.995, 1000.005]
        #   999.99999 ∈ [999.995, 1000.005] → MATCH
        #

        ions = [
            make_frag(prec_id = UInt32(1), mz = 1000.0f0),
        ]
        masses      = Vector{Union{Missing,Float32}}([1000.01f0])
        intensities = Vector{Union{Missing,Float32}}([1f5])

        matches   = [FragmentMatch{Float32}() for _ in 1:100]
        unmatched = [UnmatchedIon() for _ in 1:100]

        # Without offset: should NOT match (10 ppm gap > 5 ppm tolerance)
        tight_no_offset = MassErrorModel(0.0f0, (5.0f0, 5.0f0))
        nmatched, nunmatched = matchPeaks!(
            matches, unmatched, ions, length(ions),
            masses, intensities, tight_no_offset, typemax(Float32))

        @test nmatched   == 0
        @test nunmatched == 1

        # With +10 ppm offset: corrected peak ≈ 1000.0, now within ±5 ppm
        tight_with_offset = MassErrorModel(10.0f0, (5.0f0, 5.0f0))
        matches   = [FragmentMatch{Float32}() for _ in 1:100]
        unmatched = [UnmatchedIon() for _ in 1:100]

        nmatched, nunmatched = matchPeaks!(
            matches, unmatched, ions, length(ions),
            masses, intensities, tight_with_offset, typemax(Float32))

        @test nmatched   == 1
        @test nunmatched == 0
        @test Pioneer.getPeakInd(matches[1]) == UInt32(1)
    end

end  # matchPeaks! testset

# ═══════════════════════════════════════════════════════════════════════════
# buildDesignMatrix! tests
# ═══════════════════════════════════════════════════════════════════════════

@testset "buildDesignMatrix!" begin

    # -------------------------------------------------------------------
    @testset "Test 1: Basic 2-precursor matrix" begin
        #
        # Matches (sorted by peak_ind, prec_id):
        #   peak 1 → prec 1 (pred=2.0, obs=100), prec 2 (pred=3.0, obs=100)
        #   peak 2 → prec 1 (pred=5.0, obs=200)
        #   peak 3 → prec 2 (pred=7.0, obs=300)
        #
        # Misses:
        #   prec 1 (pred=1.0)
        #
        # Expected dense nzval matrix (rows=peaks+misses, cols=precursors):
        #         col1  col2
        # row 1 [ 2.0   3.0 ]   ← peak 1
        # row 2 [ 5.0   0   ]   ← peak 2
        # row 3 [ 0     7.0 ]   ← peak 3
        # row 4 [ 1.0   0   ]   ← miss for prec 1
        #
        test_matches = sort([
            make_match(predicted_intensity=2.0f0, intensity=100.0f0, peak_ind=UInt32(1), prec_id=UInt32(1)),
            make_match(predicted_intensity=3.0f0, intensity=100.0f0, peak_ind=UInt32(1), prec_id=UInt32(2)),
            make_match(predicted_intensity=5.0f0, intensity=200.0f0, peak_ind=UInt32(2), prec_id=UInt32(1)),
            make_match(predicted_intensity=7.0f0, intensity=300.0f0, peak_ind=UInt32(3), prec_id=UInt32(2)),
        ], lt = Pioneer.ion_match_lt)

        test_misses = [
            make_miss(prec_id=UInt32(1), predicted_intensity=1.0f0),
        ]

        H = SA(UInt32(50))
        id_to_col = ArrayDict(UInt32, UInt16, 10)

        buildDesignMatrix!(H, test_matches, test_misses,
                           length(test_matches), length(test_misses), id_to_col)

        @test H.n == 2   # 2 precursors
        @test H.m == 4   # 3 peaks + 1 miss = 4 rows

        expected_nzval = Float32[
            2.0  3.0;
            5.0  0.0;
            0.0  7.0;
            1.0  0.0;
        ]
        @test to_dense(H) ≈ expected_nzval

        # Verify x values: observed intensities for matches, 0 for misses
        expected_x = Float32[
            100.0  100.0;
            200.0    0.0;
              0.0  300.0;
              0.0    0.0;
        ]
        @test to_dense_x(H) ≈ expected_x
    end

    # -------------------------------------------------------------------
    @testset "Test 2: Same precursor, same peak — nzval accumulation" begin
        #
        # Two fragments from prec 1 both match peak 1:
        #   frag A: predicted=3.0, frag B: predicted=7.0
        #
        # Expected: single entry at (row=1, col=1) with nzval = 3.0 + 7.0 = 10.0
        #
        test_matches = [
            make_match(predicted_intensity=3.0f0, intensity=500.0f0, peak_ind=UInt32(1), prec_id=UInt32(1)),
            make_match(predicted_intensity=7.0f0, intensity=500.0f0, peak_ind=UInt32(1), prec_id=UInt32(1)),
        ]
        # Already sorted since same peak_ind and prec_id

        H = SA(UInt32(50))
        id_to_col = ArrayDict(UInt32, UInt16, 10)

        buildDesignMatrix!(H, test_matches, UnmatchedIon[],
                           length(test_matches), 0, id_to_col)

        @test H.n_vals == 1  # single entry (accumulated)
        @test H.n == 1       # 1 precursor
        @test H.m == 1       # 1 peak

        @test H.nzval[1] ≈ 10.0f0   # 3.0 + 7.0
        @test H.x[1] ≈ 500.0f0      # observed intensity from the peak
    end

    # -------------------------------------------------------------------
    @testset "Test 3: Misses only for precursors with matches" begin
        #
        # Prec 1: has a match at peak 1.
        # Prec 2: only has a miss (no matches).
        #
        # Expected: prec 2's miss is NOT in the matrix.
        #
        test_matches = [
            make_match(predicted_intensity=5.0f0, intensity=100.0f0, peak_ind=UInt32(1), prec_id=UInt32(1)),
        ]

        test_misses = [
            make_miss(prec_id=UInt32(1), predicted_intensity=2.0f0),    # prec 1 has matches → included
            make_miss(prec_id=UInt32(2), predicted_intensity=3.0f0),    # prec 2 has NO matches → excluded
        ]

        H = SA(UInt32(50))
        id_to_col = ArrayDict(UInt32, UInt16, 10)

        buildDesignMatrix!(H, test_matches, test_misses,
                           length(test_matches), length(test_misses), id_to_col)

        @test H.n == 1       # only prec 1 (prec 2 excluded)
        @test H.n_vals == 2  # 1 match + 1 miss (for prec 1 only)

        expected = Float32[5.0; 2.0;;]   # 2×1 matrix
        @test to_dense(H) ≈ expected
    end

    # -------------------------------------------------------------------
    @testset "Test 4: Observed intensities (x values)" begin
        #
        # 2 matches with different observed intensities, 1 miss (x=0).
        #
        test_matches = sort([
            make_match(predicted_intensity=1.0f0, intensity=999.0f0, peak_ind=UInt32(1), prec_id=UInt32(1)),
            make_match(predicted_intensity=1.0f0, intensity=777.0f0, peak_ind=UInt32(2), prec_id=UInt32(1)),
        ], lt = Pioneer.ion_match_lt)

        test_misses = [
            make_miss(prec_id=UInt32(1), predicted_intensity=1.0f0),
        ]

        H = SA(UInt32(50))
        id_to_col = ArrayDict(UInt32, UInt16, 10)

        buildDesignMatrix!(H, test_matches, test_misses,
                           length(test_matches), length(test_misses), id_to_col)

        @test H.n_vals == 3

        # x values: matches carry observed intensity, miss has 0
        x_vals = H.x[1:H.n_vals]
        matched_flags = H.matched[1:H.n_vals]

        for i in 1:H.n_vals
            if matched_flags[i]
                @test x_vals[i] > 0.0f0   # matches have observed intensity
            else
                @test x_vals[i] == 0.0f0   # misses have x=0
            end
        end
    end

    # -------------------------------------------------------------------
    @testset "Test 5: colptr correctness after sortSparse! (off-by-one)" begin
        #
        # Build a 3-precursor matrix and verify that every entry is reachable
        # through the colptr ranges after sortSparse! (called inside buildDesignMatrix!).
        #
        test_matches = sort([
            make_match(predicted_intensity=1.0f0, intensity=100.0f0, peak_ind=UInt32(1), prec_id=UInt32(1)),
            make_match(predicted_intensity=2.0f0, intensity=100.0f0, peak_ind=UInt32(1), prec_id=UInt32(2)),
            make_match(predicted_intensity=3.0f0, intensity=200.0f0, peak_ind=UInt32(2), prec_id=UInt32(3)),
            make_match(predicted_intensity=4.0f0, intensity=300.0f0, peak_ind=UInt32(3), prec_id=UInt32(1)),
            make_match(predicted_intensity=5.0f0, intensity=300.0f0, peak_ind=UInt32(3), prec_id=UInt32(3)),
        ], lt = Pioneer.ion_match_lt)

        test_misses = [
            make_miss(prec_id=UInt32(2), predicted_intensity=0.5f0),
        ]

        H = SA(UInt32(50))
        id_to_col = ArrayDict(UInt32, UInt16, 10)

        buildDesignMatrix!(H, test_matches, test_misses,
                           length(test_matches), length(test_misses), id_to_col)

        # colptr sentinel value must be n_vals + 1
        @test H.colptr[H.n + 1] == H.n_vals + 1

        # Every entry must be reachable through colptr iteration
        entries = collect_via_colptr(H)
        @test length(entries) == H.n_vals

        # Verify all nzval values are present (order doesn't matter)
        found_nzvals = sort([e[3] for e in entries])
        direct_nzvals = sort(H.nzval[1:H.n_vals])
        @test found_nzvals ≈ direct_nzvals

        # Verify initResiduals! works correctly with this matrix
        w = ones(Float32, H.n)
        r = zeros(Float32, H.m)
        Pioneer.initResiduals!(r, H, w)

        # r[row] = sum(nzval for entries in that row) - x[row]
        # Every row should have a finite residual
        for i in 1:H.m
            @test isfinite(r[i])
        end
    end

end  # buildDesignMatrix! testset

# ═══════════════════════════════════════════════════════════════════════════
# sortSparse! tests (complements sortSparse_colptr_bug.jl)
# ═══════════════════════════════════════════════════════════════════════════

@testset "sortSparse! (additional)" begin

    @testset "Entries sorted by column after sortSparse!" begin
        sa = SA(UInt32(20))
        # Entries in reverse column order
        sa.n_vals = 4
        sa.colval[1] = UInt16(3); sa.rowval[1] = UInt32(1); sa.nzval[1] = 30.0f0
        sa.colval[2] = UInt16(1); sa.rowval[2] = UInt32(2); sa.nzval[2] = 10.0f0
        sa.colval[3] = UInt16(2); sa.rowval[3] = UInt32(1); sa.nzval[3] = 20.0f0
        sa.colval[4] = UInt16(1); sa.rowval[4] = UInt32(3); sa.nzval[4] = 11.0f0
        for i in 1:4; sa.x[i] = Float32(i); sa.matched[i] = true; sa.isotope[i] = UInt8(0); end

        sortSparse!(sa)

        # After sort, colval must be non-decreasing
        for i in 1:sa.n_vals-1
            @test sa.colval[i] <= sa.colval[i+1]
        end

        # x values must stay with their original entries (sidecar sort)
        # col 1 had nzvals [10, 11], col 2 had [20], col 3 had [30]
        entries = collect_via_colptr(sa)
        col1 = sort([e[3] for e in entries if e[1] == UInt16(1)])
        col2 = sort([e[3] for e in entries if e[1] == UInt16(2)])
        col3 = sort([e[3] for e in entries if e[1] == UInt16(3)])
        @test col1 == [10.0f0, 11.0f0]
        @test col2 == [20.0f0]
        @test col3 == [30.0f0]
    end

end  # sortSparse! testset

println("\nAll tests passed!")
