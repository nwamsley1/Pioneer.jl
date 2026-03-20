# Unit tests for getDistanceMetrics (FirstPass variant)
#
# Target: Pioneer.getDistanceMetrics (spectralDistanceMetrics.jl:326-414)
# Single-pass spectral scoring. Computes 6 metrics per precursor:
# fitted_spectral_contrast, gof, max_matched_residual, max_unmatched_residual,
# fitted_manhattan_distance, percent_theoretical_ignored.
#
# Run: julia --project=. test/UnitTests/test_distance_metrics.jl

using Test
using Pioneer

# ---------------------------------------------------------------------------
# Aliases
# ---------------------------------------------------------------------------
const SA                     = Pioneer.SparseArray
const sortSparse!            = Pioneer.sortSparse!
const SpectralScoresFirstPass = Pioneer.SpectralScoresFirstPass
const getDistanceMetrics     = Pioneer.getDistanceMetrics

# ---------------------------------------------------------------------------
# Helper: build SparseArray from (col, row, nzval, x, matched, isotope) tuples
# ---------------------------------------------------------------------------
function build_test_sparse(entries)
    n = length(entries)
    sa = SA(UInt32(max(n + 10, 50)))
    sa.n_vals = n
    for (i, (col, row, nz, xv, m, iso)) in enumerate(entries)
        sa.colval[i] = UInt16(col)
        sa.rowval[i] = UInt32(row)
        sa.nzval[i]  = Float32(nz)
        sa.x[i]      = Float32(xv)
        sa.matched[i] = m
        sa.isotope[i] = UInt8(iso)
    end
    sortSparse!(sa)
    for i in 1:sa.n_vals
        sa.m = max(sa.m, Int64(sa.rowval[i]))
    end
    return sa
end

"""Run getDistanceMetrics and return the scores vector."""
function compute_metrics(H, w)
    r = zeros(Float32, H.m)
    scores = [SpectralScoresFirstPass(
        zero(Float32), zero(Float32), zero(Float32),
        zero(Float32), zero(Float32), zero(Float32),
        zero(Float32)
    ) for _ in 1:H.n]
    getDistanceMetrics(w, r, H, scores)
    return scores
end

# ═══════════════════════════════════════════════════════════════════════════
# Tests
# ═══════════════════════════════════════════════════════════════════════════

@testset "getDistanceMetrics (FirstPass)" begin

    # -------------------------------------------------------------------
    @testset "4.1 Single-precursor imperfect fit" begin
        # H: 1 col, 2 rows (all matched). nzval=[1, 2], x=[10, 22], w=[10]
        # Fitted: [10, 20]. Residuals: r=[0, -2].
        entries = [
            (1, 1, 1.0, 10.0, true, 0),
            (1, 2, 2.0, 22.0, true, 0),
        ]
        H = build_test_sparse(entries)
        w = Float32[10.0]
        scores = compute_metrics(H, w)
        s = scores[1]

        # Compute expected values in Float64 for precision
        nzval = [1.0, 2.0]; x_obs = [10.0, 22.0]; wv = 10.0
        fitted = wv .* nzval            # [10, 20]
        r_vals = wv .* nzval .- x_obs   # [0, -2]
        shadow = fitted .- r_vals       # [10, 22]

        x_sum = sum(x_obs)                          # 32
        manhattan = sum(abs.(fitted .- x_obs))       # 2
        sum_residuals = sum(abs.(r_vals))             # 2
        sum_fitted_matched = sum(fitted)              # 30
        fitted_dotp = sum(shadow .* fitted)           # 540
        fitted_dotp_norm1 = sum(shadow.^2)            # 584
        sum_fitted_sq = sum(fitted.^2)                # 500

        denom = sqrt(fitted_dotp_norm1) * sqrt(sum_fitted_sq)
        exp_sc = fitted_dotp / denom
        exp_gof = -log2(sum_residuals / sum_fitted_matched)
        exp_mmr = -log2(maximum(abs.(r_vals)) / sum_fitted_matched)
        exp_mur = -log2(0.0 / sum_fitted_matched + 1e-10)  # no unmatched
        exp_fmd = -log2(manhattan / x_sum + 1e-10)

        @test s.fitted_spectral_contrast ≈ Float32(exp_sc)  atol=0.01
        @test s.gof                      ≈ Float32(exp_gof) atol=0.1
        @test s.max_matched_residual     ≈ Float32(exp_mmr) atol=0.1
        @test s.max_unmatched_residual   ≈ Float32(exp_mur) atol=0.5
        @test s.fitted_manhattan_distance ≈ Float32(exp_fmd) atol=0.1
        @test s.percent_theoretical_ignored == 0.0f0
    end

    # -------------------------------------------------------------------
    @testset "4.2 Zero-weight precursor → all scores zero" begin
        entries = [
            (1, 1, 1.0, 10.0, true, 0),
            (1, 2, 2.0, 20.0, true, 0),
        ]
        H = build_test_sparse(entries)
        w = Float32[0.0]  # zero weight
        scores = compute_metrics(H, w)
        s = scores[1]

        @test s.fitted_spectral_contrast   == 0.0f0
        @test s.gof                        == 0.0f0
        @test s.max_matched_residual       == 0.0f0
        @test s.max_unmatched_residual     == 0.0f0
        @test s.fitted_manhattan_distance  == 0.0f0
        @test s.percent_theoretical_ignored == 0.0f0
    end

    # -------------------------------------------------------------------
    @testset "4.3 Two precursors with shared row" begin
        # H: 2 cols, 3 rows. col1: rows[1,2] nzval=[2,3], col2: rows[2,3] nzval=[1,4]
        # x = [5, 9, 6] (row 2 shared, same x=9 in both entries)
        # w = [3, 2]
        entries = [
            (1, 1, 2.0, 5.0, true, 0),   # col 1, row 1
            (1, 2, 3.0, 9.0, true, 0),   # col 1, row 2
            (2, 2, 1.0, 9.0, true, 0),   # col 2, row 2 (shared)
            (2, 3, 4.0, 6.0, true, 0),   # col 2, row 3
        ]
        H = build_test_sparse(entries)
        w = Float32[3.0, 2.0]
        scores = compute_metrics(H, w)

        # --- Expected residuals ---
        # r[1] = -5 + 3*2 = 1
        # r[2] = -9 + 3*3 + 2*1 = 2
        # r[3] = -6 + 2*4 = 2

        # --- Col 1 (w=3) ---
        # Entry (row 1): x=5, fitted=6, shadow=6-1=5
        #   manhattan=|6-5|=1, r_abs=1
        #   matched: sum_fitted_matched=6, fitted_dotp=30, fitted_dotp_norm1=25
        # Entry (row 2): x=9, fitted=9, shadow=9-2=7
        #   manhattan+=|9-9|=0→1, r_abs=2
        #   matched: sum_fitted_matched=15, fitted_dotp=93, fitted_dotp_norm1=74
        s1 = scores[1]
        x_sum1 = Float32(5 + 9)
        manh1  = Float32(1 + 0)
        sum_res1 = Float32(1 + 2)
        sfm1 = Float32(15); fdp1 = Float32(93); fdn1 = Float32(74); sfsq1 = Float32(36 + 81)
        denom1 = sqrt(fdn1) * sqrt(sfsq1)

        @test s1.fitted_spectral_contrast ≈ fdp1 / denom1            atol=0.01
        @test s1.gof                      ≈ -log2(sum_res1 / sfm1)   atol=0.1
        @test s1.max_matched_residual     ≈ -log2(Float32(2) / sfm1) atol=0.1
        @test s1.percent_theoretical_ignored == 0.0f0

        # --- Col 2 (w=2) ---
        # Entry (row 2): x=9, fitted=2, shadow=2-2=0
        #   manhattan=|2-9|=7, r_abs=2
        # Entry (row 3): x=6, fitted=8, shadow=8-2=6
        #   manhattan+=|8-6|=2→9, r_abs=2
        s2 = scores[2]
        sfm2 = Float32(2 + 8); sum_res2 = Float32(2 + 2)
        fdp2 = Float32(0*2 + 6*8); fdn2 = Float32(0 + 36); sfsq2 = Float32(4 + 64)
        denom2 = sqrt(fdn2) * sqrt(sfsq2)

        @test s2.fitted_spectral_contrast ≈ fdp2 / denom2            atol=0.01
        @test s2.gof                      ≈ -log2(sum_res2 / sfm2)   atol=0.1
        @test s2.percent_theoretical_ignored == 0.0f0
    end

    # -------------------------------------------------------------------
    @testset "4.4 Mixed matched/unmatched entries" begin
        # 1 precursor, 3 entries: 2 matched + 1 unmatched (miss, x=0)
        # w = [4]
        entries = [
            (1, 1, 2.0, 10.0, true,  0),  # matched
            (1, 2, 3.0, 15.0, true,  0),  # matched
            (1, 3, 1.0,  0.0, false, 0),  # unmatched (miss)
        ]
        H = build_test_sparse(entries)
        w = Float32[4.0]
        scores = compute_metrics(H, w)
        s = scores[1]

        # Residuals: r=[−10+8, −15+12, 0+4] = [−2, −3, 4]
        # Entry 1 (matched): fitted=8, shadow=8−(−2)=10, r_abs=2
        # Entry 2 (matched): fitted=12, shadow=12−(−3)=15, r_abs=3
        # Entry 3 (UNmatched): fitted=4, shadow=4−4=0, r_abs=4

        # Verify max_unmatched_residual captures the unmatched entry's residual
        sum_fitted = Float32(8 + 12 + 4)  # 24
        exp_mur = -log2(Float32(4) / sum_fitted + Float32(1e-10))
        @test s.max_unmatched_residual ≈ exp_mur atol=0.2

        # Verify max_matched_residual only considers matched entries (max=3)
        sum_fitted_matched = Float32(8 + 12)  # 20
        exp_mmr = -log2(Float32(3) / sum_fitted_matched)
        @test s.max_matched_residual ≈ exp_mmr atol=0.2

        # Verify unmatched does NOT contribute to fitted_dotp (cosine numerator)
        # fitted_dotp = 10*8 + 15*12 = 260 (only matched entries)
        # fitted_dotp_norm1 = 100 + 225 = 325 (only matched shadow^2)
        # sum_fitted_sq = 64 + 144 + 16 = 224 (ALL entries)
        exp_sc = Float32(260) / (sqrt(Float32(325)) * sqrt(Float32(224)))
        @test s.fitted_spectral_contrast ≈ exp_sc atol=0.01

        # Contrast with a hypothetical "all matched" version:
        # If the unmatched entry were matched, fitted_dotp would include 0*4=0,
        # so numerator stays 260 but denom changes. The key point is unmatched
        # entries push sum_fitted_sq up without increasing fitted_dotp.
        @test s.fitted_spectral_contrast < 1.0f0

        @test s.percent_theoretical_ignored == 0.0f0
    end

    # -------------------------------------------------------------------
    @testset "4.5 percent_theoretical_ignored always zero for FirstPass" begin
        # Verified implicitly in all tests above, but explicitly test
        # with a multi-precursor case
        entries = [
            (1, 1, 1.0, 5.0, true, 0),
            (1, 2, 2.0, 8.0, true, 0),
            (2, 3, 3.0, 9.0, true, 0),
            (2, 4, 1.0, 3.0, true, 0),
        ]
        H = build_test_sparse(entries)
        w = Float32[2.0, 3.0]
        scores = compute_metrics(H, w)

        for col in 1:H.n
            @test scores[col].percent_theoretical_ignored == 0.0f0
        end
    end

end  # getDistanceMetrics testset

println("\nAll getDistanceMetrics tests passed!")
