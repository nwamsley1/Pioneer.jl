# Stress-tests for the ML feature pipeline guarding against Inf/NaN.
#
# Background: a 1-fragment match where the matched peak's centroid landed
# on the predicted m/z to Float32 precision produced `error == 0` →
# `log2(0) = -Inf` in the `error` feature, which propagated into the
# probit IRLS solve and tripped LAPACK's chkfinite. The fix added explicit
# floors. These tests stress the floors directly so future edits don't
# silently regress them.

using Test
using Pioneer
using LinearAlgebra

const _SparseArrayFused        = Pioneer.SparseArrayFused
const _SpectralScoresMainSearch = Pioneer.SpectralScoresMainSearch
const _getDistanceMetrics      = Pioneer.getDistanceMetrics
const _psm_getPoisson          = Pioneer.psm_getPoisson
const _psm_logfac              = Pioneer.psm_logfac
const _pack_meta               = Pioneer.pack_meta

# ---------------------------------------------------------------------------
# Build a 1-column SparseArrayFused with `peaks` rows. `matched` is a
# Bool vector; `nzval[i]` is the predicted (theoretical) intensity for
# row i; `x[i]` is the observed intensity. Used by getDistanceMetrics
# stress tests below.
# ---------------------------------------------------------------------------
function _make_single_col_H(nzval::Vector{Float32}, x::Vector{Float32},
                            matched::Vector{Bool};
                            ranks::Vector{UInt8}=zeros(UInt8, length(nzval)))
    n = length(nzval)
    @assert length(x) == n && length(matched) == n && length(ranks) == n
    H = _SparseArrayFused(UInt32(max(n, 1)))
    H.n_vals = n
    H.m = n
    H.n = 1
    @inbounds for i in 1:n
        H.rowval[i] = UInt32(i)
        H.nzval[i]  = nzval[i]
        H.x[i]      = x[i]
        H.meta[i]   = _pack_meta(matched[i], UInt8(0))
        H.rank[i]   = ranks[i]
    end
    # colptr layout: [start_of_col_1, start_of_col_2, ..., n_vals+1]
    resize!(H.colptr, 2)
    H.colptr[1] = 1
    H.colptr[2] = n + 1
    H
end

function _score_single(w_val::Float32, nzval, x, matched;
                       ranks=zeros(UInt8, length(nzval)))
    H = _make_single_col_H(nzval, x, matched; ranks=ranks)
    w = Float32[w_val]
    r = zeros(Float32, H.m)
    spectral_scores = Vector{_SpectralScoresMainSearch{Float16, Float32}}(undef, 1)
    _getDistanceMetrics(w, r, H, spectral_scores)
    spectral_scores[1]
end

@testset "Feature finiteness — psm_getPoisson" begin
    # Hypothesis from production: zero-error or zero-fragment edge cases
    # should never produce NaN or Inf in any feature.

    # Sanity: typical case
    p = _psm_getPoisson(2.5, 3)
    @test isfinite(p)

    # observed == 0 (logfac(0) used to be 0*log(0) = NaN)
    @test isfinite(_psm_getPoisson(2.5, 0))
    @test isfinite(_psm_getPoisson(0.0, 0))

    # lam == 0, observed > 0 (used to be log(0^k) = -Inf)
    @test isfinite(_psm_getPoisson(0.0, 3))
    @test isfinite(_psm_getPoisson(0.0f0, 1))

    # Very small λ
    @test isfinite(_psm_getPoisson(1e-30, 5))

    # Large k — Stirling form should still be finite
    @test isfinite(_psm_getPoisson(50.0, 100))

    # Float32 input still produces a finite scalar (call site converts to Float16)
    @test isfinite(_psm_getPoisson(2.5f0, 3))

    # logfac(0) = 0 by convention
    @test _psm_logfac(0) == 0.0
    @test _psm_logfac(-5) == 0.0
    @test isfinite(_psm_logfac(1))
    @test isfinite(_psm_logfac(255))
end

@testset "Feature finiteness — getDistanceMetrics edge cases" begin
    # Case A — perfect single-fragment match: residual = 0 across the board.
    # Used to give gof = -log2(0) = +Inf and max_matched_residual = +Inf.
    s = _score_single(1.0f0,
                      Float32[1.0],   # nzval (theoretical)
                      Float32[1.0],   # x (observed)
                      Bool[true])
    @test isfinite(s.gof)
    @test isfinite(s.max_matched_residual)
    @test isfinite(s.max_unmatched_residual)
    @test isfinite(s.fitted_manhattan_distance)
    @test isfinite(s.fitted_hellinger)

    # Case B — perfect multi-fragment match
    s = _score_single(1.0f0,
                      Float32[1.0, 2.0, 3.0],
                      Float32[1.0, 2.0, 3.0],
                      Bool[true, true, true])
    for v in (s.gof, s.max_matched_residual, s.max_unmatched_residual,
              s.fitted_manhattan_distance, s.fitted_hellinger)
        @test isfinite(v)
    end

    # Case C — zero-weight column (early exit path, all zeros)
    s = _score_single(0.0f0,
                      Float32[1.0, 2.0],
                      Float32[1.0, 2.0],
                      Bool[true, true])
    for v in (s.gof, s.max_matched_residual, s.max_unmatched_residual,
              s.fitted_manhattan_distance, s.fitted_hellinger)
        @test isfinite(v)
        @test v == zero(Float16)
    end

    # Case D — single matched peak with residual exactly 0; the unmatched
    # max-residual branch is also degenerate (no unmatched fragments).
    s = _score_single(1.0f0,
                      Float32[5.0],
                      Float32[5.0],
                      Bool[true])
    @test isfinite(s.gof)
    @test isfinite(s.max_matched_residual)
    @test isfinite(s.max_unmatched_residual)

    # Case E — fitted peaks but observed all zero (extreme noise scenario)
    s = _score_single(1.0f0,
                      Float32[1.0, 1.0],
                      Float32[0.0, 0.0],
                      Bool[true, true])
    for v in (s.gof, s.max_matched_residual, s.max_unmatched_residual,
              s.fitted_manhattan_distance, s.fitted_hellinger)
        @test isfinite(v)
    end

    # Case F — purely unmatched fragments (sum_of_fitted_peaks_matched = 0
    # → branch returns zero, must not produce NaN)
    s = _score_single(1.0f0,
                      Float32[1.0, 1.0],
                      Float32[1.0, 1.0],
                      Bool[false, false])
    for v in (s.gof, s.max_matched_residual, s.max_unmatched_residual,
              s.fitted_manhattan_distance, s.fitted_hellinger)
        @test isfinite(v)
    end

    # Case G — random stress, 200 trials with arbitrary intensities and
    # match flags. Floors must hold for every shape of input.
    rng_state = 1234
    for trial in 1:200
        n = mod(trial, 5) + 1
        # deterministic pseudo-random without rand()
        rng_state = (rng_state * 1103515245 + 12345) & 0x7fffffff
        nzval = Float32[abs(sin(Float32(trial * i))) for i in 1:n]
        x     = Float32[abs(cos(Float32(trial * i + 7))) for i in 1:n]
        matched = Bool[(i + trial) % 2 == 0 for i in 1:n]
        # Sprinkle in zeros to hit floor branches
        if trial % 7 == 0
            x[1] = 0.0f0
        end
        if trial % 11 == 0
            nzval[end] = 0.0f0
        end
        s = _score_single(Float32(0.5 + (trial % 3) * 0.25), nzval, x, matched)
        @test isfinite(s.gof)
        @test isfinite(s.max_matched_residual)
        @test isfinite(s.max_unmatched_residual)
        @test isfinite(s.fitted_manhattan_distance)
        @test isfinite(s.fitted_hellinger)
    end

    s = _score_single(1.0f0,
                      Float32[2.0, 5.0, 7.0],
                      Float32[1.0, 6.0, 3.0],
                      Bool[true, true, true];
                      ranks=UInt8[1, 2, 1])
    @test s.fitted_frag1_int == 9.0f0
    @test s.fitted_frag2_int == 5.0f0
    @test s.shadow_frag1_int == 4.0f0
    @test s.shadow_frag2_int == 6.0f0
    @test s.fitted_frag3_int == 0.0f0
    @test s.shadow_frag3_int == 0.0f0
end

@testset "Feature finiteness — log2 floors at the PSM emission site" begin
    # These mirror the exact transforms applied in ScoredPSMs.jl:
    #   error feature  = log2(max(error_sum, 1e-20))
    #   intensity feat = log2(max(b_int+y_int, 1e-20) / max(spectrum_int, 1e-20))
    err_floor(x)    = Float16(log2(max(Float32(x), Float32(1e-20))))
    intensity(b_y, spec) = Float16(log2(max(Float32(b_y), Float32(1e-20)) /
                                        max(Float32(spec), Float32(1e-20))))

    # Zero error → finite
    @test isfinite(err_floor(0.0))
    @test isfinite(err_floor(0.0f0))
    @test isfinite(err_floor(1e-25))

    # Typical
    @test isfinite(err_floor(2.5))
    @test isfinite(err_floor(0.001))

    # Zero spectrum intensity (degenerate empty scan)
    @test isfinite(intensity(100.0, 0.0))
    @test isfinite(intensity(0.0, 0.0))
    @test isfinite(intensity(0.0, 1000.0))

    # Typical
    @test isfinite(intensity(100.0, 1e6))
end
