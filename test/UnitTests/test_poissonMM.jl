# Unit tests for solvePoissonMM_fast! — bisection fallback vs Float64 accumulators
#
# Tests two fixes for the "stuck weight" bug where columns with tiny a_ij values
# and a large warm-start have L2 < ε (curvature too small for Newton), causing
# the solver to leave the weight unchanged.
#
# Uses serialized real-world problems from chromatogram deconvolution:
#   - plateau: warm-start 2.87e9, should converge to ~5e6
#   - good: warm-start 924K, should converge to ~1.06e6

using Test
using Serialization
using Pioneer: SparseArray, initMu!

const DIAG_DIR = "/Users/nathanwamsley/Data/For_Figures/OlsenAstralThreeProteome200ng/OlsenAstralThreeProteome200ng_explorislike_newlib_PMM"
const PLATEAU_FILE = joinpath(DIAG_DIR, "pmm_diag_plateau_rt29.646_iters5.dat")
const GOOD_FILE = joinpath(DIAG_DIR, "pmm_diag_good_rt29.925_iters5.dat")

# ─────────────────────────────────────────────────────────────
# Helper: reconstruct problem from serialized Dict
# ─────────────────────────────────────────────────────────────
function load_problem(path::String)
    p = open(deserialize, path)
    n_rows = p[:n_rows]; n_cols = p[:n_cols]; n_vals = p[:n_vals]

    Hs = SparseArray(n_vals + 100)
    Hs.m = n_rows; Hs.n = n_cols; Hs.n_vals = n_vals
    Hs.colptr[1:n_cols+1] .= p[:colptr]
    Hs.rowval[1:n_vals] .= p[:rowval]
    Hs.nzval[1:n_vals] .= p[:nzval]

    return Hs, p
end

function setup_arrays(Hs, p)
    n_rows = p[:n_rows]; n_cols = p[:n_cols]
    w_pre = p[:w_before_all]

    X = zeros(Float32, max(n_cols + 1, length(Hs.x)))
    X[1:n_cols] .= w_pre
    y = zeros(Float32, max(n_rows + 1, length(Hs.x)))
    y[1:n_rows] .= p[:y]
    mu = zeros(Float32, max(n_rows + 1, length(Hs.x)))
    initMu!(mu, Hs, X)

    return X, y, mu
end

# ─────────────────────────────────────────────────────────────
# Approach 1: Float64 accumulators for L1/L2
# ─────────────────────────────────────────────────────────────
function solvePoissonMM_float64accum!(Hs::SparseArray{Ti, T},
                                       μ::Vector{T},
                                       y::Vector{T},
                                       X₁::Vector{T},
                                       max_iter_outer::Int64,
                                       relative_convergence_threshold::T;
                                       max_inner_iter::Int64 = Int64(5)) where {Ti<:Integer, T<:AbstractFloat}

    # Y-scaling
    y_scale = T(0)
    @inbounds for i in 1:Hs.m
        y[i] > y_scale && (y_scale = y[i])
    end
    if y_scale > T(1)
        @inbounds for i in 1:Hs.m; y[i] /= y_scale; end
        @inbounds for j in 1:Hs.n; X₁[j] /= y_scale; end
        initMu!(μ, Hs, X₁)
    end

    colptr = Hs.colptr; rowval = Hs.rowval; nzval = Hs.nzval; ncols = Hs.n
    max_weight = T(0)
    # Float64 epsilon — much tighter than Float32's 1e-6
    ε_f64 = 1e-15
    ε_mu = T(1e-6)  # floor for μ clamping (stays Float32)

    iter = 0; converged = false
    while iter < max_iter_outer
        _diff = T(0)
        weight_floor = iter >= 5 ? max_weight * T(1e-4) : T(0)
        max_weight = T(0)

        for col in 1:ncols
            col_start = colptr[col]
            col_end   = colptr[col + 1] - 1
            X_before  = X₁[col]

            # ── Float64 accumulation for L1, L2 ──
            L1 = 0.0  # Float64
            L2 = 0.0  # Float64
            @inbounds @fastmath for i in col_start:col_end
                row   = rowval[i]
                a_ij  = Float64(nzval[i])
                inv_μ = 1.0 / max(Float64(μ[row]), Float64(ε_mu))
                y_i   = Float64(y[row])
                L1   += a_ij * (1.0 - y_i * inv_μ)
                L2   += a_ij * a_ij * y_i * inv_μ * inv_μ
            end
            if L2 < ε_f64
                @inbounds @fastmath for i in col_start:col_end
                    a_ij  = Float64(nzval[i])
                    inv_μ = 1.0 / max(Float64(μ[rowval[i]]), Float64(ε_mu))
                    L2   += a_ij * a_ij * inv_μ
                end
            end

            # ── Inner Newton (same structure, but L1/L2 in Float64) ──
            for _k in 1:max_inner_iter
                (L2 <= ε_f64 || isnan(L1)) && break

                X0 = X₁[col]
                X₁[col] = max(T(Float64(X₁[col]) - L1 / L2), zero(T))
                delta = X₁[col] - X0

                done = iszero(X₁[col]) || abs(delta) / max(abs(X₁[col]), T(1e-10)) < T(1e-3)

                if !done && _k < max_inner_iter
                    L1 = 0.0; L2 = 0.0
                    @inbounds @fastmath for i in col_start:col_end
                        row     = rowval[i]
                        a_ij    = Float64(nzval[i])
                        μ[row] += T(a_ij * Float64(delta))
                        inv_μ   = 1.0 / max(Float64(μ[row]), Float64(ε_mu))
                        y_i     = Float64(y[row])
                        L1     += a_ij * (1.0 - y_i * inv_μ)
                        L2     += a_ij * a_ij * y_i * inv_μ * inv_μ
                    end
                    if L2 < ε_f64
                        @inbounds @fastmath for i in col_start:col_end
                            a_ij  = Float64(nzval[i])
                            inv_μ = 1.0 / max(Float64(μ[rowval[i]]), Float64(ε_mu))
                            L2   += a_ij * a_ij * inv_μ
                        end
                    end
                else
                    @inbounds @fastmath for i in col_start:col_end
                        μ[rowval[i]] += nzval[i] * delta
                    end
                    break
                end
            end

            δx = abs(X₁[col] - X_before)
            if X₁[col] > max_weight; max_weight = X₁[col]; end
            if X₁[col] > weight_floor
                rc = δx / max(abs(X₁[col]), T(1e-10))
                rc > _diff && (_diff = rc)
            end
        end

        _diff < relative_convergence_threshold && (converged = true; break)
        iter += 1
    end

    if y_scale > T(1)
        @inbounds for i in 1:Hs.m; y[i] *= y_scale; end
        @inbounds for j in 1:Hs.n; X₁[j] *= y_scale; end
        initMu!(μ, Hs, X₁)
    end

    return (converged, iter)
end

# ─────────────────────────────────────────────────────────────
# Approach 2: Bisection fallback (current implementation in repo)
# ─────────────────────────────────────────────────────────────
# Uses the production solvePoissonMM_fast! which already has bisection.
using Pioneer: solvePoissonMM_fast!

# ─────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────
@testset "PoissonMM solver — stuck weight fixes" begin
    Hs_p, p_plateau = load_problem(PLATEAU_FILE)
    Hs_g, p_good = load_problem(GOOD_FILE)

    target_plateau = p_plateau[:target_col]
    target_good = p_good[:target_col]

    @testset "Bisection fallback (production)" begin
        # Plateau case: weight should move from ~2.87e9 to something reasonable
        X, y, mu = setup_arrays(Hs_p, p_plateau)
        w_before = X[target_plateau]
        conv, iters = solvePoissonMM_fast!(Hs_p, mu, y, X, p_plateau[:max_iter], p_plateau[:max_diff])
        w_after = X[target_plateau]

        println("  Bisection PLATEAU: $(w_before) → $(w_after)  conv=$conv iters=$iters")
        @test conv
        @test w_after != w_before  # must not be stuck
        @test w_after < w_before / 10  # should decrease substantially
        @test w_after < 1e8  # should be in reasonable range

        # Good case: should still converge correctly
        X, y, mu = setup_arrays(Hs_g, p_good)
        w_before = X[target_good]
        conv, iters = solvePoissonMM_fast!(Hs_g, mu, y, X, p_good[:max_iter], p_good[:max_diff])
        w_after = X[target_good]

        println("  Bisection GOOD: $(w_before) → $(w_after)  conv=$conv iters=$iters")
        @test conv
        @test w_after != w_before
        @test 1e5 < w_after < 1e8  # should be in reasonable range
    end

    @testset "Float64 accumulators" begin
        # Plateau case
        X, y, mu = setup_arrays(Hs_p, p_plateau)
        w_before = X[target_plateau]
        conv, iters = solvePoissonMM_float64accum!(Hs_p, mu, y, X, p_plateau[:max_iter], p_plateau[:max_diff])
        w_after = X[target_plateau]

        println("  Float64 PLATEAU: $(w_before) → $(w_after)  conv=$conv iters=$iters")
        @test conv
        @test w_after != w_before
        @test w_after < w_before / 10
        @test w_after < 1e8

        # Good case
        X, y, mu = setup_arrays(Hs_g, p_good)
        w_before = X[target_good]
        conv, iters = solvePoissonMM_float64accum!(Hs_g, mu, y, X, p_good[:max_iter], p_good[:max_diff])
        w_after = X[target_good]

        println("  Float64 GOOD: $(w_before) → $(w_after)  conv=$conv iters=$iters")
        @test conv
        @test w_after != w_before
        @test 1e5 < w_after < 1e8
    end

    @testset "Both approaches agree" begin
        # Run both on plateau, compare results
        X_bis, y_bis, mu_bis = setup_arrays(Hs_p, p_plateau)
        X_f64, y_f64, mu_f64 = setup_arrays(Hs_p, p_plateau)

        solvePoissonMM_fast!(Hs_p, mu_bis, y_bis, X_bis, p_plateau[:max_iter], p_plateau[:max_diff])
        solvePoissonMM_float64accum!(Hs_p, mu_f64, y_f64, X_f64, p_plateau[:max_iter], p_plateau[:max_diff])

        w_bis = X_bis[target_plateau]
        w_f64 = X_f64[target_plateau]

        println("  Agreement PLATEAU: bisection=$(w_bis)  float64=$(w_f64)  ratio=$(w_bis/w_f64)")
        # They should be in the same ballpark (within 2x)
        @test 0.5 < w_bis / w_f64 < 2.0

        # Run both on good case
        X_bis, y_bis, mu_bis = setup_arrays(Hs_g, p_good)
        X_f64, y_f64, mu_f64 = setup_arrays(Hs_g, p_good)

        solvePoissonMM_fast!(Hs_g, mu_bis, y_bis, X_bis, p_good[:max_iter], p_good[:max_diff])
        solvePoissonMM_float64accum!(Hs_g, mu_f64, y_f64, X_f64, p_good[:max_iter], p_good[:max_diff])

        w_bis = X_bis[target_good]
        w_f64 = X_f64[target_good]

        println("  Agreement GOOD: bisection=$(w_bis)  float64=$(w_f64)  ratio=$(w_bis/w_f64)")
        @test 0.5 < w_bis / w_f64 < 2.0
    end
end
