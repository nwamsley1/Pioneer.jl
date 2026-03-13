# CLAUDE.md — Poisson Regression Solver Development

## Current State (2026-03-13)

### What exists

Four coordinate-descent solvers for non-negative Poisson GLM spectral deconvolution:

- **OLS** (`solveOLS!` in `spectralLinearRegression_reference.jl`): Fastest. One dot-product + update per coordinate. Minimizes SSE.
- **Huber** (`solveHuber!` in `spectralLinearRegression_reference.jl`): Newton-bisection per coordinate. ~5x slower than OLS. Produces near-identical solutions to OLS (cor≥0.98) because the Huber δ is large enough that all residuals stay in the quadratic zone.
- **PMM** (`solvePoissonMM!` in `spectralPoissonRegression.jl`): Observed-Hessian coordinate descent maximizing Poisson log-likelihood. Uses y-scaling (divide y by max(y)) to normalize step sizes. Supports inner Newton iterations (K steps per coordinate per sweep).
- **PMM fast** (`solvePoissonMM_fast!` in `spectralPoissonRegression.jl`): Optimized PMM with identical algorithm. Inlines derivative computation and μ update, then fuses the μ-update from inner iteration K with the derivative computation for iteration K+1 into a single pass over each column's nonzeros. This reduces memory passes from 2K to K+1 for K inner iterations. Zero allocations, same signature as `solvePoissonMM!`.

All solvers share:
- `SparseArray` CSC structure (`SparseArray.jl`)
- Weight-floor heuristic: after 5 iterations, columns below `max_weight * threshold` are excluded from convergence check (OLS/Huber use 1e-7, PMM uses 1e-4)

### PMM fast optimization details (added 2026-03-13)

`solvePoissonMM_fast!` applies one optimization to the inner loop of `solvePoissonMM!`:

**Fused μ-update + derivative computation.** In the original solver, each inner Newton iteration requires two passes over the column's nonzeros: one to compute derivatives (L1, L2), one to update μ. The fast solver fuses the μ-update from iteration K with the derivative computation for iteration K+1 into a single pass. For K=5 inner iterations this reduces passes from 10 to 6 (40% fewer memory accesses in the hot loop).

Other optimizations considered and rejected:
- **Precomputed a_ij²**: Would save one multiply per nonzero in L2 but requires allocating a `Vector{T}(n_vals)` buffer. Rejected — no allocations allowed in the solver hot path.
- **Column skipping** (skip zero columns after 2+ sweeps): Would save derivative computation for stable zero columns but requires allocating a `Vector{Int8}(n_cols)` tracker. Rejected for the same reason.
- **Fast reciprocal** (Quake-style bit trick for 1/μ): Profiling showed division was only ~5% of derivative time (~19/370 samples). Expected gain ~3% total, not worth reduced accuracy.

### Benchmark results (100 real problems, 11 runs each, median timing)

| Solver | Median (μs) | Mean (μs) | vs OLS | vs PMM orig |
|--------|------------|-----------|--------|-------------|
| OLS | 2.3 | 11.0 | 1.0x | — |
| PMM fast | 9.4 | 28.4 | 4.5x | 1.37x faster |
| PMM original | 12.8 | 37.7 | 6.1x | 1.0x |

Solution agreement between PMM original and PMM fast: correlation = 1.000000 on all 100 problems, relative LL difference < 2.5e-8 (Float32 rounding only).

### Inner Newton iterations (added 2026-03-13)

All three PMM solvers (`solvePoissonMM!`, `solvePoissonMM_v2!`, `solvePoissonMM_opt!`) now accept `max_inner_iter` kwarg (default 5). Per coordinate per outer sweep, up to K observed-Hessian Newton steps are taken, with early exit when the weight hits zero or the relative change < 1e-3.

**Key finding: total Newton steps are constant across K.** K only controls how many outer sweeps those steps are packed into. Median total Newton steps ≈ 270-290 regardless of K=1,3,5,10.

| K | Outer iters (median) | Total Newton (median) | Wall time (total, 100 problems) |
|---|---|---|---|
| 1 | 22 | 288 | 0.012s |
| 3 | 9 | 299 | 0.004s |
| 5 | 7 | 285 | 0.004s |
| 10 | 5 | 291 | 0.004s |

**Recommended: K=5** (default). Cuts outer iterations ~3x vs K=1, wall time ~2.7x faster due to reduced per-sweep overhead. No LL penalty.

### Convergence analysis

#### Weight-based vs LL-based convergence

Tested both approaches across 100 real problems:
- **Weight-based**: `max(relative_weight_change) < rel_conv` per sweep
- **LL-based**: `|LL_new - LL_old| / (|LL_new| + 1) < ll_tol` per sweep (costs one extra `poissonLogLikelihood` call)

LL-based convergence at `ll_tol=1e-6` with K=5 takes median 6 outer iterations (vs 7 for weight-based at `rel_conv=0.001`). Both produce equivalent solutions.

#### Weight floor: 1e-4 (verified, no change needed)

Tested weight floors 1e-2, 1e-3, 1e-4. All three produce identical results: same iteration counts (median 7), same Newton steps (median 285), same LL, same wall time. Only 4/100 problems differ by 1-2 iterations. The floor only affects when the convergence *check* triggers — the actual weights converge regardless. Kept at 1e-4 (most conservative, ensures even small weights pass the convergence check, zero cost).

#### EXT_UNSTABLE was a test artifact (CRITICAL FINDING)

The "extended run" stability test (run +200 iters from converged solution, check if weights move >1%) was reporting 31/100 problems as unstable. Deep investigation revealed:

1. **The solver converges fully in 6-14 outer iterations** (median 8) for ALL 24 "unstable" problems. By iteration 15-30, `diff` reaches Float32 machine epsilon (~1e-7) and stays there through iteration 2000.

2. **The instability was caused by the test harness**, not the solver. The `extended_pmm_inner_run!` function re-applies y-scaling from scratch on the already-converged solution. This double scale/unscale round-trip in Float32 introduces ~1-3% numerical jitter on small weights. The solver itself is perfectly stable.

3. **L2 regularization does not help** (and actively hurts). Even λ=1e-6 in y-scaled space degrades LL by 2.4% median (99/100 problems worse). λ=1e-4 loses 28% of LL. The problem was never insufficient curvature — it was a measurement artifact.

### Production recommendations

Based on the full benchmark suite:

1. **Use K=5 inner iterations** (default, already set)
2. **Use weight floor = 1e-4** (original value, verified correct)
3. **Weight-based convergence with `rel_conv=0.001`** is sufficient and cheaper than LL-based (no extra data pass). LL-based at 1e-6 is equivalent but adds ~15% overhead.
4. **No L2 regularization** — it degrades solution quality with no convergence benefit.

### Key files

| File | Purpose |
|------|---------|
| `spectralPoissonRegression.jl` | PMM solvers + Poisson LL + helpers (production code) |
| `spectralLinearRegression_reference.jl` | OLS, Huber solvers (reference implementations) |
| `SparseArray.jl` | Minimal CSC sparse matrix struct |
| `benchmark_solvers.jl` | 100-problem timing + solution comparison (OLS, PMM, Huber) |
| `benchmark_pmm_optimization.jl` | PMM original vs PMM fast A/B benchmark |
| `benchmark_inner_iterations.jl` | K=1,3,5,10 inner iteration benchmark |
| `benchmark_ll_convergence.jl` | K × LL_tol grid benchmark |
| `benchmark_l2_reg.jl` | L2 regularization benchmark (conclusion: don't use) |
| `diagnose_convergence.jl` | Original convergence diagnostic |
| `diagnose_ext_unstable.jl` | Per-weight analysis of "unstable" weights |
| `diagnose_slow_convergence.jl` | 2000-iter trace proving full convergence |
| `profile_solvers.jl` | CPU profiling with flame graphs |
| `FINDINGS.md` | Detailed results and design decisions |

### Problem data

100 serialized `.jls` files in `/Users/n.t.wamsley/Desktop/solveHuber_problems/`. Each contains a Dict with fields: `:x`, `:rowval`, `:colval`, `:nzval`, `:colptr`, `:m`, `:n`, `:n_vals`, `:max_iter_outer`, `:max_diff`, `:huber_delta`, `:lambda`, `:max_iter_newton`, `:max_iter_bisection`, `:scan_idx`.

**Important `load_sa` detail:** `data[:x]` is per-row (length m), but `SparseArray.x` is per-nonzero (length n_vals). Must expand: `x_nz[k] = x_row[rowval[k]]`. All scripts include this fix.
