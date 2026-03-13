# CLAUDE.md — Poisson Regression Solver Development

## Current State (2026-03-13)

### What exists

Three coordinate-descent solvers for non-negative Poisson GLM spectral deconvolution:

- **OLS** (`solveOLS!` in `spectralLinearRegression_reference.jl`): Fastest. One dot-product + update per coordinate. Minimizes SSE.
- **Huber** (`solveHuber!` in `spectralLinearRegression_reference.jl`): Newton-bisection per coordinate. ~5x slower than OLS. Produces near-identical solutions to OLS (cor≥0.98) because the Huber δ is large enough that all residuals stay in the quadratic zone.
- **PMM** (`solvePoissonMM!` in `spectralPoissonRegression.jl`): Observed-Hessian coordinate descent maximizing Poisson log-likelihood. Uses y-scaling (divide y by max(y)) to normalize step sizes. ~3x more iterations than OLS/Huber.

All solvers share:
- `SparseArray` CSC structure (`SparseArray.jl`)
- Same convergence criterion: `max(relative_weight_change) < rel_conv` per sweep
- Weight-floor heuristic: after 5 iterations, columns below `max_weight * threshold` are excluded from convergence check (OLS/Huber use 1e-7, PMM uses 1e-4)

### What we just did

#### 1. 100-problem benchmark (`benchmark_solvers.jl`)
Ran OLS, PMM, Huber on 100 real-world problems from `/Users/nathanwamsley/Desktop/solveHuber_problems/`. Results in `FINDINGS.md` under "100-Problem Benchmark". Key findings:
- OLS fastest, Huber adds no value over OLS
- PMM produces fundamentally different solutions (cor with Huber mean=0.68, min=-0.56)
- PMM always has higher SSE but higher Poisson LL (different objective)

#### 2. Convergence diagnostic (`diagnose_convergence.jl`)
Wrote and ran a diagnostic to determine whether PMM's poor Huber correlation is convergence failure or a genuine objective difference. For each problem it reports:
- Final `_diff` value and convergence status for all 3 solvers
- Poisson LL for all 3 solvers (OLS and Huber too, not just PMM)
- Extended PMM run: +200 iterations from "converged" PMM solution
- Flags: PMM_NOCONV, PMM_LL<OLS, PMM_ANTI (cor<0), EXT_UNSTABLE (extended run changes >1%)

**Results:**
- **All 100 problems converge** for all 3 solvers (none hit max_iter)
- **PMM always has highest Poisson LL** (100/100) — NOT convergence failure
- **31/100 problems are EXT_UNSTABLE**: PMM declares convergence (passes `_diff < 0.01`) but running 200 more iterations changes weights by up to 35%
- Only 1 problem has PMM anti-correlated with Huber

**Diagnosis:** The `rel_conv = 0.01` threshold is too loose for PMM. A single "quiet" sweep can pass the threshold transiently while the solution is still far from the true optimum. OLS doesn't have this problem because each coordinate step is exact (quadratic objective).

#### 3. Tightened PMM convergence (current change)
In `diagnose_convergence.jl`, changed PMM's convergence threshold from `0.01` (loaded from problem data) to `0.001` via `PMM_REL_CONV = Float32(0.001)`. OLS and Huber still use the original `0.01`. This has NOT been run yet — needs re-execution to see impact on iteration count and EXT_UNSTABLE count.

### What to do next

1. **Run the updated diagnostic** to see how tightening PMM rel_conv to 0.001 affects:
   - PMM iteration counts (expect increase from ~21 median)
   - EXT_UNSTABLE count (expect decrease from 31/100)
   - Whether PMM still always has highest Poisson LL

   ```bash
   julia --project=. src/utils/ML/poissonRegression/diagnose_convergence.jl
   ```

2. **If 0.001 is still insufficient**, consider:
   - Even tighter threshold (0.0001)
   - Objective-based convergence: stop when `|LL_new - LL_old| / (|LL_new| + 1) < tol` instead of max weight change. Already tested on scan342335 (see FINDINGS.md "Objective-Based Convergence" section) — cleaner but costs one extra data pass per sweep.

3. **Production integration decision**: Once convergence is solid, decide whether to:
   - Replace Huber with PMM in production (different solutions — needs end-to-end peptide quantification evaluation)
   - Use OLS warm-start → PMM refinement (5 OLS iters + PMM converges ~2-3x faster, see FINDINGS.md)
   - Keep Huber as default, offer PMM as option

### Key files

| File | Purpose |
|------|---------|
| `spectralPoissonRegression.jl` | PMM solver + Poisson LL + helper functions |
| `spectralLinearRegression_reference.jl` | OLS, Huber solvers (reference implementations) |
| `SparseArray.jl` | Minimal CSC sparse matrix struct |
| `benchmark_solvers.jl` | 100-problem timing + solution comparison |
| `diagnose_convergence.jl` | Convergence diagnostic (current focus) |
| `profile_solvers.jl` | CPU profiling with flame graphs |
| `FINDINGS.md` | Detailed results and design decisions |

### Problem data

100 serialized `.jls` files in `/Users/nathanwamsley/Desktop/solveHuber_problems/`. Each contains a Dict with fields: `:x`, `:rowval`, `:colval`, `:nzval`, `:colptr`, `:m`, `:n`, `:n_vals`, `:max_iter_outer`, `:max_diff`, `:huber_delta`, `:lambda`, `:max_iter_newton`, `:max_iter_bisection`, `:scan_idx`.

**Important `load_sa` detail:** `data[:x]` is per-row (length m), but `SparseArray.x` is per-nonzero (length n_vals). Must expand: `x_nz[k] = x_row[rowval[k]]`. All scripts include this fix.
