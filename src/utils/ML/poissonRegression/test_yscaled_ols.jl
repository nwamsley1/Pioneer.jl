## Test Y-scaled OLS vs Original OLS — synthetic problems matching real data
##
## Key insight: real spectral deconvolution problems have:
##   - nzval (library spectra) ~ O(1e-3), range [1e-7, 5e-3]
##   - x (observed intensities) ~ O(1e4–1e5)
##   - This makes colnorm2 = Σ A_ij² ~ O(1e-6), so Newton steps = gradient/1e-6 → huge
##   - With w=ones cold start, the back-and-forth across many columns overflows Float32
##
## Run:  julia src/utils/ML/poissonRegression/test_yscaled_ols.jl

using Printf, Statistics, Random

include(joinpath(@__DIR__, "SparseArray.jl"))
include(joinpath(@__DIR__, "solveOLS.jl"))

# ── Helpers ──────────────────────────────────────────────────────

function sse(r::Vector{Float32}, m::Int)
    s = 0.0; @inbounds for i in 1:m; s += Float64(r[i])^2; end; s
end

function compute_sse(sa, w)
    r = zeros(Float32, sa.m)
    @inbounds for n in 1:sa.n_vals
        if iszero(r[sa.rowval[n]]); r[sa.rowval[n]] = -sa.x[n]; end
    end
    @inbounds for col in 1:sa.n
        for n in sa.colptr[col]:(sa.colptr[col+1]-1)
            r[sa.rowval[n]] += w[col] * sa.nzval[n]
        end
    end
    return sse(r, sa.m)
end

"""
    make_spectral_problem(rng, m, n, nnz_per_col, nzval_scale, y_max)

Build a synthetic problem matching real spectral deconvolution characteristics:
  - nzval entries ~ Uniform(0, nzval_scale)  (typically 1e-3)
  - observed x ~ Uniform(0, y_max)           (typically 1e4–1e5)
  - Each column has ~nnz_per_col nonzeros
"""
function make_spectral_problem(rng, m::Int, n::Int, nnz_per_col::Int, nzval_scale::Float32, y_max::Float32)
    # Build sparse structure
    rows = Int64[]
    cols = UInt16[]
    vals = Float32[]
    obs_per_row = zeros(Float32, m)  # track observed value per row
    colptr = Int64[1]

    # Generate observed values
    y = Float32[rand(rng) * y_max for _ in 1:m]

    for j in 1:n
        # Pick random rows for this column
        row_set = sort(randperm(rng, m)[1:min(nnz_per_col, m)])
        for i in row_set
            push!(rows, i)
            push!(cols, UInt16(j))
            push!(vals, Float32(rand(rng)) * nzval_scale)
        end
        push!(colptr, length(rows) + 1)
    end

    # Build x array: observed value for each nonzero entry
    x_arr = Float32[y[rows[i]] for i in 1:length(rows)]

    n_vals = length(rows)
    sa = SparseArray{Int64, Float32}(
        n_vals, m, n,
        rows, cols, vals,
        ones(Bool, n_vals),
        zeros(UInt8, n_vals),
        x_arr,
        colptr
    )
    return sa
end

# ── Test Problems ────────────────────────────────────────────────
# Modeled after real data:
#   nzval ~ [1e-7, 5e-3], x ~ [0, 5e5], nnz_per_col ~ 10-15

problems = [
    # (m,   n,   nnz_per_col, nzval_scale, y_max,       description)
    (50,    7,   12, 5f-3, 1f3,      "7 cols, y~1e3, nzval~5e-3"),
    (50,    7,   12, 5f-3, 1f4,      "7 cols, y~1e4, nzval~5e-3"),
    (50,    7,   12, 5f-3, 1f5,      "7 cols, y~1e5, nzval~5e-3 (real-like)"),
    (50,    7,   12, 5f-3, 1f6,      "7 cols, y~1e6, nzval~5e-3"),
    (100,  12,   12, 5f-3, 3f4,      "12 cols, y~3e4 (real-like small)"),
    (200,  30,   12, 5f-3, 1f5,      "30 cols, y~1e5 (real-like medium)"),
    (500, 100,   14, 3f-3, 5f5,      "100 cols, y~5e5 (real-like large)"),
    (1000, 200,  14, 3f-3, 5f5,      "200 cols, y~5e5 (stress test)"),
    (2600, 304,  12, 5f-3, 6f4,      "304 cols, y~6e4 (matches prob 108289)"),
    # Edge cases
    (50,    7,   12, 5f-3, 1f8,      "7 cols, y~1e8 (extreme)"),
    (200,  30,   12, 1f-3, 1f5,      "30 cols, y~1e5, nzval~1e-3 (smaller A)"),
    (200,  30,   12, 1f-2, 1f5,      "30 cols, y~1e5, nzval~1e-2 (larger A)"),
    # Easy baseline (nzval ~ O(1), like my broken synthetic test)
    (200,  30,   12, 1f0,  1f5,      "30 cols, y~1e5, nzval~1 (WRONG scale)"),
]

max_outer = Int64(1000)
rel_conv  = Float32(0.01)

# Warmup
sa_w = make_spectral_problem(MersenneTwister(0), 20, 5, 8, 5f-3, 1f4)
r_w = zeros(Float32, sa_w.m); w_w = ones(Float32, sa_w.n); cn_w = zeros(Float32, sa_w.n)
initResiduals!(r_w, sa_w, w_w); solveOLS!(sa_w, r_w, w_w, cn_w, Int64(10), Float32(0.01))
w_w .= 1f0; initResiduals!(r_w, sa_w, w_w); solveOLS_yscaled!(sa_w, r_w, w_w, cn_w, Int64(10), Float32(0.01))

# ── Section 1: Cold start w=ones ──────────────────────────────

println("=" ^ 140)
println("  SECTION 1: OLS vs Y-SCALED OLS — Cold start w=ones")
println("  (This is where overflow happens: w=ones → huge initial gradient / tiny L2)")
println("=" ^ 140)
@printf("\n  %-45s %4s %5s │ %12s %12s │ %12s %12s │ %7s\n",
        "Problem", "Cols", "Rows",
        "OLS_SSE", "OLS_maxW",
        "Yscl_SSE", "Yscl_maxW",
        "Status")
println("  " * "─" ^ 130)

n_ols_nan = 0; n_yscl_nan = 0; n_fixed = 0

for (idx, (m, n, nnz_pc, nzs, ym, desc)) in enumerate(problems)
    global n_ols_nan, n_yscl_nan, n_fixed
    sa = make_spectral_problem(MersenneTwister(42 + idx), m, n, nnz_pc, nzs, ym)

    # Original OLS
    r1 = zeros(Float32, sa.m); w1 = ones(Float32, sa.n); cn1 = zeros(Float32, sa.n)
    initResiduals!(r1, sa, w1)
    solveOLS!(sa, r1, w1, cn1, max_outer, rel_conv)
    ols_sse = compute_sse(sa, w1)
    ols_maxw = maximum(abs.(w1[1:sa.n]))

    # Y-scaled OLS
    r2 = zeros(Float32, sa.m); w2 = ones(Float32, sa.n); cn2 = zeros(Float32, sa.n)
    initResiduals!(r2, sa, w2)
    solveOLS_yscaled!(sa, r2, w2, cn2, max_outer, rel_conv)
    ys_sse = compute_sse(sa, w2)
    ys_maxw = maximum(abs.(w2[1:sa.n]))

    o_nan = isnan(ols_sse) || isinf(ols_sse)
    y_nan = isnan(ys_sse) || isinf(ys_sse)
    o_nan && (n_ols_nan += 1)
    y_nan && (n_yscl_nan += 1)
    (o_nan && !y_nan) && (n_fixed += 1)

    status = if o_nan && !y_nan; "FIXED"
    elseif o_nan && y_nan; "BOTH_BAD"
    elseif !o_nan && !y_nan; "OK"
    else; "REGRESS"
    end

    @printf("  %-45s %4d %5d │ %12.4e %12.4e │ %12.4e %12.4e │ %7s\n",
            desc, sa.n, sa.m, ols_sse, ols_maxw, ys_sse, ys_maxw, status)
end

println("\n  " * "─" ^ 100)
@printf("  NaN/Inf SSE:  OLS=%d   Y-scaled=%d   (of %d problems)\n", n_ols_nan, n_yscl_nan, length(problems))
@printf("  FIXED:        %d problems (OLS failed → Y-scaled succeeded)\n", n_fixed)

# ── Section 2: Cold start w=zeros (production) ─────────────────

println("\n\n" * "=" ^ 140)
println("  SECTION 2: OLS vs Y-SCALED OLS — Cold start w=zeros (production init)")
println("=" ^ 140)
@printf("\n  %-45s %4s %5s │ %12s %12s │ %12s %12s │ %7s\n",
        "Problem", "Cols", "Rows",
        "OLS_SSE", "OLS_maxW",
        "Yscl_SSE", "Yscl_maxW",
        "Match?")
println("  " * "─" ^ 130)

for (idx, (m, n, nnz_pc, nzs, ym, desc)) in enumerate(problems)
    sa = make_spectral_problem(MersenneTwister(42 + idx), m, n, nnz_pc, nzs, ym)

    # OLS w=zeros
    r1 = zeros(Float32, sa.m); w1 = zeros(Float32, sa.n); cn1 = zeros(Float32, sa.n)
    initResiduals!(r1, sa, w1)
    solveOLS!(sa, r1, w1, cn1, max_outer, rel_conv)
    ols_sse = compute_sse(sa, w1)
    ols_maxw = maximum(abs.(w1[1:sa.n]))

    # Y-scaled w=zeros
    r2 = zeros(Float32, sa.m); w2 = zeros(Float32, sa.n); cn2 = zeros(Float32, sa.n)
    initResiduals!(r2, sa, w2)
    solveOLS_yscaled!(sa, r2, w2, cn2, max_outer, rel_conv)
    ys_sse = compute_sse(sa, w2)
    ys_maxw = maximum(abs.(w2[1:sa.n]))

    o_nan = isnan(ols_sse) || isinf(ols_sse)
    y_nan = isnan(ys_sse) || isinf(ys_sse)
    # Check if solutions match
    match_str = if o_nan || y_nan
        "NaN"
    else
        rel_diff = abs(ols_sse - ys_sse) / max(ols_sse, 1.0)
        rel_diff < 1e-3 ? "YES" : @sprintf("%.1e", rel_diff)
    end

    @printf("  %-45s %4d %5d │ %12.4e %12.4e │ %12.4e %12.4e │ %7s\n",
            desc, sa.n, sa.m, ols_sse, ols_maxw, ys_sse, ys_maxw, match_str)
end

# ── Section 3: Detailed diagnosis of overflow mechanism ─────────

println("\n\n" * "=" ^ 140)
println("  SECTION 3: Overflow mechanism — per-iteration diagnostics")
println("  (One representative problem that overflows with OLS, comparing both solvers)")
println("=" ^ 140)

# Use a problem that should overflow: 100 cols, y~5e5, nzval~3e-3
sa_diag = make_spectral_problem(MersenneTwister(99), 500, 100, 14, 3f-3, 5f5)
println("\n  Problem: 100 cols, 500 rows, nzval~3e-3, y_max~5e5")

# Check L2 values
cn_diag = zeros(Float32, sa_diag.n)
for col in 1:sa_diag.n
    s = 0f0
    for i in sa_diag.colptr[col]:(sa_diag.colptr[col+1]-1)
        s += sa_diag.nzval[i] * sa_diag.nzval[i]
    end
    cn_diag[col] = s
end
@printf("\n  colnorm2 (L2 = Σ A_ij²): min=%.4e  max=%.4e  mean=%.4e\n",
        minimum(cn_diag[1:sa_diag.n]), maximum(cn_diag[1:sa_diag.n]), mean(cn_diag[1:sa_diag.n]))
println("  This is the denominator in Newton step. When it's ~1e-5, a gradient of 1 → step of 1e5.")

# Run OLS with per-iteration logging
println("\n  --- OLS (unscaled) iteration log ---")
r_d = zeros(Float32, sa_diag.m); w_d = ones(Float32, sa_diag.n)
initResiduals!(r_d, sa_diag, w_d)

for iter in 1:10
    max_diff = 0f0
    max_w = 0f0
    for col in 1:sa_diag.n
        L2 = cn_diag[col]
        iszero(L2) && continue
        L1 = 0f0
        for k in sa_diag.colptr[col]:(sa_diag.colptr[col+1]-1)
            L1 += sa_diag.nzval[k] * r_d[sa_diag.rowval[k]]
        end
        X0 = w_d[col]
        w_d[col] = max(w_d[col] - L1 / L2, 0f0)
        updateResiduals!(sa_diag, r_d, col, w_d[col], X0)
        δ = abs(w_d[col] - X0)
        w_d[col] > max_w && (max_w = w_d[col])
        if w_d[col] > 0
            rc = δ / abs(w_d[col])
            rc > max_diff && (max_diff = rc)
        end
    end
    s = sse(r_d, sa_diag.m)
    @printf("  iter %2d: max_diff=%.4e  max_w=%.4e  SSE=%.4e  any_nan_r=%s  any_nan_w=%s\n",
            iter, max_diff, max_w, s,
            any(isnan, r_d[1:sa_diag.m]) ? "YES" : "no",
            any(isnan, w_d[1:sa_diag.n]) ? "YES" : "no")
    (isnan(max_diff) || max_diff < rel_conv) && break
end

# Now run y-scaled version
println("\n  --- Y-scaled OLS iteration log ---")
r_d2 = zeros(Float32, sa_diag.m); w_d2 = ones(Float32, sa_diag.n)
initResiduals!(r_d2, sa_diag, w_d2)

# Scale
y_scale = 0f0
for n in 1:sa_diag.n_vals; sa_diag.x[n] > y_scale && (y_scale = sa_diag.x[n]); end
@printf("  y_scale = %.4e\n", y_scale)

inv_s = 1f0 / y_scale
for n in 1:sa_diag.n_vals; sa_diag.x[n] *= inv_s; end
for j in 1:sa_diag.n; w_d2[j] *= inv_s; end
initResiduals!(r_d2, sa_diag, w_d2)

for iter in 1:10
    max_diff = 0f0
    max_w = 0f0
    for col in 1:sa_diag.n
        L2 = cn_diag[col]
        iszero(L2) && continue
        L1 = 0f0
        for k in sa_diag.colptr[col]:(sa_diag.colptr[col+1]-1)
            L1 += sa_diag.nzval[k] * r_d2[sa_diag.rowval[k]]
        end
        X0 = w_d2[col]
        w_d2[col] = max(w_d2[col] - L1 / L2, 0f0)
        updateResiduals!(sa_diag, r_d2, col, w_d2[col], X0)
        δ = abs(w_d2[col] - X0)
        w_d2[col] > max_w && (max_w = w_d2[col])
        if w_d2[col] > 0
            rc = δ / abs(w_d2[col])
            rc > max_diff && (max_diff = rc)
        end
    end
    s = sse(r_d2, sa_diag.m)
    @printf("  iter %2d: max_diff=%.4e  max_w=%.4e  SSE=%.4e  any_nan_r=%s  any_nan_w=%s\n",
            iter, max_diff, max_w, s,
            any(isnan, r_d2[1:sa_diag.m]) ? "YES" : "no",
            any(isnan, w_d2[1:sa_diag.n]) ? "YES" : "no")
    (isnan(max_diff) || max_diff < rel_conv) && break
end

# Unscale
for n in 1:sa_diag.n_vals; sa_diag.x[n] *= y_scale; end
for j in 1:sa_diag.n; w_d2[j] *= y_scale; end

ols_final = compute_sse(sa_diag, w_d)
yscl_final = compute_sse(sa_diag, w_d2)
@printf("\n  Final SSE:  OLS=%.4e  Y-scaled=%.4e\n", ols_final, yscl_final)

println("\n" * "=" ^ 140)
println("  DONE")
println("=" ^ 140)
