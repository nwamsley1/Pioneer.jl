#=
test_real_data.jl — Validate new pentadiagonal solver against real chromatogram data.

Tests against 200 real WH problems captured from an Astral 3-file run.
Each file contains: x, y, w (inputs), n, λ, and z (reference output from CHOLMOD).

Run:  julia --project=. --threads 10 scripts/whittaker_henderson/test_real_data.jl
=#

using JLD2, Printf

# ============================================================================
# New solver (copy from test_solvers.jl — self-contained)
# ============================================================================

mutable struct WHWorkspace
    d0::Vector{Float64}
    d1::Vector{Float64}
    d2::Vector{Float64}
    ld1::Vector{Float64}
    ld2::Vector{Float64}
    diag::Vector{Float64}
    rhs::Vector{Float64}
    z_f64::Vector{Float64}
    z::Vector{Float32}
    h::Vector{Float64}
    n_max::Int
end

function WHWorkspace(n_max::Int)
    WHWorkspace(
        zeros(Float64, n_max),
        zeros(Float64, n_max - 1),
        zeros(Float64, n_max - 2),
        zeros(Float64, n_max - 1),
        zeros(Float64, n_max - 2),
        zeros(Float64, n_max),
        zeros(Float64, n_max),
        zeros(Float64, n_max),
        zeros(Float32, n_max),
        zeros(Float64, n_max - 1),
        n_max
    )
end

function build_DtD_bands!(ws::WHWorkspace, x::AbstractVector{Float32}, n::Int, λ::Float32)
    h = ws.h
    d0 = ws.d0; d1 = ws.d1; d2 = ws.d2
    λ64 = Float64(λ)

    @inbounds for i in 1:n;     d0[i] = 0.0; end
    @inbounds for i in 1:n-1;   d1[i] = 0.0; h[i] = Float64(x[i+1]) - Float64(x[i]); end
    @inbounds for i in 1:n-2;   d2[i] = 0.0; end

    @inbounds for j in 1:n-2
        hj = h[j]; hj1 = h[j+1]; s = hj + hj1
        a = 1.0 / (hj * s)
        b = -1.0 / (hj * hj1)
        c = 1.0 / (hj1 * s)

        d0[j]   += λ64 * a * a
        d0[j+1] += λ64 * b * b
        d0[j+2] += λ64 * c * c
        d1[j]   += λ64 * a * b
        d1[j+1] += λ64 * b * c
        d2[j]   += λ64 * a * c
    end
end

function whitsmddw_new!(ws::WHWorkspace,
                         x::AbstractVector{Float32},
                         y::AbstractVector{Float32},
                         w::AbstractVector{Float32},
                         n::Int,
                         λ::Float32)
    n >= 3 || error("Need at least 3 points")
    n <= ws.n_max || error("n=$n exceeds workspace $(ws.n_max)")

    d0 = ws.d0; d1 = ws.d1; d2 = ws.d2
    ld1 = ws.ld1; ld2 = ws.ld2; diag = ws.diag
    rhs = ws.rhs; z64 = ws.z_f64; z = ws.z

    build_DtD_bands!(ws, x, n, λ)

    @inbounds for i in 1:n
        d0[i] += Float64(w[i])
        rhs[i] = Float64(w[i]) * Float64(y[i])
    end

    # LDL' factorization
    @inbounds diag[1] = d0[1]
    if n >= 2
        @inbounds begin
            ld1[1] = d1[1] / diag[1]
            diag[2] = d0[2] - ld1[1] * d1[1]
        end
    end
    @inbounds for i in 3:n
        ld2[i-2] = d2[i-2] / diag[i-2]
        ld1[i-1] = (d1[i-1] - ld2[i-2] * ld1[i-2] * diag[i-2]) / diag[i-1]
        diag[i] = d0[i] - ld1[i-1]^2 * diag[i-1] - ld2[i-2]^2 * diag[i-2]
    end

    # Forward solve
    @inbounds z64[1] = rhs[1]
    if n >= 2
        @inbounds z64[2] = rhs[2] - ld1[1] * z64[1]
    end
    @inbounds for i in 3:n
        z64[i] = rhs[i] - ld1[i-1] * z64[i-1] - ld2[i-2] * z64[i-2]
    end

    # Diagonal solve
    @inbounds for i in 1:n
        z64[i] = z64[i] / diag[i]
    end

    # Back solve
    @inbounds for i in n-1:-1:1
        z64[i] -= ld1[i] * z64[i+1]
        if i + 2 <= n
            z64[i] -= ld2[i] * z64[i+2]
        end
    end

    # Convert to Float32 and clamp
    @inbounds for i in 1:n
        z[i] = max(Float32(z64[i]), 0f0)
    end

    return @view z[1:n]
end

# ============================================================================
# Test against real data
# ============================================================================

function run_real_data_tests()
    test_dir = "/Users/nathanwamsley/Desktop/wh_test_cases"
    files = filter(f -> endswith(f, ".jld2"), readdir(test_dir))
    sort!(files)

    n_files = length(files)
    println("=" ^ 70)
    println("Validating pentadiagonal solver against $n_files real chromatogram problems")
    println("=" ^ 70)

    # Pre-allocate workspace for the largest problem
    max_n = 0
    for f in files
        d = load(joinpath(test_dir, f))
        max_n = max(max_n, d["n"])
    end
    ws = WHWorkspace(max_n + 10)
    println("Max problem size: n=$max_n, workspace allocated for $(ws.n_max)")

    # Track statistics
    max_abs_errors = Float64[]
    max_rel_errors = Float64[]
    problem_sizes = Int[]
    n_pass = 0
    n_fail = 0
    worst_file = ""
    worst_err = 0.0

    for f in files
        d = load(joinpath(test_dir, f))
        x = d["x"]::Vector{Float32}
        y = d["y"]::Vector{Float32}
        w = d["w"]::Vector{Float32}
        n = d["n"]::Int
        λ = d["λ"]::Float32
        z_ref = d["z"]::Vector{Float32}

        # Skip problems too small for d=2
        if n < 3
            continue
        end

        z_new = collect(whitsmddw_new!(ws, x, y, w, n, λ))

        max_abs = maximum(abs.(z_ref .- z_new))
        scale = maximum(abs.(z_ref)) + 1f-10
        max_rel = max_abs / scale

        push!(max_abs_errors, max_abs)
        push!(max_rel_errors, max_rel)
        push!(problem_sizes, n)

        # Use relative error — real data has huge intensities (1e6+) so absolute thresholds fail
        if max_rel < 1f-5  # 0.001% relative precision
            n_pass += 1
        else
            n_fail += 1
            @printf("  FAIL %-30s n=%3d  max_abs=%.2e  rel=%.2e\n", f, n, max_abs, max_rel)
        end

        if max_abs > worst_err
            worst_err = max_abs
            worst_file = f
        end
    end

    # Summary
    println("\n" * "-" ^ 70)
    @printf("Results: %d / %d passed (threshold: relative error < 1e-5)\n", n_pass, n_pass + n_fail)
    println("-" ^ 70)

    # Error distribution
    sorted_abs = sort(max_abs_errors)
    println("\nMax absolute error distribution across $(length(sorted_abs)) problems:")
    for (label, pct) in [("min", 0.0), ("p25", 0.25), ("p50", 0.5), ("p75", 0.75),
                          ("p90", 0.90), ("p95", 0.95), ("p99", 0.99), ("max", 1.0)]
        idx = max(1, ceil(Int, pct * length(sorted_abs)))
        idx = min(idx, length(sorted_abs))
        @printf("  %-4s  %.2e\n", label, sorted_abs[idx])
    end

    sorted_rel = sort(max_rel_errors)
    println("\nMax relative error distribution:")
    for (label, pct) in [("min", 0.0), ("p50", 0.5), ("p95", 0.95), ("p99", 0.99), ("max", 1.0)]
        idx = clamp(ceil(Int, pct * length(sorted_rel)), 1, length(sorted_rel))
        @printf("  %-4s  %.2e\n", label, sorted_rel[idx])
    end

    @printf("\nWorst case: %s (max_abs=%.2e, rel=%.2e)\n", worst_file, worst_err,
            maximum(max_rel_errors))

    # Problem size distribution
    sorted_n = sort(problem_sizes)
    @printf("\nProblem size distribution: min=%d  median=%d  max=%d\n",
            sorted_n[1], sorted_n[length(sorted_n)÷2], sorted_n[end])

    # Benchmark
    println("\n" * "-" ^ 70)
    println("Benchmark: solve all $n_files problems")
    println("-" ^ 70)

    # Warm up
    d = load(joinpath(test_dir, files[1]))
    whitsmddw_new!(ws, d["x"], d["y"], d["w"], d["n"], d["λ"])

    alloc = @allocated begin
        t = @elapsed for f in files
            d = load(joinpath(test_dir, f))
            x = d["x"]::Vector{Float32}
            y = d["y"]::Vector{Float32}
            w = d["w"]::Vector{Float32}
            n = d["n"]::Int
            λ = d["λ"]::Float32
            whitsmddw_new!(ws, x, y, w, n, λ)
        end
    end

    # Time just the solver (excluding JLD2 load)
    # Pre-load all data
    all_data = [(load(joinpath(test_dir, f))) for f in files]

    alloc_solve = @allocated begin
        t_solve = @elapsed for d in all_data
            whitsmddw_new!(ws, d["x"], d["y"], d["w"], d["n"], d["λ"])
        end
    end

    @printf("  Solver only: %.3f ms for %d problems (%.1f μs/call)\n",
            t_solve * 1000, n_files, t_solve * 1e6 / n_files)
    @printf("  Solver alloc: %.1f KB total (%.1f bytes/call)\n",
            alloc_solve / 1e3, alloc_solve / n_files)
    @printf("  Extrapolated to 170k calls: %.2f s, %.1f MB\n",
            t_solve / n_files * 170_000, alloc_solve / 1e6 / n_files * 170_000)

    println("\n" * "=" ^ 70)
    if n_fail == 0
        println("ALL $(n_pass) REAL-DATA TESTS PASSED")
    else
        println("$n_fail / $(n_pass + n_fail) TESTS FAILED")
    end
    println("=" ^ 70)

    return n_fail == 0
end

run_real_data_tests()
