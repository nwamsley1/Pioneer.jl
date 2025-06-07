using Pioneer
using Random

# Create a minimal test case
function test_minimal()
    Random.seed!(123)
    
    # Small problem
    n_cols = 10
    n_rows = 20
    
    # Create sparse matrix
    row_indices = Int32[]
    col_indices = Int32[]
    values = Float32[]
    
    for col in 1:n_cols
        for _ in 1:5  # 5 entries per column
            push!(row_indices, rand(1:n_rows))
            push!(col_indices, col)
            push!(values, rand(Float32))
        end
    end
    
    Hs = Pioneer.SparseArray(row_indices, col_indices, values, n_rows, n_cols)
    
    # Test case 1: All zeros
    println("Test 1: All zero initial weights")
    X₁ = zeros(Float32, n_cols)
    r = rand(Float32, n_rows)
    
    iters = Pioneer.solveHuber!(Hs, r, X₁, 
        Float32(100), Float32(0.001), 50, 100, 10,
        Float32(0.01), Float32(0.01), Float32(0.01), Float32(0.05),
        Pioneer.NoNorm())
    
    println("  Iterations: $iters")
    println("  Non-zero weights: $(sum(X₁ .> 0))")
    println("  Max weight: $(maximum(X₁))")
    
    # Test case 2: Mixed weights
    println("\nTest 2: Mixed initial weights")
    X₁ = zeros(Float32, n_cols)
    X₁[1:5] .= [1e-3, 1e-1, 1.0, 1e2, 1e4]
    r = rand(Float32, n_rows)
    
    iters = Pioneer.solveHuber!(Hs, r, X₁, 
        Float32(100), Float32(0.001), 50, 100, 10,
        Float32(0.01), Float32(0.01), Float32(0.01), Float32(0.05),
        Pioneer.NoNorm())
    
    println("  Iterations: $iters")
    println("  Weight range: $(minimum(X₁[X₁ .> 0])) to $(maximum(X₁))")
end

# Test with timing
function test_with_timing()
    Random.seed!(123)
    
    # Medium problem
    n_cols = 100
    n_rows = 200
    
    # Create sparse matrix
    row_indices = Int32[]
    col_indices = Int32[]
    values = Float32[]
    
    for col in 1:n_cols
        for _ in 1:10
            push!(row_indices, rand(1:n_rows))
            push!(col_indices, col)
            push!(values, rand(Float32))
        end
    end
    
    Hs = Pioneer.SparseArray(row_indices, col_indices, values, n_rows, n_cols)
    X₁ = rand(Float32, n_cols) .* 100
    r = rand(Float32, n_rows)
    
    println("\nTiming test:")
    println("Problem size: $n_cols columns, $n_rows rows")
    
    # Time the solve
    start_time = time()
    iters = Pioneer.solveHuber!(Hs, r, X₁, 
        Float32(100), Float32(0.001), 50, 100, 100,
        Float32(0.01), Float32(0.01), Float32(0.01), Float32(0.05),
        Pioneer.NoNorm())
    elapsed = time() - start_time
    
    println("  Time: $(round(elapsed, digits=3)) seconds")
    println("  Iterations: $iters")
    println("  Time per iteration: $(round(elapsed/iters * 1000, digits=2)) ms")
end

# Check if the issue is with the relative tolerance calculation
function test_tolerance_calculation()
    println("\nTesting tolerance calculation overhead:")
    
    # Simulate the inner loop calculations
    n = 1_000_000
    weights = 10.0f0 .^ (rand(Float32, n) * 8 .- 2)  # 1e-2 to 1e6
    max_weight = maximum(weights)
    min_threshold = max_weight / 1e4f0
    weight_range = max_weight - min_threshold
    base_tol = 0.01f0
    tol_scale_factor = 0.99f0 * base_tol
    
    # Version 1: Current implementation
    function calc_tol_current(w, min_threshold, weight_range, base_tol, tol_scale_factor)
        if w < min_threshold
            return base_tol
        else
            weight_scale = (w - min_threshold) / weight_range
            weight_scale = clamp(weight_scale, 0.0f0, 1.0f0)
            return base_tol - tol_scale_factor * weight_scale
        end
    end
    
    # Version 2: Without clamp
    function calc_tol_no_clamp(w, min_threshold, weight_range, base_tol, tol_scale_factor)
        if w < min_threshold
            return base_tol
        else
            weight_scale = (w - min_threshold) / weight_range
            return base_tol - tol_scale_factor * weight_scale
        end
    end
    
    # Version 3: Simple fixed tolerance
    function calc_tol_fixed(w, min_threshold, weight_range, base_tol, tol_scale_factor)
        return base_tol
    end
    
    # Benchmark
    println("  Current implementation:")
    @time for w in weights
        calc_tol_current(w, min_threshold, weight_range, base_tol, tol_scale_factor)
    end
    
    println("  Without clamp:")
    @time for w in weights
        calc_tol_no_clamp(w, min_threshold, weight_range, base_tol, tol_scale_factor)
    end
    
    println("  Fixed tolerance:")
    @time for w in weights
        calc_tol_fixed(w, min_threshold, weight_range, base_tol, tol_scale_factor)
    end
end

# Run tests
println("=== Debugging Spectral Linear Regression ===")
test_minimal()
test_with_timing()
test_tolerance_calculation()