using Pioneer
using BenchmarkTools
using SparseArrays
using Random

# Create test data
function create_test_problem(n_cols=1000, n_rows=5000, density=0.1)
    Random.seed!(123)
    
    # Create sparse matrix
    Hs = Pioneer.SparseArray(
        rand(1:n_rows, Int(n_cols * n_rows * density)),
        rand(1:n_cols, Int(n_cols * n_rows * density)),
        rand(Float32, Int(n_cols * n_rows * density)),
        n_rows,
        n_cols
    )
    
    # Create initial weights (mix of zeros and non-zeros)
    X₁ = zeros(Float32, n_cols)
    for i in 1:n_cols÷2
        X₁[rand(1:n_cols)] = 10.0f0^(rand() * 8 - 2)  # Range from 1e-2 to 1e6
    end
    
    # Create residuals
    r = rand(Float32, n_rows) * 100
    
    return Hs, r, X₁
end

# Test different configurations
function benchmark_configurations()
    println("Creating test problem...")
    Hs, r, X₁ = create_test_problem(1000, 5000, 0.1)
    
    # Original version (comment out the new code and uncomment this for comparison)
    # function solve_original!(Hs, r, X₁)
    #     Pioneer.solveHuber!(Hs, r, copy(X₁), 
    #         Float32(100), Float32(0.001), 50, 100, 1000,
    #         Float32(0.01), Float32(0.01), Float32(0.01), Float32(0.05),
    #         Pioneer.NoNorm())
    # end
    
    # Current version with adaptive tolerance
    function solve_adaptive!(Hs, r, X₁)
        Pioneer.solveHuber!(Hs, r, copy(X₁), 
            Float32(100), Float32(0.001), 50, 100, 1000,
            Float32(0.01), Float32(0.01), Float32(0.01), Float32(0.05),
            Pioneer.NoNorm())
    end
    
    # Version without adaptive tolerance (first iteration only)
    function solve_simple!(Hs, r, X₁)
        X₁_copy = copy(X₁)
        # Manually set relative tolerance in newton_bisection calls
        # This would require modifying the code to accept a fixed tolerance
        Pioneer.solveHuber!(Hs, r, X₁_copy, 
            Float32(100), Float32(0.001), 50, 100, 1,  # Only 1 outer iteration
            Float32(0.01), Float32(0.01), Float32(0.01), Float32(0.05),
            Pioneer.NoNorm())
    end
    
    println("\nBenchmarking adaptive tolerance version:")
    @btime solve_adaptive!($Hs, $r, $X₁)
    
    println("\nBenchmarking single iteration:")
    @btime solve_simple!($Hs, $r, $X₁)
    
    # Profile the adaptive version
    println("\nProfiling adaptive tolerance version...")
    X₁_copy = copy(X₁)
    @profview Pioneer.solveHuber!(Hs, r, X₁_copy, 
        Float32(100), Float32(0.001), 50, 100, 1000,
        Float32(0.01), Float32(0.01), Float32(0.01), Float32(0.05),
        Pioneer.NoNorm())
end

# Test specific bottlenecks
function test_bottlenecks()
    println("Testing potential bottlenecks...")
    
    # Test 1: Clamp function
    x = rand(Float32, 1000000)
    println("\nClamp benchmark:")
    @btime clamp.($x, 0.0f0, 1.0f0)
    
    # Test 2: Division vs multiplication
    a = rand(Float32, 1000000)
    b = 2.5f0
    println("\nDivision:")
    @btime $a ./ $b
    println("Multiplication by inverse:")
    @btime $a .* (1.0f0 / $b)
    
    # Test 3: Conditional branches
    weights = 10.0f0 .^ (rand(Float32, 1000000) * 8 .- 2)
    threshold = 1.0f0
    println("\nConditional branches:")
    @btime for w in $weights
        if w < $threshold
            x = 0.1f0
        else
            x = 0.001f0
        end
    end
end

# Analyze convergence behavior
function analyze_convergence()
    println("\nAnalyzing convergence behavior...")
    Hs, r, X₁ = create_test_problem(100, 500, 0.1)
    
    # Track iterations
    iterations = Pioneer.solveHuber!(Hs, r, X₁, 
        Float32(100), Float32(0.001), 50, 100, 1000,
        Float32(0.01), Float32(0.01), Float32(0.01), Float32(0.05),
        Pioneer.NoNorm())
    
    println("Total outer iterations: $iterations")
    
    # Check weight distribution
    non_zero_weights = X₁[X₁ .> 0]
    if length(non_zero_weights) > 0
        println("Weight statistics:")
        println("  Min: $(minimum(non_zero_weights))")
        println("  Max: $(maximum(non_zero_weights))")
        println("  Dynamic range: $(maximum(non_zero_weights) / minimum(non_zero_weights))")
    end
end

# Main execution
println("=== Testing Spectral Linear Regression Performance ===")
benchmark_configurations()
test_bottlenecks()
analyze_convergence()