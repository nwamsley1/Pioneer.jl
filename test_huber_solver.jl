#!/usr/bin/env julia

"""
Standalone test script for profiling the Huber solver implementation.
Loads a saved problem from Desktop and runs various performance tests.
"""

using JLD2
using BenchmarkTools
using Profile
using PProf
using LinearAlgebra
using SparseArrays

# Include the source file
include("src/utils/ML/spectralLinearRegression.jl")
include("src/structs/SparseArray.jl")

function load_test_problem(filepath::String)
    println("Loading test problem from: $filepath")
    data = load(filepath)
    
    # Reconstruct the SparseArray
    Hs = SparseArray(
        data["Hs_rowval"],
        data["Hs_colval"],
        data["Hs_nzval"],
        data["Hs_colptr"],
        data["Hs_n_vals"],
        data["Hs_n"],
        data["Hs_m"],
        data["Hs_x"]
    )
    
    r = data["r_initial"]
    X₁ = data["X1_initial"]
    δ = data["delta"]
    λ = data["lambda"]
    accuracy_newton = data["accuracy_newton"]
    accuracy_bisection = data["accuracy_bisection"]
    max_diff = data["max_diff"]
    
    # Reconstruct regularization type
    reg_type_str = data["regularization_type"]
    reg_type = if reg_type_str == "NoNorm"
        NoNorm()
    elseif reg_type_str == "L1Norm"
        L1Norm()
    elseif reg_type_str == "L2Norm"
        L2Norm()
    else
        NoNorm()
    end
    
    println("Problem size: $(Hs.n) variables, $(Hs.m) constraints")
    println("Non-zeros: $(Hs.n_vals)")
    println("Delta: $δ, Lambda: $λ")
    
    return Hs, r, X₁, δ, λ, accuracy_newton, accuracy_bisection, max_diff, reg_type
end

function benchmark_solver(Hs, r, X₁, δ, λ, accuracy_newton, accuracy_bisection, max_diff, reg_type)
    println("\n=== Running Benchmarks ===")
    
    # Create copies for multiple runs
    r_copy = copy(r)
    X₁_copy = copy(X₁)
    
    # Warmup run
    println("Warmup run...")
    solveHuber!(Hs, copy(r), copy(X₁), δ, λ, 100, 100, 1000, 
                accuracy_newton, accuracy_bisection, 1e-6, max_diff, reg_type;
                debug_capture = false)
    
    # Benchmark
    println("\nBenchmarking (5 samples)...")
    result = @benchmark solveHuber!(
        $Hs, r_bench, X_bench, $δ, $λ, 100, 100, 1000,
        $accuracy_newton, $accuracy_bisection, 1e-6, $max_diff, $reg_type;
        debug_capture = false
    ) setup=(r_bench=copy($r_copy); X_bench=copy($X₁_copy)) samples=5 evals=1
    
    display(result)
    
    return result
end

function profile_solver(Hs, r, X₁, δ, λ, accuracy_newton, accuracy_bisection, max_diff, reg_type)
    println("\n=== Profiling Solver ===")
    
    # Clear previous profile data
    Profile.clear()
    
    # Profile the solver
    @profile begin
        for i in 1:3  # Run 3 times to get more samples
            r_prof = copy(r)
            X_prof = copy(X₁)
            solveHuber!(Hs, r_prof, X_prof, δ, λ, 100, 100, 1000,
                        accuracy_newton, accuracy_bisection, 1e-6, max_diff, reg_type;
                        debug_capture = false)
        end
    end
    
    # Generate profile report
    Profile.print(format=:flat, sortedby=:count, mincount=10)
    
    # Save pprof output
    pprof_file = "/Users/nathanwamsley/Desktop/huber_solver_profile.pb.gz"
    pprof(; out=pprof_file)
    println("\nProfile saved to: $pprof_file")
    println("View with: pprof -http=localhost:8080 $pprof_file")
end

function test_convergence_behavior(Hs, r, X₁, δ, λ, accuracy_newton, accuracy_bisection, max_diff, reg_type)
    println("\n=== Testing Convergence Behavior ===")
    
    # Test with different dynamic ranges
    dynamic_ranges = [1e2, 1e3, 1e4, 1e5]
    
    for dr in dynamic_ranges
        println("\nTesting with dynamic range: $dr")
        
        # Modify solveHuber to use this dynamic range
        # For now, just run with default and measure iterations
        r_test = copy(r)
        X_test = copy(X₁)
        
        start_time = time()
        iters = solveHuber!(Hs, r_test, X_test, δ, λ, 100, 100, 1000,
                           accuracy_newton, accuracy_bisection, 1e-6, max_diff, reg_type;
                           debug_capture = false)
        elapsed = time() - start_time
        
        # Calculate some statistics
        n_nonzero = sum(abs.(X_test) .> 1e-10)
        max_weight = maximum(abs.(X_test))
        min_nonzero_weight = minimum(abs.(X_test[abs.(X_test) .> 1e-10]))
        
        println("  Iterations: $iters")
        println("  Time: $(round(elapsed, digits=3))s")
        println("  Non-zero weights: $n_nonzero / $(length(X_test))")
        println("  Weight range: $min_nonzero_weight to $max_weight")
    end
end

function main()
    problem_file = "/Users/nathanwamsley/Desktop/huber_test_problem.jld2"
    
    if !isfile(problem_file)
        println("ERROR: Test problem not found at $problem_file")
        println("Please run SearchDIA first to generate the test problem.")
        return
    end
    
    # Load the problem
    Hs, r, X₁, δ, λ, accuracy_newton, accuracy_bisection, max_diff, reg_type = load_test_problem(problem_file)
    
    # Run benchmarks
    benchmark_result = benchmark_solver(Hs, r, X₁, δ, λ, accuracy_newton, accuracy_bisection, max_diff, reg_type)
    
    # Profile the solver
    profile_solver(Hs, r, X₁, δ, λ, accuracy_newton, accuracy_bisection, max_diff, reg_type)
    
    # Test convergence behavior
    test_convergence_behavior(Hs, r, X₁, δ, λ, accuracy_newton, accuracy_bisection, max_diff, reg_type)
    
    println("\n=== Test Complete ===")
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end