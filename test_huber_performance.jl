#!/usr/bin/env julia

"""
Simple performance test for the Huber solver.
Run this after capturing a problem with SearchDIA.
"""

using JLD2
using Statistics
using Printf

# Add Pioneer source to load path
push!(LOAD_PATH, dirname(@__DIR__))

# Load required modules
using Pioneer

function test_huber_performance(problem_file::String="/Users/nathanwamsley/Desktop/huber_test_problem.jld2")
    if !isfile(problem_file)
        println("ERROR: Test problem not found at $problem_file")
        println("Please run SearchDIA first to generate the test problem.")
        return
    end
    
    println("Loading test problem...")
    data = load(problem_file)
    
    # Reconstruct the SparseArray
    # Create an empty SparseArray and fill its fields
    N = length(data["Hs_rowval"])
    Hs = Pioneer.SparseArray(N)
    Hs.n_vals = data["Hs_n_vals"]
    Hs.m = data["Hs_m"]
    Hs.n = data["Hs_n"]
    Hs.rowval = data["Hs_rowval"]
    Hs.colval = data["Hs_colval"]
    Hs.nzval = data["Hs_nzval"]
    Hs.colptr = data["Hs_colptr"]
    Hs.x = data["Hs_x"]
    # Load matched and isotope if available, otherwise use defaults
    if haskey(data, "Hs_matched")
        Hs.matched = data["Hs_matched"]
    end
    if haskey(data, "Hs_isotope")
        Hs.isotope = data["Hs_isotope"]
    end
    
    # Load other parameters
    r_original = data["r_initial"]
    X₁_original = data["X1_initial"]
    δ = data["delta"]
    λ = data["lambda"]
    accuracy_newton = data["accuracy_newton"]
    accuracy_bisection = data["accuracy_bisection"]
    max_diff = data["max_diff"]
    
    # Reconstruct regularization type
    reg_type = if data["regularization_type"] == "L1Norm"
        Pioneer.L1Norm()
    elseif data["regularization_type"] == "L2Norm"
        Pioneer.L2Norm()
    else
        Pioneer.NoNorm()
    end
    
    println("\nProblem info:")
    println("  Variables: $(Hs.n)")
    println("  Constraints: $(Hs.m)")
    println("  Non-zeros: $(Hs.n_vals)")
    println("  Sparsity: $(round(100*Hs.n_vals/(Hs.n*Hs.m), digits=2))%")
    
    # Check if this was a slow-converging problem
    if haskey(data, "iterations")
        println("  Original iterations: $(data["iterations"])")
        println("  Original max_x: $(data["final_max_x"])")
    end
    
    # Run timing tests
    println("\nRunning performance tests...")
    times = Float64[]
    iterations = Int[]
    
    n_runs = 5
    for i in 1:n_runs
        r = copy(r_original)
        X₁ = copy(X₁_original)
        
        # Time the solver (with debug capture disabled)
        t_start = time()
        iters = Pioneer.solveHuber!(
            Hs, r, X₁, δ, λ, 
            100, 100, 1000,
            accuracy_newton, accuracy_bisection, 
            1e-6, max_diff, reg_type;
            debug_capture = false
        )
        t_elapsed = time() - t_start
        
        push!(times, t_elapsed)
        push!(iterations, iters)
        
        # Quick stats on the solution
        n_nonzero = sum(abs.(X₁) .> 1e-10)
        max_weight = maximum(abs.(X₁))
        
        @printf("  Run %d: %.3fs, %d iterations, %d non-zero weights, max=%.2e\n", 
                i, t_elapsed, iters, n_nonzero, max_weight)
    end
    
    # Summary statistics
    println("\nSummary:")
    @printf("  Average time: %.3f ± %.3fs\n", mean(times), std(times))
    @printf("  Average iterations: %.1f ± %.1f\n", mean(iterations), std(iterations))
    @printf("  Time per iteration: %.3fms\n", 1000*mean(times)/mean(iterations))
    
    # Test with modified parameters
    println("\nTesting parameter sensitivity...")
    
    # Test with tighter tolerance
    println("\n  With tighter max_diff ($(max_diff/10)):")
    r = copy(r_original)
    X₁ = copy(X₁_original)
    t_start = time()
    iters = Pioneer.solveHuber!(
        Hs, r, X₁, δ, λ, 
        100, 100, 1000,
        accuracy_newton, accuracy_bisection, 
        1e-6, max_diff/10, reg_type;
        debug_capture = false
    )
    t_elapsed = time() - t_start
    @printf("    Time: %.3fs, Iterations: %d\n", t_elapsed, iters)
    
    # Test with looser tolerance
    println("\n  With looser max_diff ($(max_diff*10)):")
    r = copy(r_original)
    X₁ = copy(X₁_original)
    t_start = time()
    iters = Pioneer.solveHuber!(
        Hs, r, X₁, δ, λ, 
        100, 100, 1000,
        accuracy_newton, accuracy_bisection, 
        1e-6, max_diff*10, reg_type;
        debug_capture = false
    )
    t_elapsed = time() - t_start
    @printf("    Time: %.3fs, Iterations: %d\n", t_elapsed, iters)
end

# Run the test
if abspath(PROGRAM_FILE) == @__FILE__
    test_huber_performance()
end