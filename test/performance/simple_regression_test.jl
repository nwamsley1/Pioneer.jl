using Pioneer

# Test if the zero-weight issue is fixed
function test_zero_weight_convergence()
    println("Testing newton_bisection with zero initial weight...")
    
    # Create minimal sparse array
    sa = Pioneer.SparseArray(1000)
    
    # Add some entries manually
    sa.n_vals = 100
    sa.m = 50
    sa.n = 10
    
    # Fill with test data
    for i in 1:100
        sa.rowval[i] = rand(1:50)
        sa.colval[i] = UInt16(rand(1:10))
        sa.nzval[i] = rand(Float32)
        sa.x[i] = rand(Float32)
    end
    
    # Sort and build column pointers
    Pioneer.sortSparse!(sa)
    
    # Test newton_bisection with zero initial weight
    r = rand(Float32, 50)
    X₁ = zeros(Float32, 10)  # All zeros!
    
    # Call newton_bisection for first column
    col = 1
    δ = 100.0f0
    λ = 0.001f0
    
    println("Initial X₁[1] = ", X₁[1])
    
    # This should now work even with X₁[1] = 0
    change = Pioneer.newton_bisection!(
        sa, r, X₁, col, δ, λ,
        50,      # max_iter_newton
        100,     # max_iter_bisection  
        0.01f0,  # accuracy_newton
        0.01f0,  # accuracy_bisection
        Pioneer.NoNorm(),
        0.1f0    # rel_tol
    )
    
    println("Change after newton_bisection: ", change)
    println("Final X₁[1] = ", X₁[1])
    println("Success: ", X₁[1] != 0.0 ? "✓ Weight changed from zero" : "✗ Weight still zero")
end

# Run the test
test_zero_weight_convergence()