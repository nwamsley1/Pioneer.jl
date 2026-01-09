#!/usr/bin/env julia
# Test script to verify Polynomial coefficient ordering

using Polynomials

println("=" ^ 70)
println("Testing Julia Polynomial coefficient storage order")
println("=" ^ 70)

# Test 1: Explicit construction
println("\n### Test 1: Explicit Polynomial Construction ###")
p1 = Polynomial([1.0, 2.0, 3.0, 4.0])
println("Polynomial([1.0, 2.0, 3.0, 4.0]):")
println("  Represents: 1.0 + 2.0*x + 3.0*x² + 4.0*x³")
println("  p1.coeffs = ", p1.coeffs)
println("  Evaluated at x=2: ", p1(2.0))
println("  Manual calculation: 1.0 + 2.0*2 + 3.0*4 + 4.0*8 = ", 1.0 + 2.0*2 + 3.0*4 + 4.0*8)

# Test 2: Trailing zeros are dropped
println("\n### Test 2: Trailing Zeros Are Dropped ###")
p2_full = Polynomial([1.0, 2.0, 3.0, 0.0])
println("Polynomial([1.0, 2.0, 3.0, 0.0]):")
println("  Represents: 1.0 + 2.0*x + 3.0*x² + 0.0*x³")
println("  p2_full.coeffs = ", p2_full.coeffs)
println("  Length: ", length(p2_full.coeffs), " (Expected: 3, not 4!)")

# Test 3: Which coefficient is missing?
println("\n### Test 3: Identifying Missing Coefficient ###")
p3_expected = Polynomial([5.0, -3.0, 2.0, 0.0])  # Cubic term is zero
p3_actual = Polynomial([5.0, -3.0, 2.0])
println("Expected: 5.0 - 3.0*x + 2.0*x² + 0.0*x³")
println("  p3_expected.coeffs = ", p3_expected.coeffs)
println("Created as: Polynomial([5.0, -3.0, 2.0])")
println("  p3_actual.coeffs = ", p3_actual.coeffs)
println("  Are they equal? ", p3_expected == p3_actual)
println("  Evaluated at x=3: ", p3_actual(3.0), " vs ", 5.0 - 3.0*3 + 2.0*9)

# Test 4: Padding restoration
println("\n### Test 4: Restoring Missing Coefficients ###")
p4_short = Polynomial([1.0, 2.0, 3.0])  # Missing cubic term
println("Short polynomial: ", p4_short.coeffs)
expected_length = 4
if length(p4_short.coeffs) < expected_length
    padded_coeffs = vcat(p4_short.coeffs, zeros(expected_length - length(p4_short.coeffs)))
    println("After padding: ", padded_coeffs)
    println("Padding added at: END (highest degree terms)")
    println("This gives us: $(padded_coeffs[1]) + $(padded_coeffs[2])*x + $(padded_coeffs[3])*x² + $(padded_coeffs[4])*x³")
end

# Test 5: B-spline basis example
println("\n### Test 5: B-Spline Basis Example ###")
b0 = Polynomial([0, 0, 0, 1])/6
b1 = Polynomial([1, 3, 3, -3])/6
println("b0 = Polynomial([0, 0, 0, 1])/6")
println("  Represents: (0 + 0*x + 0*x² + 1*x³)/6 = x³/6")
println("  b0.coeffs = ", b0.coeffs)
println("\nb1 = Polynomial([1, 3, 3, -3])/6")
println("  Represents: (1 + 3*x + 3*x² - 3*x³)/6")
println("  b1.coeffs = ", b1.coeffs)

# Test 6: Linear combination that produces zero cubic term
println("\n### Test 6: Linear Combination Creating Zero Cubic Term ###")
# Create a combination where cubic terms cancel
c = [1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # Coefficients for linear combination
result = c[1] * b0 + c[2] * b1
println("Linear combination: 1.0*b0 + (-1.0)*b1")
println("  b0 has cubic: 1/6")
println("  b1 has cubic: -3/6 = -1/2")
println("  Combined cubic: 1*(1/6) + (-1)*(-1/2) = 1/6 + 1/2 = ", 1.0*(1/6) + (-1.0)*(-1/2))
println("Result coeffs: ", result.coeffs)
println("Expected length: 4, Actual length: ", length(result.coeffs))

println("\n" * "=" * 70)
println("CONCLUSION:")
println("  - Polynomial stores coefficients as [a₀, a₁, a₂, a₃] for a₀ + a₁x + a₂x² + a₃x³")
println("  - Trailing zeros (highest degree) are AUTOMATICALLY dropped")
println("  - To restore: vcat(coeffs, zeros(missing_count)) - adds at END")
println("  - This is CORRECT because missing terms are always highest degree")
println("=" * 70)
