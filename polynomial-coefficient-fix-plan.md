# Fix Plan: Polynomial Coefficient Padding Bug

## Overview
Fix the dimension mismatch error in `uniformBasisCubicSpline.jl` by ensuring all polynomials contribute exactly `degree + 1` coefficients, regardless of trailing zeros.

## Files to Modify

### 1. src/utils/ML/uniformBasisCubicSpline.jl

#### Change 1: Fix UniformSplinePenalized (Line 665-684)

**Current Code:**
```julia
    # Build piecewise polynomials
    XPoly = buildPieceWise(knots, bin_width, spline_basis)
    piecewise_polynomials = XPoly * c
    n_total_coeffs = n_knots * (degree + 1)

    # DEBUG: Check polynomial coefficient counts
    poly_coeff_lengths = [length(poly.coeffs) for poly in piecewise_polynomials]
    expected_length = degree + 1
    problematic_polys = findall(x -> x != expected_length, poly_coeff_lengths)
    if !isempty(problematic_polys)
        @warn "UniformSplinePenalized: Some polynomials have incorrect coefficient counts!"
        @warn "  n_knots: $n_knots, degree: $degree"
        @warn "  Expected coeffs per poly: $expected_length"
        @warn "  Actual coefficient counts: $poly_coeff_lengths"
        @warn "  Problematic polynomial indices: $problematic_polys"
        @warn "  Total coeffs expected: $n_total_coeffs"
        @warn "  Total coeffs actual: $(sum(poly_coeff_lengths))"
    end

    coeffs = SVector{n_total_coeffs}(vcat([poly.coeffs for poly in piecewise_polynomials]...))
```

**New Code:**
```julia
    # Build piecewise polynomials
    XPoly = buildPieceWise(knots, bin_width, spline_basis)
    piecewise_polynomials = XPoly * c
    n_total_coeffs = n_knots * (degree + 1)

    # Extract coefficients, padding to ensure exactly (degree+1) per polynomial
    # The Polynomial type drops trailing zeros, so we must pad them back
    expected_length = degree + 1
    coeff_vecs = [
        let c = poly.coeffs
            if length(c) < expected_length
                # Pad with zeros for missing trailing coefficients
                vcat(c, zeros(T, expected_length - length(c)))
            else
                # Take exactly the expected number (shouldn't exceed, but be defensive)
                c[1:expected_length]
            end
        end
        for poly in piecewise_polynomials
    ]
    coeffs = SVector{n_total_coeffs}(vcat(coeff_vecs...))
```

**Lines to replace:** 665-684
**Explanation:** Removes debug code and replaces direct coefficient extraction with padding logic.

---

#### Change 2: Fix UniformSpline (Line 211-242)

**Current Code:**
```julia
    # insure the input is sorted so the spline fitting is stable under
    # different numbers of threads
    if issorted(t)
        sorted_t = t
        sorted_u = u
    else
        perm = sortperm(t)
        sorted_t = t[perm]
        sorted_u = u[perm]
    end

    spline_basis = getSplineBasis(degree)
    _first = minimum(sorted_t)
    _last = maximum(sorted_t)
    bin_width = (_last - _first)/(n_knots - 1)
    knots = collect(LinRange(_first, _last, n_knots))
    X = buildDesignMat(sorted_t, collect(knots), bin_width, spline_basis)
    c = X\sorted_u
    XPoly = buildPieceWise(knots, bin_width, spline_basis)
    piecewise_polynomials = XPoly*c
    n_coeffs = n_knots*(degree + 1)
    coeffs = SVector{n_coeffs}(vcat([polynomial.coeffs for polynomial in piecewise_polynomials]...))


    UniformSpline{n_coeffs, T}(
        coeffs,
        degree,
        _first,
        _last,
        bin_width
    )
```

**New Code:**
```julia
    # insure the input is sorted so the spline fitting is stable under
    # different numbers of threads
    if issorted(t)
        sorted_t = t
        sorted_u = u
    else
        perm = sortperm(t)
        sorted_t = t[perm]
        sorted_u = u[perm]
    end

    spline_basis = getSplineBasis(degree)
    _first = minimum(sorted_t)
    _last = maximum(sorted_t)
    bin_width = (_last - _first)/(n_knots - 1)
    knots = collect(LinRange(_first, _last, n_knots))
    X = buildDesignMat(sorted_t, collect(knots), bin_width, spline_basis)
    c = X\sorted_u
    XPoly = buildPieceWise(knots, bin_width, spline_basis)
    piecewise_polynomials = XPoly*c
    n_coeffs = n_knots*(degree + 1)

    # Extract coefficients, padding to ensure exactly (degree+1) per polynomial
    # The Polynomial type drops trailing zeros, so we must pad them back
    expected_length = degree + 1
    coeff_vecs = [
        let c = polynomial.coeffs
            if length(c) < expected_length
                # Pad with zeros for missing trailing coefficients
                vcat(c, zeros(T, expected_length - length(c)))
            else
                # Take exactly the expected number (shouldn't exceed, but be defensive)
                c[1:expected_length]
            end
        end
        for polynomial in piecewise_polynomials
    ]
    coeffs = SVector{n_coeffs}(vcat(coeff_vecs...))

    UniformSpline{n_coeffs, T}(
        coeffs,
        degree,
        _first,
        _last,
        bin_width
    )
```

**Lines to replace:** 211-242
**Explanation:** Applies same padding logic to prevent future failures in standard UniformSpline.

---

### 2. src/Routines/SearchDIA/CommonSearchUtils/rt_alignment_utils.jl

#### Change 3: Remove Excessive Debug Logging (Lines 167-186, 203-212)

**Current Code (Lines 167-186):**
```julia
    # DEBUG: Log input parameters
    @user_info "=== make_spline_monotonic DEBUG ==="
    @user_info "  Input array lengths: rt_data=$(length(rt_data)), irt_data=$(length(irt_data))"
    @user_info "  n_knots parameter: $n_knots"

    # 1. Sample from original spline at uniform grid
    rt_min, rt_max = extrema(rt_data)
    n_sample_points = length(rt_data) - 1
    @user_info "  Computed n_sample_points: $n_sample_points (length(rt_data) - 1)"
    @user_info "  RT range: [$rt_min, $rt_max]"

    rt_grid = collect(LinRange(Float32(rt_min), Float32(rt_max), n_sample_points))
    irt_grid = [original_spline(r) for r in rt_grid]

    @user_info "  Grid lengths: rt_grid=$(length(rt_grid)), irt_grid=$(length(irt_grid))"

    # 2. Find median split point
    median_rt = median(rt_data)
    median_idx = argmin(abs.(rt_grid .- median_rt))
    @user_info "  Median RT: $median_rt, median_idx: $median_idx"
```

**New Code:**
```julia
    # 1. Sample from original spline at uniform grid
    rt_min, rt_max = extrema(rt_data)
    n_sample_points = length(rt_data) - 1
    rt_grid = collect(LinRange(Float32(rt_min), Float32(rt_max), n_sample_points))
    irt_grid = [original_spline(r) for r in rt_grid]

    # 2. Find median split point
    median_rt = median(rt_data)
    median_idx = argmin(abs.(rt_grid .- median_rt))
```

**Lines to remove:** 167-170, 175-176, 181, 186

---

**Current Code (Lines 203-212):**
```julia
    # 5. Refit spline to filtered data with penalty for smoothness
    @user_info "  After filtering - Grid lengths: rt_grid=$(length(rt_grid)), irt_grid=$(length(irt_grid))"
    @user_info "  Calling UniformSplinePenalized with:"
    @user_info "    irt_grid length: $(length(irt_grid))"
    @user_info "    rt_grid length: $(length(rt_grid))"
    @user_info "    degree: 3"
    @user_info "    n_knots: $n_knots"
    @user_info "    lambda: 1.0"
    @user_info "    penalty_order: 2"
    @user_info "  irt_grid range: [$(minimum(irt_grid)), $(maximum(irt_grid))]"
    @user_info "  rt_grid range: [$(minimum(rt_grid)), $(maximum(rt_grid))]"
```

**New Code:**
```julia
    # 5. Refit spline to filtered data with penalty for smoothness
```

**Lines to remove:** 203-212

---

**Current Code (Lines 413-418):**
```julia
        # Apply monotonic enforcement to prevent backwards slopes at edges
        @user_info "Applying monotonic enforcement with bidirectional cumulative max filter"
        @user_info "DEBUG: Before calling make_spline_monotonic:"
        @user_info "  valid_psms rows: $(nrow(valid_psms))"
        @user_info "  valid_psms[!, :rt] length: $(length(valid_psms[!, :rt]))"
        @user_info "  valid_psms[!, :irt_predicted] length: $(length(valid_psms[!, :irt_predicted]))"
        @user_info "  n_knots: $n_knots"
```

**New Code:**
```julia
        # Apply monotonic enforcement to prevent backwards slopes at edges
        @user_info "Applying monotonic enforcement with bidirectional cumulative max filter"
```

**Lines to remove:** 414-418

---

#### Change 4: Keep Stack Trace (Lines 420-425) - NO CHANGE

Keep the existing stack trace logging in the catch block - this is useful for debugging:

```julia
    catch e
        # Unexpected failure during spline fitting - fall back to linear model
        @user_warn "RT spline fitting failed unexpectedly ($e), falling back to linear model \n"
        @user_warn "Full stack trace:"
        for (exc, bt) in Base.catch_stack()
            showerror(stderr, exc, bt)
            println(stderr)
        end
```

**Lines:** 420-425 - KEEP AS IS

---

## Summary of Changes

| File | Location | Change Type | Reason |
|------|----------|-------------|--------|
| `uniformBasisCubicSpline.jl` | Lines 665-684 | Fix + cleanup | Fix UniformSplinePenalized, remove debug code |
| `uniformBasisCubicSpline.jl` | Lines 228-232 | Fix | Fix UniformSpline (same bug, not yet triggered) |
| `rt_alignment_utils.jl` | Lines 167-186 | Cleanup | Remove excessive debug logging |
| `rt_alignment_utils.jl` | Lines 203-212 | Cleanup | Remove excessive debug logging |
| `rt_alignment_utils.jl` | Lines 413-418 | Cleanup | Remove excessive debug logging |
| `rt_alignment_utils.jl` | Lines 420-425 | Keep | Retain useful error stack trace |
| `test/runtests.jl` | Line 38 | Add imports | Import UniformSplinePenalized + fit_spline_with_ransac |
| `test/UnitTests/uniformBassisCubicSpline.jl` | Entire file | Comprehensive tests | Add zero coefficient tests for both spline types |

## Unit Test Updates

### 3. test/runtests.jl

#### Change 5: Add UniformSplinePenalized Import (Line 38)

**Current Code:**
```julia
using Pioneer: UniformSpline  # For uniformBasisCubicSpline.jl test
```

**New Code:**
```julia
using Pioneer: UniformSpline, UniformSplinePenalized, fit_spline_with_ransac  # For uniformBasisCubicSpline.jl test
```

**Lines to replace:** 38
**Explanation:** Import additional spline functions for comprehensive testing.

---

### 4. test/UnitTests/uniformBassisCubicSpline.jl

#### Change 6: Replace Entire Test File

**Current Code:**
```julia
#Example data. Approximate a sine curve
N = 200
t = collect(LinRange(0.0, 4*π, N))
u = sin.(t)
#u .+= randn(N)./50
plot(t, u, seriestype=:scatter)
test_spline = UniformSpline(u, t, 3, 20)
#plot!(LinRange(-1, 4*π+1, 500), test_spline.(LinRange(0-1, 4*π+1, 500)), linewidth = 3, alpha = 0.5)
#using DataInterpolations
#test_old = BSplineApprox(u, t, 4, 20, :Uniform, :Uniform, extrapolate = true)
#plot!(LinRange(-1, 4*π+1, 500), test_old.(LinRange(0-1, 4*π+1, 500)), linewidth = 3, alpha = 0.5)
UniformSpline(u, t, 3, 3)

@test maximum(test_spline.(t) .- u) .< 1e-3
```

**New Code:**
```julia
@testset "UniformSpline Tests" begin
    @testset "Basic sine wave approximation" begin
        # Example data: Approximate a sine curve
        N = 200
        t = collect(LinRange(0.0, 4*π, N))
        u = sin.(t)

        test_spline = UniformSpline(u, t, 3, 20)
        @test maximum(abs.(test_spline.(t) .- u)) < 1e-3

        # Test with minimum knots
        test_spline_min = UniformSpline(u, t, 3, 3)
        @test test_spline_min isa UniformSpline
    end

    @testset "Zero coefficient handling - UniformSpline" begin
        # Test case 1: Nearly constant data (should produce zero higher-order coefficients)
        N = 100
        t = collect(LinRange(0.0f0, 10.0f0, N))
        u = 5.0f0 .+ 0.0001f0 .* randn(Float32, N)  # Nearly constant

        # Use many knots to force edge polynomials to have zero coefficients
        spline = UniformSpline(u, t, 3, 25)
        @test spline isa UniformSpline

        # Verify spline evaluates correctly at all input points
        predictions = [spline(ti) for ti in t]
        @test all(isfinite.(predictions))
        @test length(predictions) == N
    end

    @testset "Zero coefficient handling - Perfectly linear data" begin
        # Test case 2: Perfectly linear data (zero curvature)
        N = 100
        t = collect(LinRange(0.0f0, 10.0f0, N))
        u = 2.0f0 .+ 3.0f0 .* t  # y = 2 + 3x (perfectly linear)

        # Many knots with linear data should produce polynomials with zero cubic terms
        spline = UniformSpline(u, t, 3, 30)
        @test spline isa UniformSpline

        # Verify linear relationship is preserved
        predictions = [spline(ti) for ti in t]
        @test maximum(abs.(predictions .- u)) < 1e-4
    end

    @testset "Zero coefficient handling - Boundary polynomials" begin
        # Test case 3: Data that goes flat at boundaries
        N = 150
        t = collect(LinRange(0.0f0, 10.0f0, N))
        # Flat at edges, curved in middle
        u = [ti < 1.0f0 ? 1.0f0 : (ti > 9.0f0 ? 3.0f0 : 1.0f0 + sin(ti)) for ti in t]

        # High knot count to test boundary handling
        spline = UniformSpline(u, t, 3, 40)
        @test spline isa UniformSpline

        predictions = [spline(ti) for ti in t]
        @test all(isfinite.(predictions))
    end
end

@testset "UniformSplinePenalized Tests" begin
    @testset "Basic sine wave with penalty" begin
        N = 200
        t = collect(LinRange(0.0, 4*π, N))
        u = sin.(t)

        # Test with low penalty
        spline_low = UniformSplinePenalized(u, t, 3, 20, 0.01, 2)
        @test maximum(abs.(spline_low.(t) .- u)) < 1e-2

        # Test with medium penalty
        spline_med = UniformSplinePenalized(u, t, 3, 20, 0.1, 2)
        @test maximum(abs.(spline_med.(t) .- u)) < 0.1

        # Test with high penalty (should be smoother, less accurate)
        spline_high = UniformSplinePenalized(u, t, 3, 20, 1.0, 2)
        @test maximum(abs.(spline_high.(t) .- u)) < 0.5
    end

    @testset "Zero penalty equals standard spline" begin
        N = 100
        t = collect(LinRange(0.0, 10.0, N))
        u = sin.(t) .+ cos.(t ./ 2)

        spline_standard = UniformSpline(u, t, 3, 10)
        spline_no_penalty = UniformSplinePenalized(u, t, 3, 10, 0.0, 2)

        # Should give nearly identical results
        test_points = LinRange(0.0, 10.0, 50)
        for tp in test_points
            @test abs(spline_standard(tp) - spline_no_penalty(tp)) < 1e-6
        end
    end

    @testset "Zero coefficient handling - Nearly constant data with high penalty" begin
        # This is the case that triggers the bug!
        N = 100
        t = collect(LinRange(0.0f0, 10.0f0, N))
        u = 5.0f0 .+ 0.0001f0 .* randn(Float32, N)  # Nearly constant

        # High penalty with many knots should produce zero coefficients at boundaries
        spline = UniformSplinePenalized(u, t, 3, 25, 1.0f0, 2)
        @test spline isa UniformSpline

        # Verify evaluation works
        predictions = [spline(ti) for ti in t]
        @test all(isfinite.(predictions))
        @test length(predictions) == N
    end

    @testset "Zero coefficient handling - Monotonic filtered data" begin
        # Simulate the monotonic enforcement scenario from make_spline_monotonic
        N = 5821  # From actual failure case
        t = collect(LinRange(5.2123957f0, 112.212006f0, N))

        # Create smooth data with some noise
        u_base = LinRange(-0.05f0, 27.19f0, N)
        u = collect(u_base) .+ 0.1f0 .* randn(Float32, N)

        # Apply monotonic filtering (cummax-like smoothing)
        for i in 2:N
            if u[i] < u[i-1]
                u[i] = u[i-1]
            end
        end

        # This combination (smooth data + high penalty + many knots) triggers the bug
        spline = UniformSplinePenalized(u, t, 3, 59, 1.0f0, 2)
        @test spline isa UniformSpline

        # Verify evaluation at boundaries (where zero coefficients are likely)
        @test isfinite(spline(t[1]))
        @test isfinite(spline(t[end]))
        @test isfinite(spline(median(t)))
    end

    @testset "Zero coefficient handling - Different penalty orders" begin
        N = 100
        t = collect(LinRange(0.0f0, 10.0f0, N))
        u = 5.0f0 .+ 0.01f0 .* randn(Float32, N)

        # Test all supported penalty orders
        for order in [1, 2, 3]
            spline = UniformSplinePenalized(u, t, 3, 20, 0.5f0, order)
            @test spline isa UniformSpline
            @test all(isfinite.([spline(ti) for ti in t]))
        end
    end
end

@testset "fit_spline_with_ransac Tests" begin
    @testset "RANSAC with outliers" begin
        N = 200
        t = collect(LinRange(0.0, 10.0, N))
        u = sin.(t)

        # Add outliers
        u_noisy = copy(u)
        outlier_indices = rand(1:N, 20)
        u_noisy[outlier_indices] .+= randn(20) .* 2.0

        # RANSAC should handle outliers better than standard fitting
        spline_ransac = fit_spline_with_ransac(u_noisy, t, 3, 10, 0.1, 2)
        @test spline_ransac isa UniformSpline

        # Predictions should be reasonable
        predictions = [spline_ransac(ti) for ti in t]
        @test all(isfinite.(predictions))
    end

    @testset "RANSAC zero coefficient handling" begin
        # Test RANSAC with data that produces zero coefficients
        N = 100
        t = collect(LinRange(0.0f0, 10.0f0, N))
        u = 5.0f0 .+ 0.1f0 .* randn(Float32, N)

        spline = fit_spline_with_ransac(u, t, 3, 25, 1.0f0, 2,
                                        ransac_iterations=20,
                                        ransac_sample_size=30)
        @test spline isa UniformSpline
        @test all(isfinite.([spline(ti) for ti in t]))
    end
end
```

**Explanation:**
- Comprehensive test coverage for both UniformSpline and UniformSplinePenalized
- Specific tests for zero coefficient scenarios that trigger the bug
- Tests edge cases: constant data, linear data, boundary effects, high penalties
- Tests the actual failure case (N=5821, n_knots=59, λ=1.0)
- Verifies RANSAC handling as well

---

## Testing Plan

After applying fixes:

1. **Run unit tests:**
   ```bash
   julia --project=. test/runtests.jl
   ```

2. **Expected unit test results:**
   - All UniformSpline tests pass (including zero coefficient cases)
   - All UniformSplinePenalized tests pass (including high penalty + many knots)
   - All fit_spline_with_ransac tests pass
   - No DimensionMismatch errors in any test case

3. **Run integration test:**
   ```julia
   SearchDIA("/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/KYSE70_TEST/KYSE70_TEST_A.json")
   ```

4. **Expected integration result:**
   - No `DimensionMismatch` errors
   - No debug warnings about polynomial coefficients
   - RT spline fitting succeeds
   - All 4 files process successfully

5. **Verify output:**
   - Check that monotonic RT splines are being used (not falling back to linear)
   - Confirm cleaner log output without excessive debug messages

## Rollback Plan

If issues arise, the changes are isolated and can be reverted:

1. Revert coefficient padding logic → original direct extraction
2. Restore debug logging if needed for further investigation
3. All changes are in two files only, easy to git revert

## Code Review Checklist

- [ ] Both UniformSpline functions use identical padding logic
- [ ] Padding preserves original coefficient values (only adds zeros)
- [ ] Type safety maintained (zeros created with correct type T)
- [ ] Debug logging removed but error handling retained
- [ ] No changes to algorithm logic, only coefficient extraction
- [ ] Works for all polynomial lengths (< expected, = expected, > expected)
- [ ] Unit tests added for UniformSplinePenalized
- [ ] Unit tests cover zero coefficient edge cases
- [ ] Unit tests replicate actual failure scenario (N=5821, n_knots=59)
- [ ] RANSAC tests include zero coefficient handling
