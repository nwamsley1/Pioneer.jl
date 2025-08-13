# Mass Tolerance Estimation with Contaminated Distributions

## Problem Statement

The current mass tolerance estimation uses quantiles (1% and 99%) on the mass error distribution. However, this distribution is contaminated:

1. **True matches** (~70-90%): Follow an asymmetric Cauchy or Laplace distribution centered around the true mass error
2. **False matches** (~10-30%): Uniformly distributed across the entire search tolerance window

Using quantiles on this mixture leads to overestimation of tolerances because the uniform contamination inflates the tails.

## Mathematical Background

### Mixture Model
The observed distribution can be modeled as:
```
f(x) = π * f_true(x | μ, σ_left, σ_right) + (1-π) * f_uniform(x | -tol, +tol)
```
Where:
- π = proportion of true matches (unknown)
- f_true = asymmetric Cauchy/Laplace distribution
- f_uniform = uniform distribution across search window

### Why Quantiles Fail
- Quantiles at 1% and 99% capture the extremes
- With 20% contamination, these extremes are dominated by uniform noise
- Result: tolerance estimates reflect the search window, not the true error distribution

## Proposed Solutions

### Solution 1: Robust Scale Estimation with MAD
**Approach**: Use Median Absolute Deviation (MAD) instead of quantiles

```julia
function fit_mass_err_model_mad(params, fragments)
    ppm_errs = [calc_ppm_error(match.theoretical_mz, match.match_mz) for match in fragments]
    mass_err = median(ppm_errs)
    ppm_errs .-= mass_err
    
    # Use MAD for robust scale estimation
    mad_value = mad(ppm_errs, normalize=false)
    
    # Asymmetric MAD
    left_errors = filter(x -> x < 0, ppm_errs)
    right_errors = filter(x -> x > 0, ppm_errs)
    mad_left = isempty(left_errors) ? mad_value : mad(left_errors, normalize=false)
    mad_right = isempty(right_errors) ? mad_value : mad(right_errors, normalize=false)
    
    # Scale factor for desired coverage (e.g., 99% for Laplace)
    scale_factor = 4.0  # Configurable
    
    return MassErrorModel(
        Float32(mass_err),
        (Float32(mad_left * scale_factor), Float32(mad_right * scale_factor))
    ), ppm_errs
end
```

**Pros**:
- MAD is robust to up to 50% contamination
- Captures the scale of the true distribution
- Simple to implement

**Cons**:
- Assumes symmetric contamination
- Fixed scale factor may not be optimal for all cases

### Solution 2: Iterative Trimming (RANSAC-like)
**Approach**: Iteratively remove outliers and refit

```julia
function fit_mass_err_model_iterative(params, fragments; max_iters=5)
    ppm_errs = [calc_ppm_error(match.theoretical_mz, match.match_mz) for match in fragments]
    
    for iter in 1:max_iters
        mass_err = median(ppm_errs)
        centered_errs = ppm_errs .- mass_err
        
        # Estimate scale using central portion
        iqr = quantile(centered_errs, 0.75) - quantile(centered_errs, 0.25)
        threshold = 3.0 * iqr  # Tunable parameter
        
        # Keep only inliers
        mask = abs.(centered_errs) .< threshold
        if sum(mask) < length(ppm_errs) * 0.5
            break  # Stop if removing too many points
        end
        ppm_errs = ppm_errs[mask]
    end
    
    # Final fit on cleaned data
    mass_err = median(ppm_errs)
    ppm_errs .-= mass_err
    
    # Use higher quantiles on cleaned data
    l_bound = quantile(ppm_errs, 0.005)  # Tighter after cleaning
    r_bound = quantile(ppm_errs, 0.995)
    
    return MassErrorModel(
        Float32(mass_err),
        (Float32(abs(l_bound)), Float32(abs(r_bound)))
    ), ppm_errs
end
```

**Pros**:
- Actively removes contamination
- Can handle asymmetric distributions
- Adaptive to data characteristics

**Cons**:
- May remove legitimate outliers
- Requires tuning of threshold
- Computationally more expensive

### Solution 3: Mixture Model with EM Algorithm
**Approach**: Explicitly model the mixture and estimate parameters

```julia
function fit_mass_err_model_em(params, fragments; max_iters=50)
    ppm_errs = [calc_ppm_error(match.theoretical_mz, match.match_mz) for match in fragments]
    
    # Initialize with robust estimates
    μ = median(ppm_errs)
    σ = mad(ppm_errs, normalize=false)
    π = 0.8  # Initial guess: 80% true matches
    tol = maximum(abs.(extrema(ppm_errs)))
    
    for iter in 1:max_iters
        # E-step: Calculate responsibilities
        p_true = π * laplace_pdf.(ppm_errs, μ, σ)
        p_uniform = (1-π) / (2*tol)
        responsibilities = p_true ./ (p_true .+ p_uniform)
        
        # M-step: Update parameters
        π_new = mean(responsibilities)
        μ_new = sum(responsibilities .* ppm_errs) / sum(responsibilities)
        
        # Weighted MAD for scale
        centered = ppm_errs .- μ_new
        σ_new = sum(responsibilities .* abs.(centered)) / sum(responsibilities)
        
        # Check convergence
        if abs(π_new - π) < 0.01 && abs(μ_new - μ) < 0.1
            break
        end
        
        π, μ, σ = π_new, μ_new, σ_new
    end
    
    # Set tolerance based on true component
    scale_factor = -log(0.01)  # For 99% coverage of Laplace
    
    return MassErrorModel(
        Float32(μ),
        (Float32(σ * scale_factor), Float32(σ * scale_factor))
    ), ppm_errs .- μ
end
```

**Pros**:
- Principled statistical approach
- Estimates contamination fraction
- Can adapt to different distributions

**Cons**:
- Complex implementation
- May not converge
- Assumes specific distributional forms

### Solution 4: Hybrid Approach with Quality Scoring
**Approach**: Weight fragments by match quality before fitting

```julia
function fit_mass_err_model_weighted(params, fragments)
    # Calculate quality scores for each fragment
    scores = [calculate_match_quality(frag) for frag in fragments]
    
    # Convert to weights (higher score = higher weight)
    weights = scores ./ sum(scores)
    
    ppm_errs = [calc_ppm_error(match.theoretical_mz, match.match_mz) for match in fragments]
    
    # Weighted median for center
    mass_err = weighted_median(ppm_errs, weights)
    ppm_errs .-= mass_err
    
    # Weighted quantiles for bounds
    l_bound = weighted_quantile(ppm_errs, weights, 0.01)
    r_bound = weighted_quantile(ppm_errs, weights, 0.99)
    
    return MassErrorModel(
        Float32(mass_err),
        (Float32(abs(l_bound)), Float32(abs(r_bound)))
    ), ppm_errs
end

function calculate_match_quality(fragment)
    # Combine multiple quality metrics
    intensity_score = log10(fragment.intensity + 1)
    rank_score = 1.0 / fragment.rank
    delta_score = exp(-abs(fragment.delta_mz) / 5.0)  # Decay with mass error
    
    return intensity_score * rank_score * delta_score
end
```

**Pros**:
- Uses existing match information
- Smooth weighting (no hard cutoffs)
- Interpretable quality metrics

**Cons**:
- Quality metrics may be biased
- Still affected by contamination in weights

### Solution 5: Core Density Estimation
**Approach**: Focus on the high-density core of the distribution

```julia
function fit_mass_err_model_core(params, fragments; coverage=0.8)
    ppm_errs = [calc_ppm_error(match.theoretical_mz, match.match_mz) for match in fragments]
    
    # Find the densest region containing 'coverage' fraction of points
    mass_err = median(ppm_errs)
    centered = ppm_errs .- mass_err
    
    # Sort by absolute distance from median
    sorted_abs = sort(abs.(centered))
    
    # Find radius containing core fraction
    core_radius = sorted_abs[round(Int, length(sorted_abs) * coverage)]
    
    # Keep only core points
    core_mask = abs.(centered) .<= core_radius
    core_errors = centered[core_mask]
    
    # Fit asymmetric bounds on core
    left_core = filter(x -> x < 0, core_errors)
    right_core = filter(x -> x > 0, core_errors)
    
    # Use 95% quantiles of core (less extreme due to pre-filtering)
    l_bound = isempty(left_core) ? -core_radius : quantile(left_core, 0.05)
    r_bound = isempty(right_core) ? core_radius : quantile(right_core, 0.95)
    
    # Scale up for full coverage
    scale_factor = 1.5  # Account for trimming
    
    return MassErrorModel(
        Float32(mass_err),
        (Float32(abs(l_bound) * scale_factor), 
         Float32(abs(r_bound) * scale_factor))
    ), ppm_errs .- mass_err
end
```

**Pros**:
- Simple and robust
- Naturally handles asymmetry
- No distributional assumptions

**Cons**:
- Coverage parameter needs tuning
- May underestimate for heavy-tailed distributions

## Recommended Implementation Strategy

### Phase 1: Immediate Improvement
Replace current quantile method with **Solution 1 (MAD-based)** as it's:
- Simple to implement
- Significant improvement over quantiles
- Well-understood statistical properties

### Phase 2: Enhanced Robustness
Implement **Solution 5 (Core Density)** as primary method with:
- MAD-based fallback
- Configurable coverage parameter
- Diagnostic output for validation

### Phase 3: Advanced Options
Add **Solution 4 (Weighted)** as an option when fragment quality scores are available

## Configuration Parameters

Add to JSON configuration:
```json
"parameter_tuning": {
    "mass_error_estimation": {
        "method": "mad",  // Options: "quantile", "mad", "core", "weighted"
        "mad_scale_factor": 4.0,
        "core_coverage": 0.8,
        "iterative_max_iters": 5,
        "iterative_threshold_iqr": 3.0
    }
}
```

## Validation Strategy

1. **Synthetic Data Testing**:
   - Generate known mixture distributions
   - Validate tolerance estimates against ground truth

2. **Real Data Metrics**:
   - Track fraction of PSMs within estimated tolerance
   - Monitor convergence behavior
   - Compare methods on same datasets

3. **Diagnostic Plots**:
   - Histogram with fitted bounds
   - Q-Q plots against theoretical distributions
   - Contamination fraction estimates

## Expected Benefits

1. **More Accurate Tolerances**: 20-40% reduction in tolerance width
2. **Better FDR Control**: Fewer false positives from overly wide tolerances
3. **Improved Convergence**: More stable parameter estimation
4. **Robustness**: Handles varying contamination levels

## Implementation Priority

1. **High Priority**: MAD-based estimation (Solution 1)
2. **Medium Priority**: Core density method (Solution 5)
3. **Low Priority**: EM algorithm (Solution 3) - complex but powerful
4. **Optional**: Weighted approach (Solution 4) - when quality scores available

## Testing Approach

```julia
# Unit test for each method
function test_mass_tolerance_estimation()
    # Generate contaminated distribution
    true_errors = rand(Laplace(0, 5), 800)  # 80% true
    false_errors = rand(Uniform(-20, 20), 200)  # 20% false
    all_errors = vcat(true_errors, false_errors)
    
    # Test each method
    mad_model = fit_mass_err_model_mad(params, all_errors)
    core_model = fit_mass_err_model_core(params, all_errors)
    
    # Validate coverage
    true_coverage_mad = count(x -> within_tolerance(x, mad_model), true_errors) / length(true_errors)
    true_coverage_core = count(x -> within_tolerance(x, core_model), true_errors) / length(true_errors)
    
    @test true_coverage_mad > 0.95
    @test true_coverage_core > 0.95
end
```

## Migration Path

1. Add new methods alongside existing quantile method
2. Add method selection parameter to JSON
3. Default to MAD method for new analyses
4. Maintain quantile method for backward compatibility
5. Deprecate quantile method after validation period

## Risk Assessment

- **Low Risk**: MAD method is well-established
- **Medium Risk**: Core density requires parameter tuning
- **High Risk**: EM algorithm may not converge

## Success Metrics

1. Reduction in mass tolerance width: Target 20-40%
2. Maintenance of true positive rate: >95%
3. Reduction in false positive rate: >30%
4. Improved convergence rate: <3 iterations typical

## References

1. Robust Statistics: Huber, P. J. (1981)
2. MAD for outlier detection: Leys et al. (2013)
3. Mixture models in proteomics: Käll et al. (2007)
4. RANSAC: Fischler & Bolles (1981)