# MS1 Fitting and Target/Decoy Handling in Pioneer.jl

## Overview

Pioneer.jl uses a sophisticated sparse linear regression system for MS1-level quantification that can potentially handle target and decoy sequences with identical masses. This document analyzes how the system manages these scenarios, regularization settings, and the separation (or lack thereof) between target and decoy fitting.

## Key Findings Summary

- **L2 Regularization**: **OFF by default** (λ = 0.0, reg_type = "none")
- **Target/Decoy Separation**: **NOT separated** - targets and decoys are fit together in the same design matrix
- **Identical Masses**: Result in **separate columns** in the design matrix (based on unique precursor IDs)
- **Potential Issues**: Highly correlated or identical columns may cause numerical instability

## Design Matrix Construction

### Core Algorithm (`buildDesignMatrix.jl`)

The design matrix construction follows this pattern:

```julia
function buildDesignMatrix!(H::SparseArray{UInt32,Float32},
                            matches::Vector{m},
                            misses::Vector{m},
                            nmatches::Int64,
                            nmisses::Int64,
                            precID_to_col::ArrayDict{UInt32, UInt16})
```

**Key Process:**

1. **Column Assignment** (`buildDesignMatrix.jl:62-65`):
   ```julia
   if iszero(precID_to_col[getPrecID(match)])
       prec_col += one(UInt32)
       update!(precID_to_col, getPrecID(match), UInt16(prec_col))
   end
   ```

2. **Matrix Structure**:
   - **Rows**: Unique spectrum peaks across all scans
   - **Columns**: Unique precursor IDs (includes both targets and decoys)
   - **Values**: Predicted fragment intensities for matched peaks, zero for missed peaks

### Handling Identical Masses

**Critical Insight**: Even if two precursors (target and decoy) have identical masses, they have **different precursor IDs**. This means:

- Target precursor ID: e.g., `getPrecID(target) = 12345`
- Decoy precursor ID: e.g., `getPrecID(decoy) = 54321`
- **Result**: Two separate columns in design matrix

**Example Scenario**:
```
Target: PEPTIDEK (mass = 904.509)  → Column 1
Decoy:  PEPTIDEK_DECOY (mass = 904.509) → Column 2

Design Matrix H:
     Col1  Col2  Col3  ...
Peak1:  0.5   0.5   0.0
Peak2:  0.3   0.3   0.0
Peak3:  0.0   0.0   1.2
```

**Potential Problem**: Columns 1 and 2 are nearly identical (both explain the same peaks with similar intensities), which can cause:
- **Multicollinearity**: Poor numerical conditioning
- **Unstable Solutions**: Small changes in data lead to large coefficient changes
- **Arbitrary Coefficient Distribution**: Total signal may be split randomly between target and decoy

## Regularization Configuration

### Default Settings (`defaultSearchParams.json:123-134`)

```json
"deconvolution": {
    "lambda": 0.0,           // NO regularization penalty
    "reg_type": "none",      // Regularization type: "none", "L1", or "L2"
    "huber_delta": 300,
    "max_diff": 0.01
}
```

### Regularization Types (`spectralLinearRegression.jl:18-46`)

1. **NoNorm** (default):
   ```julia
   getRegL1(λ, xk, ::NoNorm) = zero(Float32)  # No L1 penalty
   getRegL2(λ, xk, ::NoNorm) = zero(Float32)  # No L2 penalty
   ```

2. **L2Norm**:
   ```julia
   getRegL1(λ, xk, ::L2Norm) = λ*Float32(2)*xk  # L1 derivative of L2 penalty
   getRegL2(λ, xk, ::L2Norm) = Float32(2)*λ     # L2 derivative of L2 penalty
   ```

3. **L1Norm**:
   ```julia
   getRegL1(λ, xk, ::L1Norm) = λ                # L1 penalty
   getRegL2(λ, xk, ::L1Norm) = zero(Float32)    # No L2 contribution
   ```

**Impact**: With default settings (λ=0, reg_type="none"), there is **no regularization** to handle multicollinearity from identical-mass precursors.

## Huber Solver Implementation

### Main Solver (`spectralLinearRegression.jl:238-290`)

```julia
function solveHuber!(Hs::SparseArray{Ti, T},
                     r::Vector{T},
                     X₁::Vector{T},
                     δ::T,           # Huber delta parameter
                     λ::T,           # Regularization lambda (default: 0.0)
                     ...
                     regularization_type::RegularizationType)  # default: NoNorm
```

**Algorithm**:
1. **Outer Loop**: Coordinate descent over all precursors
2. **Inner Loop**: Newton-Raphson + Bisection for each coefficient
3. **Convergence**: Based on relative change in coefficients (`max_diff = 0.01`)

### Key Properties

1. **Robust Loss Function**: Huber loss is less sensitive to outliers than squared loss
2. **Non-negative Constraints**: Coefficients constrained to `≥ 0` (`buildDesignMatrix.jl:148`)
3. **No Multicollinearity Handling**: Without regularization, solver may struggle with highly correlated columns

### Newton-Raphson with Huber Loss (`spectralLinearRegression.jl:109-150`)

The solver uses a sophisticated optimization that includes:

1. **Fast Inverse Square Root** (`spectralLinearRegression.jl:71-76`):
   ```julia
   # Quake's algorithm for computing 1/sqrt(1 + (r/δ)²)
   R = RS
   int32 = reinterpret(UInt32, R)
   int32 = 0x5f3759df - int32 >> 1
   R = reinterpret(Float32, int32)
   R *= 1.5f0 - RS * 0.5f0 * R^2
   ```

2. **Derivatives with Regularization** (`spectralLinearRegression.jl:83`):
   ```julia
   return Float32(L1) + getRegL1(λ, xk, regularization_type),
          Float32(L2) + getRegL2(λ, xk, regularization_type)
   ```

## Target/Decoy Separation Analysis

### Evidence for Joint Fitting

1. **Single Design Matrix**: Both targets and decoys are included in the same `buildDesignMatrix!` call
2. **Single Solver Call**: `solveHuber!` operates on the complete matrix containing both targets and decoys
3. **No Filtering**: No evidence of target/decoy separation in the quantification pipeline

### Implications

**Advantages**:
- **Unified Framework**: Consistent treatment of all precursors
- **Shared Peak Explanation**: Multiple precursors can contribute to explaining the same peak

**Disadvantages**:
- **Multicollinearity**: Identical-mass precursors create nearly-identical columns
- **Unstable Quantification**: Coefficients may vary arbitrarily between target and decoy
- **Difficult Interpretation**: Hard to determine "true" abundance when signal is split

## Potential Issues and Solutions

### Problem 1: Multicollinearity from Identical Masses

**Issue**: Target and decoy with same mass create highly correlated columns

**Current Behavior**: Without regularization, the solver may:
- Converge to arbitrary solutions
- Be sensitive to numerical precision
- Split signal unpredictably between target and decoy

**Potential Solutions**:

1. **Enable L2 Regularization**:
   ```json
   "deconvolution": {
       "lambda": 0.01,      // Small L2 penalty
       "reg_type": "L2"
   }
   ```

2. **Precursor Grouping**: Combine target/decoy with identical masses into single column

3. **Regularization with Grouping**: Use group lasso to select one representative per mass

### Problem 2: No Mass-Based Deduplication

**Issue**: Multiple precursors with identical masses all get separate columns

**Evidence**: Column assignment based on `getPrecID(match)`, not mass

**Potential Solutions**:

1. **Mass-Based Grouping**:
   ```julia
   # Group precursors by mass before matrix construction
   mass_to_precursors = group_by_mass(precursors)
   for (mass, prec_group) in mass_to_precursors
       # Assign single column per mass group
   end
   ```

2. **Post-Processing**: Combine coefficients for identical-mass precursors after solving

## Code References

### Key Files and Functions

1. **Design Matrix Construction**:
   - `src/Routines/SearchDIA/CommonSearchUtils/buildDesignMatrix.jl:18-104`
   - Function: `buildDesignMatrix!()`

2. **Huber Solver**:
   - `src/utils/ML/spectralLinearRegression.jl:238-290`
   - Function: `solveHuber!()`

3. **Regularization Implementation**:
   - `src/utils/ML/spectralLinearRegression.jl:18-46`
   - Types: `NoNorm`, `L1Norm`, `L2Norm`

4. **Configuration**:
   - `assets/example_config/defaultSearchParams.json:123-134`
   - Section: `optimization.deconvolution`

### Search Method Integration

The MS1 fitting is used across multiple search methods:

1. **QuadTuningSearch**: `src/Routines/SearchDIA/SearchMethods/QuadTuningSearch/utils.jl`
2. **SecondPassSearch**: `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl`
3. **IntegrateChromatogramsSearch**: `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl`
4. **HuberTuningSearch**: `src/Routines/SearchDIA/SearchMethods/HuberTuningSearch/utils.jl`

## Recommendations

### For Current Implementation

1. **Monitor for Multicollinearity**: Check condition number of design matrices
2. **Consider L2 Regularization**: Small λ (e.g., 0.01) may improve stability
3. **Post-Processing**: Sum coefficients for identical-mass precursors in downstream analysis

### For Future Development

1. **Mass-Based Grouping**: Group precursors by mass during matrix construction
2. **Advanced Regularization**: Implement group lasso for mass-based selection
3. **Separate Target/Decoy Fitting**: Consider separate design matrices if biological interpretation requires it

## Conclusion

Pioneer.jl's MS1 fitting system treats targets and decoys jointly without regularization by default. While this provides a unified framework, it may lead to numerical instability when target and decoy precursors have identical masses. The current implementation assigns separate columns to each unique precursor ID, regardless of mass, which can result in multicollinearity. Enabling L2 regularization or implementing mass-based grouping could improve the robustness of quantification in these scenarios.