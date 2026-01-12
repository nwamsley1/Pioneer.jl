# RT-Dependent Tolerance Implementation Plan

## Current System Analysis

### How RT Tolerance Currently Works

**1. Calculation (FirstPassSearch/utils.jl:664-706)**
```julia
function get_irt_errs(fwhms, prec_to_irt, params)
    # Peak width contribution: median_fwhm + n_std * MAD(fwhm)
    fwhms = map(x -> x[:median_fwhm] + params.fwhm_nstd * x[:mad_fwhm], fwhms)

    # Cross-run variance: sqrt(var_irt / (n-1)) for precursors in >2 runs
    variance_ = collect(skipmissing(
        map(x -> (x[:n] > 2) ? sqrt(x[:var_irt]/(x[:n] - 1)) : missing, prec_to_irt)
    ))
    irt_std = isempty(variance_) ? 0.0f0 : median(variance_) * params.irt_nstd

    # Final: peak_fwhm + cross_run_std (single value per file)
    return map(x -> Float32(x + irt_std), fwhms)
end
```

**2. Storage (SearchTypes.jl:226)**
```julia
irt_errors::Dict{Int64, Float32}  # Single value per file
```

**3. Usage Pattern (11 locations across 5 search methods)**
```julia
irt_tol = getIrtErrors(search_context)[ms_file_idx]  # Static per file
irt = getRtIrtModel(search_context, ms_file_idx)(observed_rt)
window = [irt - irt_tol, irt + irt_tol]  # Symmetric, constant width
```

### The Problem

- **Single tolerance value per file** - doesn't account for:
  - Peak width variation across gradient (often wider at start/end)
  - RT prediction accuracy varying with iRT (less calibration data at extremes)
  - Chromatographic effects (band broadening at high RT)

- **Real-world observation**: LC gradients typically show:
  - Wider peaks at gradient start (loading effects)
  - Narrower peaks in middle
  - Broader peaks at end (column equilibration)

---

## Proposed Solution: iRT-Dependent Tolerance Spline

### Approach Overview

1. **Bin PSMs by iRT** during FirstPassSearch tolerance calculation
2. **Calculate local tolerance** for each iRT bin
3. **Fit a smoothing spline** to tolerance vs iRT
4. **Store as callable model** (similar to existing RT conversion models)
5. **Evaluate at query time** to get iRT-specific tolerance

### Design Decision: Spline-Based Model

**Why spline instead of lookup table?**
- Smooth transitions between regions
- Handles sparse data gracefully
- Consistent with existing RT model infrastructure
- Natural extrapolation at edges

---

## Implementation Details

### Step 1: Create New Model Type

**File: `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl`**

```julia
# New struct for iRT-dependent tolerance
struct IrtDependentTolerance
    spline::UniformSpline    # Tolerance as function of iRT
    min_tol::Float32         # Floor tolerance (prevent too narrow)
    max_tol::Float32         # Ceiling tolerance (prevent too wide)
end

# Make it callable
function (model::IrtDependentTolerance)(irt::Float32)
    raw_tol = model.spline(irt)
    return clamp(raw_tol, model.min_tol, model.max_tol)
end

# Fallback for backward compatibility
struct ConstantTolerance
    value::Float32
end
(model::ConstantTolerance)(irt::Float32) = model.value

# Union type for storage
const IrtTolerance = Union{IrtDependentTolerance, ConstantTolerance, Float32}
```

**Update storage type:**
```julia
# Change from:
irt_errors::Dict{Int64, Float32}
# To:
irt_errors::Dict{Int64, IrtTolerance}
```

### Step 2: Modify Tolerance Calculation

**File: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`**

```julia
function get_irt_errs_dependent(
    psms::DataFrame,          # PSMs with :irt_predicted, :fwhm columns
    n_bins::Int = 10,         # Number of iRT bins
    min_psms_per_bin::Int = 50,
    params::FirstPassSearchParameters
)::Dictionary{Int64, IrtTolerance}

    result = Dictionary{Int64, IrtTolerance}()

    for file_idx in unique(psms.ms_file_idx)
        file_psms = filter(row -> row.ms_file_idx == file_idx, psms)

        if nrow(file_psms) < min_psms_per_bin * 3
            # Too few PSMs - fall back to constant tolerance
            result[file_idx] = ConstantTolerance(calculate_constant_tolerance(file_psms, params))
            continue
        end

        # Bin by iRT
        irt_range = extrema(file_psms.irt_predicted)
        bin_edges = range(irt_range[1], irt_range[2], length=n_bins+1)

        bin_tolerances = Float32[]
        bin_centers = Float32[]

        for i in 1:n_bins
            bin_mask = (file_psms.irt_predicted .>= bin_edges[i]) .&
                       (file_psms.irt_predicted .< bin_edges[i+1])
            bin_psms = file_psms[bin_mask, :]

            if nrow(bin_psms) >= min_psms_per_bin
                # Calculate local tolerance: median FWHM + MAD * scale
                local_fwhm = median(bin_psms.fwhm)
                local_mad = mad(bin_psms.fwhm, normalize=true)
                local_tol = local_fwhm + params.fwhm_nstd * local_mad

                push!(bin_tolerances, local_tol)
                push!(bin_centers, (bin_edges[i] + bin_edges[i+1]) / 2)
            end
        end

        if length(bin_tolerances) >= 3
            # Fit smoothing spline
            tol_spline = UniformSpline(
                Float32.(bin_tolerances),
                Float32.(bin_centers),
                3,  # cubic
                min(length(bin_tolerances) - 1, 5)  # limited knots for smoothness
            )

            result[file_idx] = IrtDependentTolerance(
                tol_spline,
                minimum(bin_tolerances) * 0.5f0,  # min floor
                maximum(bin_tolerances) * 2.0f0   # max ceiling
            )
        else
            # Fall back to constant
            result[file_idx] = ConstantTolerance(median(bin_tolerances))
        end
    end

    return result
end
```

### Step 3: Update Usage Locations (11 places)

**Pattern change:**

```julia
# OLD (constant tolerance):
irt_tol = getIrtErrors(search_context)[ms_file_idx]

# NEW (evaluate at specific iRT):
irt_tol_model = getIrtErrors(search_context)[ms_file_idx]
irt_tol = irt_tol_model(irt)  # Evaluate spline at current iRT
```

**Files to update:**
1. `SecondPassSearch/utils.jl` - Lines 169, 374
2. `HuberTuningSearch/utils.jl` - Line 218
3. `IntegrateChromatogramsSearch/utils.jl` - Lines 244, 478
4. `LibrarySearch.jl` - Lines 242, 263
5. `ParameterTuningSearch.jl` - Lines 166, 855

### Step 4: Backward Compatibility

```julia
# Getter that handles both old Float32 and new model types
function get_irt_tolerance(search_context, file_idx, irt::Float32)
    tol = getIrtErrors(search_context)[file_idx]
    if tol isa Float32
        return tol
    else
        return tol(irt)
    end
end
```

---

## Alternative Approaches Considered

### Option A: Percentile-Based Scaling
```julia
# Scale tolerance by RT percentile
irt_percentile = (irt - irt_min) / (irt_max - irt_min)
scale_factor = 1.0 + 0.5 * abs(irt_percentile - 0.5)  # Higher at edges
adjusted_tol = base_tol * scale_factor
```
**Pros:** Simple, no spline fitting
**Cons:** Assumes U-shaped pattern, not data-driven

### Option B: Lookup Table with Interpolation
```julia
# Store tolerance array, interpolate at query time
struct TolLookup
    irt_points::Vector{Float32}
    tol_values::Vector{Float32}
end
```
**Pros:** Fast lookup
**Cons:** Discrete jumps, more storage

### Option C: Piecewise Linear (Selected for simplicity alternative)
```julia
# Three-region model: early, middle, late
struct PiecewiseTolerance
    early_tol::Float32
    mid_tol::Float32
    late_tol::Float32
    breakpoint_1::Float32
    breakpoint_2::Float32
end
```
**Pros:** Very simple, interpretable
**Cons:** Sharp transitions, may not fit all gradients

---

## Recommended Implementation Order

1. **Phase 1: Data Collection** (minimal code change)
   - Add FWHM binning to FirstPassSearch diagnostic output
   - Log tolerance vs iRT relationship to verify assumption

2. **Phase 2: Model Implementation**
   - Add `IrtDependentTolerance` type
   - Implement spline fitting in `get_irt_errs_dependent`
   - Update storage type

3. **Phase 3: Integration**
   - Update all 11 usage locations
   - Add backward compatibility layer
   - Test with existing config files

4. **Phase 4: Configuration**
   - Add parameter to enable/disable RT-dependent tolerance
   - Add `min_psms_for_dependent` threshold

---

## Verification Plan

1. **Unit Tests:**
   - Test `IrtDependentTolerance` evaluation
   - Test fallback to constant tolerance
   - Test spline fitting with synthetic data

2. **Integration Tests:**
   - Run SearchDIA with RT-dependent tolerance enabled
   - Compare PSM counts vs constant tolerance
   - Check that early/late RT regions have appropriate tolerances

3. **Visual Verification:**
   - Plot tolerance vs iRT curve for each file
   - Compare to FWHM distribution across gradient

---

## Key Files to Modify

| File | Changes |
|------|---------|
| `SearchTypes.jl` | Add `IrtDependentTolerance`, update storage type |
| `FirstPassSearch/utils.jl` | Add `get_irt_errs_dependent()` |
| `SecondPassSearch/utils.jl` | Update 2 locations to evaluate at iRT |
| `HuberTuningSearch/utils.jl` | Update 1 location |
| `IntegrateChromatogramsSearch/utils.jl` | Update 2 locations |
| `LibrarySearch.jl` | Update 2 locations |
| `ParameterTuningSearch.jl` | Update 2 locations |
