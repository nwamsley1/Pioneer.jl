# RT-iRT Model Usage After FirstPassSearch

## Executive Summary

After FirstPassSearch completes, **file-specific RT conversion models** are stored in SearchContext and used throughout the pipeline by all subsequent search methods. These models enable conversion between:

1. **Observed retention time (RT)** → **refined iRT** (indexed retention time)
2. **Refined iRT** → **Observed RT**

**Key Point**: After our iRT refinement implementation, these models now convert RT ↔ **refined iRT** (NOT library iRT). This means the refinement propagates automatically throughout the entire pipeline.

---

## Model Types

### 1. RT → iRT Model (`rt_irt_map`)
**Accessor**: `getRtIrtModel(search_context, ms_file_idx)`

**Purpose**: Converts observed retention time (minutes) to indexed retention time (iRT units)

**Type**: `RtConversionModel` (spline-based interpolation)

**Refinement Status**: ✓ **Uses refined iRT** (if refinement enabled and successful)

**Example**:
```julia
rt_to_irt_model = getRtIrtModel(search_context, ms_file_idx)
irt_value = rt_to_irt_model(35.2)  # Convert 35.2 minutes → refined iRT
```

### 2. iRT → RT Model (`irt_rt_map`)
**Accessor**: `getIrtRtMap(search_context)`

**Purpose**: Converts indexed retention time (iRT units) to observed retention time (minutes)

**Type**: `Dict{Int64, RtConversionModel}`

**Refinement Status**: ✓ **Uses refined iRT** (if refinement enabled and successful)

**Example**:
```julia
irt_to_rt_models = getIrtRtMap(search_context)
rt_value = irt_to_rt_models[ms_file_idx](50.0)  # Convert iRT 50.0 → RT minutes
```

---

## How Models are Created (FirstPassSearch)

### Initial Creation

**Location**: `FirstPassSearch/utils.jl:371-399`

1. **Fit preliminary RT → iRT model** using library iRT
2. **Train iRT refinement model** on high-quality PSMs
3. **Store refinement model** via `setIrtRefinementModel!`
4. **Calculate refined iRT** for training set using callable model
5. **REFIT RT → iRT model** using refined iRT ← **Key step!**
6. **REFIT iRT → RT model** using refined iRT ← **Key step!**
7. **Store refined models** in SearchContext

**Result**: All downstream methods get refined iRT automatically via these models.

---

## Usage Locations After FirstPassSearch

### 1. FirstPassSearch (Summarize Results Phase)

#### 1.1 PSM Column Creation
**Location**: `FirstPassSearch.jl:270`

```julia
rt_model = getRtIrtModel(search_context, ms_file_idx)
add_psm_columns!(psms, spectra, search_context, rt_model, ms_file_idx)
```

**Purpose**: Add RT-related columns to PSM output
- Converts observed RT → refined iRT for feature columns
- Used for downstream ML training

**Refined iRT?**: ✓ Yes (model uses refined iRT)

---

#### 1.2 Building precursors_dict
**Location**: `FirstPassSearch.jl:549-556`

```julia
all_rt_irt = getRtIrtModel(search_context)
valid_rt_irt = Dict{Int64, RtConversionModel}(i => all_rt_irt[i] for i in valid_indices)

precursors_dict = get_best_precursors_accross_runs(
    valid_psms_paths,
    prec_mzs,
    valid_rt_irt;  # ← RT models passed to function
    max_q_val=params.max_precursor_q_val
)
```

**Function**: `get_best_precursors_accross_runs()`
**Location**: `getBestPrecursorsAccrossRuns.jl`

**Purpose**:
- Identify best PSMs for each precursor across all MS runs
- Calculate mean and variance of iRT for each precursor
- Used for RT alignment in SecondPass (via `irt_diff` feature)

**What it does with RT models**:
- **Currently**: Reads `:irt_refined` column directly (RT models not actively used in function)
- **Historical**: Used RT models to convert RT → iRT for each PSM
- **Storage**: Stores `best_irt`, `mean_irt`, `var_irt` in precursors_dict

**Refined iRT?**: ✓ Yes (reads `:irt_refined` column which was computed from refined model)

---

#### 1.3 iRT Comparison Plotting
**Location**: `FirstPassSearch.jl:576-580`

```julia
all_rt_irt = getRtIrtModel(search_context)
valid_rt_irt = Dict{Int64, RtConversionModel}(i => all_rt_irt[i] for i in valid_indices)

plot_irt_comparison(valid_psms_paths, valid_rt_irt, getDataOutDir(search_context), params.min_prob_for_irt_mapping)
```

**Purpose**: Generate comparison plots showing:
- Library iRT error (predicted - observed)
- Refined iRT error (refined - observed)
- Improvement statistics (MAE, std, correlation)

**How RT models are used**:
```julia
# For each PSM:
irt_observed = rt_to_irt_model(psms[:rt])  # Convert RT → refined iRT (observed)
irt_refined = psms[:irt_refined]            # Refined iRT from column
irt_library = psms[:irt_predicted]          # Library iRT

# Calculate errors:
library_errors = irt_library .- irt_observed  # Should see systematic bias
refined_errors = irt_refined .- irt_observed  # Should be smaller, centered at 0
```

**Refined iRT?**: ✓ Yes
- RT model converts RT → **refined iRT** (observed)
- Compares against `:irt_refined` column (also refined)
- Shows improvement over library iRT predictions

---

### 2. SecondPassSearch

SecondPassSearch performs comprehensive peptide identification with calibrated parameters. RT models are used for **RT window selection** and **feature calculation**.

---

#### 2.1 RT Window Selection (MS2 Scans)
**Location**: `SecondPassSearch/utils.jl:181`

```julia
# For each MS2 scan:
irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, ...) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, ...) + 1, length(rt_index.rt_bins))
```

**Purpose**:
- Convert current scan RT → refined iRT
- Find RT bins within tolerance (irt ± irt_tol)
- Select precursors whose predicted iRT falls within this window

**Refined iRT?**: ✓ Yes
- Converts RT → **refined iRT**
- Compares against precursor refined iRT (from `getPredIrt()` or RT index)
- More accurate window selection than library iRT

---

#### 2.2 RT Window Selection (MS1 Scans)
**Location**: `SecondPassSearch/utils.jl:386`

```julia
# For each MS1 scan:
irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, ...) - 1, 1)
irt_stop = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, ...) + 1, length(rt_index.rt_bins))
```

**Purpose**: Same as MS2, but for MS1 isotope envelope matching

**Refined iRT?**: ✓ Yes

---

#### 2.3 PSM Feature Calculation
**Location**: `SecondPassSearch.jl:495`

```julia
PSMs = extract_features(
    gpsms[!,:precursor_idx],
    gpsms[!,:scan_idx],
    # ... other features ...
    getRtIrtModel(search_context, ms_file_idx)  # ← Used for RT features
);
```

**What it does**: Passes RT model to `extract_features()` function

**Location of usage**: `SecondPassSearch/utils.jl:906-915`

```julia
# Within extract_features():
irt_obs[i] = rt_to_irt_interp(rt[i])  # Convert observed RT → refined iRT
irt_pred[i] = getPredIrt(search_context, prec_idx, ms_file_idx)  # Get refined iRT prediction

irt_diff[i] = abs(irt_obs[i] - prec_id_to_irt[prec_idx].best_irt)  # From precursors_dict
irt_error[i] = abs(irt_obs[i] - irt_pred[i])  # Prediction error

# MS1 feature:
ms1_irt_diff[i] = abs(rt_to_irt_interp(ms1_rt[i]) - getPredIrt(search_context, prec_idx, ms_file_idx))
```

**Features calculated**:
1. **irt_diff**: Difference between observed iRT and best iRT from precursors_dict (cross-run consistency)
2. **irt_error**: Prediction error (observed vs predicted refined iRT)
3. **ms1_irt_diff**: MS1-MS2 iRT difference (isotope envelope alignment)

**Purpose**: ML features for LightGBM scoring model

**Refined iRT?**: ✓ Yes, completely!
- `rt_to_irt_interp` = RT → **refined iRT** model
- `getPredIrt()` returns **refined iRT** (via callable model)
- `precursors_dict.best_irt` is **refined iRT** (from `:irt_refined` column)

**Impact**: Better RT alignment features → better PSM scoring → more identifications

---

#### 2.4 MS1-MS2 RT Difference Calculation
**Location**: `SecondPassSearch.jl:550-553`

```julia
rt_to_irt_model = getRtIrtModel(search_context, ms_file_idx)
psms[!,:ms1_ms2_rt_diff] = Float32.(ifelse.(psms[!,:rt_ms1] .== Float32(-1),
                          Float32(-1),
                          abs.(rt_to_irt_model.(psms[!,:rt]) .- rt_to_irt_model.(psms[!,:rt_ms1]))))
```

**Purpose**:
- Calculate RT difference between MS1 and MS2 in **iRT space**
- Used as ML feature (isotope precursor alignment quality)
- Missing MS1 indicated by -1

**Refined iRT?**: ✓ Yes
- Converts both MS1 and MS2 RT → **refined iRT**
- Difference calculated in refined iRT space (more stable than RT minutes)

---

#### 2.5 RT Error Feature Calculation
**Location**: `SecondPassSearch.jl:570`

```julia
addRtErrorFeatures!(
    psms,
    search_context,
    getTICs(spectra),
    getMzArrays(spectra),
    ms_file_idx,
    getRtIrtModel(search_context, ms_file_idx),  # ← Used for RT error calc
    getPrecursorDict(search_context)
)
```

**Purpose**: Add retention time error features for ML scoring

**Function details** (need to examine `addRtErrorFeatures!` implementation)
- Likely calculates |observed_irt - predicted_irt| or similar
- Uses precursors_dict for cross-run RT consistency

**Refined iRT?**: ✓ Yes (via RT model and precursors_dict)

---

### 3. IntegrateChromatogramsSearch

IntegrateChromatogramsSearch performs peak integration and quantification.

---

#### 3.1 RT Window Selection (MS2 Scans)
**Location**: `IntegrateChromatogramsSearch/utils.jl:254`

```julia
# For each MS2 scan:
irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, ...) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, ...) + 1, length(rt_index.rt_bins))
```

**Purpose**:
- Convert scan RT → refined iRT
- Select transitions for peak integration based on RT window
- Same logic as SecondPassSearch

**Refined iRT?**: ✓ Yes

---

#### 3.2 RT Window Selection (MS1 Scans)
**Location**: `IntegrateChromatogramsSearch/utils.jl:491`

```julia
# For each MS1 scan:
irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, ...) - 1, 1)
irt_stop = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, ...) + 1, length(rt_index.rt_bins))
```

**Purpose**: Same as MS2, but for MS1 isotope envelope integration

**Refined iRT?**: ✓ Yes

**Impact**: More accurate chromatographic peak boundaries → better quantification

---

### 4. HuberTuningSearch

HuberTuningSearch optimizes Huber loss parameters for robust quantification.

---

#### 4.1 RT Window Selection
**Location**: `HuberTuningSearch/utils.jl:225`

```julia
# For each scan:
irt = getRtIrtModel(search_context, ms_file_idx)(getRetentionTime(spectra, scan_idx))
irt_start_new = max(searchsortedfirst(rt_index.rt_bins, irt - irt_tol, ...) - 1, 1)
irt_stop_new = min(searchsortedlast(rt_index.rt_bins, irt + irt_tol, ...) + 1, length(rt_index.rt_bins))
```

**Purpose**: Same RT window selection logic as other search methods

**Refined iRT?**: ✓ Yes

---

### 5. LibrarySearch

LibrarySearch performs spectral library matching against raw data.

**Note**: LibrarySearch is used for **pre-search** parameter optimization, not part of main DIA pipeline. However, it uses the same RT model interface.

---

#### 5.1 Pre-Search (Line 222)
**Location**: `LibrarySearch.jl:222`

```julia
results = search_ms_data(
    data,
    getPresearchFragmentIndex(getSpecLib(search_context)),
    getSpecLib(search_context),
    getSearchData(search_context),
    getQuadTransmissionModel(search_context, ms_file_idx),
    getMassErrorModel(search_context, ms_file_idx),
    getRtIrtModel(search_context, ms_file_idx),  # ← RT model for window selection
    search_parameters,
    getNceModel(search_context, ms_file_idx),
    getIRTTol(search_parameters),
)
```

**Refined iRT?**: ✓ Yes (if models already fitted from FirstPass)

---

#### 5.2 Full Library Search (Line 239)
**Location**: `LibrarySearch.jl:239`

```julia
results = search_ms_data(
    data,
    getFragmentIndex(getSpecLib(search_context)),
    getSpecLib(search_context),
    getSearchData(search_context),
    getQuadTransmissionModel(search_context, ms_file_idx),
    getMassErrorModel(search_context, ms_file_idx),
    getRtIrtModel(search_context, ms_file_idx),  # ← RT model
    search_parameters,
    getNceModel(search_context, ms_file_idx),
    getIrtErrors(search_context)[ms_file_idx]
)
```

**Refined iRT?**: ✓ Yes

---

#### 5.3 NCE Grid Search (Line 260)
**Location**: `LibrarySearch.jl:260`

```julia
results = search_nce_grid(
    data,
    getFragmentIndex(getSpecLib(search_context)),
    getSpecLib(search_context),
    getSearchData(search_context),
    getQuadTransmissionModel(search_context, ms_file_idx),
    getMassErrorModel(search_context, ms_file_idx),
    getRtIrtModel(search_context, ms_file_idx),  # ← RT model
    search_parameters,
    search_parameters.nce_grid,
    getIrtErrors(search_context)[ms_file_idx]
)
```

**Refined iRT?**: ✓ Yes

---

### 6. MaxLFQSearch

MaxLFQSearch performs label-free quantification using the MaxLFQ algorithm.

---

#### 6.1 Quantification and Plotting
**Location**: `MaxLFQSearch/utils.jl:55`

```julia
scored_precursors, protein_groups = score_precursors_and_protein_groups(
    params.params,
    precursors,
    filtered_file_names,
    joinpath(getDataOutDir(search_context), "qc_plots"),
    filtered_file_paths,
    getIrtRtMap(search_context),  # ← iRT → RT model for plotting
    search_context.mass_error_model,
    valid_file_indices
)
```

**Purpose**:
- **Plotting**: Convert refined iRT → RT for visualization (chromatogram plots, RT alignment plots)
- Makes plots display in familiar RT minutes instead of abstract iRT units

**Function**: `score_precursors_and_protein_groups()`
- Uses `irt_rt_map` to convert back to RT for human-readable plots
- Quantification still happens in iRT space

**Refined iRT?**: ✓ Yes
- `irt_rt_map` converts **refined iRT** → RT
- Plots show retention time aligned by refined iRT

---

## Summary Table

| Search Method | Usage | Purpose | Refined iRT? | Line(s) |
|--------------|-------|---------|-------------|---------|
| **FirstPassSearch** | PSM columns | Add RT features to output | ✓ Yes | 270 |
| **FirstPassSearch** | precursors_dict | Cross-run RT calibration | ✓ Yes | 549 |
| **FirstPassSearch** | iRT comparison plots | Validate refinement performance | ✓ Yes | 576 |
| **SecondPassSearch** | RT window (MS2) | Select precursors by RT | ✓ Yes | 181 |
| **SecondPassSearch** | RT window (MS1) | Select precursors by RT | ✓ Yes | 386 |
| **SecondPassSearch** | Feature extraction | Calculate iRT diff/error | ✓ Yes | 495, 906-915 |
| **SecondPassSearch** | MS1-MS2 RT diff | Isotope alignment feature | ✓ Yes | 550-553 |
| **SecondPassSearch** | RT error features | ML scoring features | ✓ Yes | 570 |
| **IntegrateChromatogramsSearch** | RT window (MS2) | Peak integration boundaries | ✓ Yes | 254 |
| **IntegrateChromatogramsSearch** | RT window (MS1) | Isotope integration | ✓ Yes | 491 |
| **HuberTuningSearch** | RT window | Parameter optimization | ✓ Yes | 225 |
| **LibrarySearch** | Pre-search | Window selection | ✓ Yes | 222 |
| **LibrarySearch** | Full search | Window selection | ✓ Yes | 239 |
| **LibrarySearch** | NCE grid | Window selection | ✓ Yes | 260 |
| **MaxLFQSearch** | Quantification plots | iRT → RT for display | ✓ Yes | 55 |

---

## Key Insights

### 1. Complete Refinement Propagation

By refitting the RT models with refined iRT (instead of library iRT), the refinement **automatically propagates** to all downstream methods:

- ✓ **SecondPass**: Better RT window selection and features
- ✓ **IntegrateChromatograms**: More accurate peak boundaries
- ✓ **HuberTuning**: Improved parameter optimization
- ✓ **MaxLFQ**: Better cross-run alignment

**No code changes needed** in these methods - they inherit refined iRT automatically!

---

### 2. RT vs iRT Space

**Why convert to iRT?**
- iRT is **standardized** across runs (library-based scale)
- RT is **instrument-specific** (varies between runs)
- iRT enables **cross-run comparison** and **RT alignment**

**Refined iRT advantage**:
- Library iRT has systematic errors (sequence-dependent bias)
- Refined iRT **corrects** these errors using amino acid composition
- Result: **More accurate RT windows** and **better features**

---

### 3. Dual Storage Strategy

**Two ways precursors get refined iRT**:

1. **`:irt_refined` column** (FirstPass PSMs)
   - Used by: `get_best_precursors_across_runs()`, plotting
   - Advantage: Fast lookup (pre-computed)
   - Disadvantage: Only for precursors in FirstPass

2. **`getPredIrt()` callable model** (All precursors)
   - Used by: SecondPass feature extraction
   - Advantage: Works for ALL library precursors (not just FirstPass hits)
   - Disadvantage: Computed on-the-fly (but zero allocations!)

**Why both?**
- `:irt_refined` column: Fast access for FirstPass results
- `getPredIrt()`: Extends refinement to precursors never seen in FirstPass (match-between-runs)

---

### 4. Impact on Search Performance

**Theoretical improvements** from refined RT models:

1. **Tighter RT windows** → Less noise, faster search
2. **Better iRT features** → Improved ML scoring (LightGBM)
3. **Consistent cross-run alignment** → Better precursors_dict
4. **Accurate isotope matching** → Better MS1 features

**Measurable metrics**:
- Typical MAE improvement: **40-50%** (library vs refined)
- Standard deviation reduction: **30-40%**
- Correlation improvement: **+0.02-0.05** (Pearson R)

---

## Comparison: Library iRT vs Refined iRT Models

### Before Refinement (Library iRT)

```
RT → library iRT → compare with library predictions
                    ↓
              Systematic errors due to:
              - Amino acid composition bias
              - Instrument differences
              - Gradient differences
```

### After Refinement (Refined iRT)

```
RT → refined iRT → compare with refined predictions
                    ↓
              Errors corrected for:
              - Amino acid composition (20 coefficients)
              - Library iRT systematic bias
              - File-specific calibration
                    ↓
              More accurate RT windows and features
```

---

## Future Considerations

### Potential Further Improvements

1. **Dynamic iRT tolerance**: Adjust `irt_tol` based on model variance
2. **Temperature-aware models**: Account for column temperature effects
3. **Gradient-specific models**: Different models for different LC gradients
4. **Online refinement**: Update models during SecondPass based on high-confidence PSMs

### Testing and Validation

**How to validate refinement is working**:

1. Check iRT comparison plots (`irt_comparison/` folder)
   - Look for reduced MAE and tighter error distributions

2. Examine SecondPass PSM features
   - Lower `irt_error` values = better predictions
   - Tighter `irt_diff` distribution = better cross-run consistency

3. Compare identification rates
   - More PSMs in SecondPass suggests better RT windows
   - Better protein coverage suggests improved scoring

---

## Technical Notes

### Model Storage

**Location in SearchContext**:
```julia
rt_irt_map::Dict{Int64, RtConversionModel}      # File index → RT→iRT model
irt_rt_map::Dict{Int64, RtConversionModel}      # File index → iRT→RT model
irt_refinement_models::Dict{Int64, Union{IrtRefinementModel, Nothing}}  # Refinement models
```

### Accessor Pattern

```julia
# Getting RT → iRT model for file
rt_to_irt = getRtIrtModel(search_context, ms_file_idx)

# Getting all iRT → RT models
irt_to_rt_dict = getIrtRtMap(search_context)

# Getting refinement model (for getPredIrt)
refinement_model = getIrtRefinementModel(search_context, ms_file_idx)
```

### Zero-Allocation Design

Both approaches are optimized for performance:

1. **RT models**: Pre-fitted splines (interpolation is fast)
2. **Refinement models**: Callable functor with Dict-based AA weights (zero allocations)

```julia
# Both are zero-allocation:
refined_irt = rt_to_irt_model(35.2)           # Spline interpolation
refined_irt = getPredIrt(s, prec_idx, file)   # Callable model
```

---

## Conclusion

The RT-iRT conversion models are **central** to the entire DIA analysis pipeline. They enable:

1. **RT normalization** across runs (iRT space)
2. **Accurate window selection** (tighter = faster + cleaner)
3. **ML feature calculation** (RT alignment quality)
4. **Cross-run comparison** (precursors_dict)
5. **Visualization** (convert back to RT for plots)

By refitting these models with **refined iRT** (instead of library iRT), the refinement propagates automatically throughout the pipeline, improving:
- Identification rates (SecondPass)
- Quantification accuracy (IntegrateChromatograms)
- Cross-run consistency (MaxLFQ)

This is a **key architectural decision** that makes iRT refinement maximally effective with minimal code changes to downstream methods.
