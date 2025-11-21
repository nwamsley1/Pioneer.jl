# iRT Refinement Implementation Plan - FINAL VERSION

## Overview
Implement file-specific iRT refinement using amino acid composition to improve retention time prediction accuracy. After FirstPassSearch, train a linear model to predict the error between library iRT and observed iRT based on amino acid counts and library iRT. Use this model to refine iRT values for all subsequent searches.

---

## Key Concepts

### Terminology
- **library_irt** (`:irt_predicted` column): Original iRT from spectral library
- **observed_irt**: RT converted to iRT via `rt_to_irt` spline (library-based)
- **refined_irt** (`:refined_irt` column): Library iRT corrected using refinement model
- **observed_refined_irt** (`:irt` column in SecondPassSearch+): RT converted via `rt_to_refined_irt` spline

### Spline/Model Storage in SearchContext

| Field Name | Purpose | When Created | Used By |
|------------|---------|--------------|---------|
| `rt_irt_map` | RT → library_irt | FirstPassSearch | FirstPassSearch only |
| `irt_rt_map` | library_irt → RT (inverse) | FirstPassSearch | FirstPassSearch only |
| `rt_to_refined_irt_map` | RT → refined_irt | FirstPassSearch (if refinement succeeds) | SecondPassSearch+ |
| `refined_irt_to_rt_map` | refined_irt → RT (inverse) | FirstPassSearch (if refinement succeeds) | Plotting/debugging |
| `irt_refinement_models` | Per-sequence refinement: library_irt → refined_irt | FirstPassSearch | All searches (apply to library values) |

### Data Flow
```
FirstPassSearch Workflow:
  1. Fit rt_to_irt spline (RT → library_irt), store in rt_irt_map
  2. Fit irt_to_rt spline (library_irt → RT), store in irt_rt_map
  3. Calculate observed_irt = rt_to_irt(rt) for top PSMs
  4. Train IrtRefinementModel: predict (library_irt - observed_irt) from AA counts + library_irt
  5. If refinement improves MAE:
     a. Apply model to get refined_irt for top PSMs
     b. Fit rt_to_refined_irt spline (RT → refined_irt), store in rt_to_refined_irt_map
     c. Fit refined_irt_to_rt spline (refined_irt → RT), store in refined_irt_to_rt_map
  6. Add :refined_irt column to PSM files (applies model to all PSMs)
  7. Build RT indices using refined_irt values

SecondPassSearch+ (all subsequent searches):
  - Get RT → refined_irt spline from rt_to_refined_irt_map
  - Convert scan RT → observed_refined_irt using rt_to_refined_irt spline
  - Get IrtRefinementModel and apply to library_irt → refined_irt for search windows
  - Calculate irt_error = |observed_refined_irt - refined_irt|
```

### Column Names Across Searches

| Column | FirstPassSearch | SecondPassSearch+ |
|--------|----------------|-------------------|
| `:rt` | Observed RT | Observed RT |
| `:irt_predicted` | Library iRT | Library iRT |
| `:irt` | Observed library iRT (from rt_irt_map) | Observed refined iRT (from rt_to_refined_irt_map) |
| `:refined_irt` | Added post-search | Already present |

**Key insight**: The `:irt` column semantics change between FirstPassSearch (uses library spline) and later searches (uses refined spline). This is acceptable as searches access different splines from SearchContext.

---

## Implementation Steps

### 1. Create IrtRefinementModel Struct

**File:** `src/structs/RetentionTimeConversionModel.jl`

**Location:** After line 58 (after `RtConversionModel() = IdentityModel()`)

**ADD THIS CODE:**

```julia
"""
    IrtRefinementModel

File-specific model to refine library iRT predictions using amino acid composition.

# Fields
- `use_refinement::Bool`: Whether refinement improves validation MAE
- `aa_coefficients::Dict{Char, Float32}`: Per-AA weights (20 standard AAs)
- `intercept::Float32`: Model intercept
- `irt_coefficient::Float32`: Weight for library_irt feature
- `mae_original::Float32`: Validation MAE without refinement
- `mae_refined::Float32`: Validation MAE with refinement
- `r2_train::Float32`: Training R²
- `r2_val::Float32`: Validation R²

# Callable Interface
Model is callable: `refined_irt = model(sequence::String, library_irt::Float32)`

# Algorithm
Predicts error = library_irt - observed_irt, then:
refined_irt = library_irt - predicted_error

# Example
```julia
model = IrtRefinementModel(true, aa_weights, 0.5f0, 0.1f0, ...)
refined = model("PEPTIDE", 50.0f0)  # Returns refined iRT
```
"""
struct IrtRefinementModel
    use_refinement::Bool
    aa_coefficients::Dict{Char, Float32}
    intercept::Float32
    irt_coefficient::Float32
    mae_original::Float32
    mae_refined::Float32
    r2_train::Float32
    r2_val::Float32
end

"""
    (model::IrtRefinementModel)(sequence::String, library_irt::Float32) -> Float32

Apply iRT refinement to a sequence. Zero-allocation via Dict lookup.
"""
function (model::IrtRefinementModel)(sequence::String, library_irt::Float32)::Float32
    if !model.use_refinement
        return library_irt
    end

    # Calculate predicted error
    error_pred = model.intercept + model.irt_coefficient * library_irt

    # Add AA contributions
    for aa in sequence
        if haskey(model.aa_coefficients, aa)
            error_pred += model.aa_coefficients[aa]
        end
    end

    # Return refined iRT (subtract predicted error)
    return library_irt - error_pred
end

# 20 standard amino acids for model training
const STANDARD_AAS = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
                       'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
```

---

### 2. Update SearchContext

**File:** `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl`

#### Change 2a: Add fields to struct (after line 223)

**CURRENT CODE:**
```julia
    irt_rt_map::Dict{Int64, RtConversionModel}
    rt_irt_map::Dict{Int64, RtConversionModel}
    precursor_dict::Base.Ref{Dictionary}
```

**CHANGED CODE:**
```julia
    irt_rt_map::Dict{Int64, RtConversionModel}
    rt_irt_map::Dict{Int64, RtConversionModel}
    rt_to_refined_irt_map::Dict{Int64, RtConversionModel}              # NEW FIELD
    refined_irt_to_rt_map::Dict{Int64, RtConversionModel}              # NEW FIELD
    irt_refinement_models::Dict{Int64, Union{IrtRefinementModel, Nothing}}  # NEW FIELD
    precursor_dict::Base.Ref{Dictionary}
```

#### Change 2b: Initialize in constructor (lines 267-269)

**CURRENT CODE:**
```julia
            Dict{Int64, RtConversionModel}(),  # line 267: irt_rt_map
            Dict{Int64, RtConversionModel}(),  # line 268: rt_irt_map
            Ref{Dictionary}(),                 # line 269: precursor_dict
```

**CHANGED CODE:**
```julia
            Dict{Int64, RtConversionModel}(),                          # line 267: irt_rt_map
            Dict{Int64, RtConversionModel}(),                          # line 268: rt_irt_map
            Dict{Int64, RtConversionModel}(),                          # NEW LINE: rt_to_refined_irt_map
            Dict{Int64, RtConversionModel}(),                          # NEW LINE: refined_irt_to_rt_map
            Dict{Int64, Union{IrtRefinementModel, Nothing}}(),         # NEW LINE: irt_refinement_models
            Ref{Dictionary}(),                                         # precursor_dict (now later line)
```

#### Change 2c: Add getters/setters (after line 524)

**ADD AFTER setRtIrtMap!:**
```julia
"""
    getRtToRefinedIrtMap(s::SearchContext)

Get dictionary of RT → refined_iRT conversion models for all files.
"""
getRtToRefinedIrtMap(s::SearchContext) = s.rt_to_refined_irt_map

"""
    getRefinedIrtToRtMap(s::SearchContext)

Get dictionary of refined_iRT → RT conversion models for all files.
"""
getRefinedIrtToRtMap(s::SearchContext) = s.refined_irt_to_rt_map

"""
    setRtToRefinedIrtMap!(s::SearchContext, rcm::RtConversionModel, index::Integer)

Store RT → refined_iRT conversion model for MS file index.
"""
function setRtToRefinedIrtMap!(s::SearchContext, rcm::RtConversionModel, index::I) where {I<:Integer}
    s.rt_to_refined_irt_map[index] = rcm
end

"""
    setRefinedIrtToRtMap!(s::SearchContext, rcm::RtConversionModel, index::Integer)

Store refined_iRT → RT conversion model for MS file index.
"""
function setRefinedIrtToRtMap!(s::SearchContext, rcm::RtConversionModel, index::I) where {I<:Integer}
    s.refined_irt_to_rt_map[index] = rcm
end

"""
    getIrtRefinementModel(s::SearchContext, index::Integer)

Get iRT refinement model for MS file index. Returns nothing if not found.
"""
function getIrtRefinementModel(s::SearchContext, index::I) where {I<:Integer}
    return get(s.irt_refinement_models, index, nothing)
end

"""
    setIrtRefinementModel!(s::SearchContext, model::Union{IrtRefinementModel, Nothing}, index::Integer)

Store iRT refinement model for MS file index.
"""
function setIrtRefinementModel!(s::SearchContext, model::Union{IrtRefinementModel, Nothing}, index::I) where {I<:Integer}
    s.irt_refinement_models[index] = model
end
```

#### Change 2d: Update setPrecursorDict! signature (line 525)

**CURRENT CODE:**
```julia
function setPrecursorDict!(s::SearchContext, dict::Dictionary{UInt32, @NamedTuple{best_prob::Float32, best_ms_file_idx::UInt32, best_scan_idx::UInt32, best_irt::Float32, mean_irt::Union{Missing, Float32}, var_irt::Union{Missing, Float32}, n::Union{Missing, UInt16}, mz::Float32}})
    s.precursor_dict[] = dict
end
```

**CHANGED CODE:**
```julia
function setPrecursorDict!(s::SearchContext, dict::Dictionary{UInt32, @NamedTuple{best_prob::Float32, best_ms_file_idx::UInt32, best_scan_idx::UInt32, best_refined_irt::Float32, mean_refined_irt::Union{Missing, Float32}, var_refined_irt::Union{Missing, Float32}, n::Union{Missing, UInt16}, mz::Float32}})
    s.precursor_dict[] = dict
end
```

---

### 3. Create Helper Functions for iRT Refinement

**File:** `src/Routines/SearchDIA/CommonSearchUtils/irt_refinement_utils.jl` (NEW FILE)

[Full file contents - same as before, includes prepare_features_dataframe, fit_irt_refinement_model, add_refined_irt_column!]

---

### 4. Modify map_retention_times!() in FirstPassSearch

**File:** `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

**CURRENT CODE (lines 345-473):** [existing map_retention_times! function]

**REPLACE WITH:**

```julia
function map_retention_times!(
    search_context::SearchContext,
    results::FirstPassSearchResults,
    params::FirstPassSearchParameters
)

    valid_files = get_valid_file_indices(search_context)
    all_psms_paths = getFirstPassPsms(getMSData(search_context))

    # Get sequences from spectral library (needed for refinement)
    precursors = getPrecursors(getSpecLib(search_context))
    sequences = getSequence(precursors)

    for ms_file_idx in valid_files
        if is_file_failed(search_context, ms_file_idx)
            continue
        end

        psms_path = all_psms_paths[ms_file_idx]
        psms = Arrow.Table(psms_path)
        best_hits = psms[:prob].>params.min_prob_for_irt_mapping

        @user_info "File $ms_file_idx: Fitting RT alignment models..."

        try
            if params.use_robust_fitting
                # === STEP 1: Fit RT → library_iRT spline ===
                @user_info "  Step 1: Fitting RT → library_iRT spline..."
                best_psms_df = DataFrame(
                    rt = psms[:rt][best_hits],
                    irt_predicted = psms[:irt_predicted][best_hits]
                )

                rt_to_library_irt, valid_rt, valid_library_irt, irt_mad = Pioneer.fit_irt_model(
                    best_psms_df;
                    lambda_penalty = Float32(0.1),
                    ransac_threshold = 1000,
                    min_psms = 10,
                    spline_degree = 3,
                    max_knots = 7,
                    outlier_threshold = Float32(5.0)
                )

                # Store in rt_irt_map (RT → library_irt)
                setRtIrtMap!(search_context, rt_to_library_irt, ms_file_idx)

                # === STEP 2: Fit library_iRT → RT inverse spline ===
                @user_info "  Step 2: Fitting library_iRT → RT inverse spline..."
                irt_to_rt_df = DataFrame(
                    rt = valid_library_irt,
                    irt_predicted = valid_rt
                )
                library_irt_to_rt, _, _, _ = Pioneer.fit_irt_model(
                    irt_to_rt_df;
                    lambda_penalty = Float32(0.1),
                    ransac_threshold = 1000,
                    min_psms = 10,
                    spline_degree = 3,
                    max_knots = 7,
                    outlier_threshold = Float32(5.0)
                )
                # Store in irt_rt_map (library_irt → RT)
                setIrtRtMap!(search_context, library_irt_to_rt, ms_file_idx)

                # === STEP 3: Train iRT refinement model ===
                @user_info "  Step 3: Training iRT refinement model..."

                # Calculate observed_irt for best PSMs using library spline
                observed_irt = [rt_to_library_irt(rt) for rt in Float32.(psms[:rt][best_hits])]

                # Get sequences for best PSMs
                best_sequences = [sequences[idx] for idx in psms[:precursor_idx][best_hits]]

                # Train refinement model
                refinement_model = fit_irt_refinement_model(
                    best_sequences,
                    Float32.(psms[:irt_predicted][best_hits]),
                    observed_irt,
                    ms_file_idx=ms_file_idx,
                    min_psms=20,
                    train_fraction=0.67
                )

                # Store refinement model
                setIrtRefinementModel!(search_context, refinement_model, ms_file_idx)

                # === STEP 4: Fit RT → refined_iRT spline (if refinement succeeds) ===
                if !isnothing(refinement_model) && refinement_model.use_refinement
                    @user_info "  Step 4: Fitting RT → refined_iRT spline..."

                    # Apply refinement model to get refined_irt for best PSMs
                    refined_irt_values = [refinement_model(seq, lib_irt)
                                         for (seq, lib_irt) in zip(best_sequences,
                                                                    Float32.(psms[:irt_predicted][best_hits]))]

                    # Fit RT → refined_iRT spline
                    refined_psms_df = DataFrame(
                        rt = Float32.(psms[:rt][best_hits]),
                        irt_predicted = refined_irt_values
                    )

                    rt_to_refined_irt, valid_rt_refined, valid_refined_irt, _ = Pioneer.fit_irt_model(
                        refined_psms_df;
                        lambda_penalty = Float32(0.1),
                        ransac_threshold = 1000,
                        min_psms = 10,
                        spline_degree = 3,
                        max_knots = 7,
                        outlier_threshold = Float32(5.0)
                    )

                    # Store in rt_to_refined_irt_map (NEW FIELD - DO NOT overwrite rt_irt_map!)
                    setRtToRefinedIrtMap!(search_context, rt_to_refined_irt, ms_file_idx)

                    # === STEP 5: Fit refined_iRT → RT inverse spline ===
                    @user_info "  Step 5: Fitting refined_iRT → RT inverse spline..."
                    refined_irt_to_rt_df = DataFrame(
                        rt = valid_refined_irt,
                        irt_predicted = valid_rt_refined
                    )
                    refined_irt_to_rt, _, _, _ = Pioneer.fit_irt_model(
                        refined_irt_to_rt_df;
                        lambda_penalty = Float32(0.1),
                        ransac_threshold = 1000,
                        min_psms = 10,
                        spline_degree = 3,
                        max_knots = 7,
                        outlier_threshold = Float32(5.0)
                    )
                    # Store in refined_irt_to_rt_map (NEW FIELD)
                    setRefinedIrtToRtMap!(search_context, refined_irt_to_rt, ms_file_idx)
                else
                    @user_info "  No refinement applied, RT → refined_iRT splines not created"
                    # Store identity models for refined splines (fallback)
                    setRtToRefinedIrtMap!(search_context, IdentityModel(), ms_file_idx)
                    setRefinedIrtToRtMap!(search_context, IdentityModel(), ms_file_idx)
                end

                # === STEP 6: Add :refined_irt column to PSMs file ===
                @user_info "  Step 6: Adding :refined_irt column to PSMs table..."
                add_refined_irt_column!(psms_path, refinement_model, search_context)

                # Generate plots if requested
                if params.plot_rt_alignment
                    plot_rt_alignment_firstpass(
                        valid_rt, valid_library_irt, rt_to_library_irt, ms_file_idx, getDataOutDir(search_context)
                    )
                end

            else
                # Legacy path (no robust fitting)
                best_rts = psms[:rt][best_hits]
                best_irts = psms[:irt_predicted][best_hits]

                irt_to_rt_spline = UniformSpline(best_rts, best_irts, 3, 5)
                rt_to_irt_spline = UniformSpline(best_irts, best_rts, 3, 5)

                setRtIrtMap!(search_context, SplineRtConversionModel(rt_to_irt_spline), ms_file_idx)
                setIrtRtMap!(search_context, SplineRtConversionModel(irt_to_rt_spline), ms_file_idx)
                setRtToRefinedIrtMap!(search_context, IdentityModel(), ms_file_idx)
                setRefinedIrtToRtMap!(search_context, IdentityModel(), ms_file_idx)
                setIrtRefinementModel!(search_context, nothing, ms_file_idx)

                # Add refined_irt column (copy of irt_predicted since no refinement)
                add_refined_irt_column!(psms_path, nothing, search_context)
            end

        catch e
            file_name = try
                getFileIdToName(getMSData(search_context), ms_file_idx)
            catch
                "file_$ms_file_idx"
            end

            n_good_psms = try
                sum(best_hits)
            catch
                0
            end

            @user_warn "RT mapping failed for MS data file: $file_name ($n_good_psms good PSMs found, need >100 for spline). Using identity RT model."

            # Set identity models for all spline types
            identity_model = IdentityModel()
            setRtIrtMap!(search_context, identity_model, ms_file_idx)
            setIrtRtMap!(search_context, identity_model, ms_file_idx)
            setRtToRefinedIrtMap!(search_context, identity_model, ms_file_idx)
            setRefinedIrtToRtMap!(search_context, identity_model, ms_file_idx)
            setIrtRefinementModel!(search_context, nothing, ms_file_idx)
        end
    end

    # Set identity models for failed files
    ms_data = getMassSpecData(search_context)
    for failed_idx in 1:length(ms_data.file_paths)
        if getFailedIndicator(ms_data, failed_idx)
            file_name = try
                getFileIdToName(getMSData(search_context), failed_idx)
            catch
                "file_$failed_idx"
            end
            @user_warn "Setting identity RT models for failed file: $file_name"
            setRtIrtMap!(search_context, IdentityModel(), failed_idx)
            setIrtRtMap!(search_context, IdentityModel(), failed_idx)
            setRtToRefinedIrtMap!(search_context, IdentityModel(), failed_idx)
            setRefinedIrtToRtMap!(search_context, IdentityModel(), failed_idx)
            setIrtRefinementModel!(search_context, nothing, failed_idx)
        end
    end

    return nothing
end
```

**Key changes:**
- Line ~460: Store in `rt_irt_map` (DO NOT overwrite later)
- Line ~474: Store in `irt_rt_map` (DO NOT overwrite later)
- Line ~523: Store RT → refined_irt in `rt_to_refined_irt_map` (NEW)
- Line ~542: Store refined_irt → RT in `refined_irt_to_rt_map` (NEW)

---

### 5. Update Precursor Dictionary Field Names

[Same as before - change :best_irt to :best_refined_irt, etc.]

---

### 6. Update get_best_precursors_accross_runs

[Same as before - update all field names to use refined_irt]

---

### 7. Update create_rt_indices!()

[Same as before - use :best_refined_irt, pass prec_to_refined_irt]

---

### 8. Update buildRTIndex.jl Comments

[Same as before - clarify comments about refined_irt usage]

---

### 9. Update SecondPassSearch to Use Refined Splines

**File:** `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/SecondPassSearch.jl`

**CHANGES NEEDED:**

SecondPassSearch (and all subsequent searches) must use `rt_to_refined_irt_map` instead of `rt_irt_map`.

**Find where SecondPassSearch accesses RT → iRT conversion:**

Likely in `processChunk!()` or similar, there will be code like:
```julia
rt_to_irt = getRtIrtMap(search_context)[ms_file_idx]
observed_irt = rt_to_irt(observed_rt)
```

**CHANGE TO:**
```julia
rt_to_refined_irt = getRtToRefinedIrtMap(search_context)[ms_file_idx]
observed_refined_irt = rt_to_refined_irt(observed_rt)
```

**Search for these patterns:**
- `getRtIrtMap(search_context)` → change to `getRtToRefinedIrtMap(search_context)`
- Variable names like `rt_to_irt` → change to `rt_to_refined_irt`
- Variable names like `observed_irt` → change to `observed_refined_irt`

---

### 10. Update add_main_search_columns!() Documentation

**File:** `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

**NO CODE CHANGES**, but add documentation comment:

```julia
"""
    add_main_search_columns!(psms::DataFrame, rt_irt::T, ...)

Adds essential columns to PSM DataFrame for scoring and analysis.

# Important Note on rt_irt Parameter
The `rt_irt` parameter contains different splines depending on the search stage:
- FirstPassSearch: Contains RT → library_iRT spline (from rt_irt_map)
- SecondPassSearch+: Should contain RT → refined_iRT spline (from rt_to_refined_irt_map)

The `:irt` column created will contain:
- FirstPassSearch: observed library iRT values
- SecondPassSearch+: observed refined iRT values

...
"""
function add_main_search_columns!(...)
```

---

### 11. Update Imports and Exports

[Same as before]

---

## Summary of Spline/Model Fields

| Field in SearchContext | Populated When | Contains | Used By |
|------------------------|----------------|----------|---------|
| `rt_irt_map` | FirstPassSearch Step 1 | RT → library_irt | FirstPassSearch only |
| `irt_rt_map` | FirstPassSearch Step 2 | library_irt → RT | FirstPassSearch only |
| `irt_refinement_models` | FirstPassSearch Step 3 | IrtRefinementModel(seq, lib_irt) → refined_irt | All searches (for library values) |
| `rt_to_refined_irt_map` | FirstPassSearch Step 4 | RT → refined_irt | SecondPassSearch+ |
| `refined_irt_to_rt_map` | FirstPassSearch Step 5 | refined_irt → RT | Debugging/plotting |

---

## Testing Strategy

1. Verify all 5 fields populated correctly in SearchContext after FirstPassSearch
2. Verify SecondPassSearch uses `rt_to_refined_irt_map` not `rt_irt_map`
3. Verify `:refined_irt` column present in all PSM files
4. Verify RT indices built with refined_irt values
5. Compare MAE before/after refinement

---

## Files to Modify

| File | Changes |
|------|---------|
| `RetentionTimeConversionModel.jl` | Add IrtRefinementModel + STANDARD_AAS |
| `SearchTypes.jl` | Add 3 new fields + getters/setters |
| `irt_refinement_utils.jl` | NEW FILE |
| `FirstPassSearch/utils.jl` | Update map_retention_times!() with 6-step workflow |
| `FirstPassSearch/utils.jl` | Update PrecToIrtType field names |
| `FirstPassSearch/utils.jl` | Update create_rt_indices!() |
| `getBestPrecursorsAccrossRuns.jl` | Update all `irt` → `refined_irt` |
| `buildRTIndex.jl` | Update comments |
| `SecondPassSearch/*.jl` | Change `getRtIrtMap` → `getRtToRefinedIrtMap` |

---

## Next Steps After Approval

Ready to implement with clear, consistent naming throughout!
