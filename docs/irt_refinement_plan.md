# iRT Prediction Refinement Implementation Plan

## Overview

After FirstPassSearch fits the RT→iRT retention time alignment model for each MS file, we can improve library iRT predictions by learning file-specific correction models. These models predict the residual error between observed iRT (from RT transform) and library-predicted iRT based on peptide sequence composition.

**Key Insight**: Systematic biases in library iRT predictions correlate with amino acid composition. By modeling these biases, we can provide more accurate iRT predictions for subsequent search methods, improving precursor candidate selection and reducing search space.

---

## Algorithm Summary

Based on analysis in `/Users/nathanwamsley/Documents/irtRefine/refineIrt.ipynb`:

### 1. Feature Engineering
For each PSM with precursor_idx `pid`:
- **Amino Acid Counts**: Count occurrences of each of 20 standard amino acids in sequence
  - `count_A`, `count_C`, `count_D`, ..., `count_Y` (20 features)
- **Library iRT**: Include original `irt_predicted` value (1 feature)
- **Total Features**: 21

### 2. Target Variable
- **iRT Error (Residual)**: `irt_err = observed_irt - irt_predicted`
  - `observed_irt = rt_to_irt_model(empirical_rt)`
  - `irt_predicted` = library prediction from spectral library

### 3. Model Training (Per-File)
- **Model Type**: Linear regression (GLM.jl with identity link)
- **Training Set**: 2/3 of high-confidence PSMs (random split)
- **Validation Set**: 1/3 of high-confidence PSMs
- **Formula**: `irt_err ~ count_A + count_C + ... + count_Y + irt_predicted`

### 4. Model Validation
- Compare validation set MAE:
  - **MAE_original**: `mean(abs(observed_irt - irt_predicted))`
  - **MAE_refined**: `mean(abs(observed_irt - irt_refined))`
  - Where: `irt_refined = irt_predicted - predicted_irt_err`
- **Decision**: Use refined iRT only if `MAE_refined < MAE_original`

### 5. Final Model Retraining
- **If validation shows improvement**: Retrain model on **full dataset** (train + validation combined)
  - Validation model used only to decide whether refinement helps
  - Final model uses all available data for best parameter estimates
  - This is standard ML best practice: validation for decision, all data for final model

### 6. Application
- Apply **retrained** model to **all precursors** in library (not just observed PSMs)
- Update **SearchContext pred_irt** dictionary with refined values for all precursors
- Add **irt_refined column** to FirstPass PSM files (for transparency/debugging)
- Downstream methods automatically use refined iRT via existing `getPredIrt()` calls

---

## Feature Toggle

**Master Parameter**: `enable_irt_refinement` in `params.json`

```json
{
  "first_pass_search": {
    "enable_irt_refinement": false  // Set to true to enable refinement
  }
}
```

**Default**: `false` (conservative - maintains current behavior)

**Behavior**:
- **`false`** (disabled):
  - FirstPassSearch uses library iRT predictions unchanged
  - No model training, no SearchContext updates
  - Zero performance overhead
  - Identical to current Pioneer behavior

- **`true`** (enabled):
  - Train per-file linear models to predict iRT errors
  - Update SearchContext `pred_irt` with refined predictions
  - Add `irt_refined` column to FirstPass PSM files
  - Downstream methods automatically benefit (no code changes needed)
  - Small overhead: ~1-2 seconds per file

**Recommendation**:
- Start with `false` for initial runs
- Enable (`true`) after validating improvement on a subset of your data
- Compare search results with/without refinement to assess benefit for your specific dataset

---

## Downstream Method Integration

All downstream search methods automatically benefit from refined iRT predictions with **zero code changes** because they access iRT values through SearchContext:

### SecondPassSearch
**File**: `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl`
**How it uses iRT**:
- RT window selection for precursor candidates
- Builds RT index using `getPredIrt(search_context, precursor_idx)`
- **Automatically uses refined iRT** if FirstPassSearch updated SearchContext

### ScoringSearch
**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/`
**How it uses iRT**:
- RT-based scoring features (e.g., RT deviation from expected)
- Accesses via `getPredIrt(search_context, precursor_idx)`
- **Automatically uses refined iRT** for better RT deviation estimates

### IntegrateChromatogramsSearch
**File**: `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/`
**How it uses iRT**:
- Defines RT windows for extracting ion chromatograms (XICs)
- Converts predicted iRT → empirical RT via `getIrtRtMap(search_context, file_idx)`
- **Automatically uses refined iRT** for more accurate XIC extraction windows

**Key Point**: No modifications needed to downstream methods. They already use SearchContext accessor methods that will return refined values after FirstPassSearch updates `pred_irt`.

---

## Integration Point

### Location: FirstPassSearch
**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`
**Function**: `RT_iRT_Calibration!`
**Line**: After line 383 (`setRtIrtMap!(search_context, rt_model, ms_file_idx)`)

**Current Flow**:
```julia
# Line 372-383: Fit RT → iRT model
rt_model, valid_rt, valid_irt, irt_mad = Pioneer.fit_irt_model(...)
setRtIrtMap!(search_context, rt_model, ms_file_idx)

# NEW: Add iRT refinement here
irt_refinement_model = fit_irt_refinement_model(
    psms, rt_model, search_context, ms_file_idx
)

# Line 390-399: Fit inverse iRT → RT model
irt_model, _, _, _ = Pioneer.fit_irt_model(...)
```

---

## Detailed Implementation

### Step 1: Create Utility Module

**New File**: `src/Routines/SearchDIA/CommonSearchUtils/irt_refinement_utils.jl`

```julia
"""
    IrtRefinementModel

Stores linear model for predicting iRT residuals from amino acid composition.

# Fields
- `coefficients::Dict{Symbol, Float64}`: AA count coefficients + intercept
- `irt_coef::Float64`: Coefficient for irt_predicted feature
- `intercept::Float64`: Model intercept
- `r2_train::Float32`: Training R²
- `r2_val::Float32`: Validation R²
- `mae_original::Float32`: Validation MAE without refinement
- `mae_refined::Float32`: Validation MAE with refinement
- `use_refinement::Bool`: Whether refined predictions are better
- `n_train::Int`: Number of training PSMs
- `n_val::Int`: Number of validation PSMs
"""
struct IrtRefinementModel
    coefficients::Dict{Symbol, Float64}
    irt_coef::Float64
    intercept::Float64
    r2_train::Float32
    r2_val::Float32
    mae_original::Float32
    mae_refined::Float32
    use_refinement::Bool
    n_train::Int
    n_val::Int
end

"""
    count_amino_acids(sequence::AbstractString) -> Dict{Symbol, Int}

Count occurrences of each amino acid in a peptide sequence.

# Arguments
- `sequence`: Peptide sequence string (e.g., "PEPTIDE")

# Returns
Dictionary mapping amino acid symbols to counts (e.g., :count_A => 2)
"""
function count_amino_acids(sequence::AbstractString)
    aa_symbols = [:A, :C, :D, :E, :F, :G, :H, :I, :K, :L,
                  :M, :N, :P, :Q, :R, :S, :T, :V, :W, :Y]

    counts = Dict{Symbol, Int}()
    for aa in aa_symbols
        counts[Symbol("count_", aa)] = count(==(aa), sequence)
    end

    return counts
end

"""
    prepare_irt_refinement_features(
        psms::Arrow.Table,
        precursors::LibraryPrecursors
    ) -> DataFrame

Prepare feature matrix for iRT refinement modeling.

# Arguments
- `psms`: Arrow table of PSMs with columns :precursor_idx, :rt, :irt_predicted
- `precursors`: Library precursors (for sequence lookup)

# Returns
DataFrame with columns:
- All amino acid counts (count_A, ..., count_Y)
- irt_predicted
- rt
- precursor_idx
"""
function prepare_irt_refinement_features(
    psms::Arrow.Table,
    precursors::LibraryPrecursors
)
    n_psms = length(psms[:precursor_idx])

    # Initialize feature DataFrame
    df = DataFrame(
        precursor_idx = psms[:precursor_idx],
        rt = psms[:rt],
        irt_predicted = psms[:irt_predicted]
    )

    # Get sequences
    sequences = getSequence(precursors)

    # Add amino acid count columns
    aa_symbols = [:A, :C, :D, :E, :F, :G, :H, :I, :K, :L,
                  :M, :N, :P, :Q, :R, :S, :T, :V, :W, :Y]

    for aa in aa_symbols
        col_name = Symbol("count_", aa)
        df[!, col_name] = zeros(Int, n_psms)
    end

    # Fill in amino acid counts
    for (i, pid) in enumerate(psms[:precursor_idx])
        sequence = sequences[pid]
        for aa in aa_symbols
            col_name = Symbol("count_", aa)
            df[i, col_name] = count(==(aa), sequence)
        end
    end

    return df
end

"""
    fit_irt_refinement_model(
        psms::Arrow.Table,
        rt_model::RtConversionModel,
        search_context::SearchContext,
        ms_file_idx::Integer;
        min_psms_per_feature::Int = 20,
        train_fraction::Float64 = 0.667,
        seed::Int = 42
    ) -> Union{IrtRefinementModel, Nothing}

Fit linear model to predict iRT residuals from amino acid composition.

# Arguments
- `psms`: Arrow table of PSMs from FirstPassSearch
- `rt_model`: Fitted RT → iRT conversion model
- `search_context`: Search context for accessing precursor data
- `ms_file_idx`: Current MS file index

# Keyword Arguments
- `min_psms_per_feature`: Minimum PSMs per feature (default: 20 → 420 total)
- `train_fraction`: Fraction of data for training (default: 2/3)
- `seed`: Random seed for train/test split

# Returns
- `IrtRefinementModel` if successful
- `nothing` if insufficient data or model doesn't improve predictions
"""
function fit_irt_refinement_model(
    psms::Arrow.Table,
    rt_model::RtConversionModel,
    search_context::SearchContext,
    ms_file_idx::Integer;
    min_psms_per_feature::Int = 20,
    train_fraction::Float64 = 0.667,
    seed::Int = 42
)
    # Check minimum data requirement
    n_features = 21  # 20 AA + irt_predicted
    min_psms = n_features * min_psms_per_feature

    if length(psms[:precursor_idx]) < min_psms
        @debug_l1 "File $ms_file_idx: Insufficient PSMs for iRT refinement " *
                  "($(length(psms[:precursor_idx])) < $min_psms)"
        return nothing
    end

    # Prepare features
    precursors = getPrecursors(getSpecLib(search_context))
    features_df = prepare_irt_refinement_features(psms, precursors)

    # Calculate observed iRT and residuals
    features_df[!, :observed_irt] = [rt_model(rt) for rt in features_df.rt]
    features_df[!, :irt_err] = features_df.observed_irt .- features_df.irt_predicted

    # Train/validation split
    Random.seed!(seed)
    n_total = nrow(features_df)
    n_train = round(Int, n_total * train_fraction)
    indices = randperm(n_total)
    train_indices = indices[1:n_train]
    val_indices = indices[n_train+1:end]

    train_df = features_df[train_indices, :]
    val_df = features_df[val_indices, :]

    # Verify training set has enough samples
    if nrow(train_df) < min_psms
        @debug_l1 "File $ms_file_idx: Insufficient training PSMs " *
                  "($(nrow(train_df)) < $min_psms)"
        return nothing
    end

    # Fit linear model
    formula = @formula(irt_err ~ count_A + count_C + count_D + count_E + count_F +
                                 count_G + count_H + count_I + count_K + count_L +
                                 count_M + count_N + count_P + count_Q + count_R +
                                 count_S + count_T + count_V + count_W + count_Y +
                                 irt_predicted)

    model = lm(formula, train_df)

    # Extract coefficients
    coef_names = coefnames(model)
    coef_values = coef(model)

    coefficients = Dict{Symbol, Float64}()
    intercept = coef_values[1]  # First coefficient is intercept
    irt_coef = 0.0

    for (i, name) in enumerate(coef_names[2:end])  # Skip intercept
        if name == "irt_predicted"
            irt_coef = coef_values[i+1]
        else
            coefficients[Symbol(name)] = coef_values[i+1]
        end
    end

    # Evaluate on training set
    train_predictions = GLM.predict(model, train_df)
    train_residuals = train_df.irt_err .- train_predictions
    r2_train = r2(model)

    # Evaluate on validation set
    val_predictions = GLM.predict(model, val_df)
    val_residuals = val_df.irt_err .- val_predictions
    ss_res = sum(val_residuals.^2)
    ss_tot = sum((val_df.irt_err .- mean(val_df.irt_err)).^2)
    r2_val = 1 - (ss_res / ss_tot)

    # Compare MAE: original vs refined
    mae_original = mean(abs.(val_df.observed_irt .- val_df.irt_predicted))

    # Refined iRT = original - predicted_error
    val_df[!, :irt_refined] = val_df.irt_predicted .- val_predictions
    mae_refined = mean(abs.(val_df.observed_irt .- val_df.irt_refined))

    # Decision: use refinement if it improves validation MAE
    use_refinement = mae_refined < mae_original

    @user_info "File $ms_file_idx iRT Refinement:" *
              " Training R²=$(round(r2_train, digits=4))," *
              " Validation R²=$(round(r2_val, digits=4))," *
              " MAE Original=$(round(mae_original, digits=4))," *
              " MAE Refined=$(round(mae_refined, digits=4))," *
              " Use=$(use_refinement)\n"

    # If refinement improves MAE, retrain on FULL dataset for final model
    if use_refinement
        @debug_l1 "File $ms_file_idx: Retraining on full dataset (train + validation)"

        # Retrain on combined train + validation data
        final_model = lm(formula, features_df)

        # Extract coefficients from retrained model
        final_coef_names = coefnames(final_model)
        final_coef_values = coef(final_model)

        final_coefficients = Dict{Symbol, Float64}()
        final_intercept = final_coef_values[1]
        final_irt_coef = 0.0

        for (i, name) in enumerate(final_coef_names[2:end])
            if name == "irt_predicted"
                final_irt_coef = final_coef_values[i+1]
            else
                final_coefficients[Symbol(name)] = final_coef_values[i+1]
            end
        end

        # Use retrained coefficients for final model
        coefficients = final_coefficients
        intercept = final_intercept
        irt_coef = final_irt_coef
    end

    return IrtRefinementModel(
        coefficients,
        irt_coef,
        intercept,
        Float32(r2_train),
        Float32(r2_val),
        Float32(mae_original),
        Float32(mae_refined),
        use_refinement,
        nrow(train_df),
        nrow(val_df)
    )
end

"""
    apply_irt_refinement!(
        psms::Arrow.Table,
        model::IrtRefinementModel,
        precursors::LibraryPrecursors
    ) -> Vector{Float32}

Apply iRT refinement model to all PSMs.

# Arguments
- `psms`: Arrow table of PSMs
- `model`: Fitted IrtRefinementModel
- `precursors`: Library precursors for sequence lookup

# Returns
Vector of refined iRT predictions (same length as psms)
"""
function apply_irt_refinement!(
    psms::Arrow.Table,
    model::IrtRefinementModel,
    precursors::LibraryPrecursors
)
    n_psms = length(psms[:precursor_idx])
    irt_refined = Vector{Float32}(undef, n_psms)
    sequences = getSequence(precursors)

    for i in 1:n_psms
        pid = psms[:precursor_idx][i]
        irt_pred = psms[:irt_predicted][i]
        sequence = sequences[pid]

        # Calculate predicted error
        predicted_err = model.intercept
        predicted_err += model.irt_coef * irt_pred

        for aa_char in sequence
            col_name = Symbol("count_", aa_char)
            if haskey(model.coefficients, col_name)
                predicted_err += model.coefficients[col_name]
            end
        end

        # Refined iRT = original - predicted_error
        irt_refined[i] = Float32(irt_pred - predicted_err)
    end

    return irt_refined
end
```

### Step 2: Integration into FirstPassSearch

**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

**Modify** `RT_iRT_Calibration!` function (~line 350):

```julia
function RT_iRT_Calibration!(
    search_context::SearchContext,
    params::FirstPassSearchParameters
)
    # ... existing code ...

    for ms_file_idx in valid_files
        # ... existing code to fit RT model (lines 355-383) ...

        # NEW: Fit iRT refinement model (after line 383)
        # This entire block is skipped if enable_irt_refinement = false
        if params.enable_irt_refinement
            irt_refinement_model = fit_irt_refinement_model(
                Arrow.Table(psms_path),
                rt_model,
                search_context,
                ms_file_idx;
                min_psms_per_feature = params.irt_refinement_min_psms_per_feature,
                train_fraction = params.irt_refinement_train_fraction,
                seed = params.irt_refinement_seed
            )

            # Store model if successful and improvement verified
            if !isnothing(irt_refinement_model) && irt_refinement_model.use_refinement
                setIrtRefinementModel!(search_context, irt_refinement_model, ms_file_idx)

                # Get precursors and sequences
                precursors = getPrecursors(getSpecLib(search_context))
                sequences = getSequence(precursors)

                # Apply refinement to ALL precursors in library
                # Update SearchContext pred_irt dictionary
                for precursor_idx in 1:length(sequences)
                    irt_original = getIrt(precursors)[precursor_idx]
                    sequence = sequences[precursor_idx]

                    # Calculate predicted error from model
                    predicted_err = irt_refinement_model.intercept
                    predicted_err += irt_refinement_model.irt_coef * irt_original

                    for aa_char in sequence
                        col_name = Symbol("count_", aa_char)
                        if haskey(irt_refinement_model.coefficients, col_name)
                            predicted_err += irt_refinement_model.coefficients[col_name]
                        end
                    end

                    # Refined iRT = original - predicted_error
                    irt_refined = Float32(irt_original - predicted_err)

                    # Update SearchContext (downstream methods will use this)
                    setPredIrt!(search_context, UInt32(precursor_idx), irt_refined)
                end

                # Also add irt_refined column to PSM file (for transparency/debugging)
                psms_table = Arrow.Table(psms_path)
                irt_refined_psms = Vector{Float32}(undef, length(psms_table[:precursor_idx]))

                for (i, pid) in enumerate(psms_table[:precursor_idx])
                    # Get refined value from SearchContext
                    irt_refined_psms[i] = getPredIrt(search_context, pid)
                end

                psms_df = DataFrame(psms_table)
                psms_df[!, :irt_refined] = irt_refined_psms
                writeArrow(psms_path, psms_df)
            end
        end

        # ... rest of existing code (lines 385-449) ...
    end
end
```

### Step 3: Add SearchContext Support

**File**: `src/structs/SearchContext.jl`

```julia
# Add to SearchContext struct
mutable struct SearchContext
    # ... existing fields ...
    irt_refinement_models::Vector{Union{IrtRefinementModel, Nothing}}

    function SearchContext(...)
        # ... existing initialization ...
        new(..., Vector{Union{IrtRefinementModel, Nothing}}(nothing, n_files))
    end
end

# Add getter/setter methods
function getIrtRefinementModel(sc::SearchContext, ms_file_idx::Integer)
    return sc.irt_refinement_models[ms_file_idx]
end

function setIrtRefinementModel!(
    sc::SearchContext,
    model::IrtRefinementModel,
    ms_file_idx::Integer
)
    sc.irt_refinement_models[ms_file_idx] = model
end

function hasIrtRefinementModel(sc::SearchContext, ms_file_idx::Integer)
    return !isnothing(sc.irt_refinement_models[ms_file_idx])
end
```

### Step 4: Add Configuration Parameters

**File**: `src/structs/PioneerParameters.jl`

Add to `FirstPassSearchParameters` section in `params.json`:

```json
{
  "first_pass_search": {
    "enable_irt_refinement": false,
    "irt_refinement_min_psms_per_feature": 20,
    "irt_refinement_train_fraction": 0.667,
    "irt_refinement_seed": 42
  }
}
```

**Parameter Descriptions**:

- **`enable_irt_refinement`** (Boolean, default: `false`):
  - **Master toggle** for iRT refinement feature
  - `false`: Use library iRT predictions as-is (current behavior)
  - `true`: Train per-file correction models and update predictions
  - **Recommendation**: Start with `false`, enable after validating on your data

- **`irt_refinement_min_psms_per_feature`** (Integer, default: `20`):
  - Minimum PSMs per feature required for model training
  - Total minimum PSMs = 20 × 21 features = 420 PSMs
  - Increase for stricter requirements, decrease for small datasets

- **`irt_refinement_train_fraction`** (Float, default: `0.667`):
  - Fraction of PSMs used for training (remaining used for validation)
  - 0.667 = 2/3 train, 1/3 validation

- **`irt_refinement_seed`** (Integer, default: `42`):
  - Random seed for train/validation split (for reproducibility)

**Julia struct** (in parameters file):

```julia
struct FirstPassSearchParameters
    # ... existing fields ...
    enable_irt_refinement::Bool               # Master toggle
    irt_refinement_min_psms_per_feature::Int  # Minimum data requirement
    irt_refinement_train_fraction::Float64    # Train/val split ratio
    irt_refinement_seed::Int                  # Reproducibility seed
end
```

**Behavior When Disabled** (`enable_irt_refinement = false`):
- FirstPassSearch skips all refinement logic
- SearchContext `pred_irt` uses library iRT values unchanged
- No `irt_refined` column added to PSM files
- Downstream methods use original library predictions
- Zero performance overhead

### Step 5: No Downstream Method Changes Required

**Automatic Integration**: Downstream search methods (SecondPassSearch, ScoringSearch, IntegrateChromatogramsSearch) automatically use refined iRT predictions with **zero code changes** because:

1. They access iRT via `getPredIrt(search_context, precursor_idx)`
2. FirstPassSearch updates SearchContext `pred_irt` dictionary with refined values
3. Downstream methods transparently benefit from improved predictions

**Verification**: The `irt_refined` column in FirstPass PSM files allows verification that refinement is being applied correctly, but downstream methods don't read from PSM files - they use SearchContext exclusively.

---

## File Structure

```
src/Routines/SearchDIA/
├── CommonSearchUtils/
│   ├── irt_refinement_utils.jl          # NEW: Core refinement logic
│   └── rt_alignment_utils.jl            # EXISTING: RT↔iRT conversion
├── SearchMethods/
│   ├── FirstPassSearch/
│   │   ├── FirstPassSearch.jl
│   │   └── utils.jl                     # MODIFY: Add refinement call + SearchContext update
│   ├── SecondPassSearch/
│   │   └── utils.jl                     # NO CHANGES: automatically uses refined iRT
│   ├── ScoringSearch/
│   │   └── ...                          # NO CHANGES: automatically uses refined iRT
│   ├── IntegrateChromatogramsSearch/
│   │   └── ...                          # NO CHANGES: automatically uses refined iRT
│   └── ...
└── ...

src/structs/
├── SearchContext.jl                      # MODIFY: Add refinement models field
└── PioneerParameters.jl                  # MODIFY: Add refinement parameters

docs/
└── irt_refinement_plan.md               # THIS FILE
```

---

## Testing Strategy

### Unit Tests

**File**: `test/UnitTests/test_irt_refinement.jl`

```julia
using Test, Pioneer, DataFrames, Arrow

@testset "iRT Refinement Utils" begin
    @testset "Amino Acid Counting" begin
        counts = Pioneer.count_amino_acids("PEPTIDE")
        @test counts[:count_P] == 2
        @test counts[:count_E] == 2
        @test counts[:count_T] == 1
        @test counts[:count_A] == 0
    end

    @testset "Feature Preparation" begin
        # Create mock PSMs
        psms = Arrow.Table(...)
        precursors = ...

        features_df = Pioneer.prepare_irt_refinement_features(psms, precursors)

        @test size(features_df, 2) == 23  # 20 AA + irt_predicted + rt + precursor_idx
        @test all(col -> col in names(features_df),
                 ["count_A", "count_C", ..., "irt_predicted"])
    end

    @testset "Model Fitting" begin
        # Use synthetic data with known bias
        # ... create mock PSMs with systematic error ...

        model = Pioneer.fit_irt_refinement_model(...)

        @test !isnothing(model)
        @test model.r2_val > 0.0
        @test model.mae_refined < model.mae_original
        @test model.use_refinement == true
    end

    @testset "Model Application" begin
        # Test applying model to new PSMs
        irt_refined = Pioneer.apply_irt_refinement!(...)

        @test length(irt_refined) == nrow(psms)
        @test all(isfinite, irt_refined)
    end
end
```

### Integration Test

**File**: `test/IntegrationTests/test_irt_refinement_integration.jl`

```julia
@testset "iRT Refinement Integration" begin
    # Run FirstPassSearch with refinement enabled
    params_path = "test/data/test_irt_refinement_params.json"
    context = SearchDIA(params_path, stop_at="FirstPass")

    # Check that refinement models were created
    ms_file_idx = 1
    @test hasIrtRefinementModel(context, ms_file_idx)

    model = getIrtRefinementModel(context, ms_file_idx)
    @test !isnothing(model)
    @test model.use_refinement == true

    # Check that PSM files have irt_refined column
    psms_path = getFirstPassPsms(getMSData(context))[ms_file_idx]
    psms = Arrow.Table(psms_path)
    @test :irt_refined in propertynames(psms)

    # Verify refined values differ from original
    @test !all(psms[:irt_refined] .== psms[:irt_predicted])
end
```

### Real-World Validation

**Compare search results with/without refinement**:

```julia
# Run two searches on same data
params_no_refine = load_params("params.json")
params_no_refine["first_pass_search"]["enable_irt_refinement"] = false

params_with_refine = load_params("params.json")
params_with_refine["first_pass_search"]["enable_irt_refinement"] = true

# Compare metrics
results_no_refine = SearchDIA(params_no_refine)
results_with_refine = SearchDIA(params_with_refine)

# Expect:
# - More PSMs identified (better RT windows)
# - Higher precursor scores (more confident identifications)
# - Improved protein identification (more precursors → more proteins)
```

---

## Performance Considerations

### Memory Usage
- **Feature Matrix**: ~1.5 KB per PSM (21 features × 8 bytes × 10 overhead)
  - For 10,000 PSMs: ~15 MB
  - Acceptable for typical FirstPass output
- **Model Storage**: ~500 bytes per model (coefficients + metadata)
  - For 100 files: ~50 KB total
  - Negligible

### Computational Cost
- **Feature Extraction**: O(N × L) where N = PSMs, L = avg sequence length
  - ~0.01 ms per PSM
  - For 10,000 PSMs: ~100 ms
- **Model Training**: O(N × F²) where F = 21 features
  - Linear regression with QR decomposition
  - For 10,000 PSMs: ~50 ms
- **Model Application**: O(M × L) where M = total library precursors
  - For 100,000 precursors: ~1 second

**Total Added Time per File**: ~1-2 seconds (negligible compared to search time)

### Parallelization
- **Current**: Per-file sequential (inherent in FirstPassSearch loop)
- **Future**: Could parallelize across files if needed
- **Not a bottleneck**: 1-2 sec/file × 100 files = 2-3 minutes total

---

## Edge Cases and Error Handling

### 1. Insufficient PSMs
**Scenario**: File has < 420 PSMs (20 per feature × 21 features)

**Handling**:
```julia
if length(psms) < min_psms
    @debug_l1 "File $ms_file_idx: Skipping iRT refinement (too few PSMs)"
    return nothing
end
```

**Fallback**: Use library `irt_predicted` as-is

### 2. Model Doesn't Improve MAE
**Scenario**: `MAE_refined >= MAE_original` on validation set

**Handling**:
```julia
use_refinement = mae_refined < mae_original
if !use_refinement
    @debug_l1 "File $ms_file_idx: Refinement doesn't improve MAE, using library iRT"
end
```

**Fallback**: Don't add `irt_refined` column, downstream methods use `irt_predicted`

### 3. Model Fitting Fails
**Scenario**: GLM fitting throws exception (e.g., singular matrix)

**Handling**:
```julia
try
    model = lm(formula, train_df)
catch e
    @user_warn "File $ms_file_idx: iRT refinement failed: $e"
    return nothing
end
```

**Fallback**: Use library `irt_predicted`

### 4. Extreme Predictions
**Scenario**: Model predicts unrealistic iRT values (e.g., < -100 or > 200)

**Handling**:
```julia
function apply_irt_refinement!(...)
    # Clamp refined predictions to reasonable range
    irt_refined[i] = clamp(
        Float32(irt_pred - predicted_err),
        Float32(-50.0),  # Minimum reasonable iRT
        Float32(150.0)   # Maximum reasonable iRT
    )
end
```

### 5. Missing Amino Acids
**Scenario**: Non-standard amino acids in sequence (X, B, Z, J)

**Handling**:
```julia
function count_amino_acids(sequence::AbstractString)
    # Only count standard 20 AAs
    # Non-standard AAs contribute zero to all counts
    # Model will rely on irt_predicted for these cases
end
```

---

## Dependencies

### Existing Packages (already in Project.toml)
- `DataFrames`: Feature matrix manipulation
- `Arrow`: PSM file I/O
- `GLM`: Linear model fitting
- `StatsModels`: Formula interface
- `Random`: Train/test splitting
- `Statistics`: MAE, mean, etc.

### New Imports
Add to relevant files:
```julia
using GLM, StatsModels, Random
```

---

## Documentation

### User-Facing Documentation

**File**: `docs/src/user_guide/irt_refinement.md`

```markdown
# iRT Prediction Refinement

Pioneer can refine library iRT predictions using file-specific correction models
based on peptide sequence composition. This improves precursor candidate selection
in downstream search methods.

## Configuration

Enable in your parameters file:

```json
{
  "first_pass_search": {
    "enable_irt_refinement": true,
    "irt_refinement_min_psms_per_feature": 20,
    "irt_refinement_train_fraction": 0.667
  }
}
```

## Parameters

- **enable_irt_refinement**: Enable iRT refinement (default: true)
- **irt_refinement_min_psms_per_feature**: Minimum PSMs per feature for model training
  - Default: 20 (requires 420 total PSMs: 20 × 21 features)
  - Increase for stricter requirements, decrease for small datasets
- **irt_refinement_train_fraction**: Fraction of PSMs used for training (default: 0.667)
  - Remaining fraction used for validation

## Output

When enabled, FirstPass PSM files include an additional column:
- **irt_refined**: Refined iRT predictions (used by downstream methods)
- Falls back to **irt_predicted** if refinement fails or doesn't improve accuracy

## Diagnostics

Check console output for refinement statistics per file:
```
File 1 iRT Refinement: Training R²=0.4523, Validation R²=0.4312,
MAE Original=2.34, MAE Refined=2.18, Use=true
```
```

### Developer Documentation

**File**: `docs/src/advanced/irt_refinement_algorithm.md`

- Detailed algorithm description
- Mathematical formulation
- Feature engineering rationale
- Model selection criteria
- Validation strategy

---

## Future Enhancements (Not in Initial Implementation)

### 1. Modification Features
Add structural modification counts as features:
- `mod_M_Unimod_35` (Oxidation)
- `mod_C_Unimod_4` (Carbamidomethyl)
- Requires parsing structural_mods column

### 2. XGBoost Models
Replace linear regression with gradient boosting:
- Better capture non-linear relationships
- Requires `XGBoost.jl` dependency
- Higher computational cost

### 3. Reduced Amino Acid Alphabets
Use physicochemically similar AA groupings (e.g., Murphy et al., 2000):
- **Reduced-8**: Fewer features, more robust with limited data
- Reduces feature count from 21 to 9

### 4. Cross-File Information Sharing
Train global model on multiple files:
- Useful when individual files have insufficient PSMs
- Requires careful handling of batch effects

### 5. Ensemble Models
Combine multiple models:
- Linear + XGBoost weighted average
- Requires validation set performance tracking

---

## Success Criteria

### Minimum Viable Product (MVP)
- ✅ Linear model training on 20 AA + irt_predicted
- ✅ Per-file independent models
- ✅ Validation-based decision (use if MAE improves)
- ✅ `irt_refined` column in FirstPass PSMs
- ✅ Graceful fallback when insufficient data
- ✅ Minimal logging (R², MAE)

### Acceptance Tests
1. **Functional**: Run FirstPassSearch with refinement enabled, verify:
   - PSM files have `irt_refined` column
   - Models stored in SearchContext
   - No crashes on edge cases (few PSMs, failed models)

2. **Performance**: Compare with/without refinement:
   - Validation MAE improves by ≥5% on files where model is used
   - Total runtime increase <5% of FirstPassSearch time

3. **Robustness**: Test on diverse datasets:
   - Large libraries (100k+ precursors)
   - Small files (<500 PSMs)
   - Poor RT alignment (high MAD)

---

## Implementation Checklist

### Phase 1: Core Implementation (This PR)
- [ ] Create `irt_refinement_utils.jl` with core functions
- [ ] Add `IrtRefinementModel` struct
- [ ] Implement `count_amino_acids`
- [ ] Implement `prepare_irt_refinement_features`
- [ ] Implement `fit_irt_refinement_model` with:
  - [ ] Train/validation split
  - [ ] Validation-based decision
  - [ ] **Retrain on full dataset** if validation succeeds
- [ ] Implement `apply_irt_refinement!`
- [ ] Modify `FirstPassSearch/utils.jl` to:
  - [ ] Call refinement function
  - [ ] **Update SearchContext pred_irt** for ALL precursors
  - [ ] Add irt_refined column to PSM files
- [ ] Add SearchContext support (getter/setter methods)
- [ ] Add configuration parameters
- [ ] Write unit tests
- [ ] Write integration test
- [ ] Update documentation

### Phase 2: Validation & Benchmarking (Follow-up PR)
- [ ] End-to-end testing on diverse datasets
- [ ] Compare search results with/without refinement
- [ ] Performance profiling
- [ ] Verify downstream methods automatically use refined iRT

### Phase 3: Enhancements (Future PRs)
- [ ] Add modification features
- [ ] Implement XGBoost models
- [ ] Add reduced alphabet options
- [ ] Cross-file information sharing
- [ ] Ensemble models

---

## References

1. **Notebook Analysis**: `/Users/nathanwamsley/Documents/irtRefine/refineIrt.ipynb`
2. **RT Alignment Utils**: `src/Routines/SearchDIA/CommonSearchUtils/rt_alignment_utils.jl`
3. **FirstPassSearch**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/`
4. **GLM.jl Documentation**: https://juliastats.org/GLM.jl/stable/
5. **Murphy et al. (2000)**: Simplified amino acid alphabets for protein structure prediction

---

## Questions and Clarifications

### Resolved:
1. ✅ **File-specific vs global**: Per-file independent
2. ✅ **Storage location**: Update SearchContext pred_irt AND add `irt_refined` column to FirstPass PSMs
3. ✅ **Minimum PSMs**: 20 per feature (420 total), validate on holdout
4. ✅ **Diagnostics**: Minimal logging only
5. ✅ **Final model**: Retrain on full dataset (train + validation) after validation confirms improvement
6. ✅ **Downstream integration**: Automatic via SearchContext - no code changes needed

### Open:
- None currently

---

## Contact

For questions about this implementation:
- See notebook: `/Users/nathanwamsley/Documents/irtRefine/refineIrt.ipynb`
- Review existing RT alignment code: `rt_alignment_utils.jl`
- Check FirstPassSearch integration: `FirstPassSearch/utils.jl`
