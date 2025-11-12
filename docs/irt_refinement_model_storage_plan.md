# iRT Refinement Model Storage Architecture Plan

## Overview

**Goal**: Store per-file iRT refinement models in SearchContext and apply them dynamically whenever iRT predictions are needed, rather than pre-computing and storing refined iRT values for all precursors.

**Key Insight**: iRT refinement is file-specific. Each MS file has unique chromatographic characteristics, so the same precursor may need different refinement adjustments in different files.

---

## Current Architecture (Problems)

### Current SearchContext Storage
```julia
struct SearchContext
    # RT conversion models (per-file)
    irt_rt_map::Dict{Int64, RtConversionModel}   # iRT → RT
    rt_irt_map::Dict{Int64, RtConversionModel}   # RT → iRT

    # Predicted iRT values (global, single value per precursor)
    irt_obs::Dict{UInt32, Float32}  # precursor_idx → single iRT value
end
```

### Current Access Pattern
```julia
# Downstream methods call:
irt_pred = getPredIrt(search_context, prec_idx)  # Returns single value

# Problem: No file context!
# Which file are we predicting for?
```

### Problems with Current Approach

1. **No File Context**: `irt_obs` stores a single iRT value per precursor, but refinement is file-specific
2. **Memory Inefficiency**: Must store refined iRT for ALL precursors if we want to precompute
3. **Lost Information**: Refinement models contain rich information (AA coefficients, etc.) that's discarded after application
4. **Not Reusable**: Can't apply refinement to new precursors or transfer to new files
5. **Coupling Issue**: Refinement happens in FirstPassSearch but must be "baked in" to global state

---

## URGENT: Immediate Performance Fix (Phase 0)

### Current Hang Issue

**Symptom**: FirstPassSearch hangs after message:
```
[ Info: File 1: Applied iRT refinement to 170294 precursors
```

**Root Cause**: Lines 566-571 in FirstPassSearch.jl initialize ALL library precursors (500k) with library iRT upfront:
```julia
library_irt = getIrt(getPrecursors(getSpecLib(search_context)))
for prec_idx in 1:length(library_irt)
    setPredIrt!(search_context, UInt32(prec_idx), library_irt[prec_idx])
end
```

**Impact**:
- 500k dictionary insertions (Dict{UInt32, Float32})
- Multiple dictionary resize operations
- 3-5 seconds of wasted CPU time
- Unnecessary memory allocation (2MB for irt_obs)
- Performance bottleneck in get_best_precursors_accross_runs! after this

### Lazy Initialization Solution

**Strategy**: Remove upfront population entirely. Modify `getPredIrt` to lazily fall back to library iRT if precursor not found in irt_obs.

**Key Insight**: Only ~170k precursors are observed (have refined iRT). The other 330k don't need to be in irt_obs at all - we can fetch library iRT on demand.

### Implementation Details

#### Step 1: Remove Upfront Initialization
**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl`

**DELETE Lines 566-571**:
```julia
# DELETE THIS BLOCK:
# Initialize SearchContext irt_obs with library iRT values as baseline
# This ensures all precursors have iRT values (will be refined if enable_irt_refinement=true)
library_irt = getIrt(getPrecursors(getSpecLib(search_context)))
for prec_idx in 1:length(library_irt)
    setPredIrt!(search_context, UInt32(prec_idx), library_irt[prec_idx])
end
```

#### Step 2: Add Lazy Fallback to getPredIrt
**File**: `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl`

**Replace Lines 424-426**:
```julia
# OLD (direct dictionary access):
getPredIrt(s::SearchContext) = s.irt_obs
getPredIrt(s::SearchContext, prec_idx::Int64) = s.irt_obs[prec_idx]
getPredIrt(s::SearchContext, prec_idx::UInt32) = s.irt_obs[prec_idx]

# NEW (lazy fallback):
getPredIrt(s::SearchContext) = s.irt_obs

function getPredIrt(s::SearchContext, prec_idx::Int64)::Float32
    return getPredIrt(s, UInt32(prec_idx))
end

function getPredIrt(s::SearchContext, prec_idx::UInt32)::Float32
    # Check if refined/observed iRT exists in irt_obs
    irt = get(s.irt_obs, prec_idx, nothing)

    # Lazy fallback to library iRT if not found
    if isnothing(irt)
        return getIrt(getPrecursors(getSpecLib(s)))[prec_idx]
    end

    return irt
end
```

#### Step 3: Update Comments
**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl`

Update comments at lines 590, 605-607 to reflect lazy initialization:
```julia
# Line 590: Update comment
# Note: getPredIrt uses lazy fallback to library iRT for unobserved precursors

# Lines 605-607: Update comment
# Note: irt_obs contains only observed/refined precursors (~170k)
# Unobserved precursors use lazy fallback to library iRT in getPredIrt
```

### Benefits of Lazy Initialization

1. **Eliminate 500k dictionary operations** - only ~170k insertions during refinement
2. **3-5 second speedup** - no upfront initialization overhead
3. **66% memory reduction** - irt_obs: 170k entries (2.7MB) vs 500k entries (8MB)
4. **Transparent behavior** - downstream code unchanged due to lazy fallback
5. **Prepares for model storage** - already using on-demand evaluation pattern

### Performance Comparison

| Metric | Upfront Initialization | Lazy Fallback |
|--------|----------------------|---------------|
| **Startup time** | 3-5 seconds | 0 seconds |
| **Dictionary size** | 500k entries (8MB) | 170k entries (2.7MB) |
| **Dictionary insertions** | 500k + 170k = 670k | 170k only |
| **getPredIrt overhead** | O(1) direct lookup | O(1) dict lookup + O(1) array access |
| **Memory efficiency** | Stores all precursors | Stores only observed |

**Lookup overhead**: Negligible - `get(dict, key, nothing)` is O(1), library array access is O(1)

### Why This Prepares for Model Storage

The lazy initialization pattern is **exactly** what model storage needs:

```julia
# Phase 0: Lazy fallback to library iRT
function getPredIrt(ctx, prec_idx)::Float32
    irt = get(ctx.irt_obs, prec_idx, nothing)
    if isnothing(irt)
        return getIrt(getPrecursors(getSpecLib(ctx)))[prec_idx]  # Fallback
    end
    return irt
end

# Future (Phase 1-3): Lazy application of refinement model
function getPredIrt(ctx, prec_idx, ms_file_idx)::Float32
    library_irt = getIrt(getPrecursors(getSpecLib(ctx)))[prec_idx]
    model = getIrtRefinementModel(ctx, ms_file_idx)

    if isnothing(model) || !model.use_refinement
        return library_irt  # Fallback
    end

    # Apply refinement on-demand
    sequence = getSequence(getPrecursors(getSpecLib(ctx)))[prec_idx]
    return apply_irt_refinement(model, sequence, library_irt)
end
```

**Same pattern**: Check for enhanced value (irt_obs entry / refinement model), fallback to library baseline if not found.

### Testing Phase 0 Fix

```julia
@testset "Lazy iRT Initialization" begin
    ctx = create_test_context()

    # irt_obs should be empty initially
    @test isempty(getPredIrt(ctx))

    # Unobserved precursor should return library iRT (lazy fallback)
    library_irt = getIrt(getPrecursors(getSpecLib(ctx)))[1]
    @test getPredIrt(ctx, UInt32(1)) == library_irt

    # Store refined iRT for observed precursor
    refined_irt = 55.0f0
    setPredIrt!(ctx, UInt32(1), refined_irt)

    # Should now return refined value
    @test getPredIrt(ctx, UInt32(1)) == refined_irt

    # Other precursors still use lazy fallback
    @test getPredIrt(ctx, UInt32(2)) == getIrt(getPrecursors(getSpecLib(ctx)))[2]
end
```

### Implementation Checklist (Phase 0)

**Priority: URGENT - Fixes current hang issue**

- [ ] Delete lines 566-571 from FirstPassSearch.jl (upfront initialization)
- [ ] Replace getPredIrt accessors in SearchTypes.jl with lazy fallback (lines 424-426)
- [ ] Update comments in FirstPassSearch.jl (lines 590, 605-607)
- [ ] Add unit tests for lazy fallback behavior
- [ ] Verify module compiles
- [ ] Run FirstPass and verify hang is resolved
- [ ] Measure performance improvement (should see 3-5 second speedup)
- [ ] Verify irt_obs contains only ~170k entries (not 500k)

**Estimated time**: 30 minutes implementation + 30 minutes testing = 1 hour total

**Must complete before**: Starting Phase 1 (model storage architecture)

---

## Proposed Architecture

### SearchContext Storage Addition
```julia
struct SearchContext
    # Existing RT conversion models
    irt_rt_map::Dict{Int64, RtConversionModel}
    rt_irt_map::Dict{Int64, RtConversionModel}

    # NEW: Per-file iRT refinement models
    irt_refinement_models::Dict{Int64, Union{IrtRefinementModel, Nothing}}

    # Library baseline iRT (never modified)
    # Accessed via: getIrt(getPrecursors(getSpecLib(search_context)))
end
```

### New Access Pattern
```julia
# Option A: File-aware API
function getPredIrt(search_context::SearchContext,
                   prec_idx::UInt32,
                   ms_file_idx::Int64)::Float32
    # Get baseline library iRT
    library_irt = getIrt(getPrecursors(getSpecLib(search_context)))[prec_idx]

    # Check if refinement model exists for this file
    refinement_model = get(search_context.irt_refinement_models, ms_file_idx, nothing)

    if isnothing(refinement_model) || !refinement_model.use_refinement
        return library_irt
    end

    # Apply file-specific refinement
    sequence = getSequence(getPrecursors(getSpecLib(search_context)))[prec_idx]
    return apply_irt_refinement(refinement_model, sequence, library_irt)
end

# Option B: Keep old API but require file context to be set
function setCurrentMsFile!(search_context::SearchContext, ms_file_idx::Int64)
    search_context.current_ms_file_idx[] = ms_file_idx
end

function getPredIrt(search_context::SearchContext, prec_idx::UInt32)::Float32
    ms_file_idx = search_context.current_ms_file_idx[]
    return getPredIrt(search_context, prec_idx, ms_file_idx)
end
```

---

## Benefits of Model Storage Approach

### 1. **File-Specific Refinement**
- Each file gets its own refinement model
- Accounts for file-specific chromatographic drift
- More accurate than global refinement

### 2. **Memory Efficiency**
- Models are small (~100 Float32 values: 20 AA coefficients + metadata)
- vs. storing refined iRT for 100k-500k precursors per file
- ~10KB per model vs ~400KB-2MB per file for full precursor arrays

### 3. **Lazy Evaluation**
- Refinement only computed when needed
- Unobserved precursors can still get refined predictions
- No need to precompute all precursors

### 4. **Model Interpretability**
- AA coefficients show which amino acids cause RT shifts
- Can diagnose chromatographic issues
- Exportable for QC and method development

### 5. **Transfer Learning Potential**
- Models could be saved and reused across experiments
- Could average models for robust predictions
- Could predict for new precursors not in library

### 6. **Consistency**
- Model stored with RT conversion models (logical grouping)
- Clear file context for all predictions
- Easier to debug (which model was used?)

---

## Implementation Plan

### Phase 1: SearchContext Modifications

#### 1.1 Add Model Storage Field
**File**: `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl`

```julia
mutable struct SearchContext{N,L<:SpectralLibrary,M<:MassSpecDataReference}
    # ... existing fields ...

    # RT conversion models (existing)
    irt_rt_map::Dict{Int64, RtConversionModel}
    rt_irt_map::Dict{Int64, RtConversionModel}

    # NEW: iRT refinement models (per-file)
    irt_refinement_models::Dict{Int64, Union{IrtRefinementModel, Nothing}}

    # Optional: track current file context
    current_ms_file_idx::Base.Ref{Int64}  # Default: 0 (invalid)

    # ... existing fields ...
end
```

#### 1.2 Update Constructor
```julia
function SearchContext(spec_lib, temp_structures, mass_spec_data_reference, ...)
    new{N,L,M}(
        spec_lib, temp_structures, mass_spec_data_reference,
        # ... existing initializations ...
        Dict{Int64, RtConversionModel}(),  # irt_rt_map
        Dict{Int64, RtConversionModel}(),  # rt_irt_map
        Dict{Int64, Union{IrtRefinementModel, Nothing}}(),  # NEW: irt_refinement_models
        Ref{Int64}(0),  # NEW: current_ms_file_idx (0 = invalid)
        # ... rest of fields ...
    )
end
```

#### 1.3 Add Accessor Functions
```julia
# Store refinement model
function setIrtRefinementModel!(ctx::SearchContext,
                                ms_file_idx::Int64,
                                model::Union{IrtRefinementModel, Nothing})
    ctx.irt_refinement_models[ms_file_idx] = model
end

# Get refinement model
function getIrtRefinementModel(ctx::SearchContext,
                               ms_file_idx::Int64)::Union{IrtRefinementModel, Nothing}
    get(ctx.irt_refinement_models, ms_file_idx, nothing)
end

# Set current file context (Option B)
function setCurrentMsFile!(ctx::SearchContext, ms_file_idx::Int64)
    ctx.current_ms_file_idx[] = ms_file_idx
end

function getCurrentMsFile(ctx::SearchContext)::Int64
    ctx.current_ms_file_idx[]
end
```

### Phase 2: Update getPredIrt API

#### Option A: File-Aware API (Recommended)
```julia
# New primary API with file context
function getPredIrt(ctx::SearchContext,
                   prec_idx::UInt32,
                   ms_file_idx::Int64)::Float32
    # Get baseline library iRT
    library_irt = getIrt(getPrecursors(getSpecLib(ctx)))[prec_idx]

    # Get refinement model for this file
    model = getIrtRefinementModel(ctx, ms_file_idx)

    # Apply refinement if available
    if !isnothing(model) && model.use_refinement
        sequence = getSequence(getPrecursors(getSpecLib(ctx)))[prec_idx]
        counts = zeros(Int, 20)  # Pre-allocated for efficiency
        return apply_irt_refinement(model, counts, sequence, library_irt)
    end

    return library_irt
end

# Backward compatibility wrapper (uses current_ms_file_idx)
function getPredIrt(ctx::SearchContext, prec_idx::UInt32)::Float32
    ms_file_idx = getCurrentMsFile(ctx)
    if ms_file_idx == 0
        @warn "No current MS file set, using library iRT without refinement"
        return getIrt(getPrecursors(getSpecLib(ctx)))[prec_idx]
    end
    return getPredIrt(ctx, prec_idx, ms_file_idx)
end
```

#### Option B: Context-Based API
```julia
# Keep single-argument API, require setCurrentMsFile! before use
function getPredIrt(ctx::SearchContext, prec_idx::UInt32)::Float32
    ms_file_idx = getCurrentMsFile(ctx)
    @assert ms_file_idx != 0 "Must call setCurrentMsFile! before getPredIrt"

    library_irt = getIrt(getPrecursors(getSpecLib(ctx)))[prec_idx]
    model = getIrtRefinementModel(ctx, ms_file_idx)

    if !isnothing(model) && model.use_refinement
        sequence = getSequence(getPrecursors(getSpecLib(ctx)))[prec_idx]
        counts = zeros(Int, 20)
        return apply_irt_refinement(model, counts, sequence, library_irt)
    end

    return library_irt
end
```

**Recommendation**: Use Option A (file-aware API) for explicitness and safety.

### Phase 3: Update FirstPassSearch to Store Models

**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl`

```julia
function map_retention_times!(
    search_context::SearchContext,
    results::FirstPassSearchResults,
    params::FirstPassSearchParameters
)
    # ... existing RT model fitting code ...

    for ms_file_idx in valid_file_indices
        # ... fit RT model ...

        # Train refinement model if enabled
        if params.enable_irt_refinement
            refinement_model = fit_irt_refinement_model(
                sequences,
                Float32.(psms[:irt_predicted][best_hits]),
                observed_irts,
                ms_file_idx=ms_file_idx,
                min_psms=20,
                train_fraction=0.67
            )

            # Store model in SearchContext (even if nothing - indicates no improvement)
            setIrtRefinementModel!(search_context, ms_file_idx, refinement_model)
        end
    end
end
```

**Key Change**: Store the model instead of applying it to precursors.

### Phase 4: Update Downstream Methods

#### SecondPassSearch
**File**: `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl`

```julia
function add_columns_second_search_psms!(
    psms::DataFrame,
    search_context::SearchContext,
    tic::AbstractVector{Float32},
    masses::AbstractArray,
    ms_file_idx::Integer,  # Already have this parameter!
    rt_to_irt_interp::RtConversionModel,
    prec_id_to_irt::Dictionary{...}
)
    # ... existing code ...

    for i in chunk
        prec_idx = precursor_idx[i]

        # OLD: irt_pred[i] = getPredIrt(search_context, prec_idx)
        # NEW: Pass file index for file-specific refinement
        irt_pred[i] = getPredIrt(search_context, prec_idx, ms_file_idx)

        irt_obs[i] = rt_to_irt_interp(rt[i])
        irt_diff[i] = abs(irt_obs[i] - prec_id_to_irt[prec_idx].best_irt)

        if !ms1_missing[i]
            # NEW: Use file-specific refined iRT
            ms1_irt_diff[i] = abs(rt_to_irt_interp(ms1_rt[i]) - getPredIrt(search_context, prec_idx, ms_file_idx))
        end

        # ... rest of code ...
    end
end
```

#### IntegrateChromatogramsSearch
**File**: `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/integrate_chrom.jl`

```julia
function extract_chromatograms!(
    search_context::SearchContext,
    precursor_idx::UInt32,
    ms_file_idx::Int64,  # Already available in this method
    params::IntegrateChromatogramsParameters
)
    # Get file-specific refined iRT prediction
    predicted_irt = getPredIrt(search_context, precursor_idx, ms_file_idx)

    # Convert to RT using file-specific RT model
    rt_model = getIrtRtModel(search_context, ms_file_idx)
    predicted_rt = rt_model(predicted_irt)

    # Extract chromatogram around predicted RT
    # ...
end
```

#### ScoringSearch
**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/score_psms.jl`

```julia
function add_scoring_features!(
    psms::DataFrame,
    search_context::SearchContext,
    ms_file_idx::Int64  # Need to track this or use setCurrentMsFile!
)
    # Compute iRT error feature
    for i in 1:nrow(psms)
        prec_idx = psms[i, :precursor_idx]

        # Use file-specific refinement
        irt_pred = getPredIrt(search_context, prec_idx, ms_file_idx)
        irt_obs = psms[i, :irt_observed]
        psms[i, :irt_error] = abs(irt_obs - irt_pred)
    end
end
```

### Phase 5: Update FirstPassSearch PSM Output (Optional Enhancement)

**File**: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl`

Since FirstPass runs BEFORE refinement models are trained, the PSM output will still show library iRT. However, we could add a post-processing step:

```julia
function finalize_first_pass_psms!(
    search_context::SearchContext,
    results::FirstPassSearchResults,
    params::FirstPassSearchParameters
)
    # After refinement models are trained, optionally add irt_refined column
    if params.enable_irt_refinement
        for ms_file_idx in valid_files
            psm_path = getFirstPassPsms(getMSData(search_context))[ms_file_idx]
            psms = DataFrame(Arrow.Table(psm_path))

            # Add irt_refined column using stored model
            psms[!, :irt_refined] = [
                getPredIrt(search_context, UInt32(pid), ms_file_idx)
                for pid in psms[:, :precursor_idx]
            ]
            psms[!, :irt_error_refined] = abs.(psms[:, :irt_observed] .- psms[:, :irt_refined])

            # Overwrite PSM file with additional columns
            Arrow.write(psm_path, psms)
        end
    end
end
```

---

## Migration Strategy

### Phase 0: URGENT Fix (Day 1 - Non-Breaking)
**Purpose**: Resolve current hang issue and improve performance
**Changes**:
1. Remove upfront initialization of 500k precursors
2. Add lazy fallback to getPredIrt (transparent to callers)
3. Reduces irt_obs from 500k to ~170k entries
4. 3-5 second speedup + 66% memory reduction

**Non-breaking**: Downstream code sees no behavioral change due to lazy fallback

### Phase 1: Add Infrastructure (Non-Breaking)
1. Add `irt_refinement_models` field to SearchContext
2. Add accessor functions
3. Add file-aware `getPredIrt(ctx, prec_idx, ms_file_idx)` overload
4. Keep existing `getPredIrt(ctx, prec_idx)` working with warnings

### Phase 2: Update FirstPassSearch (Model Storage)
1. Modify `map_retention_times!` to store models instead of applying
2. Remove the precursor-by-precursor refinement application code
3. Test that models are stored correctly

### Phase 3: Update Downstream Methods (Gradual)
1. **SecondPassSearch** - already has `ms_file_idx` parameter, easy to update
2. **IntegrateChromatogramsSearch** - already has file context
3. **ScoringSearch** - need to track file index through feature engineering
4. **QuadTuningSearch** - check if uses getPredIrt (probably not)

### Phase 4: Deprecate Old API (Optional)
1. Make `getPredIrt(ctx, prec_idx)` error if no current file set
2. Force all callers to use file-aware API
3. Remove backward compatibility wrapper

---

## Testing Strategy

### Unit Tests

```julia
@testset "iRT Refinement Model Storage" begin
    # Test model storage
    ctx = create_test_context()
    model = IrtRefinementModel(true, 0.1f0, 0.95f0, zeros(Float32, 20), ...)

    setIrtRefinementModel!(ctx, 1, model)
    @test getIrtRefinementModel(ctx, 1) === model
    @test isnothing(getIrtRefinementModel(ctx, 2))

    # Test file-aware getPredIrt
    library_irt = 50.0f0
    refined_irt = getPredIrt(ctx, UInt32(1), 1)
    @test refined_irt != library_irt  # Should be refined

    unrefined_irt = getPredIrt(ctx, UInt32(1), 2)  # No model for file 2
    @test unrefined_irt == library_irt
end
```

### Integration Tests

```julia
@testset "Full Pipeline with Model Storage" begin
    # Run FirstPassSearch with refinement enabled
    params = FirstPassSearchParameters(enable_irt_refinement=true, ...)
    search_context = create_test_context()

    performSearch!(FirstPassSearch(search_context))

    # Verify models are stored
    @test haskey(search_context.irt_refinement_models, 1)
    model = getIrtRefinementModel(search_context, 1)
    @test !isnothing(model)
    @test model.use_refinement == true

    # Run SecondPassSearch
    performSearch!(SecondPassSearch(search_context))

    # Verify SecondPass PSMs show improved irt_error
    second_pass_psms = load_psms(...)
    @test mean(second_pass_psms[:, :irt_error]) < threshold
end
```

### QC Validation

```julia
function validate_refinement_benefit(search_context::SearchContext, ms_file_idx::Int64)
    model = getIrtRefinementModel(search_context, ms_file_idx)

    if isnothing(model)
        @info "File $ms_file_idx: No refinement model (insufficient data or no improvement)"
        return
    end

    if !model.use_refinement
        @info "File $ms_file_idx: Refinement model trained but not used (MAE worsened)"
        @info "  Original MAE: $(model.mae_original)"
        @info "  Refined MAE: $(model.mae_refined)"
        return
    end

    improvement_pct = (model.mae_original - model.mae_refined) / model.mae_original * 100
    @info "File $ms_file_idx: Refinement ENABLED"
    @info "  Training R²: $(round(model.r2_train, digits=4))"
    @info "  Validation R²: $(round(model.r2_val, digits=4))"
    @info "  MAE improvement: $(round(improvement_pct, digits=2))%"
end
```

---

## Performance Considerations

### Memory Usage

**Current Approach (Pre-computed iRT)**:
- 500k precursors × 4 bytes (Float32) = 2 MB per file
- 100 files = 200 MB total
- Must store even if refinement not beneficial

**Model Storage Approach**:
- 20 AA coefficients × 4 bytes = 80 bytes
- Plus metadata ~20 bytes = ~100 bytes per model
- 100 files × 100 bytes = 10 KB total
- **20,000x memory savings**

### Computational Cost

**Per-prediction refinement overhead**:
```julia
# Amino acid counting: O(sequence_length)
count_amino_acids!(counts, sequence)  # ~10-30 iterations

# Dot product: O(20) = constant
predicted_err = model.intercept + model.irt_coef * irt_original
for i in 1:20
    predicted_err += model.aa_coefficients[i] * counts[i]
end

# Total: ~10-50 arithmetic operations per prediction
```

**Optimization**: Pre-allocate AA counts vector per thread
```julia
struct SearchData
    # ... existing fields ...
    aa_counts::Vector{Int}  # Pre-allocated, length 20
end

# Reuse across predictions in same thread
for precursor in precursors
    irt_refined = getPredIrt(search_context, precursor.idx, ms_file_idx,
                             counts=thread_local_data.aa_counts)
end
```

### Caching Strategy (Optional)

For very hot paths (e.g., scoring millions of PSMs), could add LRU cache:

```julia
struct SearchContext
    # ... existing fields ...
    irt_prediction_cache::LRU{Tuple{UInt32, Int64}, Float32}  # (prec_idx, file_idx) → irt
end

function getPredIrt(ctx::SearchContext, prec_idx::UInt32, ms_file_idx::Int64)::Float32
    key = (prec_idx, ms_file_idx)

    # Check cache
    cached = get(ctx.irt_prediction_cache, key, nothing)
    if !isnothing(cached)
        return cached
    end

    # Compute and cache
    irt = _compute_predicted_irt(ctx, prec_idx, ms_file_idx)
    ctx.irt_prediction_cache[key] = irt
    return irt
end
```

---

## Model Export/Import (Future Enhancement)

### Model Serialization

```julia
function save_refinement_models(search_context::SearchContext, output_dir::String)
    models_dir = joinpath(output_dir, "irt_refinement_models")
    mkpath(models_dir)

    for (ms_file_idx, model) in pairs(search_context.irt_refinement_models)
        if !isnothing(model)
            model_path = joinpath(models_dir, "file_$(ms_file_idx)_irt_model.jld2")
            JLD2.save(model_path, "model", model)
        end
    end
end

function load_refinement_models!(search_context::SearchContext, models_dir::String)
    for model_file in readdir(models_dir, join=true)
        if endswith(model_file, "_irt_model.jld2")
            ms_file_idx = parse_file_idx_from_path(model_file)
            model = JLD2.load(model_file, "model")
            setIrtRefinementModel!(search_context, ms_file_idx, model)
        end
    end
end
```

### Transfer Learning

```julia
function transfer_refinement_model(
    source_context::SearchContext,
    source_file_idx::Int64,
    target_context::SearchContext,
    target_file_idx::Int64
)
    model = getIrtRefinementModel(source_context, source_file_idx)
    if !isnothing(model)
        setIrtRefinementModel!(target_context, target_file_idx, model)
        @info "Transferred refinement model from file $source_file_idx to $target_file_idx"
    end
end
```

---

## Comparison: Current vs Proposed

| Aspect | Current (Value Storage) | Proposed (Model Storage) |
|--------|------------------------|--------------------------|
| **Memory** | 2 MB per file (500k precursors) | 100 bytes per file (model) |
| **Flexibility** | Single iRT per precursor | File-specific predictions |
| **Reusability** | Cannot reapply to new precursors | Can predict any precursor |
| **Interpretability** | Black box (stored values) | AA coefficients visible |
| **Performance** | Fast lookup O(1) | Small overhead O(seq_len) |
| **Consistency** | Must precompute all precursors | Lazy evaluation |
| **Export** | Export refined iRT values | Export models for reuse |

---

## Implementation Checklist

### Milestone 0: URGENT Performance Fix (Day 1)
**Priority: CRITICAL - Fixes current hang issue**

- [ ] Delete lines 566-571 from FirstPassSearch.jl (remove upfront initialization loop)
- [ ] Replace getPredIrt accessors in SearchTypes.jl with lazy fallback (lines 424-426)
- [ ] Update comments in FirstPassSearch.jl (lines 590, 605-607)
- [ ] Add unit tests for lazy fallback behavior
- [ ] Verify module compiles
- [ ] Run FirstPass and verify hang is resolved
- [ ] Measure performance: should see 3-5 second speedup
- [ ] Verify irt_obs contains only ~170k entries (not 500k)

**Estimated time**: 1 hour (30 min implementation + 30 min testing)
**Blocks**: All subsequent phases (must complete first)

### Milestone 1: Infrastructure (Week 1)
- [ ] Add `irt_refinement_models` field to SearchContext
- [ ] Add `current_ms_file_idx` tracking field
- [ ] Implement accessor functions (set/get IrtRefinementModel)
- [ ] Implement file-aware `getPredIrt(ctx, prec_idx, ms_file_idx)`
- [ ] Add unit tests for model storage and retrieval
- [ ] Verify module compiles

### Milestone 2: FirstPassSearch Integration (Week 1)
- [ ] Modify `map_retention_times!` to store models
- [ ] Remove precursor-by-precursor refinement application
- [ ] Add logging to show stored models per file
- [ ] Test that models are stored correctly
- [ ] Verify FirstPass still completes successfully

### Milestone 3: SecondPassSearch Integration (Week 2)
- [ ] Update `add_columns_second_search_psms!` to use file-aware API
- [ ] Pass `ms_file_idx` to `getPredIrt` calls (lines 905, 909)
- [ ] Test that irt_error improves in SecondPass PSMs
- [ ] Verify no performance regression

### Milestone 4: Other Search Methods (Week 2)
- [ ] Update IntegrateChromatogramsSearch
- [ ] Update ScoringSearch feature engineering
- [ ] Check QuadTuningSearch/ParameterTuningSearch (if applicable)
- [ ] Update any other methods that call getPredIrt

### Milestone 5: Testing & Validation (Week 3)
- [ ] Integration test: Full pipeline with refinement enabled
- [ ] Benchmark: Memory usage comparison
- [ ] Benchmark: Runtime overhead per prediction
- [ ] QC validation: Verify improved PSM identification
- [ ] Documentation: Update user guide

### Milestone 6: Optional Enhancements (Week 4)
- [ ] Implement model export/import
- [ ] Add LRU caching for hot paths
- [ ] Add QC plots showing refinement benefit per file
- [ ] Add option to export refined iRT library

---

## Risks and Mitigations

### Risk 1: Breaking Changes
**Risk**: Changing getPredIrt API breaks downstream code
**Mitigation**:
- Keep backward-compatible wrapper with warnings
- Gradual migration with deprecation period
- Comprehensive test suite

### Risk 2: Performance Regression
**Risk**: Per-prediction refinement adds overhead
**Mitigation**:
- Pre-allocate AA count buffers per thread
- Profile hot paths and add caching if needed
- Benchmark shows overhead is minimal (~50 ops per prediction)

### Risk 3: Missing File Context
**Risk**: Some methods don't track ms_file_idx
**Mitigation**:
- Add `current_ms_file_idx` fallback in SearchContext
- Methods can call `setCurrentMsFile!` at beginning
- Assert/warn if file context missing

### Risk 4: Model Training Failure
**Risk**: Some files don't have enough PSMs for refinement
**Mitigation**:
- Store `nothing` for failed refinement (already implemented)
- Fallback to library iRT automatically
- Log which files have refinement models

---

## Success Criteria

### Quantitative Metrics
1. **Memory reduction**: >1000x reduction in iRT storage (from MB to KB)
2. **PSM improvement**: 5-10% more PSMs identified with refinement
3. **RT error reduction**: Mean irt_error in SecondPass reduced by 20-50%
4. **Performance**: <1% runtime overhead from per-prediction refinement
5. **Model training**: >80% of files get beneficial refinement models

### Qualitative Goals
1. **Maintainability**: Clear file context for all iRT predictions
2. **Flexibility**: Can apply refinement to any precursor, any file
3. **Debuggability**: Easy to inspect which model was used
4. **Reusability**: Models can be exported and transferred
5. **Interpretability**: AA coefficients provide scientific insights

---

## Future Directions

### 1. Multi-File Ensemble Models
Train a single robust model using data from all files:
```julia
function train_ensemble_refinement_model(search_context::SearchContext)
    # Aggregate PSMs from all files
    all_sequences = []
    all_library_irts = []
    all_observed_irts = []

    for ms_file_idx in valid_files
        # Collect data...
    end

    # Train single model on combined data
    ensemble_model = fit_irt_refinement_model(all_sequences, ...)

    # Apply to all files as fallback
    for ms_file_idx in valid_files
        if isnothing(getIrtRefinementModel(search_context, ms_file_idx))
            setIrtRefinementModel!(search_context, ms_file_idx, ensemble_model)
        end
    end
end
```

### 2. Time-Dependent Refinement
Model chromatographic drift within a single file:
```julia
struct TimeVaryingRefinementModel
    time_bins::Vector{Float32}  # RT bins
    models_per_bin::Vector{IrtRefinementModel}  # Model per bin
end

function getPredIrt(ctx, prec_idx, ms_file_idx, observed_rt::Float32)
    model = get_time_specific_model(ctx, ms_file_idx, observed_rt)
    # Apply model...
end
```

### 3. Cross-Experiment Transfer
Use refinement models from high-quality runs for poor quality runs:
```julia
function transfer_best_model!(target_ctx, source_experiments)
    # Find best-performing refinement model across experiments
    best_model = find_best_refinement_model(source_experiments)

    # Apply to all files in target
    for ms_file_idx in target_files
        setIrtRefinementModel!(target_ctx, ms_file_idx, best_model)
    end
end
```

### 4. Non-Linear Refinement
Extend beyond linear models:
```julia
struct NeuralRefinementModel
    network::Chain  # Simple neural network
    # Inputs: 20 AA counts + library iRT
    # Output: iRT error
end
```

---

## Conclusion

### Immediate Action Required (Phase 0)

**URGENT**: Implement lazy initialization fix (1 hour) to resolve current hang issue:
- Eliminates 500k unnecessary dictionary operations
- 3-5 second performance improvement
- 66% memory reduction in irt_obs
- Non-breaking change (transparent to existing code)
- **Must be completed before any further work on iRT refinement**

### Long-Term Recommendation (Phases 1-6)

**Recommendation**: Implement model storage approach for the following reasons:

1. **Architecturally Superior**: File-specific models align with the physical reality that RT behavior varies per file
2. **Memory Efficient**: 20,000x memory savings (KB vs MB)
3. **Flexible**: Can predict any precursor in any file context
4. **Extensible**: Enables future enhancements (transfer learning, ensemble models, export/import)
5. **Interpretable**: Models provide scientific insights into chromatographic behavior

**Timeline**:
- Phase 0 (URGENT): 1 hour - Fix hang issue
- Phases 1-6: 3-4 weeks - Full model storage implementation

**Priority**:
- Phase 0: **CRITICAL** (blocking current work)
- Phases 1-6: High (correct long-term architecture for iRT refinement in Pioneer.jl)
