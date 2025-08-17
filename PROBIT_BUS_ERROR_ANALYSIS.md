# Probit Regression Bus Error Analysis

## Problem Summary
When using probit regression instead of XGBoost for PSM scoring (<100k PSMs), a bus error occurs during Step 2 when writing merged DataFrames back to individual Arrow files. The XGBoost path works fine, but the probit path consistently crashes.

## Error Location
```
[7740] signal 10 (1): Bus error: 10
in expression starting at REPL[3]:1
getindex at ./essentials.jl:917 [inlined]
getindex at Arrow/owCB0/src/arraytypes/primitive.jl:48
```

The crash occurs in ScoringSearch.jl at line 281:
```julia
for (idx, ref) in enumerate(second_pass_refs)
    sub_df = merged_df[merged_df.ms_file_idx .== idx, :]
    write_arrow_file(ref, sub_df)  # <-- CRASH HERE
end
```

## Key Observations

### What Works
- Probit regression completes successfully
- Probabilities are properly calculated (range: 0.016 to 0.9999999)
- DataFrame has correct dimensions after probit (95 columns)
- XGBoost path with same data works perfectly

### What's Different with Probit
1. **Columns removed**: `intercept` and `cv_fold` (correctly removed)
2. **Column added**: `best_psm` (added to match XGBoost)
3. **No models returned**: Probit returns `nothing` while XGBoost returns models
4. **In-memory operation**: All PSMs stay in memory during probit

## Possible Root Causes

### 1. **Memory Corruption from In-Place Modifications**
**Hypothesis**: Probit modifies the DataFrame in-place, potentially corrupting memory references that Arrow needs.

**Evidence**:
- Bus error (not segfault) suggests misaligned memory access
- Occurs when Arrow tries to access array data
- Only happens with probit path

**Test**: 
- Create a deep copy of the DataFrame before probit regression
- Or reload PSMs from files after scoring

### 2. **Column Type Mismatch**
**Hypothesis**: A column has an incompatible type after probit that Arrow can't serialize.

**Potential problematic columns**:
- `accession_numbers` (Vector{String} - complex type)
- `best_psm` (Bool column we added manually)
- Any column that got modified during probit

**Test**:
- Remove `accession_numbers` column before writing
- Check all column types match between XGBoost and probit paths

### 3. **DataFrame Metadata Corruption**
**Hypothesis**: The DataFrame's internal metadata is corrupted after multiple in-place operations.

**Operations that modify DataFrame**:
1. Adding `q_value`, `decoy` columns
2. Probit adds `intercept`, `cv_fold` 
3. We remove `intercept`, `cv_fold`
4. We add `best_psm`
5. Multiple `select!` operations

**Test**:
- Reconstruct DataFrame from scratch after probit
- Use non-mutating operations (`select` instead of `select!`)

### 4. **Reference vs Value Issue**
**Hypothesis**: The DataFrame references data that's been garbage collected or moved.

**Evidence**:
- `best_psms = nothing; GC.gc()` happens after scoring
- But the actual data is in files that get re-read
- The merge and transform operations might create views instead of copies

**Test**:
- Force materialization of all columns: `df = DataFrame(df)`
- Avoid column views and references

### 5. **Threading/Concurrency Issue**
**Hypothesis**: Probit's threaded operations leave the DataFrame in an inconsistent state.

**Evidence**:
- Probit uses `Threads.@spawn` for parallel processing
- Arrow writing also uses threading
- Race condition or incomplete writes

**Test**:
- Run with single thread
- Add synchronization barriers

## Immediate Workarounds

### Option 1: Skip In-Memory Modification
Instead of modifying `best_psms` in place, write it to disk and reload:
```julia
# After probit scoring
temp_path = "temp_scored.arrow"
Arrow.write(temp_path, best_psms)
best_psms = DataFrame(Arrow.Table(temp_path))
```

### Option 2: Force Copy After Scoring
```julia
# After probit_regression_scoring_cv!
best_psms = DataFrame(best_psms)  # Force full copy
```

### Option 3: Remove Complex Columns
```julia
# Before writing to Arrow
if :accession_numbers in names(merged_df)
    select!(merged_df, Not(:accession_numbers))
end
```

### Option 4: Use XGBoost Path Structure
Make probit follow the exact same data flow as XGBoost, including writing intermediate files.

## Recommended Debugging Steps

1. **Add type checking before Arrow write**:
```julia
for col in names(sub_df)
    @info "Column $col: type=$(eltype(sub_df[!, col]))"
end
```

2. **Test with minimal DataFrame**:
```julia
# Keep only essential columns
essential_cols = [:prob, :target, :precursor_idx, :ms_file_idx]
sub_df_minimal = sub_df[:, essential_cols]
write_arrow_file(ref, sub_df_minimal)
```

3. **Compare memory addresses**:
```julia
@info "DataFrame pointer: $(pointer_from_objref(sub_df))"
@info "First column pointer: $(pointer(sub_df[!, 1]))"
```

4. **Force garbage collection before write**:
```julia
GC.gc()
sleep(0.1)  # Let GC finish
write_arrow_file(ref, sub_df)
```

## Most Likely Solution

Based on the evidence, the most likely issue is **memory corruption from in-place modifications**. The DataFrame goes through many mutations during probit processing, and Arrow's C-based serialization is sensitive to memory layout.

**Recommended fix**:
```julia
# In score_psms.jl, after probit_regression_scoring_cv! returns:
best_psms = DataFrame(best_psms)  # Force full copy
# Or even safer:
best_psms = copy(best_psms)  # Deep copy
```

This ensures Arrow gets a clean, properly aligned DataFrame to serialize.

## Alternative: Match XGBoost Exactly

Since XGBoost works, we could make probit follow the exact same pattern:
1. Don't modify `best_psms` in place
2. Write scores to separate files
3. Let the downstream merge handle everything
4. Avoid adding/removing columns during scoring