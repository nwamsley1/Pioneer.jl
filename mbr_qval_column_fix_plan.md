# MBR Q-value Column Fix - Detailed Implementation Plan

## Problem Statement

After implementing the MBR column naming fix, MaxLFQSearch and other downstream code fails because:
- When MBR is ON: We create `MBR_boosted_global_qval` and `MBR_boosted_qval` but downstream code expects `global_qval` and `qval`
- When MBR is OFF: We create `global_qval` and `qval` (correct)
- Step 10 recalculates `qval` but not `global_qval` when MBR is enabled

The error occurs in `MaxLFQSearch.jl:202` when trying to write CSV output that expects `global_qval` column.

## Desired Behavior

### When MBR is ENABLED:
After all steps complete, the final output should have:
- `MBR_boosted_global_qval`: Q-values calculated from MBR-boosted global_prob (from Step 9)
- `MBR_boosted_qval`: Q-values calculated from MBR-boosted prec_prob (from Step 9, overwritten in Step 10)
- **NO** `global_qval` column
- **NO** `qval` column

Note: Step 10 recalculates `MBR_boosted_qval` using non-MBR prec_prob scores for unbiased protein inference.

### When MBR is DISABLED:
After all steps complete, the final output should have:
- `global_qval`: Q-values calculated from global_prob (from Step 9)
- `qval`: Q-values calculated from prec_prob (from Step 9)
- **NO** `MBR_boosted_*` columns

## Root Cause Analysis

The issue is that downstream code (MaxLFQSearch, writeCSVTables, etc.) expects certain column names to always exist:
- `global_qval`
- `qval`

These are hardcoded in various places and need to be made conditional based on `match_between_runs`.

---

## Implementation Plan

### PHASE 1: Identify All References to Q-value Columns

**Search for hardcoded references to:**
- `:global_qval`
- `:qval`
- `:pep`

**Files to check:**
1. `MaxLFQSearch/MaxLFQSearch.jl`
2. `WriteOutputs/writeCSVTables.jl`
3. `WriteOutputs/writePlots.jl`
4. Any other files that read the output of ScoringSearch

**Actions:**
- Grep for `global_qval` in the codebase
- Grep for `:qval` (with colon to avoid matching variable names)
- Document each location and its context

---

### PHASE 2: Update MaxLFQSearch to Use Conditional Column Names

**File**: `src/Routines/SearchDIA/SearchMethods/MaxLFQSearch/MaxLFQSearch.jl`

#### Change 2.1: Update Column References in summarize_results!

**Current pattern** (example from line 202):
```julia
# Likely uses :global_qval or :qval directly
writePrecursorCSV(path, file_names, normalized, proteins)
```

**New pattern**:
```julia
# Pass match_between_runs flag to CSV writer
writePrecursorCSV(path, file_names, normalized, proteins, params.match_between_runs)
```

#### Change 2.2: Update Any DataFrame Operations

Look for operations like:
```julia
df.global_qval
df[!, :qval]
filter(row -> row.global_qval < threshold, df)
```

Replace with conditional column selection:
```julia
qval_col = params.match_between_runs ? :MBR_boosted_qval : :qval
global_qval_col = params.match_between_runs ? :MBR_boosted_global_qval : :global_qval

df[!, qval_col]
filter(row -> row[qval_col] < threshold, df)
```

---

### PHASE 3: Update writeCSVTables.jl for Conditional Columns

**File**: `src/Routines/SearchDIA/WriteOutputs/writeCSVTables.jl`

#### Change 3.1: Add match_between_runs Parameter to writePrecursorCSV

**Current signature** (approximate):
```julia
function writePrecursorCSV(
    long_precursors_path::String,
    file_names::Vector{String},
    normalized::Bool,
    proteins::LibraryProteins;
    write_csv::Bool = true,
    batch_size::Int64 = 10_000
)
```

**New signature**:
```julia
function writePrecursorCSV(
    long_precursors_path::String,
    file_names::Vector{String},
    normalized::Bool,
    proteins::LibraryProteins,
    match_between_runs::Bool;
    write_csv::Bool = true,
    batch_size::Int64 = 10_000
)
```

#### Change 3.2: Update Column Selection in makeWideFormat

**Location**: Line 147 (from error stack trace)

**Current pattern** (likely):
```julia
function makeWideFormat(df, ...)
    # Uses hardcoded :global_qval or :qval
    wide_df = unstack(df, [:precursor_idx, ...], :ms_file_idx, :global_qval)
end
```

**New pattern**:
```julia
function makeWideFormat(df, match_between_runs::Bool, ...)
    # Use conditional column names
    qval_col = match_between_runs ? :MBR_boosted_qval : :qval
    global_qval_col = match_between_runs ? :MBR_boosted_global_qval : :global_qval

    # Check which column exists and use it
    value_col = hasproperty(df, global_qval_col) ? global_qval_col : qval_col
    wide_df = unstack(df, [:precursor_idx, ...], :ms_file_idx, value_col)
end
```

#### Change 3.3: Update Any Column Lists or Selections

Look for places where columns are selected by name:
```julia
select(df, :precursor_idx, :global_qval, :qval, ...)
```

Replace with conditional selection:
```julia
base_cols = [:precursor_idx, ...]
qval_cols = if match_between_runs
    [:MBR_boosted_global_qval, :MBR_boosted_qval]
else
    [:global_qval, :qval]
end
select(df, vcat(base_cols, qval_cols), ...)
```

---

### PHASE 4: Update writePlots.jl for Conditional Columns

**File**: `src/Routines/SearchDIA/WriteOutputs/writePlots.jl`

#### Change 4.1: Search for Q-value References

**Actions:**
- Grep for `global_qval` in writePlots.jl
- Grep for `:qval` in writePlots.jl
- Update any plotting functions that use these columns

**Pattern**:
```julia
function plot_qvalues(df, match_between_runs::Bool)
    qval_col = match_between_runs ? :MBR_boosted_qval : :qval
    scatter(df[!, qval_col])
end
```

---

### PHASE 5: Update get_quant_necessary_columns (Again)

**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/scoring_interface.jl`

We already updated this in Phase 5 of the original fix, but we need to verify the q-value columns are correct:

**Current code** (from our previous fix):
```julia
function get_quant_necessary_columns(match_between_runs::Bool)
    base_cols = [
        :precursor_idx,
        :global_prob,
        :prec_prob,
        :trace_prob,
        :global_qval,  # ← This needs to be conditional!
        :run_specific_qval,
        ...
    ]

    if match_between_runs
        return vcat(base_cols, [
            :MBR_boosted_global_prob,
            :MBR_boosted_prec_prob,
            :MBR_boosted_trace_prob,
            :MBR_candidate,
            :MBR_transfer_q_value
        ])
    else
        return base_cols
    end
end
```

**Updated code**:
```julia
function get_quant_necessary_columns(match_between_runs::Bool)
    base_cols = [
        :precursor_idx,
        :global_prob,
        :prec_prob,
        :trace_prob,
        :run_specific_qval,
        :prec_mz,
        :pep,
        :weight,
        :target,
        :rt,
        :irt_obs,
        :missed_cleavage,
        :Mox,
        :isotopes_captured,
        :scan_idx,
        :entrapment_group_id,
        :ms_file_idx
    ]

    if match_between_runs
        # Add MBR-specific columns INCLUDING q-value columns
        return vcat(base_cols, [
            :MBR_boosted_global_prob,
            :MBR_boosted_prec_prob,
            :MBR_boosted_trace_prob,
            :MBR_boosted_global_qval,
            :MBR_boosted_qval,
            :MBR_candidate,
            :MBR_transfer_q_value
        ])
    else
        # Add standard q-value columns
        return vcat(base_cols, [
            :global_qval,
            :qval
        ])
    end
end
```

---

### PHASE 6: Update Step 10 to Recalculate global_qval Too

**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`

**Current Step 10** (from our previous fix):
```julia
step10_time = @elapsed begin
    if params.match_between_runs
        # Sort by non-MBR prec_prob
        sort_file_by_keys!(passing_refs, :prec_prob, :target; reverse=[true,true])

        # Merge by non-MBR prec_prob
        stream_sorted_merge(passing_refs, results.merged_quant_path, :prec_prob, :target;
                            batch_size=10_000_000, reverse=[true,true])

        # Calculate q-value spline using non-MBR scores
        qval_interp_nonMBR = get_qvalue_spline(
            results.merged_quant_path, :prec_prob, false;
            min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin,
            fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
        )

        # Add non-MBR qval column
        recalculate_qvalue_pipeline = TransformPipeline() |>
            add_interpolated_column(:qval, :prec_prob, qval_interp_nonMBR)

        passing_refs = apply_pipeline_batch(
            passing_refs,
            recalculate_qvalue_pipeline,
            passing_psms_folder
        )
    else
        # No Step 10 needed when MBR is disabled
    end
end
```

**Problem**: We recalculate `:qval` but not `:global_qval`, and we're using the wrong target column names.

**Updated Step 10**:
```julia
step10_time = @elapsed begin
    if params.match_between_runs
        # Sort by non-MBR prec_prob
        sort_file_by_keys!(passing_refs, :prec_prob, :target; reverse=[true,true])

        # Merge by non-MBR prec_prob
        stream_sorted_merge(passing_refs, results.merged_quant_path, :prec_prob, :target;
                            batch_size=10_000_000, reverse=[true,true])

        # Calculate q-value spline using non-MBR precursor scores
        qval_interp_nonMBR = get_qvalue_spline(
            results.merged_quant_path, :prec_prob, false;
            min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin,
            fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
        )

        # Recalculate MBR_boosted_qval column with non-MBR scores
        # (overwrites the MBR-boosted qval with unbiased qval for protein inference)
        recalculate_qvalue_pipeline = TransformPipeline() |>
            add_interpolated_column(:MBR_boosted_qval, :prec_prob, qval_interp_nonMBR)

        passing_refs = apply_pipeline_batch(
            passing_refs,
            recalculate_qvalue_pipeline,
            passing_psms_folder
        )

        # Now sort and merge by non-MBR global_prob for global q-values
        sort_file_by_keys!(passing_refs, :global_prob, :target; reverse=[true,true])
        stream_sorted_merge(passing_refs, results.merged_quant_path, :global_prob, :target;
                            batch_size=10_000_000, reverse=[true,true])

        # Calculate global q-value spline using non-MBR global scores
        global_qval_interp_nonMBR = get_qvalue_spline(
            results.merged_quant_path, :global_prob, true;
            min_pep_points_per_bin = params.precursor_q_value_interpolation_points_per_bin,
            fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
        )

        # Recalculate MBR_boosted_global_qval column with non-MBR scores
        recalculate_global_qvalue_pipeline = TransformPipeline() |>
            add_interpolated_column(:MBR_boosted_global_qval, :global_prob, global_qval_interp_nonMBR)

        passing_refs = apply_pipeline_batch(
            passing_refs,
            recalculate_global_qvalue_pipeline,
            passing_psms_folder
        )
    else
        # No Step 10 needed when MBR is disabled
    end
end
```

---

### PHASE 7: Create Helper Function for Column Name Resolution

**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/scoring_interface.jl`

Add a helper function to make conditional column selection easier:

```julia
"""
    get_qval_column_names(match_between_runs::Bool) -> NamedTuple

Returns the appropriate q-value column names based on whether MBR is enabled.

# Returns
- `qval`: Precursor q-value column name
- `global_qval`: Global precursor q-value column name
- `pep`: PEP column name
"""
function get_qval_column_names(match_between_runs::Bool)
    if match_between_runs
        return (
            qval = :MBR_boosted_qval,
            global_qval = :MBR_boosted_global_qval,
            pep = :pep  # PEP column name doesn't change
        )
    else
        return (
            qval = :qval,
            global_qval = :global_qval,
            pep = :pep
        )
    end
end
```

**Usage example**:
```julia
qval_cols = get_qval_column_names(params.match_between_runs)
df_filtered = filter(row -> row[qval_cols.qval] < 0.01, df)
```

---

### PHASE 8: Update Any Other References

**Search for additional references:**

```bash
# In Pioneer.jl root directory
grep -r "global_qval" src/ --include="*.jl"
grep -r ":qval" src/ --include="*.jl"
```

**Common patterns to look for:**
1. DataFrame column access: `df.global_qval`, `df[!, :qval]`
2. Column selection: `select(df, :qval, ...)`
3. Filtering: `filter(row -> row.qval < threshold, df)`
4. Plotting: `scatter(df.qval, ...)`
5. CSV writing: column name lists

**For each occurrence:**
- Determine if it needs to be conditional
- If yes, use the helper function or conditional logic
- If no (e.g., internal to a method that already knows the mode), document why

---

### PHASE 9: Testing

#### Test 9.1: Run with MBR Enabled
```julia
# Test with match_between_runs = true
SearchDIA("params_with_mbr.json")
```

**Verify:**
1. ScoringSearch completes without errors
2. Output files contain `MBR_boosted_qval` and `MBR_boosted_global_qval` columns
3. Output files do NOT contain `qval` or `global_qval` columns
4. MaxLFQSearch completes without column errors
5. CSV outputs are generated successfully
6. Plots are generated successfully

#### Test 9.2: Run with MBR Disabled
```julia
# Test with match_between_runs = false
SearchDIA("params_without_mbr.json")
```

**Verify:**
1. ScoringSearch completes without errors
2. Output files contain `qval` and `global_qval` columns
3. Output files do NOT contain `MBR_boosted_*` columns
4. MaxLFQSearch completes without errors
5. CSV outputs match expected format

#### Test 9.3: Inspect Arrow Files
```julia
using Arrow, DataFrames

# Check MBR mode output
df_mbr = DataFrame(Arrow.Table("path/to/mbr_output.arrow"))
@show names(df_mbr)  # Should include MBR_boosted_qval, MBR_boosted_global_qval
@assert :MBR_boosted_qval in names(df_mbr)
@assert :qval ∉ names(df_mbr)

# Check non-MBR mode output
df_no_mbr = DataFrame(Arrow.Table("path/to/no_mbr_output.arrow"))
@show names(df_no_mbr)  # Should include qval, global_qval
@assert :qval in names(df_no_mbr)
@assert :MBR_boosted_qval ∉ names(df_no_mbr)
```

---

## Summary of Changes

### Files to Modify:

1. **ScoringSearch.jl**
   - Update Step 10 to recalculate both `MBR_boosted_qval` and `MBR_boosted_global_qval`
   - Ensure correct target column names in Step 10

2. **scoring_interface.jl**
   - Fix `get_quant_necessary_columns()` to conditionally include correct q-val columns
   - Add helper function `get_qval_column_names()`

3. **MaxLFQSearch.jl**
   - Update all references to `:global_qval` and `:qval` to be conditional
   - Pass `match_between_runs` flag to downstream functions

4. **writeCSVTables.jl**
   - Add `match_between_runs` parameter to `writePrecursorCSV`
   - Update `makeWideFormat` to use conditional column names
   - Update any column selections/filters

5. **writePlots.jl** (if needed)
   - Update any plotting functions that reference q-value columns

### Key Principles:

1. **Conditional Everywhere**: Every hardcoded reference to `:qval` or `:global_qval` must check `match_between_runs`
2. **No Mixed Modes**: A file should have EITHER standard OR MBR-boosted column names, never both
3. **Helper Functions**: Use `get_qval_column_names()` for consistency
4. **Step 10 Overwrites**: When MBR is enabled, Step 10 overwrites the MBR-boosted q-values with unbiased values
5. **Documentation**: Comment each change to explain why conditional logic is needed

---

## Implementation Order

1. **Phase 1**: Search and document all q-value column references
2. **Phase 7**: Create helper function (makes subsequent steps easier)
3. **Phase 5**: Fix `get_quant_necessary_columns()` (foundation)
4. **Phase 6**: Update Step 10 to recalculate both q-value columns
5. **Phase 2**: Update MaxLFQSearch
6. **Phase 3**: Update writeCSVTables
7. **Phase 4**: Update writePlots (if needed)
8. **Phase 8**: Search for and fix any remaining references
9. **Phase 9**: Comprehensive testing

---

## Notes

- This fix ensures column names are consistent and predictable based on MBR mode
- Downstream code must be aware of MBR mode to use correct column names
- The `match_between_runs` parameter needs to be threaded through the call stack
- Step 10 is crucial: it overwrites MBR-boosted q-values with unbiased values for protein inference
- Consider adding validation checks to ensure expected columns exist before accessing them
