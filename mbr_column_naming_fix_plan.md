# MBR Column Naming Fix - Detailed Implementation Plan

## Problem Statement

Current column naming is confusing and inconsistent:
- When MBR is ON: `prob` contains MBR-boosted values, `nonMBR_prob` contains unboosted values
- When MBR is OFF: Both `prob` and `nonMBR_prob` contain the same values
- The `nonMBR_prec_prob` column is created even when MBR is disabled
- Steps 3-10 use MBR-boosted scores when they should explicitly indicate this
- Protein inference (steps 11+) needs unboosted scores but column names don't make this clear

## Desired Behavior

### When MBR is ENABLED:
- **Trace-level columns (output from percolatorSortOf.jl):**
  - `trace_prob`: Base probability (non-MBR-enhanced)
  - `MBR_boosted_trace_prob`: MBR-enhanced probability

- **Precursor-level columns (created in Step 2):**
  - `prec_prob`: Aggregated from `trace_prob` (non-MBR)
  - `MBR_boosted_prec_prob`: Aggregated from `MBR_boosted_trace_prob`
  - `global_prob`: Log-odds of `prec_prob` (non-MBR)
  - `MBR_boosted_global_prob`: Log-odds of `MBR_boosted_prec_prob`

- **Q-value columns:**
  - `qval`, `global_qval`: Calculated from MBR-boosted scores (for PSM filtering in Steps 5-9)
  - Recalculated from non-MBR scores in Step 10 for protein inference

### When MBR is DISABLED:
- **Trace-level columns (output from percolatorSortOf.jl):**
  - `trace_prob`: Base probability only
  - (No `MBR_boosted_trace_prob` column)

- **Precursor-level columns (created in Step 2):**
  - `prec_prob`: Aggregated from `trace_prob`
  - `global_prob`: Log-odds of `prec_prob`
  - (No `MBR_boosted_*` columns)

- **Q-value columns:**
  - `qval`, `global_qval`: Calculated from `trace_prob`/`prec_prob`/`global_prob`

---

## Implementation Plan

### PHASE 1: Fix percolatorSortOf.jl Output
**File**: `src/utils/ML/percolatorSortOf.jl`

**Key Change**: Output columns with correct final names from the start - no renaming needed later!

#### Change 1.1: `sort_of_percolator_in_memory!` (Lines 423-431)
**Current code:**
```julia
if match_between_runs
    # Use the final MBR probabilities for all precursors
    psms[!, :prob] = MBR_estimates
    # Store nonMBR estimates for later q-value recalculation
    psms[!, :nonMBR_prob] = nonMBR_estimates
else
    psms[!, :prob] = prob_test
    # When MBR is disabled, nonMBR_prob equals the final prob
    psms[!, :nonMBR_prob] = prob_test
end
```

**New code:**
```julia
if match_between_runs
    # Store base trace probabilities (non-MBR)
    psms[!, :trace_prob] = nonMBR_estimates
    # Store MBR-enhanced trace probabilities
    psms[!, :MBR_boosted_trace_prob] = MBR_estimates
else
    # Only store base trace probabilities
    psms[!, :trace_prob] = prob_test
    # Do NOT create MBR_boosted_trace_prob column
end
```

#### Change 1.2: `sort_of_percolator_out_of_memory!` (Similar pattern)
**Location**: Lines ~680-690 (find the equivalent section)
**Apply same logic** as Change 1.1

---

### PHASE 2: Fix ScoringSearch.jl Step 2 (Precursor Probability Calculation)
**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`

#### Change 2.1: Conditional Probability Aggregation (Lines 304-329)

**Key Insight**: Since percolatorSortOf.jl now outputs `trace_prob` and `MBR_boosted_trace_prob`, we aggregate these directly to create precursor-level columns. `transform!()` creates NEW columns without removing the source columns. So after Step 2:
- Trace-level columns still exist: `trace_prob`, `MBR_boosted_trace_prob` (if MBR enabled)
- NEW precursor-level columns created: `prec_prob`, `MBR_boosted_prec_prob`, `global_prob`, `MBR_boosted_global_prob`
- No renaming needed in Step 4!

**Current code:**
```julia
if params.match_between_runs
    prob_col = apply_mbr_filter!(merged_df, params)
else
    prob_col = :prob
end

transform!(groupby(merged_df, [:precursor_idx, :ms_file_idx]),
           prob_col => (p -> begin
               prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
               prob = clamp(prob, eps(Float32), 1.0f0 - eps(Float32))
               Float32(prob)
           end) => :prec_prob)

transform!(groupby(merged_df, :precursor_idx),
           :prec_prob => (p -> logodds(p, sqrt_n_runs)) => :global_prob)

# Also create nonMBR_prec_prob for unbiased protein scoring
transform!(groupby(merged_df, [:precursor_idx, :ms_file_idx]),
           :nonMBR_prob => (p -> begin
               prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
               prob = clamp(prob, eps(Float32), 1.0f0 - eps(Float32))
               Float32(prob)
           end) => :nonMBR_prec_prob)
```

**New code:**
```julia
if params.match_between_runs
    # Apply MBR filter to MBR_boosted_trace_prob column (modifies in place)
    apply_mbr_filter!(merged_df, params)  # Now modifies :MBR_boosted_trace_prob directly

    # Aggregate MBR-boosted scores to precursor level
    transform!(groupby(merged_df, [:precursor_idx, :ms_file_idx]),
               :MBR_boosted_trace_prob => (p -> begin
                   prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
                   prob = clamp(prob, eps(Float32), 1.0f0 - eps(Float32))
                   Float32(prob)
               end) => :MBR_boosted_prec_prob)

    transform!(groupby(merged_df, :precursor_idx),
               :MBR_boosted_prec_prob => (p -> logodds(p, sqrt_n_runs)) => :MBR_boosted_global_prob)

    # Also aggregate non-MBR scores for protein inference (steps 11+)
    transform!(groupby(merged_df, [:precursor_idx, :ms_file_idx]),
               :trace_prob => (p -> begin
                   prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
                   prob = clamp(prob, eps(Float32), 1.0f0 - eps(Float32))
                   Float32(prob)
               end) => :prec_prob)

    transform!(groupby(merged_df, :precursor_idx),
               :prec_prob => (p -> logodds(p, sqrt_n_runs)) => :global_prob)
else
    # No MBR: only aggregate base probabilities
    transform!(groupby(merged_df, [:precursor_idx, :ms_file_idx]),
               :trace_prob => (p -> begin
                   prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
                   prob = clamp(prob, eps(Float32), 1.0f0 - eps(Float32))
                   Float32(prob)
               end) => :prec_prob)

    transform!(groupby(merged_df, :precursor_idx),
               :prec_prob => (p -> logodds(p, sqrt_n_runs)) => :global_prob)
end
```

**After this step, the DataFrame contains:**
- MBR enabled: `trace_prob`, `MBR_boosted_trace_prob`, `prec_prob`, `MBR_boosted_prec_prob`, `global_prob`, `MBR_boosted_global_prob`
- MBR disabled: `trace_prob`, `prec_prob`, `global_prob`

---

### PHASE 3: Fix apply_mbr_filter! to Modify MBR_boosted_trace_prob In-Place
**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/scoring_interface.jl`

#### Change 3.1: Wrapper Function (Lines 50-70)
**Current code:**
```julia
function apply_mbr_filter!(merged_df::DataFrame, params)
    # 1) identify transfer candidates
    candidate_mask = merged_df.MBR_transfer_candidate

    # 2) identify bad transfers
    is_bad_transfer = candidate_mask .& (
         (merged_df.target .& coalesce.(merged_df.MBR_is_best_decoy, false)) .|
         (merged_df.decoy .& .!coalesce.(merged_df.MBR_is_best_decoy, false))
    )

    # 3) Apply the main filtering function
    filtered_probs = apply_mbr_filter!(merged_df, candidate_mask, is_bad_transfer, params)

    # 4) Add filtered probabilities as a new column and return column name
    merged_df[!, :MBR_filtered_prob] = filtered_probs
    return :MBR_filtered_prob
end
```

**New code:**
```julia
function apply_mbr_filter!(merged_df::DataFrame, params)
    # 1) identify transfer candidates based on MBR_boosted_trace_prob
    candidate_mask = merged_df.MBR_transfer_candidate

    # 2) identify bad transfers
    is_bad_transfer = candidate_mask .& (
         (merged_df.target .& coalesce.(merged_df.MBR_is_best_decoy, false)) .|
         (merged_df.decoy .& .!coalesce.(merged_df.MBR_is_best_decoy, false))
    )

    # 3) Apply filtering and get filtered probabilities
    filtered_probs = apply_mbr_filter!(merged_df, candidate_mask, is_bad_transfer, params)

    # 4) Modify MBR_boosted_trace_prob column IN-PLACE
    merged_df[!, :MBR_boosted_trace_prob] = filtered_probs

    # 5) No return value needed (modifies column in place)
    return nothing
end
```

#### Change 3.2: Core Filter Function (Lines 78-150)
**Current code** uses `merged_df.prob` for input scores

**New code** should use `merged_df.MBR_boosted_trace_prob` for input scores
- Line ~93: Change `return merged_df.prob` to `return merged_df.MBR_boosted_trace_prob`
- Any other references to `.prob` should be `.MBR_boosted_trace_prob`

---

### PHASE 4: Update Step 4 Column Selection (NO RENAMING NEEDED!)
**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`

#### Change 4.1: Conditional Sorting in Step 4 (Lines 355-366)

**Purpose**: Select columns and sort by the appropriate score column (MBR-boosted if available, otherwise non-MBR).

**Current code:**
```julia
quant_processing_pipeline = TransformPipeline() |>
    add_best_trace_indicator(params.isotope_tracetype, best_traces) |>
    rename_column(:prob, :trace_prob) |>
    select_columns(vcat(necessary_cols, :best_trace)) |>
    filter_rows(row -> row.best_trace; desc="keep_only_best_traces") |>
    remove_columns(:best_trace) |>
    sort_by([:global_prob, :target], rev=[true, true])
```

**New code:**
```julia
if params.match_between_runs
    quant_processing_pipeline = TransformPipeline() |>
        add_best_trace_indicator(params.isotope_tracetype, best_traces) |>
        select_columns(vcat(necessary_cols, :best_trace)) |>
        filter_rows(row -> row.best_trace; desc="keep_only_best_traces") |>
        remove_columns(:best_trace) |>
        sort_by([:MBR_boosted_global_prob, :target], rev=[true, true])
else
    quant_processing_pipeline = TransformPipeline() |>
        add_best_trace_indicator(params.isotope_tracetype, best_traces) |>
        select_columns(vcat(necessary_cols, :best_trace)) |>
        filter_rows(row -> row.best_trace; desc="keep_only_best_traces") |>
        remove_columns(:best_trace) |>
        sort_by([:global_prob, :target], rev=[true, true])
end
```

**Note**: No renaming needed! Columns already have correct names from percolatorSortOf.jl

---

### PHASE 5: Update get_quant_necessary_columns()
**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/scoring_interface.jl`

#### Change 5.1: Column List (Lines 478-503)
**Current code:**
```julia
function get_quant_necessary_columns()
    return [
        :precursor_idx,
        :global_prob,
        :prec_prob,
        :nonMBR_prob,
        :nonMBR_prec_prob,
        :trace_prob,
        ...
    ]
end
```

**New code - make it conditional:**
```julia
function get_quant_necessary_columns(match_between_runs::Bool)
    base_cols = [
        :precursor_idx,
        :global_prob,
        :prec_prob,
        :trace_prob,
        :global_qval,
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
        # Add MBR-specific columns
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

**Update caller** (Line 353 in ScoringSearch.jl):
```julia
necessary_cols = get_quant_necessary_columns(params.match_between_runs)
```

---

### PHASE 6: Update Steps 5-10 to Use Correct Columns
**File**: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`

#### Change 6.1: Step 5 - Merge by MBR_boosted_global_prob if MBR (Lines 368-373)
**Current code:**
```julia
stream_sorted_merge(filtered_refs, results.merged_quant_path, :global_prob, :target;
                   batch_size=10_000_000, reverse=[true,true])
```

**New code:**
```julia
if params.match_between_runs
    stream_sorted_merge(filtered_refs, results.merged_quant_path, :MBR_boosted_global_prob, :target;
                       batch_size=10_000_000, reverse=[true,true])
else
    stream_sorted_merge(filtered_refs, results.merged_quant_path, :global_prob, :target;
                       batch_size=10_000_000, reverse=[true,true])
end
```

#### Change 6.2: Step 6 - Use MBR_boosted_global_prob for q-values if MBR (Lines 376-381)
Pass `match_between_runs` parameter to `get_precursor_global_qval_spline()` so it uses the correct column

#### Change 6.3: Step 7 - Merge by MBR_boosted_prec_prob if MBR (Lines 383-389)
**Similar pattern to 6.1**

#### Change 6.4: Step 8 - Use MBR_boosted_prec_prob for q-values if MBR (Lines 392-399)
**Similar pattern to 6.2**

#### Change 6.5: Step 9 - Filter using MBR q-values (Lines 401-429)
**Current code:**
```julia
qvalue_filter_pipeline = TransformPipeline() |>
    add_interpolated_column(:global_qval, :global_prob, results.precursor_global_qval_interp[]) |>
    add_interpolated_column(:qval, :prec_prob, results.precursor_qval_interp[]) |>
    add_interpolated_column(:pep, :prec_prob, results.precursor_pep_interp[]) |>
    filter_by_multiple_thresholds([
        (:global_qval, params.q_value_threshold),
        (:qval, params.q_value_threshold)
    ])
```

**New code:**
```julia
if params.match_between_runs
    qvalue_filter_pipeline = TransformPipeline() |>
        add_interpolated_column(:MBR_boosted_global_qval, :MBR_boosted_global_prob, results.precursor_global_qval_interp[]) |>
        add_interpolated_column(:MBR_boosted_qval, :MBR_boosted_prec_prob, results.precursor_qval_interp[]) |>
        add_interpolated_column(:pep, :MBR_boosted_prec_prob, results.precursor_pep_interp[]) |>
        filter_by_multiple_thresholds([
            (:MBR_boosted_global_qval, params.q_value_threshold),
            (:MBR_boosted_qval, params.q_value_threshold)
        ])
else
    qvalue_filter_pipeline = TransformPipeline() |>
        add_interpolated_column(:global_qval, :global_prob, results.precursor_global_qval_interp[]) |>
        add_interpolated_column(:qval, :prec_prob, results.precursor_qval_interp[]) |>
        add_interpolated_column(:pep, :prec_prob, results.precursor_pep_interp[]) |>
        filter_by_multiple_thresholds([
            (:global_qval, params.q_value_threshold),
            (:qval, params.q_value_threshold)
        ])
end
```

#### Change 6.6: REMOVE Step 9 Column Swap (Lines 420-428)
**DELETE this entire block** - it's no longer needed since we now have separate columns:
```julia
# DELETE THIS:
swap_to_nonMBR_pipeline = TransformPipeline() |>
    rename_column(:trace_prob, :MBR_trace_prob) |>
    rename_column(:prec_prob, :MBR_prec_prob) |>
    rename_column(:nonMBR_prob, :trace_prob) |>
    rename_column(:nonMBR_prec_prob, :prec_prob)
apply_pipeline!(passing_refs, swap_to_nonMBR_pipeline)
```

#### Change 6.7: SIMPLIFY Step 10 (Lines 433-460)
**Current**: Re-calculates q-values using swapped nonMBR scores
**New**: Only needed if MBR is enabled to recalculate using non-MBR scores

**New code:**
```julia
if params.match_between_runs
    #@debug_l1 "Step 10: Re-calculate q-values using non-MBR scores..."
    step10_time = @elapsed begin
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
    end
else
    # No Step 10 needed when MBR is disabled
    step10_time = 0.0
end
```

---

### PHASE 7: Update Helper Functions That Reference Column Names

#### Change 7.1: `get_precursor_global_qval_spline()` - Add column parameter
**Location**: Find this function and add a `score_col` parameter

#### Change 7.2: `get_precursor_qval_spline()` - Add column parameter
**Location**: Find this function and add a `score_col` parameter

#### Change 7.3: `get_precursor_pep_interpolation()` - Add column parameter
**Location**: Find this function and add a `score_col` parameter

---

### PHASE 8: Testing and Validation

#### Test 8.1: With MBR Enabled
1. Run ScoringSearch with `match_between_runs = true`
2. Verify columns exist:
   - `trace_prob`, `MBR_boosted_trace_prob` (trace level)
   - `prec_prob`, `MBR_boosted_prec_prob` (precursor level)
   - `global_prob`, `MBR_boosted_global_prob` (global level)
   - `qval` (from non-MBR after Step 10), `MBR_boosted_qval` (from MBR after Step 9)
3. Verify protein inference uses `prec_prob`/`global_prob` (non-MBR)
4. Verify PSM filtering uses `MBR_boosted_qval` (Steps 5-9)

#### Test 8.2: With MBR Disabled
1. Run ScoringSearch with `match_between_runs = false`
2. Verify columns exist:
   - `trace_prob`, `prec_prob`, `global_prob`
3. Verify NO `MBR_boosted_*` or `nonMBR_*` columns
4. Verify protein inference and PSM filtering use same columns

---

## Summary of Column Names and When They Exist

| MBR Status | After percolatorSortOf | After Step 2 (Aggregation) | Q-value Columns |
|-----------|------------------------|----------------------------|-----------------|
| **Enabled** | `trace_prob`<br>`MBR_boosted_trace_prob` | `trace_prob`, `MBR_boosted_trace_prob`<br>`prec_prob`, `MBR_boosted_prec_prob`<br>`global_prob`, `MBR_boosted_global_prob` | `qval` (non-MBR, Step 10)<br>`MBR_boosted_qval` (MBR, Step 9) |
| **Disabled** | `trace_prob` | `trace_prob`, `prec_prob`, `global_prob` | `qval` |

---

## Files to Modify (Complete List)

1. **src/utils/ML/percolatorSortOf.jl**
   - Lines 423-431: Fix column assignment in `sort_of_percolator_in_memory!` - output `trace_prob` and `MBR_boosted_trace_prob`
   - Lines ~680-690: Fix column assignment in `sort_of_percolator_out_of_memory!` - output `trace_prob` and `MBR_boosted_trace_prob`

2. **src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl**
   - Lines 304-329: Step 2 - Conditional aggregation using `trace_prob` and `MBR_boosted_trace_prob`
   - Line 353: Update `get_quant_necessary_columns()` call
   - Lines 355-366: Step 4 - Remove renaming, add conditional sorting
   - Lines 368-373: Step 5 - Conditional merge key
   - Lines 376-381: Step 6 - Conditional score column
   - Lines 383-389: Step 7 - Conditional merge key
   - Lines 392-399: Step 8 - Conditional score column
   - Lines 401-418: Step 9 - Conditional filtering
   - Lines 420-428: REMOVE swap pipeline
   - Lines 433-460: Step 10 - Make conditional on MBR

3. **src/Routines/SearchDIA/SearchMethods/ScoringSearch/scoring_interface.jl**
   - Lines 50-70: Fix `apply_mbr_filter!` wrapper to modify `MBR_boosted_trace_prob`
   - Lines 78-150: Fix core filter to use `MBR_boosted_trace_prob`
   - Lines 473-504: Update `get_quant_necessary_columns()`

4. **Helper functions** (find and update):
   - `get_precursor_global_qval_spline()`
   - `get_precursor_qval_spline()`
   - `get_precursor_pep_interpolation()`

---

## Implementation Order

1. **Phase 1**: Fix percolatorSortOf.jl output (foundation) - output `trace_prob` and `MBR_boosted_trace_prob` from the start
2. **Phase 3**: Fix apply_mbr_filter! to work with `MBR_boosted_trace_prob`
3. **Phase 2**: Fix Step 2 aggregation to use `trace_prob` and `MBR_boosted_trace_prob`
4. **Phase 5**: Update get_quant_necessary_columns()
5. **Phase 4**: Update Step 4 - remove renaming, add conditional sorting
6. **Phase 7**: Update helper functions
7. **Phase 6**: Update Steps 5-10
8. **Phase 8**: Testing

---

## Notes

- This plan maintains backward compatibility for non-MBR mode
- Protein inference (Steps 11+) always uses non-MBR scores (`prec_prob`, `global_prob`)
- PSM filtering (Steps 3-9) uses MBR-boosted scores when available (`MBR_boosted_prec_prob`, `MBR_boosted_global_prob`)
- Step 10 recalculates q-values using non-MBR scores for unbiased protein scoring
- All column names are now explicit about what they contain
- **No renaming needed** - columns have correct names from the start in percolatorSortOf.jl output
- Transform operations create NEW columns, so trace-level and precursor-level columns coexist after Step 2
