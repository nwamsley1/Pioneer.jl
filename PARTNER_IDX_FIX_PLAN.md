# Fix Plan: Partner Precursor Index Assignment Order

## Problem Analysis

### Root Cause
The `partner_precursor_idx` values are assigned **before** filtering by m/z range, causing invalid row indices after filtering removes precursors.

### Current Broken Workflow
1. **`build_fasta_df()` called** (chronologer_prep.jl:198-204)
   - Inside function: `add_charge_specific_partner_columns!()` called (line 583)
   - Assigns `partner_precursor_idx` based on current row positions (e.g., rows 1-126,424)

2. **After `build_fasta_df()` returns:**
   - Calculate m/z values (lines 206-211)
   - Add length and missed cleavages (lines 213-220)
   - **FILTERING HAPPENS** (line 223): `filter!(x -> (x.mz >= prec_mz_min) & (x.mz <= prec_mz_max), fasta_df)`
   - Removes ~42,240 rows, leaving 84,184 rows

3. **Result:**
   - `partner_precursor_idx` values still point to old row positions (up to 126,424)
   - SearchDIA tries to access these invalid indices → BoundsError

### Evidence from Test Library
```julia
julia> Int64(maximum(ptable[!,:partner_precursor_idx]))
126424

julia> size(ptable, 1)
84184
```
Partner indices exceed table size by ~50%!

## Solution: Move Pairing After ALL Filtering

### Primary Solution (Recommended)

**Step 1: Remove pairing from `build_fasta_df()`**
- File: `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`
- Line 583: Comment out `seq_df = add_charge_specific_partner_columns!(seq_df)`

**Step 2: Move pairing after filtering in `prepare_chronologer_input()`**
- File: `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`
- After line 223 (the `filter!` call), add:
```julia
# Apply charge-specific target-decoy pairing AFTER all filtering
fasta_df = add_charge_specific_partner_columns!(fasta_df)
```

**Step 3: Verify no downstream filtering**
- Ensure no other filter operations happen after pairing assignment
- Partner indices must remain valid through final library output

### Alternative Solutions (if primary fails)

**Option B: Update indices after filtering**
```julia
# After filtering, remap partner indices to new row positions
old_to_new_idx = Dict{Int, Int}()
for (new_idx, row) in enumerate(eachrow(filtered_df))
    # Build mapping based on some unique identifier
    old_to_new_idx[row.original_position] = new_idx
end

# Remap partner indices
filtered_df.partner_precursor_idx = [
    haskey(old_to_new_idx, idx) ? old_to_new_idx[idx] : missing 
    for idx in filtered_df.partner_precursor_idx
]
```

**Option C: Filter before pairing entirely**
- Calculate m/z values earlier in pipeline
- Apply all filters before calling `build_fasta_df()`
- More complex restructuring required

## Implementation Details

### Files to Modify
1. `src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl`
   - Remove pairing from `build_fasta_df()` function (line 583)
   - Add pairing after filtering in `prepare_chronologer_input()` (after line 223)

### Code Changes

**Remove from `build_fasta_df()` (line 583):**
```julia
# Apply charge-specific target-decoy pairing
# seq_df = add_charge_specific_partner_columns!(seq_df)  # MOVED TO AFTER FILTERING
```

**Add to `prepare_chronologer_input()` after line 223:**
```julia
# Filter by mass range
filter!(x -> (x.mz >= prec_mz_min) & (x.mz <= prec_mz_max), fasta_df)

# Apply charge-specific target-decoy pairing AFTER all filtering is complete
fasta_df = add_charge_specific_partner_columns!(fasta_df)
```

## Expected Results

### Before Fix
- Total precursors: 84,184
- Max partner_precursor_idx: 126,424 (INVALID!)
- SearchDIA: BoundsError when accessing index 126,424+

### After Fix
- Total precursors: 84,184
- Max partner_precursor_idx: ≤ 84,184 (VALID!)
- SearchDIA: No bounds errors
- Pairing still works: `unique(pair_sizes) == [2]`

## Validation Strategy

### Test Commands
```julia
# 1. Rebuild library
BuildSpecLib("test_keap1_params.json")

# 2. Load and verify
ptable = DataFrame(Tables.columntable(Arrow.Table("/tmp/keap1_test_library/keap1_test.poin/precursors_table.arrow")))

# 3. Check partner indices are valid
max_partner_idx = Int64(maximum(skipmissing(ptable.partner_precursor_idx)))
table_size = nrow(ptable)
println("Max partner index: $max_partner_idx")
println("Table size: $table_size")
println("Valid indices: $(max_partner_idx <= table_size)")

# 4. Verify pairing still works
gptable = groupby(ptable, :pair_id)
pair_sizes = [nrow(subdf) for (key, subdf) in pairs(gptable)]
println("Pair sizes: $(unique(pair_sizes))")  # Should be [2]

# 5. Test SearchDIA
SearchDIA("path/to/search_params.json")  # Should not crash
```

### Success Criteria
1. ✅ Max partner_precursor_idx ≤ nrow(ptable)
2. ✅ All pair_ids have exactly 2 members
3. ✅ SearchDIA runs without BoundsError
4. ✅ Search results are reasonable

## Risk Assessment

### Low Risk Changes
- Moving function call location (no logic changes)
- Partner assignment logic remains identical
- Pairing algorithm unchanged

### Potential Issues
- If any downstream filtering removes rows after pairing assignment
- If pairing function has side effects we're not aware of
- If SearchDIA expects partner indices to be assigned at a specific point

### Mitigation
- Thorough testing with keap1.fasta
- Check for any other filter operations in the pipeline
- Verify SearchDIA integration works correctly

## Rollback Plan
If issues arise, simply revert the changes:
1. Uncomment pairing in `build_fasta_df()`
2. Remove pairing from `prepare_chronologer_input()`
3. Return to previous working state (with known SearchDIA error)

## Notes
- This fix addresses the immediate bounds error
- Does not change the pairing logic itself (which works correctly)
- Should be backward compatible with existing libraries
- May require SearchDIA parameter tuning if search behavior changes