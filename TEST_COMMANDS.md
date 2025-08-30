# Test Commands for Partner Precursor Index Fix

## Prerequisites
```julia
# Start Julia with proper setup
julia --threads 15 --gcthreads 7,1 --project=dev

# In Julia REPL
using Revise, Pioneer
```

## Test Steps

### Step 1: Run Validation Script
```julia
include("VALIDATE_KEAP1_FIX.jl")
```

### Step 2: Manual Verification (Alternative)
```julia
# 1. Rebuild library
BuildSpecLib("test_keap1_params.json")

# 2. Load and check
using DataFrames, Arrow, Tables
ptable = DataFrame(Tables.columntable(Arrow.Table("/tmp/keap1_test_library/keap1_test.poin/precursors_table.arrow")))

# 3. Critical test - partner indices must be valid
max_partner_idx = Int64(maximum(skipmissing(ptable.partner_precursor_idx)))
table_size = nrow(ptable)
println("Max partner index: $max_partner_idx")
println("Table size: $table_size") 
println("Valid: $(max_partner_idx <= table_size)")

# 4. Verify pairing still works
gptable = groupby(ptable, :pair_id)
pair_sizes = [nrow(subdf) for (key, subdf) in pairs(gptable)]
println("Pair sizes: $(unique(pair_sizes))")  # Should be [2]

# 5. Check balance
println("Targets: $(sum(.!ptable.is_decoy))")
println("Decoys: $(sum(ptable.is_decoy))")
```

### Step 3: Test SearchDIA (If available)
```julia
# If you have SearchDIA parameters configured
# SearchDIA("path/to/your/search_params.json")
```

## Expected Results

### Before Fix
- Max partner_precursor_idx: 126,424
- Table size: 84,184  
- Valid: **false** ❌
- SearchDIA: BoundsError

### After Fix
- Max partner_precursor_idx: ≤ 84,184
- Table size: 84,184
- Valid: **true** ✅
- SearchDIA: No bounds error

## Success Criteria
1. ✅ `max_partner_idx <= table_size` 
2. ✅ `unique(pair_sizes) == [2]`
3. ✅ Target count == Decoy count
4. ✅ No missing partner relationships
5. ✅ SearchDIA runs without BoundsError

## Quick Test Command
```julia
# One-liner to test the most critical fix
include("VALIDATE_KEAP1_FIX.jl")
```

## Rollback if Issues
If validation fails, revert with:
```bash
git checkout HEAD~1 -- src/Routines/BuildSpecLib/chronologer/chronologer_prep.jl
```