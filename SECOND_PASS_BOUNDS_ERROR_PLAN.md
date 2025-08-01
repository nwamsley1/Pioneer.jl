# Second Pass Search Bounds Error Diagnosis and Fix Plan

## Problem Summary

The Second Pass Search is failing with the exact same bounds error we fixed earlier:
```
BoundsError: attempt to access 5000-element Vector{UInt32} at index [5001]
```

Location: `rtIndexTransitionSelection.jl:92` in the assignment `precs_temp[precs_temp_size] = prec_idx`

## Root Cause Analysis

1. **Missing Fix**: The current branch (`fix-quad-tuning-empty-group-error`) was created from main before the bounds error fix was merged
2. **Already Fixed**: This exact issue was already fixed in PR #147 (commit bd1b4f4) and merged to main (commit 39252ab)
3. **Branch Divergence**: The current branch doesn't have the fix because it branched off before the fix was merged

## Solution Options

### Option 1: Cherry-pick the Fix (Quickest)
```bash
git cherry-pick bd1b4f4
```
This will apply the exact fix from the other branch to the current branch.

### Option 2: Merge Latest Main
```bash
git fetch origin
git merge origin/main
```
This will bring in all recent changes from main, including the bounds fix.

### Option 3: Rebase on Latest Main
```bash
git fetch origin
git rebase origin/main
```
This will replay your commits on top of the latest main, including the fix.

### Option 4: Manually Apply the Fix
Add the bounds check before line 92 in `rtIndexTransitionSelection.jl`:
```julia
# After line 90: precs_temp_size += 1
if precs_temp_size > length(precs_temp)
    append!(precs_temp, Vector{UInt32}(undef, length(precs_temp)))
end
# Then line 92: precs_temp[precs_temp_size] = prec_idx
```

## Recommended Approach

**Option 1 (Cherry-pick)** is the cleanest solution because:
- It brings in only the specific fix needed
- Avoids potential merge conflicts from other changes
- Maintains a clean commit history
- Quick and targeted

## Verification After Fix

The fix adds a bounds check that doubles the array size when needed:
```julia
if precs_temp_size > length(precs_temp)
    append!(precs_temp, Vector{UInt32}(undef, length(precs_temp)))
end
```

This pattern is already used in other parts of the codebase and has been tested.

## Prevention for Future

To avoid this issue in the future:
1. Always fetch and merge/rebase from main before creating new branches
2. Check for recent fixes to similar issues before starting work
3. Consider running the full test suite before major development

## Expected Outcome

After applying the fix, the Second Pass Search should complete successfully, handling cases where more than 5000 precursors are selected in an RT window.