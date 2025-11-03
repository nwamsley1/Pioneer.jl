# BuildSpecLib Test Failure Analysis and Fix Plan

## Executive Summary

All 14 BuildSpecLib test failures share the **same root cause**: test configuration files contain placeholder paths from the simplified template that conflict with the new parameter structure implemented in our recent refactoring.

**Error**: `IOError: mkdir("/path"; mode=0o777): read-only file system (EROFS)`
**Location**: `src/Routines/BuildSpecLib.jl:73` (mkpath call)

---

## Root Cause Analysis

### Problem Chain

1. **Test param files contain mixed old/new formats**
   - Line 3: `"library_path": "/path/to/output/my_library"` (new format, placeholder)
   - Line 5: `"calibration_raw_file": "/path/to/calibration/file.arrow"` (new format, placeholder)
   - Lines 106-109: OLD format with actual paths:
     ```json
     "out_dir": "/Users/.../data/test_build_spec_lib/scenario_a/.../output",
     "lib_name": "minimal_mc1",
     "new_lib_name": "minimal_mc1",
     "out_name": "minimal_mc1.tsv"
     ```

2. **Parameter validation logic in check_params.jl**
   - Lines 169-179: Processes `library_path` first
   - Creates `_lib_dir = "/path/to/output/my_library.poin"`
   - Lines 48-55: Backward compatibility code checks `if !haskey(params, "library_path")`
   - **BUT** `library_path` DOES exist (as placeholder), so backward compatibility **doesn't run**
   - Old parameters with correct paths are ignored

3. **BuildSpecLib.jl attempts directory creation**
   - Line ~73: `mkpath(params["_lib_dir"])`
   - Tries to create `/path/to/output/my_library.poin`
   - Fails trying to create root directory `/path` (read-only filesystem)

### Why This Happened

When we created the simplified template (`defaultBuildLibParamsSimplified.json`), it included placeholder paths like `/path/to/output/my_library`. These placeholders were copied into test files but never replaced with actual test paths because the tests **already had** working paths in the old `out_dir`/`lib_name` format.

The test files ended up with **both formats**, and the new format (with placeholders) took precedence over the old format (with actual paths).

---

## Affected Test Files

All test parameter files in `data/test_build_spec_lib/`:

1. ✗ `scenario_a_missed_cleavages/params_1mc.json`
2. ✗ `scenario_a_missed_cleavages/params_0mc.json`
3. ✗ `scenario_a_missed_cleavages/params_2mc.json`
4. ✗ `scenario_a_missed_cleavages/params_3mc.json`
5. ✗ `scenario_b_standard/params_altimeter.json`
6. ✗ `scenario_c_multifile/params_altimeter.json`
7. ✗ `scenario_var_mods/max_0_mods/params_var_0.json`
8. ✗ `scenario_var_mods/max_1_mods/params_var_1.json`
9. ✗ `scenario_var_mods/max_2_mods/params_var_2.json`
10. ✗ `scenario_var_mods/max_3_mods/params_var_3.json`
11. ✗ `scenario_var_mods/max_4_mods/params_var_4.json`
12. ✗ `scenario_frag_count/mult_1/params_frag_1.json`
13. ✗ `scenario_frag_count/mult_2/params_frag_2.json`
14. ✗ `scenario_frag_count/mult_3/params_frag_3.json`
15. ✗ `scenario_frag_count/mult_4/params_frag_4.json`
16. ✗ `scenario_d_comprehensive/params_altimeter.json`
17. ✗ `scenario_prosit/params_prosit.json`

---

## Detailed Error Breakdown

### Error 1-14: All BuildSpecLib Tests

**Error Message**:
```
IOError: mkdir("/path"; mode=0o777): read-only file system (EROFS)
```

**Stacktrace Key Points**:
```julia
[7] mkpath(path::String)
    @ Base.Filesystem ./file.jl:237
[8] macro expansion
    @ ~/Projects/.../BuildSpecLib.jl:73
```

**What's Happening**:
1. Test loads `params_1mc.json` (or other test param file)
2. `check_params_bsp()` sees `library_path: "/path/to/output/my_library"`
3. Line 173: `params["_lib_dir"] = "/path/to/output/my_library.poin"`
4. BuildSpecLib.jl line 73: `mkpath("/path/to/output/my_library.poin")`
5. System tries to create `/path/` directory (requires root)
6. Fails with "read-only file system" error

**Why It's Failing**:
- The placeholder path `/path/to/output/my_library` is invalid
- Attempting to create directories starting from root `/path` requires superuser permissions
- The actual valid paths in `out_dir`/`lib_name` are being ignored

---

## Fix Strategy

### Option A: Remove Placeholders (RECOMMENDED)

**Approach**: Remove placeholder `library_path` and `calibration_raw_file` from all test files, keeping only the old format parameters. The backward compatibility code will automatically convert them.

**Advantages**:
- ✓ Tests the backward compatibility code path
- ✓ Minimal changes to test files (just delete 4 lines)
- ✓ Validates that old configs can seamlessly migrate
- ✓ Less prone to path errors

**Implementation**:
1. Remove lines 2-5 from each test param file:
   ```json
   "comment_output": "=== Output library path ===",
   "library_path": "/path/to/output/my_library",
   "comment_calibration": "=== Calibration file for automatic m/z bounds detection ===",
   "calibration_raw_file": "/path/to/calibration/file.arrow",
   ```

2. Keep existing lines (currently 106-109):
   ```json
   "out_dir": "/Users/.../scenario_a_missed_cleavages/output",
   "lib_name": "minimal_mc1",
   "new_lib_name": "minimal_mc1",
   "out_name": "minimal_mc1.tsv"
   ```

3. Backward compatibility code (check_params.jl:48-55) will automatically:
   - Create `library_path = joinpath(out_dir, lib_name)`
   - Convert to new internal format with `_lib_dir`
   - Delete old parameters

**Code Changes Required**: None - this tests existing backward compatibility

---

### Option B: Replace Placeholders with Real Paths

**Approach**: Replace placeholder paths with actual test paths and remove old parameters.

**Advantages**:
- ✓ Tests the new parameter format
- ✓ More forward-looking

**Disadvantages**:
- ✗ Requires changes to all test files
- ✗ More error-prone (17 files to update with different paths)
- ✗ Doesn't test backward compatibility
- ✗ Requires defining calibration_raw_file for each test

**Implementation Example** (for params_1mc.json):
```json
{
    "library_path": "/Users/nathanwamsley/Projects/EntrapmentTests/Pioneer.jl/data/test_build_spec_lib/scenario_a_missed_cleavages/output/minimal_mc1",
    "calibration_raw_file": "/Users/nathanwamsley/Projects/EntrapmentTests/Pioneer.jl/data/test_calibration_file.arrow",
    "fasta_paths": [...],
    ...
    // REMOVE these lines:
    // "out_dir": "...",
    // "lib_name": "...",
    // "new_lib_name": "...",
    // "out_name": "..."
}
```

---

### Option C: Add Placeholder Detection (Defensive)

**Approach**: Add validation to detect and reject placeholder paths before attempting to create directories.

**Advantages**:
- ✓ Prevents cryptic errors
- ✓ Gives clear user feedback

**Implementation**:
Add to `check_params.jl` after line 170:
```julia
# Validate that paths are not placeholders
if startswith(params["library_path"], "/path/")
    throw(InvalidParametersError(
        "library_path appears to be a placeholder. Please provide an actual output path.",
        params
    ))
end
if startswith(params["calibration_raw_file"], "/path/")
    throw(InvalidParametersError(
        "calibration_raw_file appears to be a placeholder. Please provide an actual calibration file path.",
        params
    ))
end
```

**Note**: This should be combined with Option A or B, not used alone.

---

## Recommended Solution

**Use Option A** (Remove Placeholders) because:

1. **Minimal disruption**: Only requires removing 4 lines from each test file
2. **Tests backward compatibility**: Validates that old configs work seamlessly
3. **Low risk**: Uses paths that already exist and work
4. **Fast to implement**: Simple deletion, no path construction needed

**Optionally add Option C** for better user experience in production.

---

## Implementation Checklist

### Phase 1: Fix Test Files (Option A)
- [ ] Remove placeholder lines from all 17 test param files
  - [ ] Lines to remove (typically lines 2-5):
    - `"comment_output": "=== Output library path ===",`
    - `"library_path": "/path/to/output/my_library",`
    - `"comment_calibration": "=== Calibration file for automatic m/z bounds detection ===",`
    - `"calibration_raw_file": "/path/to/calibration/file.arrow",`

### Phase 2: Verify Backward Compatibility Works
- [ ] Run single test: `julia --project=. -e 'using Pkg; Pkg.test("Pioneer"; test_args=["BuildSpecLib"])'`
- [ ] Confirm tests pass
- [ ] Check that config.json outputs have new format

### Phase 3: Optional - Add Placeholder Detection (Option C)
- [ ] Add validation in check_params.jl after line 170
- [ ] Add test case for placeholder detection
- [ ] Verify error message is clear

### Phase 4: Full Test Suite
- [ ] Run complete test suite: `julia --project=. -e 'using Pkg; Pkg.test()'`
- [ ] Verify all 17 BuildSpecLib tests pass
- [ ] Commit fixes with clear message

---

## Expected Test Results After Fix

### Before Fix
```
Test Summary:                                  | Pass  Error  Total
BuildSpecLib Integration Tests                 |   58     14     72
```

### After Fix
```
Test Summary:                                  | Pass  Total
BuildSpecLib Integration Tests                 |   72     72
```

All 14 errors should convert to passes, testing:
- ✓ Backward compatibility (old → new parameter conversion)
- ✓ Path handling with actual test directories
- ✓ .poin extension handling
- ✓ All existing test scenarios (missed cleavages, mods, etc.)

---

## Prevention Strategy

### For Future Template Updates

1. **Never include placeholder paths in test files**
   - Test files should always have real, absolute paths
   - Or use relative paths that resolve correctly

2. **Separate example configs from test configs**
   - `assets/example_config/` → For users (can have placeholders)
   - `data/test_*/` → For tests (must have real paths)

3. **Template validation**
   - Add CI check to detect `/path/to/` in test config files
   - Fail build if placeholders found in test data

4. **Documentation**
   - Clearly mark example files as templates
   - Add comment to test files: `"DO NOT USE PLACEHOLDER PATHS IN TEST FILES"`

---

## Testing the Fix

### Manual Test
```bash
# From project root
cd /Users/nathanwamsley/Projects/EntrapmentTests/Pioneer.jl

# Test single scenario
julia --project=. -e 'using Pioneer; Pioneer.BuildSpecLib("data/test_build_spec_lib/scenario_a_missed_cleavages/params_1mc.json")'

# Should create:
# data/test_build_spec_lib/scenario_a_missed_cleavages/output/minimal_mc1.poin/

# Verify config.json has new format
cat data/test_build_spec_lib/scenario_a_missed_cleavages/output/minimal_mc1.poin/config.json | grep library_path
# Should show: "library_path": "/Users/.../minimal_mc1"
```

### Full Test Suite
```bash
julia --project=. -e 'using Pkg; Pkg.test("Pioneer")'
```

---

## Files to Modify

### Test Parameter Files (17 files - remove placeholder lines)
```
data/test_build_spec_lib/scenario_a_missed_cleavages/params_0mc.json
data/test_build_spec_lib/scenario_a_missed_cleavages/params_1mc.json
data/test_build_spec_lib/scenario_a_missed_cleavages/params_2mc.json
data/test_build_spec_lib/scenario_a_missed_cleavages/params_3mc.json
data/test_build_spec_lib/scenario_b_standard/params_altimeter.json
data/test_build_spec_lib/scenario_c_multifile/params_altimeter.json
data/test_build_spec_lib/scenario_var_mods/max_0_mods/params_var_0.json
data/test_build_spec_lib/scenario_var_mods/max_1_mods/params_var_1.json
data/test_build_spec_lib/scenario_var_mods/max_2_mods/params_var_2.json
data/test_build_spec_lib/scenario_var_mods/max_3_mods/params_var_3.json
data/test_build_spec_lib/scenario_var_mods/max_4_mods/params_var_4.json
data/test_build_spec_lib/scenario_frag_count/mult_1/params_frag_1.json
data/test_build_spec_lib/scenario_frag_count/mult_2/params_frag_2.json
data/test_build_spec_lib/scenario_frag_count/mult_3/params_frag_3.json
data/test_build_spec_lib/scenario_frag_count/mult_4/params_frag_4.json
data/test_build_spec_lib/scenario_d_comprehensive/params_altimeter.json
data/test_build_spec_lib/scenario_prosit/params_prosit.json
```

### Optional: Source Code (if adding placeholder detection)
```
src/Routines/BuildSpecLib/utils/check_params.jl (add validation after line 170)
```

---

## Summary

The test failures are caused by placeholder paths from the simplified template being present alongside actual test paths in the old format. The new parameter code processes the placeholders first, causing attempts to create invalid directories.

**Fix**: Remove placeholder paths from test files, allowing backward compatibility code to convert existing valid paths to the new format.

**Time Estimate**: 15-20 minutes to fix all test files + 5 minutes to run tests
**Risk Level**: Low - simple deletions, tests existing code paths
**Testing Required**: Full test suite run to verify all 17 tests pass
