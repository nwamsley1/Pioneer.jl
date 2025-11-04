# Plan: Simplify Library Prediction CLI and Improve User Experience

## Implementation Status

- ✅ **Phase 1: CLI Simplification** (Completed - Commit 58c87c73)
- ⏳ **Phase 2: Template Comments** (Pending)
- ⏳ **Phase 3: Entrapment Messages** (Pending)
- ⏳ **Phase 4: Testing and Validation** (Pending)

---

## Overview
This plan addresses three main areas of improvement for the library prediction functionality:
1. Simplify the `params-predict` command line interface
2. Improve JSON parameter template comments for better clarity
3. Conditionally display entrapment-related messages only when entrapment is enabled

---

## 1. Simplify `params-predict` Command Line Interface

### Current State
The command currently takes three positional arguments:
```bash
pioneer params-predict <library_outpath> <lib_name> <fasta_path>
```

**Example:**
```bash
pioneer params-predict yeast.poin yeast fasta/
```

### Proposed Changes
Simplify to two positional arguments by deriving the library name from the output path:
```bash
pioneer params-predict <library_outpath> <fasta_path>
```

**Example:**
```bash
pioneer params-predict yeast.poin fasta/
```

The library name will be extracted as the last component of the path (the folder/directory name).

### Implementation Details

#### Files to Modify

**1. CLI Script: `src/build/CLI/pioneer`**
- **Location:** Line 29
- **Current:**
  ```bash
  echo "  params-predict <library_outpath> <lib_name> <fasta_path> [--params-path <params_out_path>]"
  ```
- **Change to:**
  ```bash
  echo "  params-predict <library_outpath> <fasta_path> [--params-path <params_out_path>]"
  ```

- **Location:** Line 41 (example)
- **Current:**
  ```bash
  echo "  pioneer params-predict yeast.poin yeast fasta/ --params-path predict_params.json"
  ```
- **Change to:**
  ```bash
  echo "  pioneer params-predict yeast.poin fasta/ --params-path predict_params.json"
  ```

**2. Parameter Generation Entry Point: `src/Routines/GenerateParams.jl`**
- **Location:** Lines 61-99 (`main_GetBuildLibParams` function)
- **Changes needed:**
  1. Remove the `lib_name` argument from the `@add_arg_table!` block (lines 68-70)
  2. Update the `fasta_dir` argument name to `fasta_path` (line 71)
  3. Update the help text for `fasta_path` to read: "Directory containing FASTA files or path to a specific FASTA file"
  4. Modify the function call to `GetBuildLibParams` (lines 89-93):
     - Extract library name from `out_dir` using `basename(out_dir)`
     - Pass the extracted name as `lib_name` parameter

  **Pseudo-code for extraction:**
  ```julia
  # Extract library name from output path
  lib_name = basename(parsed_args[:out_dir])

  GetBuildLibParams(parsed_args[:out_dir],
                    lib_name,  # Derived from out_dir
                    parsed_args[:fasta_path];  # Renamed from fasta_dir
                    params_path=params_path,
                    simplified=simplified)
  ```

**3. Help Text for Subcommand: Update ArgParse settings**
- **Location:** `src/Routines/GenerateParams.jl`, lines 64-81
- **Current argument definitions:**
  ```julia
  "out_dir"
      help = "Output directory for library"
      arg_type = String
  "lib_name"
      help = "Name of the library"
      arg_type = String
  "fasta_dir"
      help = "Directory containing FASTA files"
      arg_type = String
  ```
- **Change to:**
  ```julia
  "out_dir"
      help = "Output directory for library"
      arg_type = String
  "fasta_path"
      help = "Directory containing FASTA files or path to a specific FASTA file"
      arg_type = String
  ```

### Testing Requirements
After implementation, test the following scenarios:
1. **Simple directory path:** `pioneer params-predict my_library.poin fasta/`
2. **Nested path:** `pioneer params-predict output/results/yeast.poin fasta/`
3. **Absolute path:** `pioneer params-predict /path/to/library.poin /path/to/fasta/`
4. **With params-path option:** `pioneer params-predict lib.poin fasta/ --params-path custom.json`
5. **Full parameter template:** `pioneer params-predict lib.poin fasta/ --full`

---

## 2. Improve JSON Parameter Template Comments

### Current Issues
- The calibration comment field uses `\n` for newlines, which don't render properly in JSON
- No comments explaining the relationship between `auto_detect_frag_bounds` and manual bounds
- Library m/z bounds section lacks clarity

### Proposed Changes

#### Files to Modify

**1. Simplified Template: `assets/example_config/defaultBuildLibParamsSimplified.json`**

**Current state (lines 5-6):**
```json
"comment_calibration": "=== Calibration file for automatic m/z bounds detection ===\n    Convert a representative RAW file to Arrow format, then provide path here.\n    Used to detect fragment and precursor m/z ranges for your instrument.",
"calibration_raw_file": "/path/to/calibration/file.arrow",
```

**Change to:**
```json
"comment_calibration": "=== Calibration file for automatic detection of fragment and precursor m/z ranges (optional) ===",
"comment_calibration": "=== Convert a representative RAW file to Arrow format and provide its path here ===",
"calibration_raw_file": "/path/to/calibration/file.arrow",
```

**Current state (line 20):**
```json
"comment_library": "=== Library m/z bounds ===\n    Note: If auto_detect_frag_bounds=true (default), these bounds are automatically\n    detected from calibration_raw_file. Set to false to use manual bounds below.",
```

**Change to:**
```json
"comment_library": "=== Library m/z bounds ===",
"comment_library": "=== If auto_detect_frag_bounds=true, bounds are determined from the calibration_raw_file; otherwise manual values are used ===",
```

**2. Full Template: `assets/example_config/defaultBuildLibParams.json`**

Make the same changes as above:
- **Lines 5-6:** Split calibration comment into two separate comment fields
- **Lines 55-56 region (in library_params section):** Add clarifying comments about auto-detection

**Before library_params section:**
```json
"comment_library": "=== Library generation parameters ===",
"library_params": {
    ...
}
```

**Change to:**
```json
"comment_library": "=== Library m/z bounds ===",
"comment_library": "=== If auto_detect_frag_bounds=true, bounds are determined from the calibration_raw_file; otherwise manual values are used ===",
"library_params": {
    ...
}
```

### Rationale
- Duplicate comment fields are allowed in JSON (same key name can appear multiple times)
- Breaking long comments into separate fields ensures each line is readable
- Makes it clearer that calibration is optional and what it's used for
- Explains the conditional logic for m/z bounds determination

---

## 3. Conditionally Display Entrapment Messages

### Current Issue
When `entrapment_r = 0` (entrapment disabled), the following messages still appear:
```
Adding entrapment target indices...
┌ Warning: Decoy column not found - proceeding without decoy filtering
┌ Warning: Available columns: [...]
[ Info: Found 623735 entrapment targets for index mapping
[ Info: Entrapment pairing Stage 2 complete (entrapment_target_idx created):
[ Info:   Mapped 623735 entries to target indices
[ Info:   Max target index: 1247464
[ Info:   Table size: 1247466 rows
[ Info:   Decoy column used: ❌ NO
   After add_entrapment_indices!: max entrapment_target_idx = 1247464
   Entrapment indices valid: ✅ YES
```

These messages are confusing when entrapment is not enabled.

### Proposed Solution
Only display entrapment-related messages when `entrapment_r > 0`.

### Implementation Details

#### Files to Modify

**1. Main Library Building: `src/Routines/BuildSpecLib.jl`**

**Location:** Lines 326-334

**Current code:**
```julia
# Add entrapment target indices if entrapment_pair_id column exists
if hasproperty(precursors_table, :entrapment_pair_id)
    println("   Adding entrapment target indices...")
    add_entrapment_indices!(precursors_table)
    ent_col = precursors_table.entrapment_target_idx
    max_entrap_target = all(ismissing, ent_col) ? 0 : Int64(maximum(skipmissing(ent_col)))
    println("   After add_entrapment_indices!: max entrapment_target_idx = $max_entrap_target")
    println("   Entrapment indices valid: $(max_entrap_target <= nrow(precursors_table) ? "✅ YES" : "❌ NO")")
end
```

**Change to:**
```julia
# Add entrapment target indices if entrapment_pair_id column exists AND entrapment is enabled
entrapment_r = get(params["fasta_digest_params"], "entrapment_r", 0)
if hasproperty(precursors_table, :entrapment_pair_id) && entrapment_r > 0
    println("   Adding entrapment target indices...")
    add_entrapment_indices!(precursors_table)
    ent_col = precursors_table.entrapment_target_idx
    max_entrap_target = all(ismissing, ent_col) ? 0 : Int64(maximum(skipmissing(ent_col)))
    println("   After add_entrapment_indices!: max entrapment_target_idx = $max_entrap_target")
    println("   Entrapment indices valid: $(max_entrap_target <= nrow(precursors_table) ? "✅ YES" : "❌ NO")")
end
```

**2. Entrapment Index Mapping: `src/Routines/BuildSpecLib/chronologer/pair_decoys.jl`**

**Function:** `add_entrapment_indices!` (starts at line 377)

**Location:** Lines 399-403 (Warning messages)

**Current code:**
```julia
# Check if decoy column exists (be robust about column names)
has_decoy_col = hasproperty(df, :decoy)
if !has_decoy_col
    @user_warn "Decoy column not found - proceeding without decoy filtering"
    @user_warn "Available columns: $(names(df))"
end
```

**Change to:**
```julia
# Check if decoy column exists (be robust about column names)
has_decoy_col = hasproperty(df, :decoy)
# Only warn about missing decoy column if we're actually processing entrapment data
# (This function is only called when entrapment_r > 0, but the warning is still confusing)
# We can silently proceed if decoy column is missing
if !has_decoy_col
    @debug_l2 "Decoy column not found - proceeding without decoy filtering"
    @debug_l2 "Available columns: $(names(df))"
end
```

**Location:** Line 426 (Info message)

**Current code:**
```julia
@user_info "Found $(length(pair_to_target)) entrapment targets for index mapping"
```

**Change to:**
```julia
# Only log if we found entrapment targets (i.e., entrapment is actually enabled)
if length(pair_to_target) > 0
    @user_info "Found $(length(pair_to_target)) entrapment targets for index mapping"
end
```

**Location:** Lines 455-465 (Final reporting)

**Current code:**
```julia
@user_info "Entrapment pairing Stage 2 complete (entrapment_target_idx created):"
@user_info "  Mapped $(n_mapped) entries to target indices"
@user_info "  Max target index: $(max_target_idx)"
@user_info "  Table size: $(n) rows"
@user_info "  Decoy column used: $(has_decoy_col ? "✅ YES" : "❌ NO")"

return df
```

**Change to:**
```julia
# Only report if we actually mapped entrapment entries
if n_mapped > 0
    @user_info "Entrapment pairing Stage 2 complete (entrapment_target_idx created):"
    @user_info "  Mapped $(n_mapped) entries to target indices"
    @user_info "  Max target index: $(max_target_idx)"
    @user_info "  Table size: $(n) rows"
    @user_info "  Decoy column used: $(has_decoy_col ? "✅ YES" : "❌ NO")"
end

return df
```

### Alternative Approach (More Conservative)
Instead of modifying the warning levels, we could check `entrapment_r` at the call site and only call `add_entrapment_indices!` when entrapment is enabled:

**In `src/Routines/BuildSpecLib.jl` (line 327):**
```julia
# Only process entrapment if it's enabled (entrapment_r > 0)
entrapment_r = get(params["fasta_digest_params"], "entrapment_r", 0)
if hasproperty(precursors_table, :entrapment_pair_id) && entrapment_r > 0
    println("   Adding entrapment target indices...")
    add_entrapment_indices!(precursors_table)
    ent_col = precursors_table.entrapment_target_idx
    max_entrap_target = all(ismissing, ent_col) ? 0 : Int64(maximum(skipmissing(ent_col)))
    println("   After add_entrapment_indices!: max entrapment_target_idx = $max_entrap_target")
    println("   Entrapment indices valid: $(max_entrap_target <= nrow(precursors_table) ? "✅ YES" : "❌ NO")")
end
```

This approach is cleaner because it prevents the function from being called at all when entrapment is disabled, avoiding all the warnings.

### Testing Requirements
Test with the following parameter configurations:

1. **Entrapment disabled** (`entrapment_r = 0`):
   - Verify NO entrapment messages appear
   - Verify library builds successfully
   - Check that `entrapment_pair_id` column is NOT created or is empty

2. **Entrapment enabled** (`entrapment_r = 1`):
   - Verify entrapment messages DO appear
   - Verify mapping is correct
   - Check that `entrapment_target_idx` is created properly

3. **Entrapment enabled with higher values** (`entrapment_r = 2, 3, 4`):
   - Verify correct number of entrapment sequences generated
   - Verify all mappings are correct
   - Check message counts match expected values

---

## Summary of Files to Modify

### Critical Files (Must Change)
1. `src/build/CLI/pioneer` - CLI help text
2. `src/Routines/GenerateParams.jl` - Argument parsing and function calls
3. `assets/example_config/defaultBuildLibParamsSimplified.json` - Simplified template comments
4. `assets/example_config/defaultBuildLibParams.json` - Full template comments
5. `src/Routines/BuildSpecLib.jl` - Conditional entrapment processing
6. `src/Routines/BuildSpecLib/chronologer/pair_decoys.jl` - Conditional message display

### Documentation Files (Should Update)
- Any user documentation referencing the `params-predict` command
- Examples in README.md or quickstart guides

---

## Implementation Order

1. **Phase 1: CLI Simplification** ✅ **COMPLETED**
   - ✅ Modify `pioneer` script help text
   - ✅ Update `GenerateParams.jl` argument parsing
   - ✅ Extract library name from output path
   - ✅ Test all CLI combinations
   - Commit: 58c87c73

2. **Phase 2: Template Comments**
   - Update simplified template
   - Update full template
   - Verify JSON parsing still works
   - Test parameter generation

3. **Phase 3: Entrapment Messages**
   - Add conditional check in `BuildSpecLib.jl`
   - Modify message levels in `pair_decoys.jl`
   - Test with `entrapment_r = 0`
   - Test with `entrapment_r > 0`

4. **Phase 4: Testing and Validation**
   - Run existing test suite
   - Add new tests for CLI changes
   - Verify backward compatibility (parameters still work)
   - Update documentation

---

## Backward Compatibility Notes

### CLI Changes
- **Breaking change:** Users must update their scripts/pipelines
- **Migration path:** Remove the explicit library name argument
- **Example migration:**
  - Old: `pioneer params-predict output/my_lib.poin my_lib fasta/`
  - New: `pioneer params-predict output/my_lib.poin fasta/`

### Parameter Files
- **Non-breaking:** Existing parameter files will continue to work
- Template changes only affect newly generated files
- Users can regenerate parameter files with improved comments

### Entrapment Messages
- **Non-breaking:** Only affects console output
- Functionality remains unchanged
- No impact on library files or search results

---

## Risk Assessment

### Low Risk
- Template comment changes (only affects readability)
- Entrapment message filtering (cosmetic improvements)

### Medium Risk
- CLI argument changes (requires user adaptation)
- Library name extraction logic (must handle edge cases)

### Mitigation Strategies
1. **Test extensively** with various path formats (relative, absolute, nested)
2. **Handle edge cases** (empty paths, special characters, trailing slashes)
3. **Provide clear error messages** if library name extraction fails
4. **Document migration** in release notes and changelog
5. **Keep GetBuildLibParams function signature** unchanged for programmatic use

---

## Future Enhancements (Out of Scope)

Consider for future versions:
1. Auto-detect whether path is a library name or full path
2. Support for custom library naming (optional `--name` flag)
3. Validate library path ends with `.poin` extension
4. Interactive mode for parameter generation
5. Configuration profiles for common use cases
