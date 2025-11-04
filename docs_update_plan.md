# Documentation Update Plan: CLI Simplification Changes

## Overview
Update all documentation to reflect the simplified `params-predict` command that was implemented in the CLI simplification changes. The library name is now automatically derived from the output path, eliminating the need for a separate `lib_name` argument.

---

## Files Requiring Updates

### 1. Windows CLI Script (CRITICAL - Missing from Phase 1!)

**File:** `src/build/CLI/pioneer.bat`

**Issue:** The Windows batch script was not updated during Phase 1 implementation. It still shows the old three-argument syntax.

**Current State (Lines 108, 120):**
```batch
echo   params-predict ^<library_outpath^> ^<lib_name^> ^<fasta_path^> [--params-path ^<params_out_path^>]
...
echo   pioneer params-predict yeast.poin yeast fasta/ --params-path predict_params.json
```

**Required Changes:**
- **Line 108:** Remove `^<lib_name^>` argument
  ```batch
  echo   params-predict ^<library_outpath^> ^<fasta_path^> [--params-path ^<params_out_path^>]
  ```
- **Line 120:** Remove `yeast` library name from example
  ```batch
  echo   pioneer params-predict yeast.poin fasta/ --params-path predict_params.json
  ```

**Priority:** HIGH (Must be included in the PR - this is a bug that we missed!)

---

### 2. Main README

**File:** `README.md`

**Current State (Line 40):**
```bash
pioneer params-predict lib_dir lib_name fasta_dir --params-path=predict_params.json
```

**Required Change:**
```bash
pioneer params-predict lib_dir fasta_dir --params-path=predict_params.json
```

**Additional Context:**
- Update the Quick Start section (lines 36-49)
- Ensure consistency with the new CLI interface
- The example uses generic names (`lib_dir`, `fasta_dir`) which is good for documentation

**Priority:** HIGH (Main entry point for users)

---

### 3. User Guide - Quick Start

**File:** `docs/src/user_guide/quickstart.md`

**Current State (Lines 38-44):**
```bash
# Traditional approach with single FASTA directory
pioneer params-predict lib_dir lib_name fasta_dir --params-path=predict_params.json

# Or with new flexible FASTA input (mixing directories and files)
# Note: CLI support for mixed inputs requires editing the JSON parameter file
pioneer params-predict lib_dir lib_name fasta_dir --params-path=predict_params.json
# Then edit predict_params.json to set fasta_paths to include specific files
```

**Required Changes:**

**New version:**
```bash
# Generate library build parameters
pioneer params-predict lib_dir fasta_dir --params-path=predict_params.json

# Edit predict_params.json to customize:
# - Use single FASTA directory (as shown above)
# - Or specify multiple directories and/or specific FASTA files by editing fasta_paths
```

**Rationale:**
- Remove redundant `lib_name` argument
- Simplify the explanation - no need to show two separate commands
- Clarify that mixed inputs are configured by editing the JSON file
- The old comment about "CLI support requires editing JSON" is now misleading since the CLI always requires JSON editing for advanced features

**Additional Changes Needed:**
- Lines 54-59: Update the "Advanced FASTA Input" tip box to clarify that the CLI now only takes one path, but the Julia API still supports arrays
  ```julia
  # Julia API still supports mixed inputs directly
  params = GetBuildLibParams(out_dir, lib_name,
      ["/path/to/uniprot/", "/custom/proteins.fasta"])
  ```
- Note that CLI users derive `lib_name` from `out_dir` automatically

**Priority:** HIGH (Primary tutorial for new users)

---

### 4. User Guide - Installation

**File:** `docs/src/user_guide/installation.md`

**Current State (Lines 60-66):**
```julia
# Option 1: Single FASTA directory (backward compatible)
params = GetBuildLibParams(out_dir, lib_name, fasta_dir)
BuildSpecLib(params)

# Option 2: Flexible input - files and/or directories
params = GetBuildLibParams(out_dir, lib_name,
    ["/path/to/dir1", "/path/to/file.fasta", "/path/to/dir2"])
BuildSpecLib(params)
```

**Assessment:** NO CHANGES NEEDED
- This shows the Julia API, not the CLI
- The `GetBuildLibParams` function signature remains unchanged
- Users calling the function programmatically still provide all three arguments

**Current State (Line 91-92):**
```
2. Generate parameter files with `pioneer params-predict` or `pioneer params-search`,
   then edit them according to [Parameter Configuration](@ref "Parameter Configuration").
```

**Assessment:** NO CHANGES NEEDED
- Generic reference without showing specific syntax
- Still accurate

**Priority:** LOW (No changes required)

---

### 5. User Guide - Parameters

**File:** `docs/src/user_guide/parameters.md`

**Action Required:** Search for any references to `params-predict` syntax

**Expected:** Likely no changes needed if it only describes parameter meanings, not CLI syntax

**Priority:** LOW (Verify only)

---

### 6. Documentation Index

**File:** `docs/src/index.md`

**Current State (Line 65):**
```julia
GetBuildLibParams
```

**Assessment:** NO CHANGES NEEDED
- This is a docstring reference for API documentation
- The function signature hasn't changed

**Priority:** LOW (No changes required)

---

## Implementation Checklist

### Phase 1: Critical Fixes
- [x] Update `src/build/CLI/pioneer.bat` (Windows script - MISSED IN ORIGINAL PR!)
  - Remove `^<lib_name^>` from help text (line 108)
  - Remove `yeast` from example (line 120)
  - Completed in commit 046eeb71

### Phase 2: User-Facing Documentation
- [x] Update `README.md`
  - Change line 40 to remove `lib_name` argument
  - Verify example makes sense in context
  - Completed in commit 83d94f94

- [x] Update `docs/src/user_guide/quickstart.md`
  - Simplify command examples (lines 38-44)
  - Update "Advanced FASTA Input" tip box to clarify CLI vs Julia API (lines 54-59)
  - Remove confusing comments about "CLI support for mixed inputs"
  - Completed in next commit

### Phase 3: Verification
- [x] Search `docs/src/user_guide/parameters.md` for any command syntax examples
  - No matches found - no changes needed
- [x] Verify all documentation is consistent
  - Grep searches confirm no remaining old syntax in .md or .bat files
- [ ] Check for any other markdown files that might have examples

---

## Testing Documentation Updates

After making changes, verify:

1. **Consistency Check:**
   ```bash
   # Search for any remaining old syntax
   grep -r "params-predict.*lib_name" .
   grep -r "params-predict.*yeast.*yeast" .
   ```

2. **Documentation Build:**
   - Ensure documentation still builds correctly
   - Check that all internal links work

3. **Example Validation:**
   - Ensure all code examples are syntactically correct
   - Verify bash examples use correct syntax
   - Verify Julia examples use correct function signatures

---

## Key Points to Communicate

When updating documentation, ensure these points are clear:

1. **CLI Simplification:**
   - Library name is automatically derived from the output path
   - Reduces redundancy and typing
   - One less argument to remember

2. **Backward Compatibility:**
   - Julia API (`GetBuildLibParams`) is unchanged
   - Existing Julia scripts continue to work
   - Only CLI command syntax changed

3. **FASTA Input Flexibility:**
   - CLI takes a single path (file or directory)
   - For mixed inputs (multiple directories/files), edit the generated JSON
   - Julia API supports array input directly

4. **Migration:**
   - Old: `pioneer params-predict output.poin mylib fasta/`
   - New: `pioneer params-predict output.poin fasta/`
   - Just remove the middle argument

---

## Example Documentation Snippets

### For README.md Quick Start Section

```bash
# Generate library build parameters
pioneer params-predict library_output.poin fasta/ --params-path=predict_params.json

# Edit predict_params.json to configure library building options
# Then build the library
pioneer predict predict_params.json

# Convert raw data
pioneer convert-raw raw_dir

# Generate search parameters
pioneer params-search library_output.poin ms_data_dir results_dir --params-path=search_params.json

# Run the search
pioneer search search_params.json
```

### For Quickstart Guide

```bash
# 1. Generate library build parameters (library name derived from output path)
pioneer params-predict my_library.poin /path/to/fasta --params-path=build_params.json

# 2. Review and edit build_params.json
#    - Adjust digestion parameters (missed cleavages, modifications, etc.)
#    - For multiple FASTA sources, edit fasta_paths array in the JSON file
#    - Set calibration file if available for accurate m/z range detection

# 3. Build the spectral library
pioneer predict build_params.json
```

---

## Commit Strategy

**Recommended approach:** Single commit for documentation updates

**Commit message template:**
```
docs: Update documentation for simplified params-predict CLI

Update all user-facing documentation to reflect the simplified
params-predict command that removes the redundant lib_name argument.

Changes:
- Fix Windows batch script (pioneer.bat) - missed in original PR
- Update README.md Quick Start section
- Simplify quickstart.md examples and clarify CLI vs Julia API
- Ensure consistency across all documentation

The library name is now automatically derived from the output path
using basename(), reducing redundancy and simplifying the user experience.

Relates to commits 58c87c73 (CLI simplification implementation)
```

---

## Risk Assessment

### Low Risk
- README.md updates (main documentation)
- Quickstart guide updates
- Clarification of CLI vs API usage

### Medium Risk
- Windows batch script (pioneer.bat)
  - This was missed in the original implementation
  - Windows users would see incorrect help text
  - Need to verify it works correctly

### Mitigation
1. Test Windows script if possible (or note in PR that Windows testing needed)
2. Ensure all grep searches show no remaining old syntax
3. Have someone review the documentation changes for clarity

---

## Future Enhancements (Out of Scope)

Consider for future documentation improvements:
1. Add animated GIFs showing CLI workflow
2. Create a migration guide document for users updating from older versions
3. Add troubleshooting section for common CLI mistakes
4. Include more real-world examples with actual file paths
5. Create a CLI cheat sheet

---

## Notes

- The Julia API documentation (docstrings) should NOT be changed as the function signature remains the same
- Only CLI command examples and user-facing guides need updates
- Ensure consistency in terminology (use "library output path" consistently)
- The term "fasta_dir" vs "fasta_path" - the new argument is called "fasta_path" which is more accurate since it can be a file or directory
