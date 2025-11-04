# Plan: Fix Old Argument Syntax in GitHub Actions Workflows

## Problem Statement

The CLI simplification implemented in commits 58c87c73, 7bfdca8f, ab6232c1, 212bdd1d, and a7d9baf5 removed the redundant `lib_name` argument from the `params-predict` command. The library name is now automatically derived from the output path using `basename()`.

**Old syntax:**
```bash
pioneer params-predict <out_dir> <lib_name> <fasta_path> [--params-path <path>]
```

**New syntax:**
```bash
pioneer params-predict <out_dir> <fasta_path> [--params-path <path>]
```

However, two GitHub Actions workflow files still use the old three-argument syntax, causing build failures with the error:
```
too many arguments
usage: <PROGRAM> [--params-path PARAMS_PATH] [--full] [out_dir]
                 [fasta_path]
```

## Affected Files

### 1. `.github/workflows/build_app_linux.yml`

**Location:** Line 142
**Current Code:**
```bash
"$pioneer_root/pioneer" params-predict "$dummy" dummy "$dummy" --params-path "$dummy/params.json"
```

**Issue:** Uses three arguments: `"$dummy"` (out_dir), `dummy` (lib_name), `"$dummy"` (fasta_path)

**Context:** This is in the "Pre-download artifacts" step that runs the params-predict command to trigger artifact downloads (IntelOpenMP, oneTBB, MKL, LightGBM) during the Linux build process.

---

### 2. `.github/workflows/build_app_windows.yml`

**Location:** Line 138
**Current Code:**
```powershell
& "$pioneerRoot\pioneer.bat" params-predict $dummy dummy $dummy --params-path "$dummy\params.json"
```

**Issue:** Uses three arguments: `$dummy` (out_dir), `dummy` (lib_name), `$dummy` (fasta_path)

**Context:** This is in the "Pre-download artifacts" step that runs the params-predict command to trigger artifact downloads during the Windows build process.

---

### 3. `.github/workflows/build_app_macos.yml`

**Status:** ✅ No changes needed - this workflow does not have a "Pre-download artifacts" step

---

## Files Verified as Clean

The following locations were checked and confirmed NOT to have old syntax:

✅ Test files (`test/**/*.jl`, `test/**/*.sh`) - No params-predict usage found
✅ Shell scripts (`**/*.sh`, `**/*.bash`) - No params-predict usage found
✅ Python scripts (`**/*.py`) - No params-predict usage found
✅ Documentation (`README.md`, `docs/**/*.md`) - Already updated in commits 83d94f94 and 84dfecbd
✅ CLI scripts (`src/build/CLI/pioneer`, `src/build/CLI/pioneer.bat`) - Already updated in commits 58c87c73 and 046eeb71

---

## Implementation Plan

### Step 1: Fix Linux Workflow

**File:** `.github/workflows/build_app_linux.yml`
**Line:** 142

**Change from:**
```bash
"$pioneer_root/pioneer" params-predict "$dummy" dummy "$dummy" --params-path "$dummy/params.json"
```

**Change to:**
```bash
"$pioneer_root/pioneer" params-predict "$dummy" "$dummy" --params-path "$dummy/params.json"
```

**Explanation:** Remove the middle `dummy` argument (lib_name). The library name will be automatically derived from `"$dummy"` using `basename()`.

---

### Step 2: Fix Windows Workflow

**File:** `.github/workflows/build_app_windows.yml`
**Line:** 138

**Change from:**
```powershell
& "$pioneerRoot\pioneer.bat" params-predict $dummy dummy $dummy --params-path "$dummy\params.json"
```

**Change to:**
```powershell
& "$pioneerRoot\pioneer.bat" params-predict $dummy $dummy --params-path "$dummy\params.json"
```

**Explanation:** Remove the middle `dummy` argument (lib_name). The library name will be automatically derived from `$dummy` using `basename()`.

---

## Testing Strategy

### Pre-Implementation Verification

1. **Confirm the error** (already done):
   - Error message: "too many arguments"
   - Source: GitHub Actions build logs for Linux and Windows

2. **Verify current behavior:**
   - Read the workflow files to confirm exact syntax
   - Confirm that macOS workflow doesn't need changes

### Post-Implementation Verification

1. **Local syntax check:**
   ```bash
   # Check that no old syntax remains in workflows
   grep -n "params-predict.*dummy.*dummy.*dummy" .github/workflows/*.yml
   # Should return no matches
   ```

2. **GitHub Actions validation:**
   - Push changes to a test branch
   - Verify that Linux build workflow completes successfully
   - Verify that Windows build workflow completes successfully
   - Confirm that artifacts are still pre-downloaded correctly

3. **Artifact download verification:**
   - Check build logs for messages like:
     ```
     Downloading artifact: IntelOpenMP
     Downloading artifact: oneTBB
     Downloading artifact: MKL
     [ Info: lib_lightgbm found in system dirs!
     ```
   - Verify no "too many arguments" error

---

## Risk Assessment

### Low Risk
- **Why:** These are straightforward one-line changes
- **Impact:** Only affects CI/CD build process, not user-facing functionality
- **Reversibility:** Easy to revert if issues arise
- **Testing:** Can be validated immediately in GitHub Actions

### Potential Issues
1. **None expected** - The changes are simple argument removal
2. The `params-predict` command should work identically with the new syntax
3. The dummy directory structure remains the same

---

## Commit Strategy

**Single commit for both workflow files:**

```
fix(ci): Update params-predict syntax in build workflows

Remove redundant lib_name argument from params-predict commands
in Linux and Windows build workflows. The library name is now
automatically derived from the output path using basename().

This fixes the "too many arguments" error in the pre-download
artifacts step.

Changes:
- build_app_linux.yml line 142: Remove middle 'dummy' argument
- build_app_windows.yml line 138: Remove middle 'dummy' argument

Fixes compatibility with CLI simplification from commits:
- 58c87c73 (CLI implementation)
- 7bfdca8f (backward compatibility removal)
- 046eeb71 (Windows batch script fix)
```

---

## Verification Checklist

After making changes, verify:

- [ ] `.github/workflows/build_app_linux.yml` line 142 updated
- [ ] `.github/workflows/build_app_windows.yml` line 138 updated
- [ ] No remaining three-argument params-predict calls:
  ```bash
  grep -r "params-predict.*dummy.*dummy.*dummy" .github/
  # Should return no matches (except in this plan.md)
  ```
- [ ] Syntax is consistent with documentation examples
- [ ] Commit message follows conventional commit format
- [ ] Push to test branch first
- [ ] Verify Linux build succeeds in GitHub Actions
- [ ] Verify Windows build succeeds in GitHub Actions
- [ ] Merge to develop branch after successful testing

---

## Related Issues

This plan addresses the root cause of the GitHub Actions failure reported by the user:
```
too many arguments
usage: <PROGRAM> [--params-path PARAMS_PATH] [--full] [out_dir]
                 [fasta_path]
```

The error occurred because the workflows were still using the old three-argument syntax after the CLI was simplified to use two arguments.

---

## Timeline

1. **Implementation:** 5 minutes (two simple edits)
2. **Local verification:** 2 minutes (grep searches)
3. **Commit:** 2 minutes
4. **Push and test:** 10-20 minutes (GitHub Actions runtime)
5. **Total:** ~30 minutes

---

## Success Criteria

✅ Linux build workflow completes without "too many arguments" error
✅ Windows build workflow completes without "too many arguments" error
✅ Artifacts are successfully pre-downloaded in both workflows
✅ No regressions in build process
✅ Consistent syntax across all documentation and workflows

---

## Notes

- The macOS workflow (`build_app_macos.yml`) does not need changes as it doesn't have a pre-download artifacts step
- The dummy directory is still created and removed correctly with the new syntax
- The `--params-path` flag remains unchanged
- This completes the CLI simplification work started in the `simplify-buildspeclib-params` branch
