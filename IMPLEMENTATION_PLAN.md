# Implementation Plan for Pioneer.jl Improvements

## Overview
This document outlines the implementation plan for two key improvements to Pioneer.jl:

1. **Enhanced FASTA Input Handling**: Support flexible input of directories and/or individual FASTA files with proper regex code mapping
2. **Reduced Koina API Warning Spam**: Convert retry warnings to debug messages and only show errors when entire batches fail

## Issue 1: Enhanced FASTA Input Handling for Library Building

### Current State
- **Location**: `src/Routines/GenerateParams.jl`, `GetBuildLibParams` function (lines 238-244)
- **Current Behavior**: Only handles directories - scans directory for `.fasta` and `.fasta.gz` files
- **Regex Handling**: Four regex vectors for different FASTA header components:
  - `fasta_header_regex_accessions`
  - `fasta_header_regex_genes` 
  - `fasta_header_regex_proteins`
  - `fasta_header_regex_organisms`

### Proposed Solution
Transform `fasta_dir` parameter to accept flexible input with proper regex code mapping:

#### Implementation Details

1. **Input Flexibility**:
   - Accept single directory, single file, or array of directories/files
   - Example inputs:
     - `"/path/to/fasta/dir"` (single directory)
     - `"/path/to/file.fasta"` (single file)  
     - `["/path/to/dir1", "/path/to/file.fasta"]` (mixed array)

2. **Regex Code Mapping Logic**:
   - **Single input + single regex set**: Apply regex to all FASTA files (whether from directory or single file)
   - **Multiple inputs + single regex set**: Apply same regex to all FASTA files from all inputs
   - **Multiple inputs + matching regex sets**: First regex applies to first input, second regex to second input, etc.
   - **Error cases**: Any mismatch between number of inputs and regex codes → error
     - More regex codes than inputs → error
     - Fewer regex codes than inputs (unless exactly 1 regex code) → error

3. **FASTA File Discovery**:
   - For directories: scan for `.fasta` and `.fasta.gz` files (existing logic)
   - For files: validate is FASTA file and use directly
   - Maintain file order within each directory for consistent behavior

4. **Code Changes**:
   - Modify `GetBuildLibParams` function parameter from `fasta_dir::String` to `fasta_inputs` (flexible type)
   - Update discovery logic to handle both files and directories
   - Implement regex code expansion/mapping logic
   - Update CLI argument parsing to accept flexible input

5. **JSON Output Structure**:
   - `fasta_paths`: All discovered FASTA file paths
   - `fasta_names`: Corresponding names (existing logic)
   - `fasta_header_regex_*`: Expanded regex arrays matching `fasta_paths` length

#### Regex Code Expansion Examples

**Example 1**: Single directory, default regex
```
Input: "/path/to/fastas/"
Regex: ["default_regex"] 
Directory contains: file1.fasta, file2.fasta
Result: Both files get "default_regex"
```

**Example 2**: Two directories, one regex code  
```
Input: ["/path/dir1/", "/path/dir2/"]
Regex: ["shared_regex"]
Result: All FASTA files from both directories get "shared_regex"
```

**Example 3**: Mixed directory/file with corresponding regex codes
```
Input: ["/path/dir1/", "/path/file.fasta"]
Regex: ["dir1_regex", "file_regex"]
Result: All files in dir1 get "dir1_regex", file.fasta gets "file_regex"
```

### Files to Modify
- `src/Routines/GenerateParams.jl` - Main implementation
- Update CLI help text and function docstrings

### Test Cases
1. **Single directory**: All FASTA files get same regex code
2. **Single file**: File gets specified regex code
3. **Multiple directories, single regex**: All files from all directories get same regex
4. **Multiple inputs, multiple regex**: Positional mapping works correctly
5. **Error cases**: Too many regex codes → proper error message
6. **Validation**: Non-existent paths, non-FASTA files → appropriate errors

---

## Issue 2: Reduce Koina API Warning Spam

### Current State
- **Location**: `src/Routines/BuildSpecLib/koina/koina_api.jl`
- **Current Behavior**: `make_koina_request` function logs warnings for every retry attempt (lines 41-44)
- **Problem**: Users see excessive warnings even when requests ultimately succeed after retries

### Proposed Solution
Convert individual retry warnings to debug messages and terminate program if entire batches fail:

#### Implementation Strategy
1. **Convert Retry Warnings to Debug Messages**:
   - Replace `@user_warn` calls in `make_koina_request` with `@debug_l2` messages
   - Users won't see retry attempts unless they enable debug level 2
   - Maintains debugging capability for developers

2. **Batch-Level Failure Handling**:
   - If an entire batch fails (all individual requests fail after max attempts), throw an error
   - Error should be informative and terminate library building process
   - No partial failure handling - either batch succeeds or program exits

#### Implementation Details

1. **Modify `make_koina_request` function**:
   - Replace lines 41-44: Change `@user_warn` to `@debug_l2`
   - Keep final error throwing behavior (line 38) unchanged
   - Individual request failures still propagate as exceptions

2. **Batch failure detection**:
   - `make_koina_batch_requests` already handles individual request failures properly
   - If a request fails completely, exception propagates and terminates batch
   - No additional batch-level logic needed - existing error handling is sufficient

3. **Error message improvement**:
   - When `make_koina_request` finally fails (line 38), provide more context
   - Include batch information and suggest possible causes (network, API limits, etc.)

#### Code Changes

**In `make_koina_request` function** (`src/Routines/BuildSpecLib/koina/koina_api.jl`):

```julia
# Replace lines 41-44:
if e isa KoinaRequestError
    @debug_l2 "Request failed (attempt $attempt): $(e.message)"
else
    @debug_l2 "Request failed (attempt $attempt): $(sprint(showerror, e))"
end
```

**Improve final error message** (line 49):
```julia
error("Koina API request failed after $max_attempts attempts. This may indicate network issues, API rate limiting, or service unavailability. Check your internet connection and try again later.")
```

### Files to Modify
- `src/Routines/BuildSpecLib/koina/koina_api.jl` - Convert warnings to debug messages

### Backward Compatibility
- Fully backward compatible - no parameter changes
- Existing error handling behavior preserved
- Only changes warning visibility level

### Test Cases
1. **All requests succeed on first try**: No output (existing behavior)
2. **Some requests succeed after retries**: No warnings visible to user, debug messages available
3. **Individual request fails completely**: Program terminates with informative error
4. **Network/API issues**: Clear error message helps user diagnose problem

---

## Implementation Approach

### Phase 1: FASTA Input Enhancement
1. Implement file vs directory detection logic
2. Add FASTA file validation
3. Update tests and documentation
4. Test with real-world use cases

### Phase 2: Koina Warning Reduction
1. Modify `make_koina_request` for optional warning suppression
2. Enhance `make_koina_batch_requests` with batch-level reporting
3. Update calling code to use new behavior
4. Test with various failure scenarios

### Phase 3: Integration and Testing
1. End-to-end testing of both improvements
2. Performance validation (ensure no regression)
3. Documentation updates
4. User acceptance testing

### Risk Mitigation
- **Backward Compatibility**: All changes use optional parameters with safe defaults
- **Testing**: Comprehensive test coverage for edge cases
- **Rollback Plan**: Both improvements are self-contained and can be reverted independently

### Success Criteria
1. **FASTA Enhancement**: Users can successfully use both directory and single-file inputs
2. **Koina Warnings**: Users see significantly fewer warnings during library building, with warnings only appearing for actual failures
3. **No Regressions**: Existing workflows continue to work unchanged
4. **Performance**: No measurable performance impact from changes

---

## Questions for User Review
1. Should we also support wildcards for FASTA file patterns (e.g., `*.fasta`)?
2. For Koina warnings, would you prefer a summary at the end of each batch or no output unless there's a true failure?
3. Are there any other file types besides `.fasta` and `.fasta.gz` that should be supported?
4. Should we add progress indicators showing retry statistics without the warning level?