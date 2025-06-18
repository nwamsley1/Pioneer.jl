# Plan to Remove Stub Implementations and Fix Conflicts

## Problem Analysis

There are serious issues with the current implementation that need to be addressed:

1. **Duplicate Function Definition**: `getProteinGroupsDict` is defined in both:
   - `FileOperations.jl` (lines 18-44) - stub implementation
   - `ScoringSearch/utils.jl` (line 775+) - real implementation

2. **Function Conflict**: Both ScoringSearch and MaxLFQSearch include FileOperations.jl, which means they're using the stub version instead of the real one.

3. **Testing Issue**: The tests are passing with stub implementations, not the real algorithms.

## Solution Plan

### Step 1: Remove Stub Implementation
- Delete the stub `getProteinGroupsDict` function from FileOperations.jl (lines 17-44)
- Update the comment to indicate we're importing from the real location

### Step 2: Fix Module Dependencies
Since FileOperations needs to call getProteinGroupsDict from ScoringSearch/utils.jl, we need to handle the circular dependency. The recommended approach is to pass the function as a parameter:

```julia
function apply_protein_inference(psm_ref::PSMFileReference, 
                               output_path::String,
                               precursors,
                               protein_inference_fn::Function;  # Pass the real function
                               min_peptides::Int=2)
```

This way, ScoringSearch can pass its own `getProteinGroupsDict` function when calling `apply_protein_inference`.

### Step 3: Update Call Sites
- Update ScoringSearch to pass the real `getProteinGroupsDict` when calling `apply_protein_inference`
- Ensure the function signature matches what's expected

### Step 4: Update Tests
- Ensure tests are using the real implementation
- Run tests to verify everything still works with real algorithms

### Step 5: Check for Other Stubs
- Search for any other stub implementations that might cause similar issues
- Ensure all algorithm wrappers call the real implementations
- **Found**: `apply_maxlfq` function has the LFQ call commented out (lines 704-707)
- Need to import and use the real LFQ function from `src/utils/maxLFQ.jl`

## Implementation Steps

1. Remove stub from FileOperations.jl
2. Add `protein_inference_fn` parameter to `apply_protein_inference`
3. Update the function body to use the passed function
4. Update any calls to `apply_protein_inference` to pass the real function
5. Run tests to ensure everything works correctly

## Benefits

- Avoids circular dependencies
- Maintains clean separation of concerns
- Makes the dependency explicit
- Allows for testing with mock functions if needed
- Ensures the real implementation is always used in production