# Convergence Logic and merge_pdfs Fix Plan

## Issue 1: Convergence Checking Logic

### Current Problem
The refactored code checks convergence after all three strategies, but the original implementation only checked after:
- Strategy 1: After adding more scans
- Strategy 3: After correcting mass bias
- NOT after Strategy 2: Expanding tolerance was always followed by bias correction

### Original Logic Flow
```
For each iteration:
  1. Add scans (doubling strategy)
  2. Try Strategy 1 (collect with current params) → CHECK CONVERGENCE
  3. If not converged: Try Strategy 2 (expand tolerance) → NO CHECK
  4. Always follow with Strategy 3 (adjust bias) → CHECK CONVERGENCE
```

### Fix Required
Modify the `execute_strategy` and main loop to:
1. Only check convergence after strategies 1 and 3
2. Strategy 2 should always be followed by Strategy 3
3. If Strategy 1 converges, skip Strategies 2 and 3

## Issue 2: merge_pdfs Error

### Current Problem
```julia
merge_pdfs(rt_plots_folder, rt_merged_path)  # WRONG - passing directory path
```

`merge_pdfs` expects a Vector of file paths, not a directory path.

### Fix Required
```julia
# Get all PDF files in the folder
pdf_files = filter(f -> endswith(f, ".pdf"), readdir(rt_plots_folder, join=true))
if !isempty(pdf_files)
    merge_pdfs(pdf_files, rt_merged_path)
end
```

## Implementation Plan

### Step 1: Fix the main convergence loop in process_file!

Change from:
```julia
# Try three strategies in sequence
for strategy in 1:3
    psms, mass_err_model, ppm_errs = execute_strategy(...)
    
    # Check convergence
    if check_and_store_convergence!(...)
        converged = true
        break
    end
end
```

To:
```julia
# Strategy 1: Try with current parameters
psms, mass_err_model, ppm_errs = execute_strategy(1, ...)
if check_and_store_convergence!(...)
    converged = true
else
    # Strategy 2: Expand tolerance (no convergence check)
    psms, mass_err_model, ppm_errs = execute_strategy(2, ...)
    
    # Strategy 3: Adjust bias (always follows Strategy 2)
    psms, mass_err_model, ppm_errs = execute_strategy(3, ...)
    if check_and_store_convergence!(...)
        converged = true
    end
end
```

### Step 2: Fix merge_pdfs calls in summarize_results!

Change from:
```julia
merge_pdfs(rt_plots_folder, rt_merged_path)
merge_pdfs(mass_error_plots_folder, mass_merged_path)
```

To:
```julia
# Get all PDF files in RT plots folder
rt_pdf_files = filter(f -> endswith(f, ".pdf"), 
                      readdir(rt_plots_folder, join=true))
if !isempty(rt_pdf_files)
    merge_pdfs(rt_pdf_files, rt_merged_path)
    @info "Merged $(length(rt_pdf_files)) RT alignment plots to $rt_merged_path"
else
    @info "No RT alignment plots to merge"
end

# Get all PDF files in mass error plots folder
mass_pdf_files = filter(f -> endswith(f, ".pdf"), 
                        readdir(mass_error_plots_folder, join=true))
if !isempty(mass_pdf_files)
    merge_pdfs(mass_pdf_files, mass_merged_path)
    @info "Merged $(length(mass_pdf_files)) mass error plots to $mass_merged_path"
else
    @info "No mass error plots to merge"
end
```

### Step 3: Simplify execute_strategy function

Since Strategy 3 now needs the mass error model from Strategy 2, we should:
1. Remove the redundant PSM collection in Strategy 3
2. Pass the mass error model from Strategy 2 to Strategy 3

Or alternatively, keep execute_strategy as is but change how we call it in the main loop.

## Benefits

1. **Correct convergence logic**: Matches the original implementation's behavior
2. **Fewer unnecessary checks**: Don't check convergence after tolerance expansion
3. **Proper PDF merging**: Correctly passes file lists to merge_pdfs
4. **Better error handling**: Check for empty file lists before merging

## Testing

After implementation:
1. Verify convergence happens only after strategies 1 and 3
2. Confirm Strategy 2 is always followed by Strategy 3
3. Test PDF merging with multiple files
4. Test with empty folders (no PDFs to merge)