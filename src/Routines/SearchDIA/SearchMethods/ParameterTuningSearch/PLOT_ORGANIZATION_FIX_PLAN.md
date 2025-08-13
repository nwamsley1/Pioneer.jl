# Fix Plan: QC Plot Organization

## Current Issue
The combined/merged PDF files are being placed in the top-level `qc_plots/` directory instead of in their respective subdirectories.

## Current Structure (Wrong)
```
qc_plots/
├── rt_alignment_plots/
│   ├── file1.pdf
│   ├── file2.pdf
│   └── ...
├── mass_error_plots/
│   ├── file1.pdf
│   ├── file2.pdf
│   └── ...
├── rt_alignment_combined.pdf    ❌ Should be in rt_alignment_plots/
└── mass_error_combined.pdf      ❌ Should be in mass_error_plots/
```

## Expected Structure (Original Behavior)
```
qc_plots/
├── rt_alignment_plots/
│   ├── file1.pdf
│   ├── file2.pdf
│   ├── ...
│   └── rt_alignment_combined.pdf    ✓ Inside subdirectory
└── mass_error_plots/
    ├── file1.pdf
    ├── file2.pdf
    ├── ...
    └── mass_error_combined.pdf      ✓ Inside subdirectory
```

## Investigation Needed

### 1. Check Original Implementation
Need to verify where the original code placed the combined PDFs. Look for:
- Path construction for merged files
- Folder references used

### 2. Current Code Analysis

Looking at `summarize_results!`:
```julia
# Current (incorrect):
qc_plots_folder = getQcPlotsFolder(results)  # This is "qc_plots/"
rt_merged_path = joinpath(qc_plots_folder, "rt_alignment_combined.pdf")  # Goes to qc_plots/
mass_merged_path = joinpath(qc_plots_folder, "mass_error_combined.pdf")   # Goes to qc_plots/
```

Should be:
```julia
# Fixed:
rt_merged_path = joinpath(rt_plots_folder, "rt_alignment_combined.pdf")     # Goes to subdirectory
mass_merged_path = joinpath(mass_error_plots_folder, "mass_error_combined.pdf")  # Goes to subdirectory
```

### 3. Verify Folder Variables
The code already has:
- `rt_plots_folder` = getRtAlignPlotFolder(search_context)
- `mass_error_plots_folder` = getMassErrPlotFolder(search_context)

These should point to the subdirectories.

## Root Cause
The merged file paths are using `qc_plots_folder` (top level) instead of the specific subdirectory folders (`rt_plots_folder` and `mass_error_plots_folder`).

## Fix Implementation

### Change Required in summarize_results!

**Current Code (lines ~750-785):**
```julia
# Merge RT alignment plots
rt_merged_path = joinpath(qc_plots_folder, "rt_alignment_combined.pdf")
rt_pdf_files = filter(f -> endswith(f, ".pdf"), readdir(rt_plots_folder, join=true))
if !isempty(rt_pdf_files)
    merge_pdfs(rt_pdf_files, rt_merged_path)
    @info "Merged $(length(rt_pdf_files)) RT alignment plots to $rt_merged_path"
else
    @info "No RT alignment plots to merge"
end

# Merge mass error plots
mass_merged_path = joinpath(qc_plots_folder, "mass_error_combined.pdf")
mass_pdf_files = filter(f -> endswith(f, ".pdf"), readdir(mass_error_plots_folder, join=true))
if !isempty(mass_pdf_files)
    merge_pdfs(mass_pdf_files, mass_merged_path)
    @info "Merged $(length(mass_pdf_files)) mass error plots to $mass_merged_path"
else
    @info "No mass error plots to merge"
end
```

**Fixed Code:**
```julia
# Merge RT alignment plots into their subdirectory
rt_merged_path = joinpath(rt_plots_folder, "rt_alignment_combined.pdf")  # Changed: use rt_plots_folder
rt_pdf_files = filter(f -> endswith(f, ".pdf"), readdir(rt_plots_folder, join=true))
# Filter out any existing combined file to avoid including it in merge
rt_pdf_files = filter(f -> !endswith(f, "_combined.pdf"), rt_pdf_files)
if !isempty(rt_pdf_files)
    merge_pdfs(rt_pdf_files, rt_merged_path)
    @info "Merged $(length(rt_pdf_files)) RT alignment plots to $rt_merged_path"
else
    @info "No RT alignment plots to merge"
end

# Merge mass error plots into their subdirectory
mass_merged_path = joinpath(mass_error_plots_folder, "mass_error_combined.pdf")  # Changed: use mass_error_plots_folder
mass_pdf_files = filter(f -> endswith(f, ".pdf"), readdir(mass_error_plots_folder, join=true))
# Filter out any existing combined file to avoid including it in merge
mass_pdf_files = filter(f -> !endswith(f, "_combined.pdf"), mass_pdf_files)
if !isempty(mass_pdf_files)
    merge_pdfs(mass_pdf_files, mass_merged_path)
    @info "Merged $(length(mass_pdf_files)) mass error plots to $mass_merged_path"
else
    @info "No mass error plots to merge"
end
```

## Additional Improvements

### 1. Filter Existing Combined Files
When re-running, we should exclude any existing `*_combined.pdf` files from the merge to avoid including old combined files in new ones.

### 2. Consistent Naming
Ensure the combined files have consistent names that are easy to identify and exclude.

## Benefits of Fix

1. **Organized Structure**: Combined plots stay with their individual plots
2. **Easier Navigation**: Users can find all related plots in one folder
3. **Cleaner Top Directory**: The qc_plots folder isn't cluttered with combined files
4. **Consistent with Original**: Matches the original implementation's behavior

## Testing

After implementation:
1. Run a parameter tuning search
2. Verify individual plots are in subdirectories
3. Verify combined plots are in the same subdirectories as individual plots
4. Check that re-running doesn't include old combined files in new merges