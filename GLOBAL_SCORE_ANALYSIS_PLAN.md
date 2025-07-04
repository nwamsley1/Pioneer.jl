# Global Score Analysis Implementation Plan

## Overview
Implement separate dataframe handling for per-file and global score analyses in the EntrapmentAnalysis module.

## Problem
The current `add_original_target_scores!` function processes scores on a per-file basis (using `ms_file_idx`), but global probability scores should be handled across all files by selecting the best score for each precursor.

## Solution
Create two separate dataframes at the beginning of the analysis:
1. `results_df` - Original data for per-file analysis (prec_prob scores)
2. `global_results_df` - Transformed data for global analysis (global_prob scores)

## Implementation Steps

### Step 1: Create transformation function
- [x] Add `create_global_results_df` function to `scoring.jl`
  - Groups by precursor_idx
  - Keeps row with maximum score for each precursor
  - Sets ms_file_idx to 0 for all rows
  - Returns new dataframe for global analysis
  - Creates deep copy to avoid aliasing issues

### Step 2: Update API function
- [x] Modify `run_efdr_analysis` in `api.jl`
  - Create both dataframes after loading data
  - Process per-file scores with original dataframe
  - Process global scores with global dataframe
  - Handle output appropriately
  - Separate comparison and calibration calculations
  - Save both dataframes if needed

### Step 3: Update examples
- [x] Update `basic_usage.jl` to demonstrate:
  - Creating global results dataframe
  - Processing both types of scores
  - Keeping separate results for per-file and global analyses
  - Saving both dataframes separately

### Step 4: Testing
- [x] Test with example data
- [x] Verify EFDR calculations are correct for both score types
- [x] Check that output is properly formatted
- [x] Created test_global_analysis.jl to verify implementation

## Code Snippets

### create_global_results_df function
```julia
function create_global_results_df(prec_results::DataFrame; score_col::Symbol=:global_prob)
    # Group by precursor_idx
    grouped = groupby(prec_results, :precursor_idx)
    
    # For each group, keep only the row with maximum score
    global_df = combine(grouped) do group
        best_idx = argmax(group[!, score_col])
        return group[best_idx:best_idx, :]
    end
    
    # Set all ms_file_idx to 0 to indicate global analysis
    global_df[!, :ms_file_idx] = 0
    
    return global_df
end
```

## Commits Plan
1. "Add create_global_results_df function for global score analysis"
2. "Update API to handle separate per-file and global analyses"
3. "Update examples to demonstrate global score analysis"
4. "Complete global score analysis implementation"

## Notes
- Setting `ms_file_idx = 0` clearly indicates global analysis
- Original functions work unchanged on their respective dataframes
- This approach maintains backward compatibility while adding new functionality