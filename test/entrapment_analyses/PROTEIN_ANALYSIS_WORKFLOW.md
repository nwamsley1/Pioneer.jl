# Protein Analysis Workflow in EntrapmentAnalysis Module

## Overview

The protein analysis workflow in the EntrapmentAnalysis module implements empirical false discovery rate (EFDR) calculations for protein-level proteomics data using entrapment sequences. This workflow differs significantly from the precursor-level analysis by focusing on protein identifications rather than individual peptide-spectrum matches.

## Key Concepts

### 1. Entrapment Proteins
- **Target proteins** (entrap_id = 0): Original proteins from the sample
- **Entrapment proteins** (entrap_id > 0): Known proteins spiked into the sample from different species
- Multiple entrapment groups allow for robust FDR estimation

### 2. Protein Pairing
Unlike precursor pairing which considers charge states and modifications, protein pairing is simpler but still considers target/decoy status:
- Proteins are paired by their name/identifier AND target/decoy status
- Target proteins are never paired with decoy proteins
- Each target protein (entrap_id=0) is paired with corresponding entrapment versions
- Round-robin assignment handles uneven group sizes
- Decoy proteins are filtered out before EFDR calculation

### 3. Scoring Architecture
The current implementation (`protein_scoring.jl`) uses a simplified approach:
- Proteins with the same name are inherently paired
- No need for explicit `entrap_pair_id` in many cases
- Original target scores are propagated to entrapment proteins for EFDR calculation
- Target/decoy status should be considered to avoid incorrect pairings

## Workflow Components

### 1. Data Structure Requirements

Protein results must contain:
```julia
- protein: String          # Protein identifier
- target: Bool            # Target/decoy status (true = target)
- entrap_id: UInt8        # Entrapment group (0 = target)
- file_name/ms_file_idx   # File identifier
- pg_score: Float32       # Protein group score
- global_pg_score: Float32 # Best score across files
- qval: Float32           # File-specific q-value
- global_qval: Float32    # Global q-value
- n_peptides: UInt32      # Peptide count
```

**Important**: Decoy proteins (target=false) are automatically filtered out before EFDR calculation.

### 2. Core Processing Steps

#### Step 1: Add Original Target Scores
```julia
add_original_target_protein_scores!(protein_results; score_col=:pg_score)
```
- For targets (entrap_id=0): Use their own score
- For entrapments (entrap_id>0): Get score from target with same protein name
- Handles both per-file and global analyses

#### Step 2: Create Global Results (if needed)
```julia
global_df = create_global_protein_results_df(protein_results; score_col=:global_pg_score)
```
- Selects best-scoring instance of each protein across all files
- Sets file identifier to indicate global analysis (0 or "global")

#### Step 3: Calculate EFDR
```julia
add_protein_efdr_columns!(protein_results; 
    method_types=[CombinedEFDR, PairedEFDR],
    score_qval_pairs=[(:pg_score, :qval)])
```
- Applies EFDR calculation methods directly to protein data
- No precursor indices or library needed

### 3. EFDR Calculation Methods

#### Combined EFDR
- Treats all entrapment groups equally
- Formula: EFDR = (n_entrapments / n_total) * r
- Where r is the ratio correction factor

#### Paired EFDR
- Uses paired target/entrapment scores
- More sophisticated handling of score distributions
- Better calibration for heterogeneous samples

### 4. Key Functions

#### `assign_protein_entrapment_pairs!`
- Groups proteins by name
- Assigns unique pair IDs using round-robin for uneven groups
- Handles multiple files per protein

#### `add_original_target_protein_scores!`
- Critical for EFDR calculation
- Maps target scores to entrapment proteins
- Supports multiple score columns

#### `create_global_protein_results_df`
- Aggregates multi-file results
- Selects best score per protein
- Maintains all metadata

#### `add_protein_efdr_columns!`
- Core EFDR calculation
- Supports multiple methods and score types
- Direct operation on protein data

#### `compare_protein_efdr_methods`
- Evaluates EFDR calibration
- Compares against q-value estimates
- Generates performance metrics

## Usage Patterns

### Basic Protein EFDR Analysis
```julia
using EntrapmentAnalysis

results = run_protein_efdr_analysis(
    "protein_groups_long.arrow";
    output_dir="protein_efdr_output",
    score_qval_pairs=[(:global_pg_score, :global_qval), (:pg_score, :qval)]
)
```

### Manual Step-by-Step Analysis
```julia
# 1. Load data
protein_df = DataFrame(Arrow.Table("protein_groups.arrow"))

# 2. Add original target scores
add_original_target_protein_scores!(protein_df; score_col=:pg_score)

# 3. Create global results
global_df = create_global_protein_results_df(protein_df; score_col=:global_pg_score)

# 4. Calculate EFDR
add_protein_efdr_columns!(global_df; 
    method_types=[CombinedEFDR, PairedEFDR],
    score_qval_pairs=[(:global_pg_score, :global_qval)])

# 5. Compare methods
comparison = compare_protein_efdr_methods(global_df, :global_qval, :global_pg_score)
```

## Key Differences from Precursor Analysis

1. **Simpler Pairing**: No charge/modification considerations
2. **Direct Data Access**: No precursor index mapping needed
3. **File Aggregation**: Built-in support for multi-file experiments
4. **Score Propagation**: Original target scores used for all calculations

## Output Structure

The analysis generates:
```
output_dir/
├── protein_efdr_analysis_report.md    # Summary report
├── protein_results_with_efdr.arrow    # Per-file results
├── global_protein_results_with_efdr.arrow  # Global results
├── efdr_comparison_pg_score.png/pdf   # Per-file plots
└── efdr_comparison_global_pg_score.png/pdf  # Global plots
```

## Important Considerations

1. **Score Direction**: Higher scores = better (unlike p-values)
2. **Missing Data**: Proteins without entrapment versions are excluded from EFDR
3. **File Handling**: Supports both single and multi-file experiments
4. **Calibration**: EFDR estimates should match actual FDR at each threshold

## Debugging Tips

1. Check entrap_id distribution:
   ```julia
   combine(groupby(df, :entrap_id), nrow)
   ```

2. Verify protein pairing:
   ```julia
   filter(x -> !ismissing(x.entrap_pair_id), df)
   ```

3. Examine score propagation:
   ```julia
   select(df, :protein, :entrap_id, :pg_score, :pg_score_original_target)
   ```

4. Monitor EFDR calibration:
   ```julia
   plot(comparison.threshold, comparison.combined_efdr, label="EFDR")
   plot!(comparison.threshold, comparison.combined_actual_fdr, label="Actual")
   ```

## Performance Optimization

- Dictionary-based lookups for O(1) score mapping
- Vectorized operations for EFDR calculations
- Minimal memory allocation through in-place updates
- Efficient groupby operations for large datasets

## Key Differences Between Protein and Precursor Analyses

### 1. Pairing Complexity

**Precursor Pairing** (`entrapment_pairing.jl`):
- Groups by: `base_pep_id`, `prec_charge`, `is_decoy`, `mod_key`
- Requires modification key extraction via `getModKey()`
- Considers charge states and modifications
- More complex grouping criteria

**Protein Pairing** (`protein_entrapment_pairing.jl`):
- Groups by: `protein` name AND `target` status
- No charge or modification considerations
- Simpler than precursor, but considers target/decoy
- Only proteins with same name AND same target status are paired
- Decoys filtered before pairing

### 2. Data Dependencies

**Precursor Analysis**:
- Requires `library_precursors` DataFrame
- Uses `precursor_idx` for mapping
- Needs `entrapment_group_id` from library
- Two-step process: library → results

**Protein Analysis**:
- Self-contained - `entrap_id` directly in results
- No library dependency for EFDR calculation
- Direct operation on protein data
- Single-step process

### 3. Score Propagation

**Precursor Scoring** (`scoring.jl`):
- Complex nested dictionary: `ms_file_idx → pair_id → target_info`
- Handles complement scores via `get_complement_score()`
- Tracks target row indices
- Uses `precursor_idx` for entrapment group lookup

**Protein Scoring** (`protein_scoring.jl`):
- Simple dictionary: `(file_id, protein_name) → score`
- Direct score mapping by protein name
- No complement score functionality needed
- Uses `entrap_id` directly from results
- Should be updated to include target status in key

### 4. Global Analysis Creation

**Precursor Global Results**:
- Groups by `precursor_idx`
- Sets `ms_file_idx = 0` for global
- Maintains precursor-level granularity

**Protein Global Results**:
- Groups by `protein` name
- Sets `ms_file_idx = 0` or `file_name = "global"`
- Protein-level aggregation

### 5. Column Requirements

**Precursor Analysis Required Columns**:
```
- precursor_idx
- ms_file_idx
- entrap_pair_id (added from library)
- score columns
```

**Protein Analysis Required Columns**:
```
- protein
- target (for filtering decoys)
- entrap_id
- file_name/ms_file_idx
- pg_score/global_pg_score
```

### 6. Implementation Complexity

**Precursor Analysis**:
- More complex due to peptide chemistry
- Requires careful handling of modifications
- Multiple grouping dimensions
- Library dependency adds complexity

**Protein Analysis**:
- Simplified by protein-level abstraction
- No chemical considerations
- Single grouping dimension (name)
- Self-contained data model

## Future Enhancements

1. Support for protein isoforms
2. Confidence intervals for EFDR estimates
3. Adaptive r parameter estimation
4. Integration with protein inference algorithms