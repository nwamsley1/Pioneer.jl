# Cleanup Plan: Remove Unused Functions After Pipeline Refactoring

## Overview
After refactoring ScoringSearch and MaxLFQSearch to use composable pipeline operations, many wrapper functions and monolithic implementations are no longer needed. This plan identifies functions to remove.

## Functions to Remove (42 total)

### 1. scoring_interface.jl (10 functions)
- [x] `process_psms_for_scoring` - Never called
- [x] `update_psms_with_protein_scores!` - Never called
- [x] `merge_protein_groups` (simple version) - Replaced by `merge_protein_groups_by_score`
- [x] `filter_protein_groups_by_qvalue` - Never called
- [x] `calculate_global_protein_scores` - Never called
- [x] `add_global_scores_to_psms!` - Never called
- [x] `calculate_protein_qvalues` - Never called
- [x] `generate_scoring_summary` - Never called
- [x] `validate_scoring_results` - Never called
- [x] `add_global_scores_to_psms_refs` - Never called

### 2. utils.jl (4 functions)
- [x] `get_psms_passing_qval` - Large commented-out function, not used
- [x] `build_protein_peptide_rows_for_file` - Never called
- [x] `plot_probit_decision_boundary` - Not called (commented out)
- [x] `calculate_qvalues_from_scores` - Never called

### 3. protein_inference_helpers.jl (8 functions - entire file)
This entire file can be removed as all functions were replaced by pipeline approach:
- [x] `build_protein_peptide_catalog` - Not used
- [x] `extract_peptide_protein_pairs` - Not used
- [x] `perform_file_protein_inference` - Not used
- [x] `update_psms_with_inference` - Not used
- [x] `create_protein_groups_from_psms` - Not used
- [x] `calculate_protein_features` - Only used in utils.jl, not in main pipeline
- [x] `process_single_file_protein_inference` - Not used
- [x] `write_protein_groups_arrow` - Only used in utils.jl, not in main pipeline

### 4. Functions mentioned in docs but don't exist
These were already removed or never implemented:
- `filter_psms_by_qvalue` - Doesn't exist in scoring_interface.jl
- `sort_and_filter_quant_tables_refs` - Doesn't exist in scoring_interface.jl

## Functions to Keep (still used)
- `merge_psm_files` - Called in ScoringSearch.jl
- `calculate_and_add_global_scores!` - Called in ScoringSearch.jl
- `update_psms_with_probit_scores_refs` - Called in ScoringSearch.jl
- `get_proteins_passing_qval_refs` - Called in ScoringSearch.jl
- `apply_probit_scores!` - Called via perform_protein_probit_regression
- `add_best_trace_indicator` - Used in pipeline operations
- `get_quant_necessary_columns` - Called in ScoringSearch.jl
- `get_best_traces` - Called in ScoringSearch.jl
- `get_qvalue_spline` - Called in ScoringSearch.jl
- `getProteinGroupsDict` - Called via apply_protein_inference
- `writeProteinGroups` - Called via apply_protein_inference
- `count_protein_peptides` - Called in ScoringSearch.jl
- `perform_protein_probit_regression` - Called in ScoringSearch.jl

## Implementation Steps

1. [x] Identify all deprecated functions
2. [x] Remove deprecated functions from scoring_interface.jl
3. [x] Remove deprecated functions from utils.jl
4. [x] Remove entire protein_inference_helpers.jl file
5. [x] Update export statements
6. [ ] Run tests to ensure nothing breaks
7. [ ] Update documentation

## Status
- [x] Plan created
- [x] Functions identified
- [x] Removal completed
- [ ] Tests passing
- [ ] Documentation updated

## Additional Work Done
- Added missing `calculate_protein_features` and `write_protein_groups_arrow` functions to utils.jl
- These functions were needed after removing protein_inference_helpers.jl
- Functions integrate with the new type system (ProteinKey, ProteinGroup, etc.)
- Replaced `get_proteins_passing_qval_refs` with pipeline operations in ScoringSearch.jl