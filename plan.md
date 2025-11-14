# Dual RT Index Second Pass Search Implementation Plan

## Overview
Modify SecondPassSearch to run twice per file with different rt_index configurations:
1. **First search**: targets-only rt_index → collect PSMs with only target precursors
2. **Second search**: targets+decoys rt_index → collect PSMs, filter to keep only decoys
3. **Combine**: Concatenate target PSMs with decoy PSMs

## Motivation
This approach separates target and decoy searches to prevent targets and decoys from competing in the same scans during deconvolution, potentially improving target identification while maintaining proper decoy statistics.

## Implementation Details

### Location
**File**: `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/SecondPassSearch.jl`
**Function**: `process_file!` (lines 297-402)

### Current Flow
```julia
# Line 312-314: Build rt_index from Arrow table
rt_index = buildRtIndex(
    DataFrame(Arrow.Table(getRtIndex(getMSData(search_context), ms_file_idx))),
    bin_rt_size = 0.1)

# Line 317-324: Perform search once
psms = perform_second_pass_search(
    spectra, rt_index, search_context, params, ms_file_idx, MS2CHROM())
```

### Proposed Flow
```julia
# 1. Load the full rt_index Arrow table
rt_df_full = DataFrame(Arrow.Table(getRtIndex(getMSData(search_context), ms_file_idx)))

# 2. Create targets-only rt_index
is_decoy = getIsDecoy(getPrecursors(getSpecLib(search_context)))
targets_mask = .!is_decoy[rt_df_full[!, :precursor_idx]]
rt_df_targets = rt_df_full[targets_mask, :]
rt_index_targets = buildRtIndex(rt_df_targets, bin_rt_size = 0.1)

# 3. First search: targets only
psms_targets = perform_second_pass_search(
    spectra, rt_index_targets, search_context, params, ms_file_idx, MS2CHROM())

# 4. Create full rt_index (targets + decoys)
rt_index_full = buildRtIndex(rt_df_full, bin_rt_size = 0.1)

# 5. Second search: all precursors
psms_all = perform_second_pass_search(
    spectra, rt_index_full, search_context, params, ms_file_idx, MS2CHROM())

# 6. Filter to keep only decoys from second search
is_decoy_lookup = getIsDecoy(getPrecursors(getSpecLib(search_context)))
decoy_mask = [is_decoy_lookup[pid] for pid in psms_all[!, :precursor_idx]]
psms_decoys = psms_all[decoy_mask, :]

# 7. Concatenate targets and decoys
psms = vcat(psms_targets, psms_decoys)
```

### MS1 Search Handling
The same dual-index approach needs to be applied to the MS1 search (lines 325-391):
- Build two rt_indices as above
- Run MS1 search twice (once with targets-only, once with full)
- Filter and concatenate MS1 PSMs the same way

## Code Changes Required

### 1. Modify `process_file!` in SecondPassSearch.jl
- Replace single rt_index building with dual rt_index approach
- Run `perform_second_pass_search` twice (once for MS2, once for MS1 if enabled)
- Filter decoys from second search
- Concatenate results

### 2. Helper Functions (if needed)
Create utility functions to avoid code duplication:
```julia
function build_dual_rt_indices(rt_df_full::DataFrame, precursors)
    is_decoy = getIsDecoy(precursors)
    targets_mask = .!is_decoy[rt_df_full[!, :precursor_idx]]
    rt_df_targets = rt_df_full[targets_mask, :]

    rt_index_targets = buildRtIndex(rt_df_targets, bin_rt_size = 0.1)
    rt_index_full = buildRtIndex(rt_df_full, bin_rt_size = 0.1)

    return rt_index_targets, rt_index_full
end

function filter_psms_by_decoy_status(psms::DataFrame, precursors, keep_decoys::Bool)
    if isempty(psms)
        return psms
    end
    is_decoy_lookup = getIsDecoy(precursors)
    mask = [is_decoy_lookup[pid] for pid in psms[!, :precursor_idx]]
    return psms[keep_decoys ? mask : .!mask, :]
end
```

## Potential Downstream Issues

### 1. Sorting and Data Order
**Issue**: The concatenated PSMs may not be sorted correctly.
**Solution**: Ensure `vcat(psms_targets, psms_decoys)` is sorted by the expected column (likely `:scan_idx` or `:rt`).
**Location**: After concatenation, add `sort!(psms, :scan_idx)` or whatever column is expected downstream.

### 2. PSM Indexing and Uniqueness
**Issue**: Downstream code may assume unique (scan_idx, precursor_idx) pairs.
**Impact**: With dual searches, the same scan may be searched twice (once with targets, once with decoys).
**Assessment**: This should be fine since we filter to keep only targets from first search and only decoys from second search, so there should be no overlap in (scan_idx, precursor_idx) pairs.

### 3. Statistical Independence
**Issue**: Are decoys truly independent if they're selected based on same RT windows?
**Assessment**: This is actually BETTER than current approach - decoys are selected based on their own RT characteristics, not artificially included in target-rich scans.

### 4. FDR Calculation
**Issue**: FDR calculations in ScoringSearch assume targets and decoys were searched under identical conditions.
**Assessment**: This change means targets and decoys MAY be searched in slightly different scans due to different rt_index composition. However, both are still searched according to their RT characteristics, so FDR should remain valid.
**Monitoring**: Check TDC (target-decoy competition) plots in ScoringSearch outputs to ensure proper decoy distribution.

### 5. Memory Usage
**Issue**: Running search twice doubles memory usage and computation time.
**Impact**: May be significant for large files.
**Mitigation**: Consider monitoring memory usage and potentially implementing batch processing if needed.

### 6. Empty DataFrame Handling
**Issue**: If first search returns no targets OR second search returns no decoys, concatenation may fail.
**Solution**: Add checks for empty DataFrames:
```julia
if isempty(psms_targets)
    psms = psms_decoys
elseif isempty(psms_decoys)
    psms = psms_targets
else
    psms = vcat(psms_targets, psms_decoys)
end
```

### 7. Column Schema Consistency
**Issue**: `vcat` requires identical column schemas. If searches return different columns, concatenation will fail.
**Assessment**: Both searches use identical parameters and functions, so schemas should match.
**Verification**: Test with small dataset first to confirm schema consistency.

### 8. Chromatogram Integration (IntegrateChromatogramSearch)
**Issue**: IntegrateChromatogramSearch may expect PSMs from both targets and decoys in a specific structure.
**Assessment**: Should be fine since we're just changing HOW PSMs are collected, not their final structure.
**Monitoring**: Verify chromatogram integration works correctly for both target and decoy PSMs.

### 9. Match Between Runs (MBR)
**Issue**: MBR uses rt_indices to transfer identifications between runs.
**Assessment**: MBR happens in IntegrateChromatogramSearch, which uses the rt_indices from FirstPassSearch (not SecondPassSearch).
**Impact**: No impact expected, but verify MBR results are consistent.

### 10. Runtime Performance
**Issue**: Doubling the number of searches per file will roughly double SecondPassSearch runtime.
**Assessment**: SecondPassSearch is already one of the most expensive steps. This will increase total pipeline time.
**Mitigation**: Consider parallelization opportunities or making this feature optional via parameters.

## Testing Strategy

### 1. Unit Testing
- Test `build_dual_rt_indices` with known target/decoy distributions
- Test `filter_psms_by_decoy_status` with sample PSM DataFrames
- Verify empty DataFrame handling

### 2. Integration Testing
- Run on small test dataset (e.g., `data/ecoli_test/`)
- Compare PSM counts: targets from search 1 + decoys from search 2 should equal total PSMs
- Verify no duplicate (scan_idx, precursor_idx) pairs
- Check that targets and decoys maintain proper ratio

### 3. Validation Checks
- Confirm all target PSMs come from targets-only search
- Confirm all decoy PSMs come from full search
- Verify RT distributions are reasonable for both targets and decoys
- Check FDR calculations in ScoringSearch are stable

### 4. Performance Testing
- Monitor memory usage during dual search
- Measure runtime increase (expect ~2x for SecondPassSearch)
- Profile for optimization opportunities

## Rollback Plan
If issues arise:
1. This implementation is isolated to `process_file!` in SecondPassSearch.jl
2. Can easily revert by checking out previous version of the file
3. No changes to data structures or interfaces
4. No impact on other search methods

## Success Criteria
1. Pipeline completes successfully on test dataset
2. PSM counts are reasonable (similar total, proper target/decoy ratio)
3. FDR curves look normal in ScoringSearch
4. Protein groups are identified correctly
5. No duplicate PSMs
6. Memory usage remains acceptable

## Alternative Approaches Considered

### Alternative 1: Single search with decoy-exclusion filter
- Modify `process_scans!` to skip decoys during target precursor quantification
- **Rejected**: Would require complex logic and may not fully separate target/decoy competition

### Alternative 2: Modify rt_index selection logic
- Build single rt_index but dynamically exclude decoys during transition selection
- **Rejected**: More complex, harder to verify correctness

### Alternative 3: Post-processing filter
- Run normal search but filter out mixed target/decoy scans
- **Rejected**: Would lose too many PSMs

## Questions for User
1. Should this be made optional via parameters (e.g., `separate_target_decoy_search = true`)?
2. Are there specific QC metrics to monitor beyond standard FDR curves?
3. Is the ~2x runtime increase for SecondPassSearch acceptable?
4. Should we implement any optimizations to reduce memory usage?
