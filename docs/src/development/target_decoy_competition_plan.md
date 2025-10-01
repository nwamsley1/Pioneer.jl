# Target/Decoy Competition Filtering Plan

## Overview

Implement target/decoy competition filtering after PSM scoring (Step 1) and before probability calculation (Step 2) in ScoringSearch. This will filter out losing peptides in target/decoy pairs, keeping only the winner based on score comparison.

## Context

### Current Pipeline Position
- **After**: `score_precursor_isotope_traces()` completes (Step 1 - Model Training & PSM Scoring)
- **Before**: Probability computation and best trace selection (Step 2)
- **Files to process**: Scored PSM files in `second_pass_folder` (e.g., `scored_psms_[file_idx].arrow`)

### Key Data Structures
- **pair_id**: UInt32 from library (`getPairIdx(precursors)`) - identifies target/decoy pairs
  - For targets and their matched decoys, `pair_id` is the same
  - Precursors without pairs have unique `pair_id` values
- **isotopes_captured**: Tuple{Int8, Int8} - identifies which isotopes were captured
  - Example: `(0, 3)` means M0 to M+3 isotopes
- **Grouping key**: `(pair_id, isotopes_captured)` uniquely identifies competing PSMs

## Requirements

### Validation Rules
1. Each `(pair_id, isotopes_captured)` group must have 1 or 2 rows:
   - **1 row**: No competition (unpaired precursor)
   - **2 rows**: Target vs decoy competition
   - **>2 rows**: ERROR - indicates data corruption

2. In 2-row groups:
   - One row must be target (`target==true`)
   - One row must be decoy (`target==false`)
   - If both are targets or both decoys: ERROR

### Winner Selection
- **Criterion**: Higher `prob` score wins
- **Actions**:
  - Winner: Keep row unchanged
  - Loser: Delete row from DataFrame

### Library pair_id Source
- **Critical**: Must use `pair_id` from the spectral library, NOT the one constructed in `score_precursor_isotope_traces`
- **Access**: Via `getPairIdx(precursors)` which returns `precursors.data[:pair_id]`
- **Mapping**: Use `precursor_idx` to map PSM rows to library pair_id values

## Implementation Plan

### Phase 1: New Function in score_psms.jl

Create `apply_target_decoy_competition!` function:

```julia
"""
    apply_target_decoy_competition!(file_paths::Vector{String},
                                   precursors::LibraryPrecursors)

Apply target/decoy competition filtering to scored PSM files.

For each file:
1. Load scored PSMs
2. Add library pair_id column
3. Group by (pair_id, isotopes_captured)
4. Validate group sizes (1 or 2 rows only)
5. For 2-row groups: keep winner (higher prob), delete loser
6. Write filtered data back to file

# Arguments
- `file_paths`: Paths to scored PSM Arrow files
- `precursors`: Library precursors for pair_id lookup

# Errors
- Throws if any group has >2 rows
- Throws if 2-row group doesn't have exactly 1 target and 1 decoy
"""
function apply_target_decoy_competition!(
    file_paths::Vector{String},
    precursors::LibraryPrecursors
)
    # Get library pair_id mapping
    pair_id_array = getPairIdx(precursors)  # Returns entire :pair_id column

    for file_path in file_paths
        # Step 1: Load file into mutable DataFrame
        # Must use columntable to materialize for modification
        df = DataFrame(Tables.columntable(Arrow.Table(file_path)))

        # Step 2: Map precursor_idx to library pair_id
        df.pair_id = [pair_id_array[pid] for pid in df.precursor_idx]

        # Step 3: Group by (pair_id, isotopes_captured)
        groups = groupby(df, [:pair_id, :isotopes_captured])

        # Step 4: Process each group
        rows_to_keep = Bool[]
        sizehint!(rows_to_keep, nrow(df))

        for group in groups
            n = nrow(group)

            # Validation
            if n > 2
                error("Group (pair_id=$(group.pair_id[1]), isotopes=$(group.isotopes_captured[1])) has $n rows (expected 1 or 2)")
            end

            if n == 1
                # No competition - keep the row
                push!(rows_to_keep, true)
            else  # n == 2
                # Validate one target, one decoy
                targets = group.target
                if count(targets) != 1
                    error("Group (pair_id=$(group.pair_id[1]), isotopes=$(group.isotopes_captured[1])) has $(count(targets)) targets (expected exactly 1)")
                end

                # Find winner (higher prob)
                probs = group.prob
                winner_idx = argmax(probs)

                # Mark winner to keep, loser to delete
                push!(rows_to_keep, 1 == winner_idx)
                push!(rows_to_keep, 2 == winner_idx)
            end
        end

        # Step 5: Filter DataFrame
        filtered_df = df[rows_to_keep, :]

        # Step 6: Remove temporary pair_id column
        select!(filtered_df, Not(:pair_id))

        # Step 7: Write back to file using FileOperations writeArrow
        writeArrow(file_path, filtered_df)

        @info "Target/decoy competition: $(nrow(df)) → $(nrow(filtered_df)) PSMs" file=basename(file_path)
    end
end
```

### Phase 2: Integration into ScoringSearch Pipeline

Modify `ScoringSearch.jl` `summarize_results!` function:

```julia
# After Step 1 (score_precursor_isotope_traces) completes
step1_time = @elapsed begin
    score_precursor_isotope_traces(...)
end

# NEW: Apply target/decoy competition filtering
if params.enable_target_decoy_competition  # Add parameter
    @info "Applying target/decoy competition filtering..."
    competition_time = @elapsed begin
        apply_target_decoy_competition!(
            valid_second_pass_psms,
            getPrecursors(getSpecLib(search_context))
        )
    end
    @info "Target/decoy competition completed in $(round(competition_time, digits=2)) seconds"
end

# Continue with Step 2 (probability calculation)
step2_time = @elapsed begin
    # Existing code for best trace selection...
end
```

### Phase 3: Parameter Addition

Add to `PioneerParameters`:

```julia
struct ScoringSearchParameters{I<:IsotopeTraceType} <: SearchParameters
    # ... existing parameters ...
    enable_target_decoy_competition::Bool
end
```

Default: `false` (opt-in feature)

## File Operations Pattern

### Following FileOperations Conventions

The implementation follows Pioneer's FileOperations patterns:

1. **Direct Arrow I/O**: Load/write using Arrow.jl directly for simple transformations
2. **In-place modification**: Filter DataFrame in memory, write back to same path
3. **writeArrow usage**: Use `writeArrow(path, df)` for Windows-compatible writing
4. **No FileReference needed**: Simple filter operation doesn't require ref tracking

### Alternative: FileReference Approach

If we want to use FileReferences (more robust):

```julia
function apply_target_decoy_competition!(
    file_paths::Vector{String},
    precursors::LibraryPrecursors
)
    pair_id_array = getPairIdx(precursors)

    for file_path in file_paths
        # Create file reference
        psm_ref = PSMFileReference(file_path)

        # Use stream_transform with competition filter
        filter_fn = df -> begin
            df.pair_id = [pair_id_array[pid] for pid in df.precursor_idx]
            # ... group and filter logic ...
            select!(filtered_df, Not(:pair_id))
            return filtered_df
        end

        # Apply transformation (writes to temp, replaces original)
        temp_path = file_path * ".tmp"
        new_ref = stream_transform(psm_ref, temp_path, filter_fn)
        mv(temp_path, file_path; force=true)
    end
end
```

## Testing Strategy

### Unit Tests

```julia
@testset "Target/Decoy Competition" begin
    # Test 1: Single precursor (no competition)
    df = DataFrame(
        precursor_idx=[1],
        isotopes_captured=[(0,3)],
        prob=[0.9],
        target=[true]
    )
    # Expected: Same DataFrame returned

    # Test 2: Target wins
    df = DataFrame(
        precursor_idx=[1, 2],
        isotopes_captured=[(0,3), (0,3)],
        prob=[0.9, 0.7],  # Target higher
        target=[true, false]
    )
    # Expected: Only target row remains

    # Test 3: Decoy wins
    df = DataFrame(
        precursor_idx=[1, 2],
        isotopes_captured=[(0,3), (0,3)],
        prob=[0.7, 0.9],  # Decoy higher
        target=[true, false]
    )
    # Expected: Only decoy row remains

    # Test 4: Multiple isotope traces (no competition between traces)
    df = DataFrame(
        precursor_idx=[1, 2, 1, 2],
        isotopes_captured=[(0,3), (0,3), (1,4), (1,4)],
        prob=[0.9, 0.7, 0.8, 0.6],
        target=[true, false, true, false]
    )
    # Expected: 2 rows (winners from each isotope group)

    # Test 5: Error case - 3 rows in group
    df = DataFrame(
        precursor_idx=[1, 1, 1],
        isotopes_captured=[(0,3), (0,3), (0,3)],
        prob=[0.9, 0.8, 0.7],
        target=[true, false, true]
    )
    # Expected: Error thrown

    # Test 6: Error case - both targets
    df = DataFrame(
        precursor_idx=[1, 2],
        isotopes_captured=[(0,3), (0,3)],
        prob=[0.9, 0.7],
        target=[true, true]
    )
    # Expected: Error thrown
end
```

### Integration Test

Run full SearchDIA pipeline with `enable_target_decoy_competition=true`:

```julia
params = GetSearchParams(...)
params.scoring.enable_target_decoy_competition = true
SearchDIA(params)
```

Verify:
- No errors during scoring
- PSM counts decrease appropriately
- Final protein groups are reasonable
- Downstream steps (integration, MaxLFQ) work correctly

## Performance Considerations

### Memory Usage
- Loads one file at a time (memory-efficient)
- DataFrame grouping is efficient for this operation
- Temporary pair_id column adds minimal overhead

### Speed
- Grouping is O(n log n) with sorting
- Should take <1 second per file for typical dataset
- Total overhead: ~5-10 seconds for 100 files

### Disk I/O
- Reads and writes each file once
- Arrow format is fast for I/O
- Could optimize with in-place file modification if needed

## Edge Cases

1. **No decoy matches**: Some targets may not have paired decoys
   - Handled by 1-row groups (no competition)

2. **Missing pair_id in library**:
   - Should not happen if library is well-formed
   - Could add validation check

3. **Tied scores**: Both have same prob
   - Current: `argmax` picks first (target if listed first)
   - Alternative: Could add tiebreaker (e.g., prefer target)

4. **Empty files**: Some files may have no PSMs
   - Safe: groupby on empty DataFrame returns empty result

5. **Different isotope traces**: Same precursor, different isotopes
   - Correct: Treated as separate competitions (no interaction)

## Migration Path

### Phase 1: Development (Current)
- Add function to score_psms.jl
- Add parameter to ScoringSearchParameters
- Add integration point in ScoringSearch.jl
- Default `enable_target_decoy_competition = false`

### Phase 2: Testing
- Unit tests for competition logic
- Integration tests with real data
- Performance benchmarking

### Phase 3: Validation
- Compare results with/without competition
- Verify protein group changes are expected
- Check FDR calibration is maintained

### Phase 4: Production
- Enable by default if validation passes
- Document in user guide
- Add to CHANGELOG

## Open Questions

1. **Parameter name**: `enable_target_decoy_competition` vs `apply_competition_filter`?
2. **Logging level**: @info or @debug for per-file stats?
3. **Tiebreaker rule**: Keep first (current) or add explicit rule?
4. **Validation strictness**: Error on issues or warning + skip?
5. **Output format**: Keep same Arrow schema or add metadata?

## Related Code References

- Grouping pattern: `src/utils/ML/percolatorSortOf.jl:182` (sort and group by pair_id, isotopes)
- Library access: `src/structs/LibraryIon.jl:785` (getPairIdx function)
- File writing: `src/utils/FileOperations/io/ArrowOperations.jl` (writeArrow pattern)
- Pipeline integration: `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl:274-289`
