# Rebase Conflict Analysis and Resolution Plan

## Date: 2025-08-28

## Branch Information
- **Feature Branch**: `improve-parameter-tuning-robustness`
- **Target Branch**: `main` 
- **Commits to Rebase**: 222 total

## Current Rebase Status
- **Progress**: On commit 11 of 222
- **Current Conflict**: Commit `8a1fcc81` - "feat(ParameterTuning): Implement graceful fallback for parameter tuning failures"

## Conflicts Encountered

### 1. Binary Files (RECURRING)
**Files**: 
- `data/.DS_Store`
- `src/.DS_Store`

**Nature**: macOS system files that track folder display preferences
**Issue**: These files appear in almost every commit (18+ times)
**Resolution**: 
  - Already in .gitignore but were tracked before
  - Need to remove from every conflicting commit
  - Created automated script to handle this
**Status**: ⚠️ ONGOING - Appearing in multiple commits

### 2. Expected Major Conflict Areas

Based on the feature branch work, the following areas are likely to have conflicts:

#### A. ParameterTuningSearch Module
**Expected Conflicts**:
- `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl`
- Main branch has 4 recent commits affecting ScoringSearch and ML utilities
- Feature branch has extensive 3-phase convergence algorithm improvements

**Resolution Strategy**:
1. Keep the enhanced 3-phase convergence algorithm from feature branch
2. Ensure compatibility with any new interfaces from main
3. Preserve diagnostic infrastructure and cross-run learning features

#### B. ScoringSearch Module  
**Expected Conflicts**:
- `src/Routines/SearchDIA/SearchMethods/ScoringSearch/ScoringSearch.jl`
- `src/Routines/SearchDIA/SearchMethods/ScoringSearch/scoring_interface.jl`
- Main branch simplified scoring_interface.jl (77 lines reduced)

**Resolution Strategy**:
1. Accept main's simplified interface
2. Ensure feature branch's protein scoring fixes are preserved
3. Verify logodds function still receives correct probability values

#### C. ML Utilities
**Expected Conflicts**:
- `src/utils/ML/ftrUtilities.jl` 
- `src/utils/ML/percolatorSortOf.jl`
- Main branch has refactored percolatorSortOf.jl (+82 lines, better structure)

**Resolution Strategy**:
1. Accept main's refactoring
2. Ensure feature branch's fixes are compatible
3. Test ML scoring pipeline after resolution

## Key Feature Branch Improvements to Preserve

### 1. ParameterTuningSearch Enhancements
- **3-Phase Convergence**: Zero bias, positive shift, negative shift phases
- **Scan Scaling**: Dynamic adjustment when insufficient PSMs found
- **Cross-Run Learning**: Apply learned parameters across files
- **Comprehensive Diagnostics**: Detailed convergence tracking
- **PDF Generation**: In-memory plot storage, combined PDFs only

### 2. Critical Bug Fixes
- **Protein Scoring**: Fixed logodds receiving log-sum instead of probabilities
  - Conversion: `p = 1 - exp(-pg_score)` when probit skipped
- **RT Bin Indexing**: Ensured 1-based indexing throughout
- **ArrowTableReference**: Fixed length method compatibility

### 3. SearchMethods Refactoring (if present)
- FileReference type hierarchy
- Algorithm wrappers for protein inference
- Generic N-key heap-based merge
- Simplified MaxLFQSearch using MSData directly

## Testing Plan After Resolution

1. **Unit Tests**:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.test()'
   ```

2. **Integration Test**:
   ```julia
   SearchDIA("./data/ecoli_test/ecoli_test_params.json")
   ```

3. **Specific Component Tests**:
   - ParameterTuningSearch convergence
   - Protein scoring with ML
   - Cross-run parameter learning
   - PDF generation pipeline

## Next Steps

1. Continue rebase with `git rebase --continue`
2. Resolve each conflict carefully:
   - For ParameterTuningSearch: Preserve feature branch enhancements
   - For ScoringSearch: Merge with main's simplifications
   - For ML utilities: Accept main's structure, preserve fixes
3. After each major conflict resolution, compile and test
4. Document any non-obvious resolutions in commit messages
5. Final comprehensive test suite run

## Resolution Commands

```bash
# To continue after resolving conflicts
git add <resolved_files>
git rebase --continue

# To skip a problematic commit (use cautiously)
git rebase --skip

# To abort if needed
git rebase --abort

# To see remaining commits
git rebase --edit-todo
```

## Notes

- The feature branch has 81 commits of improvements that must be preserved
- Main branch has important ML refactoring that should be incorporated
- Both branches have valuable changes - careful merging required
- Test thoroughly after each major conflict resolution