# Comprehensive Test Suite Plan for FileOperations.jl and FileReferences.jl

## Overview
Create a comprehensive test suite that follows Pioneer.jl's testing patterns and provides full coverage for the FileReference system and pipeline operations.

## Test Structure Analysis
Based on the existing test framework:
- Tests are in `test/UnitTests/` mirroring `src/` structure
- Main test runner is `test/runtests.jl` which includes individual test files
- Tests use standard Julia Test framework with `@testset` blocks
- Temporary files are created/cleaned in test directories using `mktempdir()`
- Tests include unit tests for individual functions and integration tests
- Toy DataFrames are created and written to Arrow files in temp folders for testing

## Plan: Create Two New Test Files

### 1. `test/UnitTests/test_file_operations.jl`
Comprehensive tests for FileOperations.jl covering:

#### A. Streaming Operations
- **`sort_file_by_keys!()`** 
  - Single key sorting (ascending/descending)
  - Multi-key sorting with different combinations
  - Sort state tracking validation
  - Large file sorting performance
  - Error handling for invalid keys

- **`stream_sorted_merge()`**
  - Heap-based merge with 2-10 files
  - Different sort key combinations
  - Memory usage validation during merge
  - Error handling for mismatched schemas
  - Batch size parameter testing

- **`apply_pipeline!()`**
  - Single operation pipelines
  - Complex multi-operation pipelines
  - Pipeline operation ordering
  - Error handling and rollback scenarios

#### B. Pipeline API
- **`TransformPipeline` Core**
  - Creation and composition with `|>` operator
  - Operation chaining validation
  - Pipeline state management

- **Column Operations**
  - `add_column()` with various compute functions
  - `rename_column()` for existing and non-existing columns
  - `select_columns()` with subset selection
  - `remove_columns()` with single and multiple columns

- **Row Operations**
  - `filter_rows()` with simple and complex predicates
  - `filter_by_multiple_thresholds()` helper function
  - Performance with large datasets

- **Sort Operations**
  - `sort_by()` with automatic state tracking
  - Multi-column sorting with mixed directions
  - Sort optimization detection

#### C. Validation Functions
- **MaxLFQ Validation**
  - `validate_maxlfq_input()` - Required columns and sort order
  - `validate_maxlfq_parameters()` - Parameter range validation
  - Integration with FileReference validation

- **Schema Validation**
  - Cross-file schema compatibility
  - Required column validation with detailed errors
  - Type compatibility checking

#### D. Memory-Efficient Operations
- **Large File Handling**
  - Batch processing with different batch sizes
  - Memory usage monitoring during operations
  - Streaming vs. materialization trade-offs

- **Performance Testing**
  - Heap-based merge scaling with file count
  - Pipeline operation overhead measurement
  - Memory leak detection

#### E. Error Handling
- **File System Errors**
  - Invalid file paths and permissions
  - Corrupted Arrow files
  - Disk space issues

- **Data Validation Errors**
  - Missing columns and schema mismatches
  - Type conversion failures
  - Sort key validation failures

- **Pipeline Errors**
  - Operation failure handling
  - Partial execution rollback
  - Error message clarity

### 2. `test/UnitTests/test_file_references_comprehensive.jl`
Enhanced version of existing `test_file_references.jl` with additional coverage:

#### A. Core FileReference Functionality
- **All FileReference Types**
  - `PSMFileReference`, `ProteinGroupFileReference`, `ProteinQuantFileReference`
  - Constructor validation with real Arrow files
  - Metadata extraction accuracy
  - File existence handling

- **Schema Management**
  - `FileSchema` immutability validation
  - `has_column()`, `get_column_or_default()` with edge cases
  - `validate_required_columns()` with complex requirements

#### B. Specialized FileReference Features
- **`ProteinQuantFileReference`**
  - `n_protein_groups()`, `n_experiments()` accuracy
  - Metadata calculation with missing values
  - Integration with MaxLFQ workflow

- **Accessor Functions**
  - `file_path()`, `schema()`, `sorted_by()`, `row_count()`, `exists()`
  - Consistency across all FileReference types
  - State mutation tracking

#### C. Sort State Management
- **Basic Sort Operations**
  - `mark_sorted!()`, `is_sorted_by()`, `ensure_sorted!()`
  - Multi-key sort state tracking
  - Sort order validation

- **Advanced Sort Features**
  - `Base.sort!()` override functionality
  - Mixed sort direction error handling
  - Sort state persistence across file operations

#### D. Paired File Operations
- **`PairedSearchFiles`**
  - Constructor validation with matching/mismatched files
  - Consistency checking across file pairs
  - Error scenarios (one file missing, corrupted, etc.)

- **File Relationship Management**
  - 1:1 pairing enforcement
  - State synchronization between paired files
  - Metadata consistency validation

#### E. Utility Functions
- **Factory Methods**
  - `create_reference()` for all FileReference types
  - `create_psm_reference()`, `create_protein_reference()`, etc.
  - Type dispatch validation

- **Diagnostic Functions**
  - `describe_reference()` output format validation
  - Reference comparison and equality
  - Debug information accuracy

## Test Data Strategy

### Toy DataFrame Creation
Create realistic test data with various schemas:

```julia
# PSM test data
psm_df = DataFrame(
    precursor_idx = UInt32[1, 2, 3, 4, 5],
    prob = Float32[0.9, 0.85, 0.8, 0.75, 0.7],
    inferred_protein_group = ["P1", "P1", "P2", "P2", "P3"],
    target = Bool[true, true, true, false, true],
    entrapment_group_id = UInt8[0, 0, 0, 0, 1],
    pg_qval = Float32[0.01, 0.02, 0.03, 0.04, 0.05],
    global_qval_pg = Float32[0.01, 0.02, 0.03, 0.04, 0.05],
    use_for_protein_quant = Bool[true, true, true, true, true],
    peak_area = Float32[1000.0, 2000.0, 1500.0, 800.0, 1200.0]
)

# Protein group test data
protein_df = DataFrame(
    protein = ["P1", "P2", "P3"],
    target = Bool[true, true, true],
    entrap_id = UInt8[0, 0, 1],
    species = ["human", "human", "human"],
    abundance = Float32[5000.0, 3000.0, 2000.0],
    n_peptides = Int64[2, 2, 1]
)

# Write to temporary Arrow files for testing
temp_dir = mktempdir()
psm_file = joinpath(temp_dir, "test_psms.arrow")
protein_file = joinpath(temp_dir, "test_proteins.arrow")
Arrow.write(psm_file, psm_df)
Arrow.write(protein_file, protein_df)
```

### Edge Case Data
- Empty DataFrames
- Single-row DataFrames  
- DataFrames with missing values
- Large DataFrames for performance testing
- Malformed data for error testing

## Integration with Existing Framework

### Update `test/runtests.jl`
Add the new test files to the main test suite:
```julia
include("./UnitTests/test_file_operations.jl")
include("./UnitTests/test_file_references_comprehensive.jl")
```

### Dependencies and Imports
Ensure proper imports in test files:
```julia
using Test
using Arrow, DataFrames, Tables
using DataStructures  # For heap operations

# Include the source files
cd(@__DIR__)
package_root = dirname(dirname(@__DIR__))
include(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods", "FileReferences.jl"))
include(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods", "FileOperations.jl"))
```

## Test Coverage Goals

### Quantitative Targets
- **FileReferences.jl**: 95%+ line coverage of all public functions
- **FileOperations.jl**: 90%+ line coverage including error paths
- **Integration**: Test interaction between references and operations
- **Performance**: Benchmark critical operations and memory usage

### Qualitative Targets
- **Error Handling**: Comprehensive validation of error scenarios
- **Edge Cases**: Empty files, large files, corrupted data
- **Type Safety**: Validation of type constraints and conversions
- **Memory Safety**: No memory leaks during large operations

## Success Metrics

1. **Functional Correctness**
   - All tests pass consistently
   - Operations maintain data integrity
   - FileReference metadata remains accurate

2. **Performance**
   - No significant performance regressions
   - Memory usage within acceptable bounds
   - Scalability with large datasets

3. **Robustness**
   - Graceful error handling with informative messages
   - Recovery from partial failures
   - Consistent behavior across edge cases

4. **Maintainability**
   - Clear test structure and documentation
   - Easy addition of new test cases
   - Integration with existing test framework

## Implementation Order

1. **Phase 1**: Create comprehensive `test_file_references_comprehensive.jl`
   - Extend existing tests with additional coverage
   - Add specialized FileReference functionality
   - Test all accessor and utility functions

2. **Phase 2**: Create `test_file_operations.jl`
   - Start with basic streaming operations
   - Add pipeline API testing
   - Include validation and error handling

3. **Phase 3**: Integration and Performance
   - Add integration tests between modules
   - Include performance benchmarks
   - Add memory usage validation

4. **Phase 4**: Update Test Framework
   - Integrate with `test/runtests.jl`
   - Ensure compatibility with existing tests
   - Add CI/CD validation

This comprehensive plan ensures thorough testing of the FileReference system while following Pioneer.jl's established patterns and maintaining high code quality standards.