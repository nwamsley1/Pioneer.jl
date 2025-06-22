"""
FileOperations Module

Comprehensive file operations system for Pioneer.jl with memory-efficient
streaming operations, type-safe file references, and composable pipelines.

## Features

- **Type-safe file references** with schema and sort state tracking
- **Memory-efficient streaming** operations for large datasets  
- **Composable pipeline API** for chaining transformations
- **High-performance merging** with heap-based algorithms
- **Algorithm integration** for protein inference and MaxLFQ
- **Cross-platform compatibility** with Windows file handling

## Architecture

```
FileOperations/
├── core/              # File references, schema, sort state
├── streaming/         # Memory-efficient operations
├── pipeline/          # Composable transformation framework
├── algorithms/        # Algorithm wrappers and validation
└── io/               # Arrow operations and file utilities
```

## Example Usage

```julia
using FileOperations

# Create file references
psm_ref = PSMFileReference("data/psms.arrow")
protein_ref = ProteinGroupFileReference("data/proteins.arrow")

# Build processing pipeline
pipeline = TransformPipeline() |>
    filter_by_threshold(:qval, 0.01) |>
    add_column(:log_score, row -> log(row.score)) |>
    sort_by([:log_score], rev=[true])

# Apply pipeline
apply_pipeline!(psm_ref, pipeline)

# Merge multiple files
merged_ref = stream_sorted_merge(
    [psm_ref1, psm_ref2], 
    "output/merged.arrow",
    :score, :target;
    reverse=[true, true]
)
```

## Memory Efficiency

All operations are designed to handle large datasets efficiently:
- Streaming operations process data in batches
- Type-stable implementations minimize allocations
- Memory usage estimation and validation
- Automatic batch size optimization

## Performance

- **4-20x speedup** in merge operations through type-stable heaps
- **Parallel processing** support for multi-file operations
- **Memory-mapped** reading where possible
- **Optimized column operations** with specialized implementations
"""

# Load all core modules in dependency order
include("core/FileSchema.jl")
include("core/FileReferences.jl")
include("core/SortStateManagement.jl")

# Load streaming operations
include("streaming/StreamOperations.jl")
include("streaming/ColumnOperations.jl")
include("streaming/MergeOperations.jl")

# Load pipeline framework
include("pipeline/PipelineAPI.jl")
include("pipeline/PipelineOperations.jl")
include("pipeline/PipelineExecution.jl")

# Load algorithm integrations
include("algorithms/ProteinInference.jl")
include("algorithms/MaxLFQOperations.jl")
include("algorithms/ValidationUtils.jl")

# Load I/O operations
include("io/ArrowOperations.jl")
include("io/FileUtilities.jl")

# Re-export all public functions from submodules
# Core types and functions
export FileReference, FileSchema, PSMFileReference, ProteinGroupFileReference, 
       ProteinQuantFileReference, PairedSearchFiles,
       file_path, schema, sorted_by, row_count, exists, n_protein_groups, n_experiments,
       has_column, get_column_or_default, validate_required_columns, validate_schema,
       mark_sorted!, is_sorted_by, ensure_sorted!, validate_exists,
       create_psm_reference, create_protein_reference, create_protein_quant_reference, 
       create_reference, describe_reference

# Streaming operations
export stream_filter, stream_transform, process_with_memory_limit, estimate_batch_size,
       add_column_to_file!, update_column_in_file!, add_column_and_sort!,
       stream_sorted_merge

# Pipeline framework
export TransformPipeline, PipelineOperation, |>,
       add_column, rename_column, select_columns, remove_columns,
       filter_rows, sort_by, add_interpolated_column, filter_by_threshold,
       filter_by_multiple_thresholds,
       apply_pipeline!, transform_and_write!, apply_pipeline_batch, write_transformed

# Algorithm integrations
export apply_protein_inference, update_psms_with_scores,
       apply_maxlfq, prepare_maxlfq_input,
       validate_maxlfq_input, validate_maxlfq_parameters, check_maxlfq_memory_requirements,
       validate_join_compatibility, validate_sort_compatibility

# I/O operations
export sort_file_by_keys!, write_arrow_file, load_dataframe, column_names, has_columns,
       safe_join_files, ensure_directory_exists, create_temp_path, cleanup_temp_files,
       safe_move_file, get_file_size_mb, estimate_processing_time