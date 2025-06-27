# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

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