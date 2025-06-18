# SearchMethods Refactoring Plan: Improving Encapsulation and Orthogonality

## Current Issues

1. **Tight Coupling**: ScoringSearch writes PSM files that MaxLFQSearch directly reads and modifies
2. **Shared State**: Both methods rely on file paths and folder structures in SearchContext
3. **Data Format Dependencies**: MaxLFQSearch expects specific columns from ScoringSearch
4. **Side Effects**: Both methods modify files in place, making testing and reuse difficult

## Proposed Improvements

### 1. Define Clear Data Interfaces

Create abstract types for data exchange:

```julia
abstract type ProteinQuantData end
abstract type PSMData end

struct ScoringSearchOutput <: PSMData
    psms::DataFrame
    protein_groups::DataFrame
    file_mappings::Dict{String,String}
end

struct MaxLFQInput <: ProteinQuantData
    psms::DataFrame
    required_columns::Vector{Symbol}
end
```

### 2. Separate Data Processing from I/O

Extract core logic into pure functions:

```julia
# Instead of reading/writing files directly:
function process_maxlfq(input::MaxLFQInput, params::MaxLFQParameters)
    # Pure computation, returns DataFrame
    return quantified_proteins
end

# Separate I/O layer:
function write_maxlfq_results(results::DataFrame, output_path::String)
    # Handle file writing
end
```

### 3. Use Result Types for Error Handling

Replace file-based communication with structured results:

```julia
struct SearchMethodResult{T}
    data::T
    metadata::Dict{String,Any}
    warnings::Vector{String}
end
```

### 4. Create Data Transformation Pipeline

Make data flow explicit:

```julia
# In SearchContext
mutable struct SearchPipeline
    results::Dict{Type{<:SearchMethod}, SearchMethodResult}
end

# Methods register their outputs:
function register_result!(pipeline::SearchPipeline, ::Type{ScoringSearch}, result)
    pipeline.results[ScoringSearch] = result
end

# Methods retrieve inputs:
function get_input(pipeline::SearchPipeline, ::Type{MaxLFQSearch})
    scoring_result = pipeline.results[ScoringSearch]
    return transform_to_maxlfq_input(scoring_result)
end
```

### 5. Decouple File Management

Use a file manager abstraction:

```julia
abstract type DataStore end

struct FileSystemStore <: DataStore
    base_path::String
end

struct InMemoryStore <: DataStore
    data::Dict{String,Any}
end

# Methods work with abstract store:
function store_results!(store::DataStore, key::String, data)
    # Implementation depends on store type
end
```

### 6. Make Dependencies Explicit

```julia
struct MaxLFQSearch <: SearchMethod
    dependencies::Vector{Type{<:SearchMethod}}
    
    MaxLFQSearch() = new([ScoringSearch])
end

# Framework checks dependencies are satisfied
function can_execute(method::SearchMethod, completed::Set{Type})
    all(dep in completed for dep in method.dependencies)
end
```

### 7. Separate Configuration from Execution

```julia
# Configuration phase
struct MaxLFQConfig
    q_value_threshold::Float32
    min_peptides::Int
    # ... other params
end

# Execution phase receives only what it needs
function execute_maxlfq(psms::DataFrame, config::MaxLFQConfig)
    # No access to global state
end
```

## Implementation Strategy

### Phase 1: Create Interfaces (Non-breaking)
- Define abstract types and result structures
- Create adapter functions to convert between old and new formats
- Add new pure function variants alongside existing ones

### Phase 2: Refactor Core Logic
- Extract pure functions from existing methods
- Create unit tests for pure functions
- Maintain backward compatibility with wrapper functions

### Phase 3: Update Search Methods
- Modify SearchMethods to use new interfaces internally
- Keep file I/O at the boundaries
- Update one method at a time to minimize risk

### Phase 4: Optimize Data Flow
- Implement in-memory data passing for small datasets
- Add streaming support for large datasets
- Profile and optimize memory usage

## Benefits

- **Testability**: Can test each component in isolation
- **Reusability**: Components can be used outside the search pipeline
- **Maintainability**: Clear interfaces make changes easier
- **Debugging**: Data flow is explicit and traceable
- **Parallelization**: Pure functions can be safely parallelized
- **Extensibility**: New search methods can be added more easily

## Example: Refactored MaxLFQ

```julia
# Pure function for MaxLFQ computation
function compute_protein_abundances(
    psms::DataFrame,
    config::MaxLFQConfig
)::DataFrame
    # Group by protein
    grouped = groupby(psms, :inferred_protein_group)
    
    # Apply MaxLFQ algorithm
    results = DataFrame()
    for group in grouped
        abundance = maxlfq_algorithm(group, config)
        push!(results, abundance)
    end
    
    return results
end

# Separate I/O concerns
function MaxLFQSearch.execute(
    input::MaxLFQInput,
    output::DataStore,
    config::MaxLFQConfig
)
    # Pure computation
    results = compute_protein_abundances(input.psms, config)
    
    # Store results
    store_results!(output, "protein_abundances", results)
    
    return SearchMethodResult(
        data = results,
        metadata = Dict("n_proteins" => nrow(results)),
        warnings = String[]
    )
end
```

## Next Steps

1. Review and approve this plan
2. Create proof-of-concept for one search method
3. Gather feedback from team
4. Create detailed implementation timeline
5. Begin phased implementation