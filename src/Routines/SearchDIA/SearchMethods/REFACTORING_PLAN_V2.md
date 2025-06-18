# Updated SearchMethods Refactoring Plan

## Current Progress

### ✅ Completed
- FileReferences.jl with basic types
- SearchResultReferences.jl 
- FileOperations.jl with streaming operations
- Comprehensive unit tests

### 🔄 Needs Revision
- Add abstract FileReference type hierarchy
- Update operations to use abstract types

## Improved Architecture with Abstract Types

### Layer 1: Abstract Reference Types

```julia
# Abstract base type for all file references
abstract type FileReference end

# Common interface all FileReferences must implement
file_path(ref::FileReference) = ref.file_path
schema(ref::FileReference) = ref.schema
sorted_by(ref::FileReference) = ref.sorted_by
row_count(ref::FileReference) = ref.row_count
exists(ref::FileReference) = ref.file_exists

# Concrete implementations
mutable struct PSMFileReference <: FileReference
    file_path::String
    schema::FileSchema
    sorted_by::Tuple{Vararg{Symbol}}
    row_count::Int64
    file_exists::Bool
end

mutable struct ProteinGroupFileReference <: FileReference
    file_path::String
    schema::FileSchema
    sorted_by::Tuple{Vararg{Symbol}}
    row_count::Int64
    file_exists::Bool
end
```

### Layer 2: Generic Operations on Abstract Types

```julia
# Generic operations that work with any FileReference
function ensure_sorted!(ref::FileReference, keys::Symbol...)
    if !is_sorted_by(ref, keys...)
        sort_file_by_keys!(ref, keys...)
    end
    return ref
end

function validate_schema(ref::FileReference, required_cols::Set{Symbol})
    validate_required_columns(schema(ref), required_cols)
end

# Type-specific operations when needed
function validate_psm_schema(ref::PSMFileReference)
    required = Set([:precursor_idx, :prob, :target, :use_for_protein_quant])
    validate_schema(ref, required)
end
```

### Layer 3: High-Level Algorithm Wrappers

```julia
# In FileOperations.jl - wraps existing algorithms
function apply_algorithm(ref::FileReference, 
                        algorithm_fn::Function, 
                        output_path::String;
                        batch_size::Int=100_000)
    validate_exists(ref)
    
    # Stream through file, applying algorithm in batches
    result_ref = stream_transform(ref, output_path, batch_size=batch_size) do batch
        algorithm_fn(batch)
    end
    
    return result_ref
end

# Specific wrappers for existing algorithms
function apply_protein_inference(psm_ref::PSMFileReference, output_path::String, args...)
    validate_psm_schema(psm_ref)
    
    apply_algorithm(psm_ref, output_path) do batch
        # Call existing protein inference from utils/proteinInference.jl
        getProteinGroupsDict(batch, args...)
    end
    
    return ProteinGroupFileReference(output_path)
end

function apply_maxlfq(psm_refs::Vector{PSMFileReference}, output_path::String, args...)
    # Validate all inputs
    for ref in psm_refs
        ensure_sorted!(ref, :inferred_protein_group, :target, 
                      :entrapment_group_id, :precursor_idx)
    end
    
    # Call existing LFQ from utils/maxLFQ.jl
    # But through our streaming interface
    ...
end
```

## Complete Implementation Plan

### Phase 1: Refactor Type Hierarchy (REVISED)
1. Add abstract `FileReference` type
2. Update PSMFileReference and ProteinGroupFileReference to inherit
3. Create common interface functions
4. Update FileOperations to use abstract types where possible
5. Add tests for polymorphic behavior

### Phase 2: Create Algorithm Wrappers
1. Add `apply_protein_inference` to FileOperations
2. Add `apply_maxlfq` to FileOperations  
3. Add `update_psms_with_scores` for streaming updates
4. Each wrapper:
   - Validates inputs through references
   - Calls existing algorithm from utils/
   - Returns appropriate reference type

### Phase 3: Update ScoringSearch
1. Create `scoring_interface.jl` with functions that:
   - Accept only reference types
   - Return reference types
   - Never directly access files
2. Update `summarize_results!` to:
   - Create references as work progresses
   - Use only reference-based operations
   - Store results in SearchContext

### Phase 4: Update MaxLFQSearch
1. Retrieve references from SearchContext
2. Validate using reference methods
3. Call wrapped algorithms
4. Return reference types

### Phase 5: Add to SearchContext
1. Add field: `method_results::Dict{Type{<:SearchMethod}, Any}`
2. Add type-safe accessors:
   ```julia
   store_results!(ctx::SearchContext, ::Type{ScoringSearch}, refs::ScoringSearchResultRefs)
   get_results(ctx::SearchContext, ::Type{ScoringSearch})::Union{Nothing, ScoringSearchResultRefs}
   ```

## Key Design Principles

1. **Abstraction Layers**:
   - References: Metadata only, no file access
   - Operations: All file access happens here
   - Interfaces: Search methods use only references

2. **Type Safety**:
   - Abstract types enable polymorphism
   - Concrete types provide specificity
   - Compiler can catch type mismatches

3. **Memory Efficiency**:
   - Streaming operations throughout
   - Never load full files unless necessary
   - Batch processing with configurable sizes

4. **Validation**:
   - Automatic at operation boundaries
   - Sort state tracked and enforced
   - Schema compatibility checked

## Example Usage Pattern

```julia
# In ScoringSearch
function score_proteins(search_context::SearchContext)
    # Get PSM references (no file access)
    psm_refs = get_psm_refs(search_context)
    
    # Validate and prepare (through references)
    for ref in psm_refs
        ensure_sorted!(ref, :precursor_idx)
        validate_psm_schema(ref)
    end
    
    # Apply algorithm (encapsulated file access)
    paired_refs = map(enumerate(psm_refs)) do (idx, psm_ref)
        pg_path = joinpath(output_dir, "proteins_$idx.arrow")
        pg_ref = apply_protein_inference(psm_ref, pg_path, precursors, ...)
        PairedSearchFiles(psm_ref, pg_ref, idx)
    end
    
    # Store results (as references)
    result_refs = ScoringSearchResultRefs(paired_refs)
    store_results!(search_context, ScoringSearch, result_refs)
end
```

## Benefits Over Current Implementation

1. **Better Abstraction**: Generic operations on abstract types
2. **Type Safety**: Can't mix up different file types
3. **Centralized Validation**: One place for each check
4. **Easier Testing**: Can mock FileReference types
5. **Clear Boundaries**: File access only in FileOperations
6. **Maintainable**: Changes to file handling in one place

## Next Immediate Steps

1. Refactor FileReferences.jl to use abstract base type
2. Update FileOperations.jl to use abstract types
3. Add algorithm wrapper functions
4. Update existing tests
5. Start modifying ScoringSearch

## Critical Clarifications

**IMPORTANT**: 
- The protein inference algorithm already exists in `src/utils/proteinInference.jl`
- The MaxLFQ algorithm already exists in `src/utils/maxLFQ.jl`
- We are NOT reimplementing these algorithms
- We are ONLY wrapping them with better abstractions and safety checks

## Success Metrics

1. All existing tests pass (no algorithm changes)
2. Better error messages when things go wrong
3. Impossible to accidentally use unsorted files
4. Clear tracking of data flow between methods
5. No performance degradation
6. Type-safe interfaces prevent common errors