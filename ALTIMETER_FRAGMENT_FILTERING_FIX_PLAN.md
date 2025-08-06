# Fix Plan for Altimeter Fragment Filtering

## Problem
Altimeter fragments from Koina are filtered for invalid ion types (precursor, immonium) AFTER ranking by intensity, causing libraries to have fewer fragments than requested.

## Proposed Solution

### Option 1: Filter During Fragment Prediction (Recommended)

**Location**: `src/Routines/BuildSpecLib/fragments/fragment_predict.jl`

Enhance the `filter_fragments!` function for `SplineCoefficientModel` to use annotation information:

```julia
function filter_fragments!(df::DataFrame, model::SplineCoefficientModel; 
                          include_precursor::Bool=false,
                          include_immonium::Bool=false,
                          ion_dictionary::Union{Dict{Int32,String}, Nothing}=nothing)
    # Basic filtering
    filter!(:mz => x -> x > 0, df)
    
    # If we have the ion dictionary, filter by ion type
    if !isnothing(ion_dictionary) && hasproperty(df, :annotation)
        filter!(row -> begin
            ion_idx = row.annotation
            if haskey(ion_dictionary, ion_idx)
                ion_name = ion_dictionary[ion_idx]
                
                # Check for precursor ions (typically "p" or containing "precursor")
                is_precursor = occursin(r"^p\d*", ion_name) || occursin("precursor", lowercase(ion_name))
                if is_precursor && !include_precursor
                    return false
                end
                
                # Check for immonium ions (single letter or "Imm" prefix)
                is_immonium = occursin(r"^[A-Z]$", ion_name) || startswith(ion_name, "Imm")
                if is_immonium && !include_immonium
                    return false
                end
            end
            return true
        end, df)
    end
end
```

### Option 2: Fix During Fragment Table Building

**Location**: `src/Routines/BuildSpecLib/fragments/fragment_parse.jl`

Add a pre-filtering step in `parse_altimeter_fragments` after parsing annotations:

```julia
function parse_altimeter_fragments(
    precursor_table::Arrow.Table,
    fragment_table::Arrow.Table,
    annotation_type::FragAnnotation,
    ion_dictionary::Dict{Int32, String},
    precursor_batch_size::Int64,
    immonium_data_path::String,
    out_dir::String,
    mods_to_sulfur_diff::Dict{String, Int8},
    iso_mod_to_mass::Dict{String, Float32},
    model_type::KoinaModelType;
    include_precursor::Bool=false,  # Add parameter
    include_immonium::Bool=false,   # Add parameter
    max_fragments_per_precursor::Int=50  # Add parameter
)
    # ... existing code ...
    
    # After parsing annotations, create a filtered fragment table
    # that excludes unwanted ion types BEFORE writing to disk
    
    # Group fragments by precursor
    grouped_fragments = groupby(fragment_table, :precursor_idx)
    
    # Filter and select top N for each precursor
    filtered_fragments = []
    for group in grouped_fragments
        # Apply ion type filters
        valid_frags = filter(row -> 
            should_include_fragment(row.annotation, 
                                   ion_annotation_to_features_dict,
                                   include_precursor,
                                   include_immonium),
            group)
        
        # Sort by intensity (should already be sorted from Koina)
        sort!(valid_frags, :intensities, rev=true)
        
        # Take top N
        top_frags = first(valid_frags, min(max_fragments_per_precursor, nrow(valid_frags)))
        push!(filtered_fragments, top_frags)
    end
    
    # Continue with filtered fragments...
```

### Option 3: Thread Parameters Through the Pipeline

**Most Complete Solution - Combines Options 1 & 2**

1. **Update `BuildSpecLib` parameters** to include ion filtering preferences
2. **Pass parameters through the prediction pipeline**
3. **Apply filtering at the earliest possible stage**

#### Step 1: Update Parameter Structure

```julia
# In parameter loading
struct FragmentFilterParams
    include_precursor::Bool
    include_immonium::Bool
    include_internal::Bool
    max_fragments_per_precursor::Int
end
```

#### Step 2: Update Fragment Prediction

```julia
# In predict_fragments
function predict_fragments(
    peptide_table_path::String,
    frags_out_path::String,
    model_type::KoinaModelType,
    instrument_type::String,
    max_koina_batches::Int,
    batch_size::Int,
    model_name::String;
    intensity_threshold::Float32 = 0.001f0,
    filter_params::FragmentFilterParams = FragmentFilterParams(false, false, true, 50)
)
```

#### Step 3: Apply Filtering Before Writing

```julia
# In predict_fragments_batch for SplineCoefficientModel
# After getting fragments from Koina but before writing
fragments_df = vcat(batch_dfs...)

# Group by precursor and filter
grouped = groupby(fragments_df, :precursor_idx)
filtered_dfs = []

for group in grouped
    # Filter by ion type first
    valid = filter_by_ion_type(group, ion_dict, filter_params)
    
    # Sort by intensity (if not already sorted)
    sort!(valid, :intensities, rev=true)
    
    # Take top N
    top_n = first(valid, min(filter_params.max_fragments_per_precursor, nrow(valid)))
    push!(filtered_dfs, top_n)
end

return vcat(filtered_dfs...)
```

## Implementation Steps

1. **Add ion dictionary loading** to the fragment prediction pipeline
2. **Create helper function** to identify ion types from Altimeter annotations
3. **Implement filtering** at the appropriate stage
4. **Update tests** to verify correct filtering
5. **Add configuration parameters** to control filtering behavior

## Testing Strategy

1. Create test case with known Altimeter output containing precursor/immonium ions
2. Verify that with `include_precursor=false, include_immonium=false`:
   - These ion types are excluded
   - The correct number of fragments remains
   - The highest-scoring valid fragments are selected
3. Compare library sizes before and after the fix
4. Ensure backward compatibility with other models

## Configuration Changes

Add to `defaultBuildLibParams.json`:

```json
{
  "fragment_filtering": {
    "filter_before_ranking": true,
    "max_fragments_per_precursor": 50,
    "include_precursor": false,
    "include_immonium": false,
    "include_internal": true
  }
}
```

## Benefits

1. **Correct fragment count**: Libraries will have the requested number of fragments
2. **Better quality**: Higher-quality b/y ions instead of precursor/immonium
3. **Consistency**: Same behavior across all prediction models
4. **User control**: Configuration options for different use cases

## Risks and Mitigation

1. **Risk**: Breaking existing libraries
   - **Mitigation**: Add backward compatibility flag
   
2. **Risk**: Performance impact from additional filtering
   - **Mitigation**: Filter early in the pipeline to minimize data processing
   
3. **Risk**: Ion dictionary might not be available
   - **Mitigation**: Fall back to current behavior if dictionary is missing

## Recommended Approach

Implement **Option 3** (complete solution) with a phased rollout:

1. **Phase 1**: Add filtering infrastructure without changing default behavior
2. **Phase 2**: Enable filtering for new libraries with configuration flag
3. **Phase 3**: Make filtering the default after validation