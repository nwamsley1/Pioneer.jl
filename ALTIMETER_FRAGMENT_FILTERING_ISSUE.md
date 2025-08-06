# Altimeter Fragment Filtering Issue

## Problem Summary

When using Altimeter with Koina, fragments are being filtered for invalid ion types (precursor, immonium) **AFTER** selecting the top X fragments by intensity, rather than **BEFORE**. This results in libraries with fewer fragments than expected.

## Current Fragment Processing Flow

### 1. Fragment Prediction (`fragment_predict.jl`)

When Altimeter returns fragments from Koina:

```julia
# Line 203-211 in predict_fragments_batch for SplineCoefficientModel
for (i, response) in enumerate(responses)
    batch_result = parse_koina_batch(model, response)
    # ...
    batch_df = batch_result.fragments
    batch_df[!, :precursor_idx] = repeat(start_idx:(start_idx + n_precursors_in_batch - one(UInt32)), 
                                          inner=batch_result.frags_per_precursor)
    filter_fragments!(batch_df, model)  # <-- PROBLEM: Minimal filtering here
    push!(batch_dfs, batch_df)
end
```

### 2. Fragment Filtering Function

The current `filter_fragments!` for Altimeter only does minimal filtering:

```julia
# Line 248-257
function filter_fragments!(df::DataFrame, model::SplineCoefficientModel)
    # Basic filtering common to all models
    filter!(:mz => x -> x > 0, df)  # Remove invalid m/z
    
    # Model-specific filtering (doesn't apply to SplineCoefficientModel)
    if model isa InstrumentSpecificModel
        filter!(row -> !occursin('i', row.annotation), df)  # Remove isotope peaks
    end
end
```

**Issue**: This doesn't filter out precursor ions, immonium ions, or other invalid ion types!

### 3. Fragment Sorting

Fragments are sorted by intensity within each precursor:

```julia
# Line 263-264
function sort_fragments!(df::DataFrame)
    sort!(df, [:precursor_idx, order(:intensities, rev=true)])
```

Note: For Altimeter, fragments come pre-sorted from Koina (line 225 comment).

### 4. Library Building (`build_poin_lib.jl`)

During library building, fragments are processed in order and filtered:

```julia
# Line 1293-1343 in getDetailedFrags for SplineCoefficientModel
for pid in ProgressBar(range(one(UInt32), n_precursors))
    rank = 1
    for frag_idx in range(frag_start_idx, frag_stop_idx)
        if !fragFilter(  # <-- Filters here, but AFTER counting rank
                frag_is_y[frag_idx],
                frag_is_b[frag_idx], 
                frag_is_p[frag_idx],  # Precursor flag
                # ... other parameters ...
                include_p,  # Whether to include precursors
                include_immonium,  # Whether to include immonium
                # ...
                )
            continue
        end
        
        # Add fragment to library
        detailed_frags[detailed_frag_idx] = SplineDetailedFrag(...)
        
        rank += 1
        if rank > max_frag_rank  # <-- Stops at max fragments
            break
        end
    end
end
```

## The Problem Illustrated

Consider a peptide where Koina returns fragments sorted by intensity:

```
1. Precursor ion (intensity: 100)  <- Most intense
2. Immonium ion (intensity: 90)
3. y8 ion (intensity: 80)
4. b7 ion (intensity: 70)
5. y6 ion (intensity: 60)
```

With `max_frag_rank = 3` and `include_p = false, include_immonium = false`:

**Current behavior:**
1. Process precursor (rank=1) → Filtered out by `fragFilter`
2. Process immonium (rank=2) → Filtered out by `fragFilter`
3. Process y8 (rank=3) → Added to library
4. Stop (rank > 3)

**Result**: Only 1 fragment in library instead of 3!

**Expected behavior:**
1. Filter out precursor and immonium FIRST
2. Take top 3 from remaining: y8, b7, y6
3. Result: 3 fragments in library

## Fragment Annotation Parsing

Altimeter uses integer indices for fragments that map to an ion dictionary:

```julia
# Line 395-405 in fragment_annotation.jl
function get_altimeter_ion_dict(ion_table_path::String)
    ion_index_to_name = Dict{Int32, String}()
    open(ion_table_path) do file
        for (i, l) in enumerate(eachline(file))
            ion_name, _, _, _ = split(l, '\t')
            ion_index_to_name[Int32(i-1)] = ion_name
        end
    end
    return ion_index_to_name
end
```

The annotations are parsed to determine ion types:

```julia
# Line 858-860 in fragment_parse.jl
# Ion type flags
frag_data.base_type == 'y',
frag_data.base_type == 'b', 
frag_data.base_type == 'p',  # Precursor
```

## Impact

This issue means that when using Altimeter:
1. Libraries may have fewer fragments than requested
2. High-intensity but unwanted ion types (precursor, immonium) consume slots
3. Lower-intensity but useful fragments (b/y ions) are excluded

## Solution Requirements

The fix must:
1. Filter invalid ion types BEFORE selecting top X fragments
2. Maintain backward compatibility with other models
3. Properly handle Altimeter's integer annotation system
4. Respect user preferences for including/excluding ion types