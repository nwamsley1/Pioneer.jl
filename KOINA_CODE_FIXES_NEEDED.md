# Koina API Code Fixes Needed

## Summary
After fixing test expectations, we have 4 remaining errors that are actual code functionality problems. **IMPORTANT**: The HTML error messages in the logs are misleading - they are secondary symptoms, not the root cause.

## Root Cause Clarification
The "HTML error page" references in test logs are **NOT** from the Koina API returning HTML instead of JSON. The real issues are:
1. Column name preservation bugs that prevent API calls from being made
2. Empty fragment results from legitimate API responses
3. The HTML errors are likely from test framework exception handling or unrelated parts of the test suite

## Remaining Errors Analysis

### Issue 1: Column Name Preservation (3 errors) 🔧
**Error**: `ArgumentError: column name :precursor_charge not found in the data frame`

**Affected Test Cases**:
- `InstrumentAgnosticModel Response Parsing`
- `AlphaPeptDeep Response Parsing` 
- `Multiple Peptide Parsing`

**Root Cause**: 
The method dispatch chain loses column names between preparation and parsing. The wrapper method for `InstrumentAgnosticModel` (line 76-82 in koina_batch_prep.jl) calls the generic method which expects `precursor_charge`, `collision_energy`, and `koina_sequence` columns, but these might not be preserved through the call chain.

**Investigation Path**:
```julia
# Call chain that needs debugging:
prepare_koina_batch(InstrumentAgnosticModel(...), data, "QE", batch_size=1)
  → wrapper method (line 76-82)
  → calls prepare_koina_batch(model, data; batch_size=batch_size)
  → generic method (line 87-114) 
  → accesses batch_data.precursor_charge (line 100)
```

**Proposed Fix**:
```julia
# In koina_batch_prep.jl, update the wrapper method:
function prepare_koina_batch(model::InstrumentAgnosticModel,
                           data::DataFrame,
                           instrument_type::String;  # ignored for agnostic models
                           batch_size::Int = 1000)::Vector{String}
    # Ensure column names match what the generic method expects
    required_cols = [:koina_sequence, :precursor_charge, :collision_energy]
    missing_cols = setdiff(required_cols, Symbol.(names(data)))
    
    if !isempty(missing_cols)
        throw(ArgumentError("Missing required columns: $missing_cols"))
    end
    
    # Call the generic method that doesn't need instrument_type
    return prepare_koina_batch(model, data; batch_size=batch_size)
end
```

### Issue 2: Empty Fragment Results (1 error) 🔧
**Error**: `BoundsError: attempt to access 0-element Vector{Int64} at index [1]`

**Affected Test Cases**:
- `End-to-End Prosit Library Building`  
- `Batch Processing Integration`

**Root Cause**:
The `frags_per_precursor` field is 0, indicating no fragments were parsed from the API response. This can happen when:
1. The API call succeeds but returns no fragments for the given peptides
2. The parsing logic fails to extract fragments from a valid response
3. The response structure doesn't match expectations

**Note**: The HTML error messages in logs are likely unrelated to this specific issue.

**Proposed Fix**:
```julia
# In koina_batch_parse.jl, add validation:
function parse_koina_batch(model::InstrumentSpecificModel,
                          response::Dict{String,Any})::KoinaBatchResult{Nothing}
    # Validate response structure
    if !haskey(response, "outputs") || isempty(response["outputs"])
        @warn "Invalid or empty Koina API response"
        return KoinaBatchResult(DataFrame(), 0, nothing)
    end
    
    # Validate shape information exists
    first_output = first(response["outputs"])
    if !haskey(first_output, "shape") || length(first_output["shape"]) < 2
        @warn "Invalid shape information in Koina response"
        return KoinaBatchResult(DataFrame(), 0, nothing)
    end
    
    df = DataFrame()
    n_precs, n_frags = first_output["shape"]
    
    # Continue with existing parsing...
    for col in response["outputs"]
        col_name = Symbol(col["name"])
        if col_name ∈ [:intensities, :mz]
            df[!, col_name] = Float32.(col["data"])
        else
            df[!, :annotation] = string.(col["data"])
        end
    end
    
    return KoinaBatchResult(df, Int64(n_frags), nothing)
end
```

### Issue 3: ~~HTML Error Page Handling~~ (REMOVED)
**Status**: This was a misdiagnosis. The HTML error messages in the logs are not from the Koina API returning HTML responses. They are secondary symptoms or from unrelated parts of the test suite. **No fix needed for this issue.**

## Priority Order

1. **High Priority**: Fix column name preservation (Issue 1)
   - Affects 3 test cases
   - Core functionality issue
   - Prevents API calls from being made

2. **Medium Priority**: Add response validation (Issue 2)
   - Prevents crashes on invalid responses
   - Improves robustness
   - Handles legitimate empty fragment responses

## Testing Strategy

After implementing fixes:
1. Run full test suite to verify fixes
2. Add specific tests for edge cases:
   - Test with missing columns
   - Test with empty API responses
   - Test with HTML error responses
3. Verify no regressions in passing tests

## Expected Outcome

After implementing these fixes:
- All 150 Koina API tests should pass
- Robust handling of API errors
- Clear error messages for debugging
- No silent failures on invalid data