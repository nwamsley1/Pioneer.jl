# Detailed Analysis: ParameterTuningSearch Refactoring and NceTuningSearch Error

## The Error Chain

### What's Happening
1. ParameterTuningSearch completes successfully
2. NceTuningSearch starts
3. NceTuningSearch calls `process_psms!` with a `Vector{Any}` instead of expected `DataFrame`
4. Method error occurs

### The Specific Error
```julia
MethodError: no method matching process_psms!(::Vector{Any}, ::BasicNonIonMobilityMassSpecData{Float32}, 
                                              ::SearchContext{10, ...}, ::NceTuningSearchParameters{...})
```

## Deep Dive: Original Implementation

### Original ParameterTuningSearch Output

Let me examine what ParameterTuningSearch originally provided to downstream methods:

1. **Mass Error Model**: Stored in `search_context.mass_error_model[ms_file_idx]`
2. **RT Model**: Stored in `search_context.rt_irt_map[ms_file_idx]`
3. **IRT Errors**: Stored in `search_context.irt_errors[ms_file_idx]`
4. **NO PSMs**: ParameterTuningSearch does NOT pass PSMs to NceTuningSearch

### How NceTuningSearch Gets Its PSMs

NceTuningSearch must collect its own PSMs! Let's trace this:

1. NceTuningSearch starts with calibrated parameters from ParameterTuningSearch
2. It performs its own library search using those parameters
3. The PSMs it collects are used for NCE model fitting

## The Real Problem

### Looking at NceTuningSearch Code

```julia
# In NceTuningSearch.jl line 181
process_psms!(psms, spectra, search_context, params)
```

The `psms` variable is being passed as `Vector{Any}` instead of `DataFrame`.

### Where Do These PSMs Come From?

Checking NceTuningSearch implementation:
1. It should call its own `collect_psms` or similar function
2. This function should return a DataFrame
3. Something is returning Vector{Any} instead

## Critical Discovery

### The Interface Between Search Methods

Each search method:
1. Uses models from SearchContext (set by previous methods)
2. Collects its own PSMs using those models
3. Does its specific processing
4. Updates SearchContext with new models/results

### What ParameterTuningSearch Should Provide

✅ **What it SHOULD store:**
- Mass error model → `search_context.mass_error_model[ms_file_idx]`
- RT conversion model → `search_context.rt_irt_map[ms_file_idx]`
- IRT error estimates → `search_context.irt_errors[ms_file_idx]`

❌ **What it should NOT store:**
- PSMs (each method collects its own)
- Fragments (temporary, method-specific)

## The Actual Problem Location

The error is NOT in ParameterTuningSearch! It's in how NceTuningSearch is collecting PSMs.

### Investigation Path

1. NceTuningSearch tries to collect PSMs
2. The collection function returns Vector{Any} instead of DataFrame
3. This happens when library_search fails or returns empty results

### Likely Root Cause

The `library_search` function in NceTuningSearch context is:
1. Not finding any PSMs (possibly due to wrong parameters)
2. Returning an empty Vector instead of empty DataFrame
3. Or hitting an error and returning wrong type

## Checking Model Storage

### Are Models Being Stored Correctly?

Let's verify the refactored ParameterTuningSearch stores models properly:

```julia
# In check_and_store_convergence!
setMassErrorModel!(search_context, ms_file_idx, mass_err_model)  ✓

# In set_rt_to_irt_model!
getIrtErrors(search_context)[ms_file_idx] = model[4] * params.irt_tol_sd  ✓

# But where is the RT model stored in SearchContext?
# This might be missing!
```

### Missing RT Model Storage!

The refactored code might not be storing the RT model in SearchContext properly!

```julia
# Original pattern:
setRtIrtMap!(search_context, rt_to_irt_model, ms_file_idx)

# Refactored - need to verify this is called
```

## The Smoking Gun

### In handle_non_convergence!
```julia
function handle_non_convergence!(results, search_context, ms_file_idx, n_attempts, warnings)
    # ...
    setMassErrorModel!(search_context, ms_file_idx, results.mass_err_model[])
    setRtIrtMap!(search_context, results.rt_to_irt_model[], ms_file_idx)  # ✓ Here for fallback
    # ...
end
```

### But in check_and_store_convergence!
```julia
function check_and_store_convergence!(...)
    # ...
    setMassErrorModel!(search_context, ms_file_idx, mass_err_model)  # ✓ Mass error stored
    
    # Store RT model
    rt_model_data = fit_irt_model(params, psms)
    set_rt_to_irt_model!(results, search_context, params, ms_file_idx, rt_model_data)  # Updates results
    # BUT DOES IT UPDATE search_context.rt_irt_map???
end
```

## Critical Finding

### The set_rt_to_irt_model! Function

```julia
function set_rt_to_irt_model!(ptsr::ParameterTuningSearchResults, search_context::SearchContext,
                              params::P, ms_file_idx::Int64, model::Tuple{...})
    ptsr.rt_to_irt_model[] = model[1]  # Updates results
    # ... resize and append ...
    getIrtErrors(search_context)[ms_file_idx] = model[4] * params.irt_tol_sd  # Updates IRT errors
    
    # MISSING: setRtIrtMap!(search_context, model[1], ms_file_idx)
end
```

**THE RT MODEL IS NOT BEING STORED IN SEARCH_CONTEXT!**

## The Complete Picture

### What's Wrong

1. **Primary Issue**: RT model not stored in SearchContext when convergence succeeds
2. **Secondary Effect**: NceTuningSearch can't find PSMs because RT model is missing
3. **Result**: Empty/wrong type returned instead of DataFrame

### Why It Works for Non-Converged Files

Non-converged files go through `handle_non_convergence!` which DOES store the RT model!

### Why It Fails for Converged Files

Converged files go through `check_and_store_convergence!` which does NOT store the RT model in SearchContext!

## Required Fixes

### Fix 1: Update set_rt_to_irt_model!

```julia
function set_rt_to_irt_model!(ptsr::ParameterTuningSearchResults, search_context::SearchContext,
                              params::P, ms_file_idx::Int64, 
                              model::Tuple{SplineRtConversionModel, Vector{Float32}, Vector{Float32}, Float32})
    ptsr.rt_to_irt_model[] = model[1]
    resize!(ptsr.irt, 0)
    resize!(ptsr.rt, 0)
    append!(ptsr.rt, model[2])
    append!(ptsr.irt, model[3])
    
    # CRITICAL FIX: Store RT model in SearchContext
    setRtIrtMap!(search_context, model[1], ms_file_idx)
    
    getIrtErrors(search_context)[ms_file_idx] = model[4] * params.irt_tol_sd
end
```

### Fix 2: Verify store_final_results!

Ensure the final results storage also updates SearchContext:

```julia
function store_final_results!(...)
    # ...
    # Make sure RT model is in SearchContext
    if hasfield(typeof(results), :rt_to_irt_model) && isassigned(results.rt_to_irt_model)
        setRtIrtMap!(search_context, results.rt_to_irt_model[], ms_file_idx)
    end
    # ...
end
```

## Summary

The refactoring broke the critical link between ParameterTuningSearch and downstream methods by:
1. Not storing the RT model in SearchContext for converged files
2. Only storing it for non-converged/fallback cases
3. This causes NceTuningSearch to fail when trying to collect PSMs

The fix is simple: ensure `setRtIrtMap!` is called whenever an RT model is fitted and accepted.