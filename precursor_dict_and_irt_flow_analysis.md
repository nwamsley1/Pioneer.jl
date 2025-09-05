# Precursor Dictionary Building and iRT Flow Analysis

## Executive Summary

This document provides a comprehensive analysis of how the `precursor_dict` is built in FirstPassSearch and how the `setPredIrt!` commands flow through the search pipeline in Pioneer.jl. Understanding this flow is critical for diagnosing issues with match-between-runs (MBR) functionality and conservative scoring behavior.

## Overview of Data Flow

```
FirstPassSearch:
├── get_best_precursors_accross_runs() → precursor_dict
├── Partner Addition (if MBR enabled) → augmented precursor_dict
├── setPredIrt!(search_context, pid, irt_value) → irt_obs[pid]
└── setPrecursorDict!(search_context, precursor_dict)

SecondPassSearch:
├── getPredIrt(search_context, prec_idx) → retrieves irt_obs[pid]
├── Uses predicted iRT for irt_error calculations
└── Adds :irt_pred column to PSM dataframes

Downstream Usage:
├── Normalization uses :irt_obs for binning
├── Quantification uses iRT differences for quality metrics
└── Protein scoring may use iRT-related features
```

## Detailed Analysis

### 1. Building the Precursor Dictionary

#### 1.1 Function: `get_best_precursors_accross_runs`

**File:** `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/getBestPrecursorsAccrossRuns.jl:49-229`

**Purpose:** Aggregates the best precursor matches across all MS runs to build a comprehensive precursor dictionary for downstream searches.

#### 1.2 Input Parameters

```julia
function get_best_precursors_accross_runs(
    psms_paths::Vector{String},          # Paths to PSM files from first pass search
    prec_mzs::AbstractVector{Float32},   # Vector of precursor m/z values
    rt_irt::Dict{Int64, RtConversionModel}; # RT-iRT conversion models per file
    max_q_val::Float32 = 0.01f0          # Q-value threshold for inclusion
)
```

#### 1.3 Data Structure Definition

The function returns a dictionary with this complex NamedTuple structure:

```julia
Dictionary{UInt32, @NamedTuple{
    best_prob::Float32,                   # Highest probability score across runs
    best_ms_file_idx::UInt32,            # File index containing best match
    best_scan_idx::UInt32,               # Scan index of best match
    best_irt::Float32,                   # iRT value of best match
    mean_irt::Union{Missing, Float32},   # Mean iRT across qualifying matches
    var_irt::Union{Missing, Float32},    # Variance in iRT across qualifying matches
    n::Union{Missing, UInt16},           # Number of qualifying matches
    mz::Float32                          # Precursor m/z value
}}
```

#### 1.4 Two-Pass Algorithm

**First Pass: Data Collection (`readPSMs!` - lines 56-129)**

For each PSM file:
1. **Read PSM data:** Extract precursor indices, q-values, probabilities, retention times, scan indices, and file indices
2. **RT-to-iRT conversion:** Convert retention times to iRT using `rt_irt[key](rts[row])`
3. **Quality filtering:** Only include PSMs with `q_value <= max_q_val` in statistics
4. **Statistical accumulation:**
   - Track best probability score and associated metadata
   - Maintain running sum of iRT values for mean calculation
   - Count number of qualifying matches per precursor

**Algorithm Logic:**
```julia
for row in eachindex(precursor_idxs)
    q_value = q_values[row]
    precursor_idx = precursor_idxs[row]
    prob = probs[row]
    irt = rt_irt(rts[row])  # Convert RT to iRT
    
    # Only count toward statistics if passes q-value threshold
    passed_q_val = (q_value <= max_q_val)
    n = passed_q_val ? one(UInt16) : zero(UInt16)
    mean_irt = passed_q_val ? irt : zero(Float32)
    
    if haskey(prec_to_best_prob, precursor_idx)
        # Update existing entry
        old_entry = prec_to_best_prob[precursor_idx]
        
        # Update best match if current is better
        if (old_entry.best_prob < prob)
            best_prob = prob
            best_irt = irt
            best_scan_idx = scan_idx
            best_ms_file_idx = ms_file_idx
        else
            # Keep existing best match info
            best_prob = old_entry.best_prob
            # ... (keep other best_ fields)
        end
        
        # Update running statistics
        mean_irt += old_entry.mean_irt
        n += old_entry.n
        
    else
        # Create new entry
        # ... (initialize with current values)
    end
end
```

**Key Insights:**
- **Best match tracking:** Maintains the highest probability match across all runs
- **Statistical computation:** Only PSMs passing q-value threshold contribute to mean iRT
- **Memory efficiency:** Updates existing entries rather than storing all observations

**Filtering Step (lines 202-211):**
```julia
# Filter to top N precursors by probability
sort!(prec_to_best_prob, by = x->x[:best_prob], alg=PartialQuickSort(1:max_precursors), rev = true)
N = 0
for key in collect(keys(prec_to_best_prob))
    N += 1
    if N > max_precursors
        delete!(prec_to_best_prob, key)
    end
end
```

This step limits the dictionary to the `max_precursors` best-scoring precursors to control memory usage and processing time.

**Second Pass: Variance Calculation (`getVariance!` - lines 130-171)**

For remaining precursors, calculate iRT variance:

```julia
for row in eachindex(precursor_idxs)
    if q_value > max_q_val
        continue  # Skip low-quality PSMs
    end
    
    if haskey(prec_to_best_prob, precursor_idx)
        # Update variance using: var += (x - mean)^2
        entry = prec_to_best_prob[precursor_idx]
        var_irt += (irt - entry.mean_irt/entry.n)^2
        # Update entry with new variance
    end
end
```

**Result:** A comprehensive dictionary containing the best precursor matches with complete statistical summaries.

### 2. Partner Addition Logic (MBR Enhancement)

**File:** `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl:502-534`

When `match_between_runs` is enabled, partner precursors are added:

```julia
if params.match_between_runs == true  # Currently hardcoded to true==true
    precursors = getPrecursors(getSpecLib(search_context))
    
    for (pid, val) in pairs(precursor_dict)
        setPredIrt!(search_context, pid, getIrt(precursors)[pid])
        partner_pid = getPartnerPrecursorIdx(precursors)[pid]
        
        if ismissing(partner_pid)
            continue  # No partner exists
        end
        
        if !haskey(precursor_dict, partner_pid)
            # Add partner with same statistics as original
            insert!(precursor_dict, partner_pid, val)
            # Give partner the iRT of the original precursor
            setPredIrt!(search_context, partner_pid, getIrt(precursors)[pid])
        else
            # Partner already exists, use its original predicted iRT
            setPredIrt!(search_context, partner_pid, getIrt(precursors)[partner_pid])
        end
    end
else
    # Non-MBR mode: just set predicted iRT to library values
    for (pid, val) in pairs(precursor_dict)
        setPredIrt!(search_context, pid, getIrt(precursors)[pid])
    end
end
```

**Critical Impact Analysis:**

1. **Dictionary Size Changes:** When MBR is enabled, the precursor dictionary can nearly double in size as partners are added
2. **iRT Assignment Strategy:**
   - **New partners:** Receive the iRT of their "discovered" counterpart
   - **Existing partners:** Keep their original library-predicted iRT
3. **Statistical Propagation:** Partner precursors inherit the statistical summary (best_prob, mean_irt, var_irt, etc.) from their counterpart

### 3. iRT Storage in SearchContext

#### 3.1 Storage Mechanism

**File:** `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl:403-404`

```julia
setPredIrt!(s::SearchContext, prec_idx::Int64, irt::Float32) = s.irt_obs[prec_idx] = irt
setPredIrt!(s::SearchContext, prec_idx::UInt32, irt::Float32) = s.irt_obs[prec_idx] = irt
```

**SearchContext Structure (lines 227, 267):**
```julia
struct SearchContext{N,L,M}
    # ... other fields ...
    irt_obs::Dict{UInt32, Float32}  # Storage for predicted iRT values
    # ... other fields ...
end

# Constructor initializes empty dictionary
irt_obs = Dict{UInt32, Float32}()
```

#### 3.2 Access Pattern

**Retrieval Functions (lines 399-401):**
```julia
getPredIrt(s::SearchContext) = s.irt_obs                     # Get entire dictionary
getPredIrt(s::SearchContext, prec_idx::Int64) = s.irt_obs[prec_idx]   # Get single value
getPredIrt(s::SearchContext, prec_idx::UInt32) = s.irt_obs[prec_idx]  # Get single value
```

### 4. Usage in Downstream Search Methods

#### 4.1 SecondPassSearch Integration

**File:** `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl:855-863`

In `add_main_search_columns!`:

```julia
function add_main_search_columns!(...)
    # ... initialization ...
    
    for i in 1:N
        # ... other column assignments ...
        
        prec_idx = precursor_idx[i]
        irt_obs[i] = rt_to_irt_interp(rt[i])                           # Convert observed RT to iRT
        irt_pred[i] = getPredIrt(search_context, prec_idx)             # Get predicted iRT from FirstPass
        irt_diff[i] = abs(irt_obs[i] - prec_id_to_irt[prec_idx].best_irt)  # Difference from empirical best
        
        if !ms1_missing[i]
            ms1_irt_diff[i] = abs(rt_to_irt_interp(ms1_rt[i]) - getPredIrt(search_context, prec_idx))
        end
        
        irt_error[i] = abs(irt_obs[i] - irt_pred[i])                   # Error from predicted iRT
    end
    
    # Add columns to PSM dataframe
    psms[!,:irt_obs] = irt_obs
    psms[!,:irt_pred] = irt_pred  # This uses the setPredIrt! values from FirstPass
    psms[!,:irt_diff] = irt_diff
    # ... other columns ...
end
```

**Key Usage Patterns:**

1. **`:irt_pred` column:** Contains the predicted iRT values set by `setPredIrt!` in FirstPassSearch
2. **`:irt_error` calculation:** Computes absolute difference between observed and predicted iRT
3. **MS1 iRT differences:** Uses predicted iRT for MS1-level retention time validation

#### 4.2 Normalization and Quantification

**File:** `src/Routines/SearchDIA/CommonSearchUtils/normalizeQuant.jl:28-42`

```julia
function get_quantification_correction(...)
    # Sort by observed iRT (which comes from :irt_obs column)
    sort!(psms, :irt_obs, alg = QuickSort)
    
    # Use iRT values for binning
    if minimum(psms[!,:irt_obs]) < min_rt
        min_rt = minimum(psms[!,:irt_obs])
    end
    
    # Create RT bins based on iRT values
    for bin_idx in 1:n_bins
        bin_start = (bin_idx - 1)*bin_size + 1
        bin_stop = bin_idx*bin_size
        median_rts[bin_idx] = median(@view(psms[bin_start:bin_stop,:irt_obs]))
        # ... quantification correction logic ...
    end
end
```

**Normalization Impact:**
- Uses `:irt_obs` (observed iRT) for retention time binning
- Predicted iRT values (`:irt_pred`) provide reference points for quality assessment
- iRT errors help identify problematic identifications

### 5. Critical Data Flow Analysis

#### 5.1 Normal Flow (Non-MBR Mode)

```
1. FirstPassSearch runs individual file searches
2. get_best_precursors_accross_runs() collects best matches
3. For each precursor in dictionary:
   setPredIrt!(search_context, pid, getIrt(precursors)[pid])
   # Sets predicted iRT to library value
4. SecondPassSearch uses getPredIrt() for :irt_pred column
5. Downstream methods use :irt_pred vs :irt_obs comparisons
```

#### 5.2 MBR Flow (Critical Differences)

```
1. FirstPassSearch runs individual file searches
2. get_best_precursors_accross_runs() collects best matches
3. Partner Addition Loop:
   For each precursor (pid, val) in original dictionary:
     setPredIrt!(search_context, pid, library_irt[pid])
     
     partner_pid = getPartner(pid)
     if partner not in dictionary:
       insert!(dictionary, partner_pid, val)  # Same statistics!
       setPredIrt!(search_context, partner_pid, library_irt[pid])  # Same iRT!
     else:
       setPredIrt!(search_context, partner_pid, library_irt[partner_pid])
4. Dictionary may double in size
5. Many partners have "borrowed" statistics and iRT values
6. SecondPassSearch processes larger precursor set with mixed predicted iRTs
```

#### 5.3 Impact on Statistical Integrity

**Issue 1: Statistical Contamination**
- Partners inherit statistical summaries (mean_irt, var_irt, best_prob) from their counterparts
- This creates artificial correlations in the data
- PSMs for partners may be scored using statistics from different precursors

**Issue 2: iRT Assignment Inconsistencies**
- Newly added partners receive the library iRT of their discoverer, not their own
- This can lead to systematic iRT errors for partner precursors
- Downstream `:irt_error` calculations may be biased

**Issue 3: Dictionary Composition Changes**
- The precursor dictionary has fundamentally different composition between MBR/non-MBR modes
- This affects all downstream processing that depends on the dictionary
- Training data for ML models becomes different even before feature selection

## Summary and Implications

### Key Findings

1. **Complex Two-Pass Algorithm:** The precursor dictionary building process is sophisticated, involving statistical aggregation across runs and variance calculations.

2. **Partner Addition Doubles Dictionary:** When MBR is enabled, the dictionary can nearly double in size through partner addition, fundamentally changing the search space.

3. **Shared Statistics:** Partners inherit statistical summaries from their counterparts, creating artificial data dependencies.

4. **Inconsistent iRT Assignment:** Partner iRT assignment strategy differs for new vs. existing partners, potentially introducing systematic errors.

5. **Pipeline-Wide Impact:** The `irt_obs` dictionary in SearchContext flows through SecondPassSearch, quantification, and normalization, amplifying any biases introduced early.

### Root Cause Analysis

The conservative results observed with MBR enabled likely stem from:

1. **Early Data Contamination:** Partners are added with "borrowed" statistics before any ML training occurs
2. **Statistical Dependencies:** Artificial correlations introduced through shared statistical summaries
3. **iRT Bias:** Inconsistent predicted iRT assignment may systematically bias error calculations
4. **Training Data Modification:** The fundamentally different precursor dictionary composition changes training data even when MBR features aren't used

### Recommendations

1. **Delay Partner Addition:** Move partner addition to after initial model training or make it conditional on MBR feature usage
2. **Independent Statistics:** Compute separate statistics for partners rather than inheriting from counterparts
3. **Consistent iRT Assignment:** Use library-predicted iRT values consistently for all precursors
4. **Data Integrity Checks:** Add validation to ensure training data consistency between MBR/non-MBR modes

This analysis reveals that the conservative behavior isn't just about feature selection—it's about fundamental changes to the data structures that flow through the entire pipeline.