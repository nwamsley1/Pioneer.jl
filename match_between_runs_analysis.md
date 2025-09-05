# Match Between Runs (MBR) Parameter Analysis

## Overview

This document analyzes all code locations affected by the `match_between_runs` parameter in Pioneer.jl and investigates why enabling MBR leads to conservative results with entrapment sequences, even when MBR features are disabled or the final MBR iteration is prevented.

## Executive Summary

The conservative behavior when MBR is enabled appears to stem from **data contamination during the training phase**, not from the MBR features themselves. Key issues identified:

1. **Early Partner Addition**: Partners are added to precursor_dict in FirstPassSearch before ML training
2. **CV Fold Contamination**: Target/decoy pairs can end up in the same CV fold
3. **Training Data Bias**: Modified training datasets affect model learning even without MBR features
4. **Feature Correlation**: MBR-specific features are initialized and may interact with other features

## Detailed Code Analysis

### 1. Parameter Definition and Default Values

**Files Affected:**
- `src/Routines/SearchDIA/ParseInputs/paramDefaults.jl:36`
- All example parameter files in `assets/example_config/` and `data/example_config/`

**Default Value:** `true`

The parameter is enabled by default across all configurations, indicating it's expected to be the standard operating mode.

### 2. Early Impact: SearchDIA Main Entry Point

**File:** `src/Routines/SearchDIA.jl:257-264`

```julia
if length(MS_TABLE_PATHS) == 1 && params.global_settings.match_between_runs
    @user_warn "Only one run detected; disabling match_between_runs (MBR)."
    params_dict = JSON.parse(params_string)
    params_dict["global"]["match_between_runs"] = false
    params_string = JSON.json(params_dict, 4)
    updated_global = (; params.global_settings..., match_between_runs=false)
    params = PioneerParameters(updated_global, ...)
end
```

**Impact:** MBR is automatically disabled for single-run experiments, suggesting the feature requires multiple runs to function properly.

### 3. Critical Impact: FirstPassSearch Partner Addition

**File:** `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl:502-534`

```julia
if params.match_between_runs==true
    #######
    #Each target has a corresponding decoy and vice versa
    #Add the complement targets/decoys to the precursor dict 
    #if the `sibling_peptide_scores` parameter is set to true
    #In the target/decoy scoring (see SearchMethods/ScoringSearch)
    #the maximum score for each target/decoy pair is shared accross runs
    #in an iterative training scheme. 
    precursors = getPrecursors(getSpecLib(search_context))
    for (pid, val) in pairs(precursor_dict)
        partner_pid = getPartnerPrecursorIdx(precursors)[pid]
        if ismissing(partner_pid)
            continue
        end

        # If the partner needs to be added, then give it the irt of the currently identified precursor
        # Otherwise if the partner was ID'ed, it should keep its original predicted iRT
        if !haskey(precursor_dict, partner_pid)
            insert!(precursor_dict, partner_pid, val)
            setPredIrt!(search_context, partner_pid, getIrt(getPrecursors(getSpecLib(search_context)))[pid])
        else
            setPredIrt!(search_context, partner_pid, getIrt(getPrecursors(getSpecLib(search_context)))[partner_pid])
        end
    end
end
```

**Critical Finding:** This code adds partner precursors (target/decoy pairs) to the precursor dictionary BEFORE any ML training occurs. This means:

1. **More precursors enter the pipeline** when MBR is enabled
2. **Partner relationships are established early** in the pipeline
3. **Training datasets are fundamentally different** between MBR and non-MBR modes

### 4. CV Fold Assignment Analysis

**File:** `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl:584-592`

```julia
fold1 = getCvFold(precursors, pid)
fold2 = getCvFold(precursors, partner_pid)

if fold1 == fold2
    same_fold_pairs += 1
else
    split_fold_pairs += 1
end
```

**CV Fold Assignment Logic (src/structs/LibraryIon.jl:689-693):**
```julia
for pid in range(1, n)
    pid_to_cv_fold[pid] = pg_to_cv_fold[accession_numbers[pid]]
end
```

**Key Finding:** CV folds are assigned based on **protein groups**, not precursor pairs. This means:
- Target and decoy precursors from the same protein can end up in the same CV fold
- This creates **data leakage** where training and test sets contain related information
- The empirical analysis in FirstPassSearch confirms that pairs can be in the same fold

### 5. Feature Initialization Impact

**File:** `src/utils/ML/percolatorSortOf.jl:616-637`

```julia
function initialize_prob_group_features!(
    psms::AbstractDataFrame,
    match_between_runs::Bool
)
    n = nrow(psms)
    psms[!, :prob]      = zeros(Float32, n)
    psms[!, :q_value]   = zeros(Float64, n)

    if match_between_runs
        psms[!, :MBR_max_pair_prob]             = zeros(Float32, n)
        psms[!, :MBR_best_irt_diff]             = zeros(Float32, n)
        psms[!, :MBR_log2_weight_ratio]         = zeros(Float32, n)
        psms[!, :MBR_log2_explained_ratio]      = zeros(Float32, n)
        psms[!, :MBR_rv_coefficient]            = zeros(Float32, n)
        psms[!, :MBR_is_best_decoy]             = trues(n)
        psms[!, :MBR_num_runs]                  = zeros(Int32, n)
        psms[!, :MBR_transfer_candidate]        = falses(n)
        psms[!, :MBR_is_missing]                = falses(n)
    end

    return psms
end
```

**Impact:** When MBR is enabled, 9 additional columns are added to the PSM dataframe, even if these features aren't used in training. This can affect:
- Memory layout and access patterns
- DataFrame operations and performance
- Potential feature correlations during model training

### 6. Training Data Selection Across Iterations

**File:** `src/utils/ML/percolatorSortOf.jl:385-396`

```julia
for (itr, num_round) in enumerate(iter_scheme)
    psms_train_itr = get_training_data_for_iteration!(psms_train, 
                                                      itr,
                                                      match_between_runs, 
                                                      max_q_value_xgboost_rescore,
                                                      max_q_value_xgboost_mbr_rescore,
                                                      min_PEP_neg_threshold_xgboost_rescore,
                                                      itr >= length(iter_scheme))
    ###################
    #Train a model on the n-1 training folds.
    train_feats = itr < length(iter_scheme) ? non_mbr_features : features
```

**Key Finding:** Feature selection differs by iteration:
- **Early iterations (1 to N-1):** Use `non_mbr_features` (excludes MBR features)
- **Final iteration (N):** Uses all `features` (includes MBR features if enabled)

**File:** `src/utils/ML/percolatorSortOf.jl:669-707`

```julia
# Also train on top scoring MBR candidates if requested
if match_between_runs && last_iter
    # Determine prob threshold for precursors passing the q-value threshold
    max_prob_threshold = minimum(
        psms_train_itr.prob[
            psms_train_itr.target .& (psms_train_itr.q_value .<= max_q_value_xgboost_rescore)
        ]
    )
    
    # ... additional MBR candidate selection logic ...
    
    # Take all decoys and targets passing q_thresh (all 0's now) or mbr_q_thresh
    psms_train_itr = subset(
        psms_train_itr,
        [:target, :q_value, :MBR_is_best_decoy, :MBR_is_missing] => ByRow((t, q, MBR_d, im) -> (!t) || (t && !im && !MBR_d && q <= max_q_value_xgboost_mbr_rescore))
    )
else
    # Take all decoys and targets passing q_thresh
    psms_train_itr = subset(
        psms_train_itr,
        [:target, :q_value] => ByRow((t,q) -> (!t) || (t && q <= max_q_value_xgboost_rescore))
    )
end
```

**Critical Impact:** Training data selection is different when MBR is enabled, even in non-final iterations:
- Different PSM subsets are selected for training
- MBR-specific filtering logic is applied
- This affects model learning from the beginning, not just the final iteration

### 7. Pair-Based Feature Computation

**File:** `src/utils/ML/percolatorSortOf.jl:534-536`

```julia
function summarize_precursors!(psms::AbstractDataFrame; q_cutoff::Float32 = 0.01f0)
    # Compute pair specific features that rely on decoys and chromatograms
    pair_groups = collect(pairs(groupby(psms, [:pair_id, :isotopes_captured])))
```

**Impact:** When MBR is enabled, precursors are processed in pair groups, which affects:
- Feature computation logic
- Memory access patterns  
- Parallel processing behavior

### 8. Score Propagation and Final Assignment

**File:** `src/utils/ML/percolatorSortOf.jl:156-171`

```julia
if match_between_runs
    # Determine which precursors failed the q-value cutoff prior to MBR
    qvals_prev = Vector{Float32}(undef, length(prob_test))
    get_qvalues!(prob_test, psms.target, qvals_prev)
    pass_mask = (qvals_prev .<= max_q_value_xgboost_rescore)
    prob_thresh = any(pass_mask) ? minimum(prob_test[pass_mask]) : typemax(Float32)
    # Label as transfer candidates only those failing the q-value cutoff but
    # whose best matched pair surpassed the passing probability threshold.
    psms[!, :MBR_transfer_candidate] .= (prob_test .< prob_thresh) .&
                                        (psms.MBR_max_pair_prob .>= prob_thresh)

    # Use the final MBR probabilities for all precursors
    psms[!, :prob] = MBR_estimates
else
    psms[!, :prob] = prob_test
end
```

**Impact:** Final probability assignment differs between MBR and non-MBR modes, with MBR using `MBR_estimates` instead of `prob_test`.

## Hypotheses for Conservative Results

### Primary Hypothesis: Training Data Contamination

**Root Cause:** The addition of partner precursors in FirstPassSearch (line 502-534) creates a fundamentally different training dataset when MBR is enabled, leading to biased model learning.

**Mechanism:**
1. **Early Partner Addition:** Partners are added before any ML training occurs
2. **Increased Dataset Size:** More precursors enter the pipeline with MBR enabled
3. **Altered Training Dynamics:** Models learn on different data distributions
4. **Conservative Bias:** Training on partner-augmented datasets may lead to more conservative scoring

### Secondary Hypothesis: CV Fold Data Leakage

**Root Cause:** CV folds are assigned based on proteins, not precursor pairs, allowing target/decoy pairs to end up in the same fold.

**Mechanism:**
1. **Same-Fold Pairs:** Target and decoy from same protein can be in same CV fold
2. **Information Leakage:** Training set contains information about test set through partner relationships
3. **Conservative Scoring:** Models become more conservative to avoid false positives from leaked information

### Tertiary Hypothesis: Feature Space Modification

**Root Cause:** MBR initialization adds 9 additional columns to the dataframe, potentially affecting model behavior even when these features aren't used.

**Mechanism:**
1. **Memory Layout Changes:** Additional columns alter DataFrame structure
2. **Implicit Correlations:** Initialized MBR columns may correlate with other features
3. **Training Artifacts:** Model learning may be subtly affected by presence of unused features

### Quaternary Hypothesis: Training Algorithm Differences

**Root Cause:** Different training data selection logic when MBR is enabled affects model learning throughout all iterations.

**Mechanism:**
1. **Different Subsets:** MBR-specific filtering creates different training subsets
2. **Early Iteration Effects:** Even non-MBR iterations are affected by MBR parameter
3. **Cascading Impact:** Early model differences compound through iterative training

## Evidence Supporting Primary Hypothesis

1. **User's Observation:** Conservative results occur even when MBR features are commented out
2. **Code Evidence:** Partner addition happens before any feature selection or training
3. **Iteration Independence:** Issue persists even when final MBR iteration is disabled
4. **Data Flow:** FirstPassSearch → SecondPassSearch → ScoringSearch pipeline means early changes affect all downstream components

## Recommended Diagnostics

### 1. Precursor Count Analysis
```julia
# Compare precursor counts between MBR and non-MBR modes
count_precursors_before_after_partner_addition()
```

### 2. CV Fold Distribution Analysis
```julia
# Analyze same-fold vs split-fold pair distributions
analyze_cv_fold_pair_distributions()
```

### 3. Training Data Composition Analysis  
```julia
# Compare training datasets between MBR and non-MBR modes
compare_training_data_composition()
```

### 4. Model Performance Analysis
```julia
# Compare model performance metrics between modes
analyze_model_performance_differences()
```

## Recommended Fixes

### 1. Delay Partner Addition
Move partner precursor addition from FirstPassSearch to after initial model training, or make it conditional on actually using MBR features.

### 2. Improve CV Fold Assignment
Ensure target/decoy pairs are always assigned to different CV folds to prevent data leakage.

### 3. Feature Space Isolation
Only initialize MBR features when they will actually be used in training.

### 4. Training Data Consistency
Ensure training data selection logic is identical between MBR and non-MBR modes for non-MBR iterations.

## Conclusion

The conservative results when MBR is enabled likely stem from **early partner addition in FirstPassSearch** rather than the MBR features themselves. This creates a fundamentally different training environment that biases the models toward conservative scoring. The issue is compounded by potential CV fold data leakage and different training data selection logic.

The most promising fix would be to delay partner addition until it's actually needed for MBR feature computation, or to make the partner addition conditional on whether MBR features will be used in the final iteration.