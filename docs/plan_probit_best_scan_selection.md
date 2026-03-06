# Plan: Probit-Based Best Scan Selection in SecondPassSearch

## Motivation

Currently in `SecondPassSearch`, the "best scan" for each precursor is selected by `argmax(weight)` — the scan with the highest deconvolution weight (i.e., the chromatographic apex). This is a single-feature heuristic. A scan at the apex may have poor spectral matching quality, while a neighboring scan with slightly lower intensity could have much better spectral contrast, error norms, or goodness-of-fit.

By training a quick probit model on a handful of discriminative features (using target/decoy labels), we can produce a per-scan probability that integrates multiple quality metrics. We then select the scan with the highest probit score as the "best scan" for each precursor, giving downstream LightGBM training better input data.

## Current Flow (SecondPassSearch.jl, lines 373-560)

```
process_search_results!()
  1. add_second_search_columns!()      # adds rt, target, charge, err_norm, etc.
  2. get_isotopes_captured!()          # adds isotopes_captured, precursor_fraction_transmitted
  3. filter! by precursor_fraction_transmitted
  4. init_summary_columns!()           # pre-allocates max_gof, smoothness, etc.
  5. groupby → get_summary_scores!()   # finds apex by argmax(weight), fills summary stats
  6. filter!(best_scan)                # keeps only apex scan per precursor
  7. add_features!()                   # adds irt_diff, prec_mz, tic, etc.
  8. initialize_prob_group_features!() # adds trace_prob, q_value placeholders
  9. write to Arrow by cv_fold
```

## Proposed Flow

Insert a probit training + scoring step between steps 3 and 4:

```
process_search_results!()
  1. add_second_search_columns!()
  2. get_isotopes_captured!()
  3. filter! by precursor_fraction_transmitted
  ─── NEW: Probit-based best scan selection ───
  3a. Train probit model on ALL scans (target/decoy labels, ~5 features)
  3b. Score all scans with probit model → new column :prescore
  3c. For each precursor group, mark best_scan = scan with max(:prescore)
  ─────────────────────────────────────────────
  4. init_summary_columns!()
  5. groupby → get_summary_scores!()   # NOW uses probit-selected apex
  6. filter!(best_scan)
  7. add_features!()
  8. initialize_prob_group_features!()
  9. write to Arrow by cv_fold
```

**Key change to `get_summary_scores!`**: The function currently does `apex_scan = argmax(weight)`. After this change, the apex scan is pre-selected via `best_scan` by the probit step, so `get_summary_scores!` should use the pre-marked `best_scan` index instead of recomputing it from weight.

## Features Available for Probit (per-scan, before summarization)

These columns exist on every scan row after step 3:

| Column | Type | Description |
|--------|------|-------------|
| `weight` | Float32 | Deconvolution weight (intensity) |
| `gof` | Float16 | Goodness of fit |
| `matched_ratio` | Float16 | Fraction of ions matched |
| `fitted_spectral_contrast` | Float16 | Spectral contrast from fitting |
| `fitted_manhattan_distance` | Float16 | Manhattan distance metric |
| `scribe` | Float16 | Spectral scoring metric |
| `err_norm` | Float16 | Normalized mass error |
| `log2_intensity_explained` | Float16 | Log2 fraction of intensity explained |
| `y_count` | UInt8 | Number of y ions matched |
| `b_count` | UInt8 | Number of b ions matched |
| `total_ions` | UInt16 | y_count + b_count + isotope_count |
| `charge` | UInt8 | Precursor charge |
| `target` | Bool | Target (true) or decoy (false) |

### Recommended Probit Feature Set (~5-7 features)

```julia
PRESCORE_FEATURES = [
    :fitted_spectral_contrast,   # how well the spectrum matches
    :err_norm,                   # mass accuracy
    :log2_intensity_explained,   # fraction of signal explained
    :gof,                        # goodness of fit
    :scribe,                     # spectral quality score
    :weight,                     # intensity (keep as one signal among many)
]
```

These are the most discriminative per-scan features. They're all available before any groupby/summarization. The probit model learns the relative importance of each, so even if `weight` is included, the model can downweight it relative to spectral quality.

## Detailed Code Changes

### File 1: `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl`

#### Change A: New function `train_and_apply_prescore!()`

Add after `init_summary_columns!` definition (~line 1215). This function:
1. Trains a probit model on all scans using target/decoy labels
2. Applies predictions to get a per-scan probability
3. Marks `best_scan = true` for the highest-scoring scan in each precursor group

```julia
"""
    train_and_apply_prescore!(psms::DataFrame)

Train a lightweight probit model on per-scan features to produce a preliminary
quality score (`prescore`) for each scan. Then for each precursor group, mark
the scan with the highest prescore as `best_scan = true`.

This replaces the default `argmax(weight)` apex selection with a multi-feature
quality-based selection.
"""
function train_and_apply_prescore!(
    psms::DataFrame;
    features::Vector{Symbol} = [
        :fitted_spectral_contrast,
        :err_norm,
        :log2_intensity_explained,
        :gof,
        :scribe,
        :weight,
    ],
    n_train_rounds::Int = 2,
    max_iter::Int = 20
)
    n = nrow(psms)
    if n < 100
        # Too few PSMs for meaningful probit training — fall back to weight-based
        return false
    end

    # Prepare labels
    targets = psms[!, :target]

    # Prepare feature matrix (convert to Float64 for probit)
    X = DataFrame(Matrix{Float64}(psms[!, features]), features)

    # Partition for parallel IRLS
    chunk_size = max(1, n ÷ (10 * Threads.nthreads()))
    data_chunks = collect(Iterators.partition(1:n, chunk_size))

    # Train probit model iteratively
    β = zeros(Float64, length(features))
    scores = zeros(Float32, n)

    for round in 1:n_train_rounds
        β = ProbitRegression(β, X, targets, data_chunks; max_iter=max_iter)

        # Score all scans
        ModelPredict!(scores, X, β, data_chunks)

        if round < n_train_rounds
            # Compute q-values for next round's training data selection
            # (only train on high-confidence targets + all decoys)
            q_vals = zeros(Float64, n)
            get_qvalues!(scores, targets, q_vals)

            # Filter training data: targets at q <= 0.01 + all decoys
            train_mask = (q_vals .<= 0.01 .&& targets) .|| .!targets
            X_train = X[train_mask, :]
            targets_train = targets[train_mask]
            chunks_train = collect(Iterators.partition(1:sum(train_mask),
                                   max(1, sum(train_mask) ÷ (10 * Threads.nthreads()))))
            β = ProbitRegression(β, X_train, targets_train, chunks_train; max_iter=max_iter)
        end
    end

    # Final scoring with probabilities
    ModelPredictProbs!(scores, X, β, data_chunks)

    # Store prescores
    psms[!, :prescore] = scores

    # Mark best_scan based on highest prescore per precursor group
    psms[!, :best_scan] .= false
    for gpsms in groupby(psms, :precursor_idx)
        best_idx = argmax(gpsms[!, :prescore])
        gpsms[best_idx, :best_scan] = true
    end

    return true
end
```

#### Change B: Modify `get_summary_scores!()` (~line 1247)

Change apex selection from `argmax(weight)` to use pre-marked `best_scan` if available, falling back to `argmax(weight)`.

**Current code (line 1247):**
```julia
apex_scan = argmax(psms[!,:weight])
```

**New code:**
```julia
# Use probit-selected best_scan if available, otherwise fall back to max weight
best_scan_col = psms[!, :best_scan]
if any(best_scan_col)
    apex_scan = findfirst(best_scan_col)
else
    apex_scan = argmax(psms[!, :weight])
end
```

The rest of `get_summary_scores!` remains unchanged — it still computes summary stats in a ±2 window around the apex and sets `best_scan[apex_scan] = true`.

### File 2: `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/SecondPassSearch.jl`

#### Change C: Call `train_and_apply_prescore!()` in `process_search_results!()`

Insert between the `filter!` by `precursor_fraction_transmitted` (line 415) and `init_summary_columns!` (line 420).

**Current code (lines 415-420):**
```julia
        filter!(row -> row.precursor_fraction_transmitted >= params.min_fraction_transmitted, psms)

        # Initialize columns for best scan selection and summary statistics
        psms[!,:best_scan] = zeros(Bool, size(psms, 1));
        init_summary_columns!(psms);
```

**New code:**
```julia
        filter!(row -> row.precursor_fraction_transmitted >= params.min_fraction_transmitted, psms)

        # Pre-score all scans with probit model for quality-based best scan selection
        psms[!,:best_scan] = zeros(Bool, size(psms, 1));
        train_and_apply_prescore!(psms)

        # Initialize columns for best scan selection and summary statistics
        init_summary_columns!(psms);
```

#### Change D: Remove prescore column before writing Arrow

The `:prescore` column is a temporary internal column. Remove it before writing to avoid schema issues across files.

Insert before the Arrow write block (~line 516):

```julia
        # Remove temporary prescore column if present
        if hasproperty(psms, :prescore)
            select!(psms, Not(:prescore))
        end
```

## What Does NOT Change

- **`add_features!()`**: Unchanged. Still operates on the filtered best_scan PSMs.
- **`initialize_prob_group_features!()`**: Unchanged. Still initializes trace_prob/q_value.
- **ScoringSearch / LightGBM training**: Unchanged. Receives the same Arrow schema, just with potentially different (better) scans selected.
- **Arrow file schema**: Unchanged. The prescore column is dropped before writing.
- **CV fold assignment**: Unchanged. cv_fold is assigned per precursor in the library.
- **Summary statistics (max_gof, smoothness, etc.)**: Still computed by `get_summary_scores!()` around the apex, just the apex may be a different scan now.

## Edge Cases

1. **Too few PSMs (<100)**: `train_and_apply_prescore!` returns `false`, `best_scan` stays all-false, and `get_summary_scores!` falls back to `argmax(weight)` as before.
2. **Single-scan precursors**: Probit selection is irrelevant — only one scan to pick. Both old and new logic select it.
3. **All-decoy files**: Probit model may be poorly trained, but `argmax(prescore)` still picks the best-scoring scan. Summary stats are still computed correctly.
4. **Feature type conversion**: Features like `gof` (Float16) and `weight` (Float32) need to be converted to Float64 for the probit regression. The DataFrame constructor `Matrix{Float64}(psms[!, features])` handles this.

## Performance Impact

- **Probit training**: ~20 IRLS iterations on ~5 features. With parallel chunking, this trains in <1 second even on 100K scans. Runs once per MS file.
- **Memory**: One additional Float32 column (`prescore`) on the pre-filtered PSM DataFrame. Dropped before Arrow write.
- **Net effect on pipeline**: Minimal wall-clock impact. The probit model is orders of magnitude cheaper than the downstream LightGBM training.

## Verification

1. Run `SearchDIA("./data/ecoli_test/ecoli_test_params.json")` — integration test should pass.
2. Compare precursor counts at 1% FDR with and without the change (use the new per-file logging).
3. Check that the selected best_scans differ from the old weight-based selection in some fraction of cases (log the count of changed selections).
4. Verify Arrow file schemas are identical before and after (no prescore column leaking through).

## Future Extensions

- **Configurable feature set**: Expose `PRESCORE_FEATURES` as a parameter.
- **Per-file vs global model**: Currently trains one probit model per file. Could train a single model across all files in a first pass, then apply per-file.
- **Feature gating**: Only use prescore selection when the probit model achieves meaningful discrimination (e.g., AUC > 0.6). Otherwise fall back to weight.
