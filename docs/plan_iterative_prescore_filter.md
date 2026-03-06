# Plan: Iterative DIA-NN-Style Prescore Filtering in SecondPassSearch

## Goal

Replace the current single-pass `process_search_results!` flow with a **two-round iterative prescore filter** that progressively eliminates low-quality (precursor_idx, scan_idx) pairs *before* the final feature collection, best-PSM reduction, and Arrow write. The key idea comes from DIA-NN: train a difference-covariance linear classifier on paired target-decoy features, then a probit model, and use the resulting q-values to discard unlikely PSMs. Repeat this filtering twice. Only the survivors go through final processing.

**Why re-run the search instead of filtering in-place:** The Huber-regularized deconvolution solves a joint linear system where all precursors' fragments in a scan contribute to the same design matrix. Removing a precursor changes the weights and residuals for all remaining precursors in that scan. Therefore, after each filter round we must re-run `perform_second_pass_search` with the reduced (precursor, scan) set so the deconvolution is re-solved with the correct competition.

## Current Flow (Single Pass)

```
process_file!
  -> perform_second_pass_search()       # fragment index search -> raw PSMs
  -> process_search_results!()
       -> add_second_search_columns!()   # RT, charge, target/decoy
       -> get_isotopes_captured!()       # isotope info
       -> filter!(weight > 0)
       -> train_and_apply_prescore!()    # probit prescore -> best_scan
       -> init_summary_columns!()
       -> get_summary_scores!()
       -> filter!(best_scan)             # one PSM per precursor
       -> MS1 join
       -> add_features!()                # iRT, aa comp, pair_id, etc.
       -> write fold arrows
```

## Proposed Flow (Iterative Filter + Re-search)

```
process_file!
  -> perform_second_pass_search()       # initial fragment index search -> ALL PSMs
  -> iterative_prescore_filter!()       # NEW: 2 rounds, each re-runs search
  |    |
  |    |  Round 1:
  |    |    -> add_second_search_columns!()
  |    |    -> get_isotopes_captured!() + basic quality filters
  |    |    -> train_and_apply_prescore!() + summary scores
  |    |    -> add_features!()  (with MS1 stubs)
  |    |    -> reduce to best PSM per precursor (by weight)
  |    |    -> form target-decoy pairs via pair_id
  |    |    -> train diff-cov classifier on 10 features using pairs
  |    |    -> apply diff-cov weights to ALL PSMs -> compute q-values
  |    |    -> train probit on targets(q<=0.05) + all decoys
  |    |    -> apply probit to ALL PSMs -> compute q-values
  |    |    -> collect surviving (precursor_idx, scan_idx) where q <= 0.50
  |    |    -> rebuild filtered scan_to_prec_idx
  |    |    -> RE-RUN perform_second_pass_search() with filtered pairs
  |    |       (re-deconvolves with reduced competition)
  |    |
  |    |  Round 2: identical pipeline on round-1 survivors
  |    |    -> same feature computation + diff-cov + probit
  |    |    -> collect surviving (precursor_idx, scan_idx) where q <= 0.50
  |    |    -> rebuild filtered scan_to_prec_idx
  |    |    -> RE-RUN perform_second_pass_search() with filtered pairs
  |    |       (final re-deconvolution)
  |    |
  |    |  Return: fresh PSMs from final re-search (raw columns only)
  |    |
  -> process_search_results!()          # UNCHANGED: final processing on survivors
       -> (existing flow: prescore, summary, best_scan, MS1, features, arrow write)
```

## Top 10 Features for Diff-Cov + Probit

Selected from LightGBM feature importance (gain):

```julia
const ITERATIVE_PRESCORE_FEATURES = [
    :max_fitted_manhattan_distance,   # 16339
    :max_matched_residual,            # 10980
    :max_gof,                         #  6162
    :gof,                             #  1004
    :max_y_ions,                      #   896
    :poisson,                         #   728
    :irt_error,                       #   683
    :y_ions_sum,                      #   555
    :max_unmatched_residual,          #   546
    :y_count,                         #   500
]
```

**Feature categories:**
- Summary features (`max_gof`, `max_y_ions`, `y_ions_sum`, `max_fitted_manhattan_distance`, `max_matched_residual`, `max_unmatched_residual`): require `get_summary_scores!` computation
- Per-scan features (`gof`, `poisson`, `y_count`): available from fragment index search
- RT features (`irt_error`): require `add_features!`

## Detailed Implementation

### Main Orchestrator: `iterative_prescore_filter!`

**File:** `SecondPassSearch/utils.jl` (new function)

```julia
"""
    iterative_prescore_filter!(psms, search_context, spectra, params, ms_file_idx)

Run two rounds of diff-cov + probit filtering to progressively eliminate
low-quality (precursor_idx, scan_idx) pairs. After each round, re-runs
perform_second_pass_search with the surviving pairs so deconvolution is
re-solved with reduced competition.

Returns a fresh PSM DataFrame from the final re-search (raw columns only,
ready for process_search_results!).
"""
function iterative_prescore_filter!(
    psms::DataFrame,
    search_context::SearchContext,
    spectra::MassSpecData,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64
)
    n_initial = nrow(psms)
    @info "=== Iterative Prescore Filter: starting with $n_initial PSMs ===\n"

    for round in 1:2
        n_round_start = nrow(psms)
        @info "--- Filter Round $round: $n_round_start PSMs entering ---\n"

        # Step 1: Add basic columns (RT, charge, target/decoy)
        add_second_search_columns!(psms,
            getRetentionTimes(spectra),
            getCharge(getPrecursors(getSpecLib(search_context))),
            getIsDecoy(getPrecursors(getSpecLib(search_context))),
            getPrecursors(getSpecLib(search_context))
        )

        # Step 2: Isotope capture info
        get_isotopes_captured!(psms,
            getIsotopeTraceType(params),
            getQuadTransmissionModel(search_context, ms_file_idx),
            getSearchData(search_context),
            psms[!, :scan_idx],
            getCharge(getPrecursors(getSpecLib(search_context))),
            getMz(getPrecursors(getSpecLib(search_context))),
            getSulfurCount(getPrecursors(getSpecLib(search_context))),
            getCenterMzs(spectra),
            getIsolationWidthMzs(spectra)
        )

        # Step 3: Basic quality filters
        filter!(row -> row.precursor_fraction_transmitted >= params.min_fraction_transmitted, psms)
        n_pre_weight = nrow(psms)
        filter!(row -> row.weight > 0.0f0, psms)
        n_post_filter = nrow(psms)
        @info "  Round $round: $n_post_filter PSMs after quality filters (removed $(n_round_start - n_post_filter))\n"

        # Step 4: Prescore for best_scan selection (used by summary scores)
        psms[!, :best_scan] = zeros(Bool, nrow(psms))
        train_and_apply_prescore!(psms)

        # Step 5: Summary scores (max_gof, max_y_ions, etc.)
        init_summary_columns!(psms)
        for (key, gpsms) in pairs(groupby(psms, getPsmGroupbyCols(getIsotopeTraceType(params))))
            get_summary_scores!(
                gpsms,
                gpsms[!, :weight],
                gpsms[!, :gof],
                gpsms[!, :matched_ratio],
                gpsms[!, :fitted_manhattan_distance],
                gpsms[!, :fitted_spectral_contrast],
                gpsms[!, :scribe],
                gpsms[!, :y_count],
                getRtIrtModel(search_context, ms_file_idx)
            )
        end

        # Step 6: ML features (irt_error, pair_id, etc.)
        # Need MS1 stubs since we haven't done the MS1 join yet
        ensure_ms1_stub_columns!(psms)
        add_features!(psms, search_context,
            getTICs(spectra), getMzArrays(spectra),
            ms_file_idx,
            getRtIrtModel(search_context, ms_file_idx),
            getPrecursorDict(search_context)
        )

        # Step 7: Reduce to best PSM per precursor (by weight) for pairing
        best_psms = get_best_psm_per_precursor(psms)
        n_unique_prec = nrow(best_psms)
        n_targets = count(best_psms[!, :target])
        n_decoys = n_unique_prec - n_targets
        @info "  Round $round: $n_unique_prec unique precursors ($n_targets T, $n_decoys D)\n"

        # Step 8: Form target-decoy pairs
        pairs_td = form_target_decoy_pairs(best_psms)
        n_pairs = length(pairs_td)
        @info "  Round $round: $n_pairs target-decoy pairs\n"

        if n_pairs < 20
            @info "  Round $round: Too few pairs ($n_pairs < 20), skipping filter\n"
            break
        end

        # Step 9: Train difference-covariance classifier on paired best PSMs
        features = filter(f -> hasproperty(best_psms, f), collect(ITERATIVE_PRESCORE_FEATURES))
        w = train_difference_covariance!(best_psms, pairs_td, features)
        if w === nothing
            @info "  Round $round: Diff-cov training failed, skipping\n"
            break
        end

        # Step 10: Apply diff-cov weights to ALL PSMs (not just best)
        apply_linear_scores!(psms, w, features)

        # Step 11: Compute q-values from diff-cov scores
        q_vals = zeros(Float64, nrow(psms))
        get_qvalues!(psms[!, :score], psms[!, :target], q_vals)
        psms[!, :q_value] = q_vals

        n_1pct_dc = count((q_vals .<= 0.01) .& psms[!, :target])
        n_5pct_dc = count((q_vals .<= 0.05) .& psms[!, :target])
        @info "  Round $round diff-cov: $n_1pct_dc targets @ 1% FDR, $n_5pct_dc @ 5% FDR\n"

        # Step 12: Select probit training data: targets at q <= 0.05 + all decoys
        targets_col = psms[!, :target]
        train_mask = BitVector(((q_vals .<= 0.05) .& targets_col) .| .!targets_col)
        n_train_t = count((q_vals .<= 0.05) .& targets_col)
        n_train_d = count(.!targets_col)
        @info "  Round $round probit training: $n_train_t targets + $n_train_d decoys = $(n_train_t + n_train_d) total\n"

        # Step 13: Train probit on same features
        beta = train_probit_on_features!(psms, train_mask, features)

        # Step 14: Apply probit to ALL PSMs -> compute q-values
        apply_probit_scores!(psms, beta, features)
        q_vals_probit = zeros(Float64, nrow(psms))
        get_qvalues!(psms[!, :probit_score], psms[!, :target], q_vals_probit)
        psms[!, :q_value] = q_vals_probit

        n_1pct_p = count((q_vals_probit .<= 0.01) .& psms[!, :target])
        n_5pct_p = count((q_vals_probit .<= 0.05) .& psms[!, :target])
        @info "  Round $round probit: $n_1pct_p targets @ 1% FDR, $n_5pct_p @ 5% FDR\n"

        # Step 15: Collect surviving (precursor_idx, scan_idx) where q <= 0.50
        keep_mask = q_vals_probit .<= 0.50
        n_keep = count(keep_mask)
        n_discard = nrow(psms) - n_keep
        @info "  Round $round: keeping $n_keep PSMs, discarding $n_discard (q > 0.50)\n"

        surviving_pairs = Set{Tuple{UInt32, UInt32}}()
        prec_col = psms[!, :precursor_idx]
        scan_col = psms[!, :scan_idx]
        for i in 1:nrow(psms)
            if keep_mask[i]
                push!(surviving_pairs, (prec_col[i], scan_col[i]))
            end
        end
        n_surviving_precs = length(unique(p[1] for p in surviving_pairs))
        @info "  Round $round: $(length(surviving_pairs)) surviving (prec, scan) pairs across $n_surviving_precs precursors\n"

        # Step 16: Re-run fragment index search with filtered pairs
        # This re-solves the deconvolution with reduced competition
        @info "  Round $round: re-running search with filtered pairs...\n"
        psms = rerun_search_with_filter(
            spectra, search_context, params, ms_file_idx, surviving_pairs
        )
        @info "  Round $round: re-search produced $(nrow(psms)) PSMs\n"
    end

    n_final = nrow(psms)
    pct_reduction = round(100.0 * (1 - n_final / n_initial), digits=1)
    @info "=== Filter complete: $n_initial -> $n_final PSMs ($pct_reduction% reduction) ===\n"

    # Return fresh PSMs from final re-search (raw columns only).
    # process_search_results! will add all derived columns from scratch.
    return psms
end
```

### Helper Function 1: `get_best_psm_per_precursor`

```julia
"""
    get_best_psm_per_precursor(psms::DataFrame) -> DataFrame

Reduce to one PSM per precursor_idx, keeping the row with highest `weight`.
Returns a copy (does not modify input).
"""
function get_best_psm_per_precursor(psms::DataFrame)
    sorted = sort(psms, [:precursor_idx, order(:weight, rev=true)])
    best = combine(groupby(sorted, :precursor_idx), first)
    return best
end
```

### Helper Function 2: `form_target_decoy_pairs`

Uses the `pair_id` column (partner's precursor_idx from the spectral library).

```julia
"""
    form_target_decoy_pairs(best_psms::DataFrame) -> Vector{Tuple{Int,Int}}

Find matched target-decoy pairs. Each pair is (target_row_idx, decoy_row_idx).
Uses the `pair_id` column which holds the partner's precursor_idx.
Only pairs where both the target and decoy precursor have a best PSM are included.
"""
function form_target_decoy_pairs(best_psms::DataFrame)
    # Build lookup: precursor_idx -> row index in best_psms
    prec_to_row = Dict{UInt32, Int}()
    for i in 1:nrow(best_psms)
        prec_to_row[best_psms[i, :precursor_idx]] = i
    end

    pairs = Vector{Tuple{Int,Int}}()
    visited = falses(nrow(best_psms))

    for i in 1:nrow(best_psms)
        visited[i] && continue
        partner_pid = best_psms[i, :pair_id]
        partner_pid == zero(UInt32) && continue

        partner_row = get(prec_to_row, partner_pid, 0)
        partner_row == 0 && continue
        visited[partner_row] && continue

        is_target_i = best_psms[i, :target]
        is_target_j = best_psms[partner_row, :target]
        # Must be one target and one decoy
        (is_target_i == is_target_j) && continue

        t_row = is_target_i ? i : partner_row
        d_row = is_target_i ? partner_row : i

        push!(pairs, (t_row, d_row))
        visited[i] = true
        visited[partner_row] = true
    end

    return pairs
end
```

### Helper Function 3: `train_difference_covariance!`

Direct port of DIA-NN's approach (`diann.cpp:9281-9422`).

```julia
"""
    train_difference_covariance!(best_psms, pairs, features) -> Union{Vector{Float64}, Nothing}

Train the DIA-NN-style difference-covariance classifier.

Math:
  delta_k = features_target[k] - features_decoy[k]   (per-pair difference)
  mean_delta = mean(delta_k)
  A = cov(delta_k)                                     (covariance of differences)
  w = (A + eps*I) \\ mean_delta                        (regularized solve)

Score for any PSM i: score_i = w' * features_i
"""
function train_difference_covariance!(
    best_psms::DataFrame,
    pairs::Vector{Tuple{Int,Int}},
    features::Vector{Symbol};
    regularization::Float64 = 1e-9
)
    N = length(pairs)
    p = length(features)

    if N < 20
        @info "    Diff-cov: too few pairs ($N), need >= 20\n"
        return nothing
    end

    # Mean feature differences (target - decoy)
    ds_mean = zeros(Float64, p)
    for (t_row, d_row) in pairs
        for j in 1:p
            ds_mean[j] += Float64(best_psms[t_row, features[j]]) -
                          Float64(best_psms[d_row, features[j]])
        end
    end
    ds_mean ./= N

    # Covariance of differences (upper triangle, then symmetrize)
    A = zeros(Float64, p, p)
    for (t_row, d_row) in pairs
        for i in 1:p
            di = (Float64(best_psms[t_row, features[i]]) -
                  Float64(best_psms[d_row, features[i]])) - ds_mean[i]
            for j in i:p
                dj = (Float64(best_psms[t_row, features[j]]) -
                      Float64(best_psms[d_row, features[j]])) - ds_mean[j]
                A[i,j] += di * dj
            end
        end
    end
    for i in 1:p, j in i:p
        A[i,j] /= (N - 1)
        A[j,i] = A[i,j]
    end
    # Tikhonov regularization
    for i in 1:p
        A[i,i] += regularization
    end

    # Solve for weights: w = (A + eps*I)^-1 * mean_delta
    w = A \ ds_mean

    # Diagnostic: print weights
    @info "    Diff-cov weights ($N pairs, $p features):\n"
    for (fname, wval) in zip(features, w)
        @info "      $(rpad(fname, 35)) $(round(wval, digits=6))\n"
    end

    return w
end
```

### Helper Function 4: `apply_linear_scores!`

```julia
"""
    apply_linear_scores!(psms, w, features)

Apply linear classifier weights to ALL PSMs: score_i = w' * features_i.
Stores result in psms[!, :score].
"""
function apply_linear_scores!(psms::DataFrame, w::Vector{Float64}, features::Vector{Symbol})
    n = nrow(psms)
    scores = zeros(Float32, n)
    for i in 1:n
        s = 0.0
        for (j, f) in enumerate(features)
            s += w[j] * Float64(psms[i, f])
        end
        scores[i] = Float32(s)
    end
    psms[!, :score] = scores
end
```

### Helper Function 5: `train_probit_on_features!`

```julia
"""
    train_probit_on_features!(psms, train_mask, features) -> Vector{Float64}

Train a probit model on the masked subset using the given features.
"""
function train_probit_on_features!(
    psms::DataFrame,
    train_mask::BitVector,
    features::Vector{Symbol};
    max_iter::Int = 20
)
    X_train = DataFrame(Matrix{Float64}(psms[train_mask, features]), features)
    targets_train = Vector{Bool}(psms[train_mask, :target])
    n_train = sum(train_mask)

    chunk_size = max(1, n_train ÷ (10 * Threads.nthreads()))
    data_chunks = Iterators.partition(1:n_train, chunk_size)

    beta = zeros(Float64, length(features))
    beta = ProbitRegression(beta, X_train, targets_train, data_chunks; max_iter=max_iter)

    @info "    Probit coefficients:\n"
    for (fname, coef) in zip(features, beta)
        @info "      $(rpad(fname, 35)) $(round(coef, digits=6))\n"
    end

    return beta
end
```

### Helper Function 6: `apply_probit_scores!`

```julia
"""
    apply_probit_scores!(psms, beta, features)

Apply probit model to ALL PSMs. Stores probability scores in :probit_score.
"""
function apply_probit_scores!(
    psms::DataFrame,
    beta::Vector{Float64},
    features::Vector{Symbol}
)
    n = nrow(psms)
    X_all = DataFrame(Matrix{Float64}(psms[!, features]), features)
    scores = zeros(Float32, n)

    chunk_size = max(1, n ÷ (10 * Threads.nthreads()))
    data_chunks = Iterators.partition(1:n, chunk_size)

    ModelPredictProbs!(scores, X_all, beta, data_chunks)
    psms[!, :probit_score] = scores
end
```

### Helper Function 7: `rerun_search_with_filter`

This is the critical function. It rebuilds the `scan_to_prec_idx` / `precursors_passed` data structures with only the surviving (precursor_idx, scan_idx) combinations, then re-runs `perform_second_pass_search`. This causes the Huber-regularized deconvolution to be re-solved with fewer precursors competing per scan.

```julia
"""
    rerun_search_with_filter(spectra, search_context, params, ms_file_idx, surviving_pairs)

Re-run the fragment index search with only the surviving (precursor_idx, scan_idx)
combinations. Rebuilds the scan_to_prec_idx mapping and re-runs
perform_second_pass_search, which re-solves the deconvolution with reduced
precursor competition per scan.
"""
function rerun_search_with_filter(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64,
    surviving_pairs::Set{Tuple{UInt32, UInt32}}
)
    # Reload original fragment index matches
    frag_match_path = getFragmentIndexMatches(getMSData(search_context), ms_file_idx)
    scan_to_prec_idx_orig, precursors_passed_orig = load_fragment_index_matches(
        frag_match_path, length(spectra)
    )

    n_original = length(precursors_passed_orig)

    # Rebuild with only surviving entries
    # scan_to_prec_idx[scan] = start:end range into precursors_passed
    new_precursors = UInt32[]
    new_scan_to_prec = Vector{Union{Missing, UnitRange{Int64}}}(missing, length(spectra))

    for scan_idx in 1:length(spectra)
        range = scan_to_prec_idx_orig[scan_idx]
        ismissing(range) && continue

        start_new = length(new_precursors) + 1
        for idx in range
            pid = precursors_passed_orig[idx]
            if (pid, UInt32(scan_idx)) in surviving_pairs
                push!(new_precursors, pid)
            end
        end
        end_new = length(new_precursors)

        if end_new >= start_new
            new_scan_to_prec[scan_idx] = start_new:end_new
        end
    end

    n_filtered = length(new_precursors)
    pct = round(100.0 * n_filtered / max(1, n_original), digits=1)
    @info "    Re-search input: $n_filtered / $n_original fragment index entries ($pct%)\n"

    # Re-run the actual search with filtered input
    psms = perform_second_pass_search(
        spectra, new_scan_to_prec, new_precursors,
        search_context, params, ms_file_idx, MS2CHROM()
    )

    return psms
end
```

### Helper Function 8: `ensure_ms1_stub_columns!`

During iterative filter rounds, `add_features!` expects `rt_ms1` and `ms1_features_missing` to exist (used for `ms1_irt_diff` computation). We add placeholder values.

```julia
"""
    ensure_ms1_stub_columns!(psms::DataFrame)

Add placeholder MS1 columns needed by add_features!() during iterative
prescore filtering, before the real MS1 join in process_search_results!.
"""
function ensure_ms1_stub_columns!(psms::DataFrame)
    n = nrow(psms)
    if !hasproperty(psms, :rt_ms1)
        psms[!, :rt_ms1] = fill(Float32(-1), n)
    end
    if !hasproperty(psms, :ms1_features_missing)
        psms[!, :ms1_features_missing] = trues(n)
    end
    if !hasproperty(psms, :ms1_ms2_rt_diff)
        psms[!, :ms1_ms2_rt_diff] = fill(Float32(-1), n)
    end
end
```

### Integration Point: `SecondPassSearch.jl` process_file!

```julia
function process_file!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:SecondPassSearchParameters}

    if check_and_skip_failed_file(search_context, ms_file_idx, "SecondPassSearch")
        return results
    end

    try
        frag_match_path = getFragmentIndexMatches(getMSData(search_context), ms_file_idx)
        scan_to_prec_idx, precursors_passed = load_fragment_index_matches(
            frag_match_path, length(spectra)
        )

        psms = perform_second_pass_search(
            spectra, scan_to_prec_idx, precursors_passed,
            search_context, params, ms_file_idx, MS2CHROM()
        )

        # NEW: iterative prescore filter (2 rounds of diff-cov + probit)
        # Returns fresh PSMs from final re-search (raw columns only)
        psms = iterative_prescore_filter!(
            psms, search_context, spectra, params, ms_file_idx
        )

        ms1_psms = DataFrame()
        results.psms[] = psms
        results.ms1_psms[] = ms1_psms
    catch e
        handle_search_error!(search_context, ms_file_idx, "SecondPassSearch", e,
                            createFallbackResults!, results)
    end

    return results
end
```

**Key point:** `iterative_prescore_filter!` returns a fresh DataFrame from the final `perform_second_pass_search` call. This DataFrame has only the raw `ComplexScoredPSM` columns — no derived columns like `target`, `irt_error`, `pair_id`, etc. So `process_search_results!` runs its normal flow without any column conflicts.

## Summary of All Changes

### New Functions (all in `SecondPassSearch/utils.jl`):
| Function | Purpose |
|----------|---------|
| `iterative_prescore_filter!` | Main orchestrator: 2 rounds of filter + re-search |
| `get_best_psm_per_precursor` | Reduce to 1 PSM per precursor (max weight) |
| `form_target_decoy_pairs` | Pair targets with decoys via `pair_id` |
| `train_difference_covariance!` | DIA-NN-style diff-cov classifier (solve A\\delta) |
| `apply_linear_scores!` | Apply linear weights: score = w'*features |
| `train_probit_on_features!` | Train probit on masked subset |
| `apply_probit_scores!` | Apply probit probabilities to all PSMs |
| `rerun_search_with_filter` | Rebuild scan_to_prec_idx and re-run search |
| `ensure_ms1_stub_columns!` | Add placeholder MS1 columns for add_features! |

### New Constant:
| Name | Value |
|------|-------|
| `ITERATIVE_PRESCORE_FEATURES` | 10 features selected by LightGBM importance |

### Modified Functions:
| Function | File | Change |
|----------|------|--------|
| `process_file!` | `SecondPassSearch.jl` | Insert `iterative_prescore_filter!` call between search and `process_search_results!` |

### NOT Modified:
- `process_search_results!` — runs unchanged on the filtered PSMs
- `ScoringSearch` — reads the same Arrow files as before

## Expected Diagnostic Output

For each MS file, you should see:
```
=== Iterative Prescore Filter: starting with 450000 PSMs ===
--- Filter Round 1: 450000 PSMs entering ---
  Round 1: 380000 PSMs after quality filters (removed 70000)
  Training probit prescore model on 380000 scans...
  Spectral pair coverage: 72000 / 380000 PSMs (18.9%) have paired target/decoy in file
  Round 1: 95000 unique precursors (55000 T, 40000 D)
  Round 1: 18000 target-decoy pairs
    Diff-cov weights (18000 pairs, 10 features):
      max_fitted_manhattan_distance       0.002341
      max_matched_residual                -0.001892
      max_gof                             0.004521
      gof                                 0.001205
      max_y_ions                          0.000891
      poisson                             -0.003102
      irt_error                           -0.005634
      y_ions_sum                          0.000445
      max_unmatched_residual              -0.001123
      y_count                             0.000678
  Round 1 diff-cov: 42000 targets @ 1% FDR, 48000 @ 5% FDR
  Round 1 probit training: 48000 targets + 190000 decoys = 238000 total
    Probit coefficients:
      max_fitted_manhattan_distance       0.1234
      max_matched_residual                -0.0891
      ...
  Round 1 probit: 44000 targets @ 1% FDR, 50000 @ 5% FDR
  Round 1: keeping 310000 PSMs, discarding 70000 (q > 0.50)
  Round 1: 310000 surviving (prec, scan) pairs across 78000 precursors
  Round 1: re-running search with filtered pairs...
    Re-search input: 310000 / 500000 fragment index entries (62.0%)
  Round 1: re-search produced 305000 PSMs
--- Filter Round 2: 305000 PSMs entering ---
  Round 2: 290000 PSMs after quality filters (removed 15000)
  ...
  Round 2: keeping 270000 PSMs, discarding 20000 (q > 0.50)
  Round 2: 270000 surviving (prec, scan) pairs across 72000 precursors
  Round 2: re-running search with filtered pairs...
    Re-search input: 270000 / 500000 fragment index entries (54.0%)
  Round 2: re-search produced 268000 PSMs
=== Filter complete: 450000 -> 268000 PSMs (40.4% reduction) ===
```

Then `process_search_results!` runs its normal flow on those 268000 fresh PSMs.

## Performance Considerations

Each round involves:
1. `add_second_search_columns!` + `get_isotopes_captured!` — fast (column ops)
2. `train_and_apply_prescore!` — moderate (probit IRLS)
3. `init_summary_columns!` + `get_summary_scores!` — fast (group-by + aggregation)
4. `add_features!` — moderate (threaded RT lookups)
5. `get_best_psm_per_precursor` — fast (sort + groupby)
6. `form_target_decoy_pairs` — fast (Dict lookup)
7. `train_difference_covariance!` — fast (small matrix solve, p=10)
8. `apply_linear_scores!` + `apply_probit_scores!` — fast (dot products)
9. `rerun_search_with_filter` — **expensive** (full deconvolution re-run)

Total: roughly **3x the current SecondPassSearch time per file** (2 extra deconvolution passes). The final `process_search_results!` is faster because it processes fewer PSMs.

## Open Questions

1. **Q-value threshold for probit training**: Currently 0.05. DIA-NN uses 0.01 for the diff-cov filter. Start with 0.05 for more training data; tune later.

2. **Should the diff-cov iterate within each round?**: The LaTeX doc describes up to 4 iterations with feature gating. Current plan does 1 diff-cov + 1 probit per round. Can add multi-iteration diff-cov later.

3. **Q-value discard threshold**: Currently 0.50. Could be more aggressive (e.g., 0.25) or more conservative (0.75). Start at 0.50 and tune based on results.
