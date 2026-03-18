#!/usr/bin/env julia
#
# Standalone LightGBM rescoring script for first-pass PSM arrow files.
# Loads arrow files, trains 2-fold CV LightGBM with iterative negative mining,
# computes q-values, and reports per-file ID counts at 1% FDR.
#
# Usage:
#   julia --project=. scripts/rescore_first_pass.jl <path_to_first_pass_psms_folder>
#
# Example:
#   julia --project=. scripts/rescore_first_pass.jl \
#       /Users/nathanwamsley/Data/Results/fp_rank10_iso1_17170a71/temp_data/first_pass_psms

using Arrow, DataFrames, Random, Printf, CSV
using LightGBM

###############################################################################
# Configuration
###############################################################################

const CANDIDATE_FEATURES = [
    :fitted_spectral_contrast, :gof, :scribe, :percent_theoretical_ignored,
    :fitted_manhattan_distance, :max_matched_residual, :max_unmatched_residual,
    :matched_ratio, :weight, :charge2, :poisson, :irt_error,
    :missed_cleavage, :Mox, :TIC, :y_count, :err_norm,
    :spectrum_peak_count, :prec_mz,
    # Multi-scan max features (present in newer arrow files)
    :max_gof, :max_fitted_manhattan_distance, :max_fitted_spectral_contrast,
    :max_scribe, :max_matched_ratio, :max_y_count,
]

const ITER_ROUNDS = [100, 200, 200]
const MAX_Q_VALUE = 0.01f0
const MIN_PEP_THRESHOLD = 0.90f0

# LightGBM hyperparameters (matching AdvancedLightGBM / LightGBMScorer defaults)
const LGBM_PARAMS = (
    max_depth       = 10,
    learning_rate   = 0.05,
    num_leaves      = 63,
    feature_fraction = 0.5,
    bagging_fraction = 0.25,
    min_data_in_leaf = 500,
    min_gain_to_split = 0.5,
)

###############################################################################
# Core helpers (inlined from Pioneer to keep script self-contained)
###############################################################################

"""Build feature matrix from DataFrame columns."""
function feature_matrix(df::AbstractDataFrame, features::Vector{Symbol})
    n = nrow(df)
    mat = Matrix{Float32}(undef, n, length(features))
    for (j, feat) in enumerate(features)
        col = df[!, feat]
        T = nonmissingtype(eltype(col))
        if eltype(col) <: Union{Missing, T}
            mat[:, j] = Float32.(coalesce.(col, zero(T)))
        else
            mat[:, j] = Float32.(col)
        end
    end
    return mat
end

"""Build a LightGBM classifier with the configured hyperparameters."""
function build_classifier(; num_iterations::Int)
    return LightGBM.LGBMClassification(
        objective        = "binary",
        metric           = ["binary_logloss"],
        learning_rate    = LGBM_PARAMS.learning_rate,
        num_iterations   = num_iterations,
        num_leaves       = LGBM_PARAMS.num_leaves,
        max_depth        = LGBM_PARAMS.max_depth,
        feature_fraction = LGBM_PARAMS.feature_fraction,
        bagging_fraction = LGBM_PARAMS.bagging_fraction,
        bagging_freq     = 1,
        min_data_in_leaf = LGBM_PARAMS.min_data_in_leaf,
        min_gain_to_split = LGBM_PARAMS.min_gain_to_split,
        num_class        = 1,
        num_threads      = Threads.nthreads(),
        verbosity        = -1,
        is_unbalance     = false,
        seed             = 1776,
        deterministic    = true,
        force_row_wise   = true,
    )
end

"""Train LightGBM and return the fitted booster (or nothing if degenerate)."""
function train_lgbm(X::Matrix{Float32}, labels::Vector{Int}; num_iterations::Int)
    if length(unique(labels)) == 1
        return nothing
    end
    model = build_classifier(; num_iterations)
    LightGBM.fit!(model, X, labels; verbosity = -1)
    return model
end

"""Predict probabilities from a fitted booster."""
function predict_lgbm(model, X::Matrix{Float32})
    if model === nothing
        return fill(0.5f0, size(X, 1))
    end
    raw = LightGBM.predict(model, X)
    return vec(Float32.(raw))
end

"""In-place q-value calculation (target-decoy approach)."""
function get_qvalues!(probs::AbstractVector{<:AbstractFloat},
                      labels::AbstractVector{Bool},
                      qvals::AbstractVector{<:AbstractFloat})
    order = sortperm(probs, rev=true, alg=QuickSort)
    targets = 0
    decoys = 0
    @inbounds for i in order
        targets += labels[i]
        decoys += (1 - labels[i])
        qvals[i] = decoys / max(targets, 1)
    end
    fdr = Inf
    @inbounds for i in reverse(order)
        if qvals[i] > fdr
            qvals[i] = fdr
        else
            fdr = qvals[i]
        end
    end
end

"""Weighted PAVA for isotonic regression (used by get_PEP!)."""
function _weighted_pava(y::Vector{Float64}, w::Vector{Float64})
    n = length(y)
    v    = Vector{Float64}(undef, n)
    wt   = Vector{Float64}(undef, n)
    lenv = Vector{Int}(undef, n)
    m = 0
    for i in 1:n
        m += 1
        @inbounds begin
            v[m]    = y[i]
            wt[m]   = w[i]
            lenv[m] = 1
        end
        while m > 1 && v[m-1] > v[m]
            @inbounds begin
                new_w   = wt[m-1] + wt[m]
                new_v   = (v[m-1]*wt[m-1] + v[m]*wt[m]) / new_w
                new_len = lenv[m-1] + lenv[m]
                m -= 1
                v[m]    = new_v
                wt[m]   = new_w
                lenv[m] = new_len
            end
        end
    end
    result = Vector{Float64}(undef, n)
    idx = 1
    for j in 1:m
        @inbounds for _ in 1:lenv[j]
            result[idx] = v[j]
            idx += 1
        end
    end
    return result
end

"""In-place PEP estimation via isotonic regression."""
function get_PEP!(scores::AbstractVector{<:AbstractFloat},
                  is_target::AbstractVector{Bool},
                  peps::AbstractVector{<:AbstractFloat};
                  doSort::Bool=true)
    N = length(scores)
    N == 0 && return
    order = doSort ? sortperm(scores, rev=true, alg=QuickSort) : collect(eachindex(scores))
    labels = Vector{Float64}(undef, N + 1)
    weights = Vector{Float64}(undef, N + 1)
    labels[1] = 0.0
    weights[1] = 1.0
    @inbounds for j in 1:N
        idx = order[j]
        labels[j+1]  = is_target[idx] ? 0.0 : 1.0
        weights[j+1] = 1.0
    end
    fitted = _weighted_pava(labels, weights)[2:end]
    pep = clamp.(fitted ./ (1.0 .- fitted), 0.0, 1.0)
    @inbounds for (j, idx) in enumerate(order)
        peps[idx] = pep[j]
    end
end

###############################################################################
# Data loading
###############################################################################

function load_psms(folder::String)
    arrow_files = filter(f -> endswith(f, ".arrow"), readdir(folder))
    isempty(arrow_files) && error("No .arrow files found in $folder")

    dfs = DataFrame[]
    for f in arrow_files
        path = joinpath(folder, f)
        df = DataFrame(Arrow.Table(path))
        println("  Loaded $f: $(nrow(df)) PSMs")
        push!(dfs, df)
    end
    psms = vcat(dfs...)
    println("Total: $(nrow(psms)) PSMs from $(length(arrow_files)) files")

    # Validate required columns
    required = [:target; CANDIDATE_FEATURES]
    missing_cols = setdiff(required, Symbol.(names(psms)))
    if !isempty(missing_cols)
        error("Missing required columns: $missing_cols")
    end
    return psms
end

###############################################################################
# Iterative 2-fold CV rescoring
###############################################################################

"""Get gain-based feature importances from a LightGBM booster."""
function get_importances(model, features::Vector{Symbol})
    model === nothing && return nothing
    try
        gains = LightGBM.gain_importance(model)
        return collect(zip(features, gains))
    catch
        return nothing
    end
end

function rescore!(psms::DataFrame, features::Vector{Symbol})
    N = nrow(psms)
    targets = psms.target
    n_targets = sum(targets)
    n_decoys = N - n_targets
    println("\nDataset: $n_targets targets, $n_decoys decoys")
    println("Features ($(length(features))): $features")
    println("Iteration rounds: $ITER_ROUNDS\n")

    # Assign CV folds (seeded)
    Random.seed!(1776)
    folds = rand(UInt8[0, 1], N)

    # Output arrays
    probs = fill(0.5f0, N)
    qvals = zeros(Float64, N)

    # Collect all trained models for feature importances
    all_models = []

    for test_fold in UInt8[0, 1]
        train_mask = folds .!= test_fold
        test_mask  = folds .== test_fold
        train_idx  = findall(train_mask)
        test_idx   = findall(test_mask)

        println("Fold $(test_fold): train=$(length(train_idx)), test=$(length(test_idx))")

        # Build feature matrices (test is fixed)
        X_test = feature_matrix(psms[test_idx, :], features)

        for (itr, num_rounds) in enumerate(ITER_ROUNDS)
            # Select training data
            if itr == 1
                sel_idx = train_idx
            else
                # Negative mining: compute PEP on training fold predictions,
                # flip worst targets to decoys, then keep q_value ≤ 0.01 targets + all decoys
                train_probs = probs[train_idx]
                train_targets = Vector{Bool}(targets[train_idx])

                # Compute PEP on training fold
                pep_vals = Vector{Float64}(undef, length(train_idx))
                get_PEP!(Float64.(train_probs), train_targets, pep_vals)

                # Flip worst targets to decoys
                effective_targets = copy(train_targets)
                pep_order = sortperm(Float64.(train_probs), rev=true)
                sorted_targets = train_targets[pep_order]
                sorted_pep = Vector{Float64}(undef, length(pep_order))
                get_PEP!(Float64.(train_probs[pep_order]), sorted_targets, sorted_pep; doSort=false)
                cutoff_pos = findfirst(x -> x >= MIN_PEP_THRESHOLD, sorted_pep)
                if cutoff_pos !== nothing
                    worst = pep_order[cutoff_pos:end]
                    effective_targets[worst] .= false
                end

                # Compute q-values with effective targets
                train_qvals = Vector{Float64}(undef, length(train_idx))
                get_qvalues!(Float64.(train_probs), effective_targets, train_qvals)

                # Keep: all decoys + targets with q ≤ threshold
                keep = BitVector([(!t) || (t && q <= MAX_Q_VALUE)
                                  for (t, q) in zip(effective_targets, train_qvals)])
                sel_idx = train_idx[keep]
            end

            # Build training matrix
            X_train = feature_matrix(psms[sel_idx, :], features)
            y_train = Int.(psms.target[sel_idx])

            n_pos = sum(y_train)
            n_neg = length(y_train) - n_pos
            println("  Iter $itr ($num_rounds rounds): train=$(length(sel_idx)) (pos=$n_pos, neg=$n_neg)")

            # Train and predict
            model = train_lgbm(X_train, y_train; num_iterations=num_rounds)
            push!(all_models, (fold=test_fold, iter=itr, model=model))
            probs[test_idx] .= predict_lgbm(model, X_test)

            # Also predict on training fold (needed for next iteration's negative mining)
            X_train_full = feature_matrix(psms[train_idx, :], features)
            probs[train_idx] .= predict_lgbm(model, X_train_full)
        end
    end

    # Final q-values and PEP
    get_qvalues!(Float64.(probs), targets, qvals)
    pep_vals = zeros(Float64, N)
    get_PEP!(Float64.(probs), targets, pep_vals)

    psms[!, :rescore_prob] = probs
    psms[!, :rescore_qvalue] = qvals
    psms[!, :rescore_PEP] = Float32.(pep_vals)

    # Print feature importances (averaged across all final-iteration models)
    print_feature_importances(all_models, features)

    return psms
end

function print_feature_importances(all_models, features::Vector{Symbol})
    # Average gain importances across the final-iteration models from each fold
    n_iters = maximum(m.iter for m in all_models)
    final_models = [m for m in all_models if m.iter == n_iters]

    gain_sums = zeros(Float64, length(features))
    n_models = 0
    for m in final_models
        imp = get_importances(m.model, features)
        imp === nothing && continue
        for (i, (_, g)) in enumerate(imp)
            gain_sums[i] += g
        end
        n_models += 1
    end
    n_models == 0 && return

    avg_gains = gain_sums ./ n_models
    order = sortperm(avg_gains, rev=true)

    println("\nFeature importances (avg gain, final iteration, $n_models models):")
    println("-"^50)
    for i in order
        @printf("  %-35s %12.1f\n", features[i], avg_gains[i])
    end
end

###############################################################################
# Reporting & output
###############################################################################

const Q_THRESHOLDS = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10]

"""Build a table of target counts per file at multiple q-value thresholds."""
function build_results_table(psms::DataFrame;
                             prob_col::Symbol=:rescore_prob,
                             label::String="rescore")
    targets = psms.target
    N = nrow(psms)

    # Compute q-values from the given probability column
    qvals = zeros(Float64, N)
    get_qvalues!(Float64.(psms[!, prob_col]), targets, qvals)

    file_ids = sort(unique(psms.ms_file_idx))

    # Build per-file, per-threshold counts
    rows = DataFrame[]
    for fid in file_ids
        mask = psms.ms_file_idx .== fid
        row = Dict{String, Any}("ms_file_idx" => Int(fid), "source" => label)
        for qt in Q_THRESHOLDS
            col = @sprintf("q<=%.1f%%", qt * 100)
            row[col] = sum(targets[mask] .& (qvals[mask] .<= qt))
        end
        push!(rows, DataFrame(row))
    end

    # Total row
    total = Dict{String, Any}("ms_file_idx" => -1, "source" => label)
    for qt in Q_THRESHOLDS
        col = @sprintf("q<=%.1f%%", qt * 100)
        total[col] = sum(targets .& (qvals .<= qt))
    end
    push!(rows, DataFrame(total))

    result = vcat(rows...)

    # Order columns nicely
    qcols = [@sprintf("q<=%.1f%%", qt * 100) for qt in Q_THRESHOLDS]
    return result[!, vcat(["ms_file_idx", "source"], qcols)]
end

function report_results(psms::DataFrame, outdir::String)
    println("\n" * "="^80)
    println("Target counts at multiple q-value thresholds")
    println("="^80)

    # Rescore results
    rescore_table = build_results_table(psms; prob_col=:rescore_prob, label="rescore")

    # Pipeline comparison if available
    has_pipeline = hasproperty(psms, :lgbm_prob)
    if has_pipeline
        pipeline_table = build_results_table(psms; prob_col=:lgbm_prob, label="pipeline")
        combined = vcat(rescore_table, pipeline_table)
        # Interleave: for each file, show rescore then pipeline
        sort!(combined, [:ms_file_idx, :source], rev=[false, true])
    else
        combined = rescore_table
    end

    # Print the table
    # Header
    qcols = names(combined)[3:end]
    header = @sprintf("%-12s %-10s", "file_idx", "source")
    for qc in qcols
        header *= @sprintf("  %10s", qc)
    end
    println(header)
    println("-"^length(header))

    for row in eachrow(combined)
        fid_str = row.ms_file_idx == -1 ? "TOTAL" : string(Int(row.ms_file_idx))
        line = @sprintf("%-12s %-10s", fid_str, row.source)
        for qc in qcols
            line *= @sprintf("  %10d", row[qc])
        end
        println(line)

        # Print blank line between file groups
        if has_pipeline && row.source == "pipeline" && row.ms_file_idx != -1
            println()
        end
    end

    # Write TSV
    outpath = joinpath(outdir, "rescore_results.tsv")
    # Replace -1 with "TOTAL" for the file
    out_df = copy(combined)
    out_df.ms_file_idx = [x == -1 ? "TOTAL" : string(Int(x)) for x in out_df.ms_file_idx]
    CSV.write(outpath, out_df; delim='\t')
    println("\nResults table written to: $outpath")

    # Also write the full rescored PSM data as arrow
    psm_outpath = joinpath(outdir, "rescored_psms.arrow")
    Arrow.write(psm_outpath, psms)
    println("Rescored PSMs written to: $psm_outpath")
end

###############################################################################
# Main
###############################################################################

function main()
    folder = if length(ARGS) >= 1
        ARGS[1]
    else
        # Default path for interactive use
        "/Users/nathanwamsley/Data/Results/fp_rank10_iso1_17170a71/temp_data/first_pass_psms"
    end

    isdir(folder) || error("Not a directory: $folder")
    println("Loading PSMs from: $folder\n")
    psms = load_psms(folder)

    # Detect available features
    available = Symbol.(names(psms))
    features = filter(f -> f in available, CANDIDATE_FEATURES)
    @assert !isempty(features) "No candidate features found in data"

    rescore!(psms, features)

    # Output to a sibling directory of the input folder
    outdir = joinpath(dirname(folder), "rescore_results")
    mkpath(outdir)
    report_results(psms, outdir)
end

main()
