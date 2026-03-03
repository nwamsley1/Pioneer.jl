#=
FRAGCORR Model Comparison — Do eigengap features improve ML scoring?
=====================================================================
Tests whether adding FRAGCORR features (multi-scan fragment correlation)
to LightGBM and probit models improves target-decoy discrimination
beyond traditional single-scan scoring features.

4 conditions: {LightGBM, Probit} × {baseline, augmented}
+ reference lines: existing pipeline `prob` and unsupervised `eigengap_raw`

Usage:
    include("scripts/fragcorr_model_comparison.jl")
=#

# Reuse all FRAGCORR scoring + q-value utilities
include(joinpath(@__DIR__, "fragcorr_unsupervised_scores.jl"))

using LightGBM
using SpecialFunctions: erf

const MC_DATA_DIR = "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/OlsenEclipse_fragcorr_test/temp_data"
const MC_PSM_DIR  = joinpath(MC_DATA_DIR, "first_pass_psms")
const MC_CORR_DIR = joinpath(MC_DATA_DIR, "first_pass_corr")

# ============================================================
# Step 1: Load & Join Data
# ============================================================

function mc_load_and_join(filename::String)
    psm_path  = joinpath(MC_PSM_DIR, filename)
    corr_path = joinpath(MC_CORR_DIR, filename)

    psm_df = DataFrame(Tables.columntable(Arrow.Table(psm_path)))
    gdf = load_fragcorr_data(corr_path)
    fragcorr_df = score_all_precursors(gdf)

    combined = innerjoin(psm_df, fragcorr_df;
                         on=[:precursor_idx, :target],
                         makeunique=true)
    println("  Joined: $(nrow(psm_df)) PSMs × $(nrow(fragcorr_df)) FRAGCORR → $(nrow(combined)) combined")
    return combined
end

function mc_load_all_files()
    files = filter(f -> endswith(f, ".arrow"), readdir(MC_PSM_DIR))
    corr_files = Set(readdir(MC_CORR_DIR))
    files = filter(f -> f in corr_files, files)

    println("\nLoading $(length(files)) file pairs...")
    dfs = DataFrame[]
    for f in files
        println("\n── $f ──")
        push!(dfs, mc_load_and_join(f))
    end

    combined = vcat(dfs...; cols=:union)
    n_target = sum(combined.target)
    n_decoy  = sum(.!combined.target)
    println("\n════════════════════════════════════════")
    println("Combined: $(nrow(combined)) rows ($n_target targets, $n_decoy decoys)")
    println("════════════════════════════════════════")
    return combined
end

# ============================================================
# Step 2: Feature Sets
# ============================================================

const BASELINE_FEATURES = [:scribe, :spectral_contrast, :city_block, :entropy_score, :y_count, :err_norm]
const FRAGCORR_FEATURES = [:eigengap_raw, :lambda1_frac_raw, :ev1_raw, :median_corr_raw, :best_frag_corr_mean]
const AUGMENTED_FEATURES = vcat(BASELINE_FEATURES, FRAGCORR_FEATURES)

# ============================================================
# Step 3a: Standalone Probit Regression (no @turbo dependency)
# ============================================================

"""
    probit_train(X, y; max_iter=30, tol=1e-2, z_bounds=(-8.0, 8.0)) → β

Simplified probit IRLS (single-threaded, no LoopVectorization).
X: n×p Matrix{Float64}, y: Vector{Bool}.
"""
function probit_train(X::Matrix{Float64}, y::Vector{Bool};
                      max_iter::Int=30, tol::Float64=1e-2,
                      z_bounds::Tuple{Float64,Float64}=(-8.0, 8.0))
    n, p = size(X)
    β = zeros(p)
    η = zeros(n)
    Z = zeros(n)
    W = zeros(n)
    β_old = zeros(p)
    old_loss = -Inf

    for _ in 1:max_iter
        copyto!(β_old, β)

        # η = X * β, clamped
        mul!(η, X, β)
        @inbounds for i in 1:n
            η[i] = clamp(η[i], z_bounds[1], z_bounds[2])
        end

        # Compute Z, W, and log-likelihood
        loss = 0.0
        @inbounds for i in 1:n
            ϕ = exp(-η[i]^2 / 2) / sqrt(2π)
            μ = (1 + erf(η[i] / sqrt(2))) / 2
            μ = clamp(μ, 1e-12, 1 - 1e-12)
            if y[i]
                Z[i] = η[i] + (1 - μ) / ϕ
                loss += log(μ)
            else
                Z[i] = η[i] - μ / ϕ
                loss += log1p(-μ)
            end
            W[i] = ϕ^2 / (μ * (1 - μ))
        end

        # Solve weighted least squares: (X'WX) \ (X'WZ)
        XW = X' * Diagonal(W)
        XWX = XW * X
        XWZ = XW * Z
        β_new = XWX \ XWZ

        # Backtracking line search
        Δβ = β_new .- β
        step = 1.0
        accepted = false
        while step > 1e-4
            β_trial = β .+ step .* Δβ
            η_trial = X * β_trial
            @inbounds for i in 1:n
                η_trial[i] = clamp(η_trial[i], z_bounds[1], z_bounds[2])
            end
            trial_loss = 0.0
            @inbounds for i in 1:n
                μ_t = clamp((1 + erf(η_trial[i] / sqrt(2))) / 2, 1e-12, 1 - 1e-12)
                trial_loss += y[i] ? log(μ_t) : log1p(-μ_t)
            end
            if trial_loss >= loss
                β .= β_trial
                loss = trial_loss
                accepted = true
                break
            end
            step *= 0.5
        end
        if !accepted
            β .+= step .* Δβ
        end

        if norm(β .- β_old) < tol || abs(loss - old_loss) < tol
            break
        end
        old_loss = loss
    end

    return β
end

"""
    probit_predict(X, β) → Vector{Float64}

Probit predicted probabilities: Φ(Xβ).
"""
function probit_predict(X::Matrix{Float64}, β::Vector{Float64})
    η = X * β
    probs = Vector{Float64}(undef, length(η))
    @inbounds for i in eachindex(η)
        probs[i] = (1 + erf(η[i] / sqrt(2))) / 2
    end
    return probs
end

# ============================================================
# Step 3b: Cross-Validated Model Training
# ============================================================

"""
    cv_fold_assignment(precursor_idx) → Vector{Int}

Deterministic 3-fold CV assignment based on precursor_idx modulo.
"""
cv_fold_assignment(precursor_idx::AbstractVector) = Int.(precursor_idx .% 3 .+ 1)

"""
    run_cv_lightgbm(df, features; n_folds=3) → (predictions, fold_models)

3-fold cross-validated LightGBM. Returns held-out predicted probabilities.
"""
function run_cv_lightgbm(df::DataFrame, features::Vector{Symbol}; n_folds::Int=3)
    folds = cv_fold_assignment(df.precursor_idx)
    preds = fill(NaN, nrow(df))
    fold_models = Vector{Any}(undef, n_folds)

    for k in 1:n_folds
        train_mask = folds .!= k
        test_mask  = folds .== k

        X_train = Matrix{Float32}(df[train_mask, features])
        y_train = Int.(df[train_mask, :target])
        X_test  = Matrix{Float32}(df[test_mask, features])

        estimator = LGBMClassification(
            objective = "binary",
            metric = ["binary_logloss"],
            num_iterations = 200,
            learning_rate = 0.1,
            max_depth = 4,
            num_leaves = 15,
            feature_fraction = 0.8,
            min_data_in_leaf = 200,
            min_gain_to_split = 0.5,
            num_threads = Threads.nthreads(),
            verbosity = -1,
            seed = 1776,
            deterministic = true,
            force_row_wise = true,
        )

        fit!(estimator, X_train, y_train; verbosity=-1)
        raw = LightGBM.predict(estimator, X_test)
        ŷ = ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
        preds[test_mask] .= Float64.(ŷ)
        fold_models[k] = estimator
    end

    return preds, fold_models
end

"""
    run_cv_probit(df, features; n_folds=3) → predictions

3-fold cross-validated probit regression. Returns held-out predicted probabilities.
"""
function run_cv_probit(df::DataFrame, features::Vector{Symbol}; n_folds::Int=3)
    folds = cv_fold_assignment(df.precursor_idx)
    preds = fill(NaN, nrow(df))

    for k in 1:n_folds
        train_mask = folds .!= k
        test_mask  = folds .== k

        X_train = Matrix{Float64}(df[train_mask, features])
        y_train = Vector{Bool}(df[train_mask, :target])
        X_test  = Matrix{Float64}(df[test_mask, features])

        β = probit_train(X_train, y_train)
        preds[test_mask] .= probit_predict(X_test, β)
    end

    return preds
end

# ============================================================
# Step 4: Evaluate & Compare
# ============================================================

function run_model_comparison(df::DataFrame;
                              thresholds=[0.001, 0.005, 0.01, 0.02, 0.05, 0.10])
    is_target = df.target
    n_targets = sum(is_target)
    n_decoys  = sum(.!is_target)

    # ── Train all 4 conditions ──
    conditions = Dict{String, Vector{Float64}}()

    println("\n── Training LightGBM (baseline) ──")
    preds_lgbm_base, _ = run_cv_lightgbm(df, BASELINE_FEATURES)
    conditions["LightGBM baseline"] = preds_lgbm_base

    println("\n── Training LightGBM (augmented) ──")
    preds_lgbm_aug, lgbm_aug_models = run_cv_lightgbm(df, AUGMENTED_FEATURES)
    conditions["LightGBM augmented"] = preds_lgbm_aug

    println("\n── Training Probit (baseline) ──")
    preds_probit_base = run_cv_probit(df, BASELINE_FEATURES)
    conditions["Probit baseline"] = preds_probit_base

    println("\n── Training Probit (augmented) ──")
    preds_probit_aug = run_cv_probit(df, AUGMENTED_FEATURES)
    conditions["Probit augmented"] = preds_probit_aug

    # ── Reference scores ──
    conditions["Pipeline prob (ref)"] = Float64.(df.prob)
    conditions["eigengap_raw (unsup)"] = Float64.(df.eigengap_raw)

    # ── Compute q-values for each condition ──
    all_qvals = Dict{String, Vector{Float64}}()
    for (name, scores) in conditions
        # All model outputs and prob are higher-is-better
        hib = true
        # eigengap_raw is also higher-is-better
        all_qvals[name] = compute_qvalues(scores, is_target; higher_is_better=hib)
    end

    # ── Print summary table ──
    thresh_strs = [string(round(t * 100, digits=1)) * "%" for t in thresholds]
    w = 32 + 10 * length(thresholds)
    println("\n\nTargets: $n_targets  |  Decoys: $n_decoys")
    println("═"^w)
    println(rpad("condition", 32), join(rpad.(thresh_strs, 10)))
    println("─"^w)

    # Print in a logical order
    ordered_names = [
        "LightGBM baseline", "LightGBM augmented",
        "Probit baseline", "Probit augmented",
        "Pipeline prob (ref)", "eigengap_raw (unsup)",
    ]
    for name in ordered_names
        qvals = all_qvals[name]
        counts = [targets_at_qvalue(qvals, is_target, t) for t in thresholds]
        println(rpad(name, 32), join(rpad.(string.(counts), 10)))
    end
    println("═"^w)

    # ── LightGBM feature importances (augmented model, last fold) ──
    println("\n── LightGBM Feature Importances (augmented, fold 3) ──")
    try
        gains = LightGBM.gain_importance(lgbm_aug_models[end])
        feat_imp = collect(zip(AUGMENTED_FEATURES, gains))
        sort!(feat_imp; by=x -> x[2], rev=true)
        for (feat, gain) in feat_imp
            println("  ", rpad(string(feat), 28), round(gain, digits=2))
        end
    catch e
        println("  Could not extract importances: $e")
    end

    # ── Plot q-value curves ──
    n_points = 200
    qv_range = range(0.0, 0.10, length=n_points)

    colors_map = Dict(
        "LightGBM baseline"    => :blue,
        "LightGBM augmented"   => :red,
        "Probit baseline"      => :cyan,
        "Probit augmented"     => :magenta,
        "Pipeline prob (ref)"  => :gray,
        "eigengap_raw (unsup)" => :orange,
    )
    styles_map = Dict(
        "LightGBM baseline"    => :solid,
        "LightGBM augmented"   => :solid,
        "Probit baseline"      => :dash,
        "Probit augmented"     => :dash,
        "Pipeline prob (ref)"  => :dot,
        "eigengap_raw (unsup)" => :dot,
    )

    p = plot(xlabel="q-value", ylabel="targets passing",
             title="FRAGCORR Model Comparison: targets at q-value threshold",
             legend=:topleft, size=(900, 550))

    for name in ordered_names
        qvals = all_qvals[name]
        counts = [targets_at_qvalue(qvals, is_target, t) for t in qv_range]
        plot!(p, collect(qv_range), counts;
              label=name, lw=2.5,
              color=colors_map[name],
              ls=styles_map[name])
    end

    display(p)

    return (plot=p, conditions=conditions, qvals=all_qvals,
            lgbm_aug_models=lgbm_aug_models)
end

# ============================================================
# Main
# ============================================================

println("\n╔══════════════════════════════════════════════════════════╗")
println("║  FRAGCORR Model Comparison: ML scoring with FRAGCORR   ║")
println("╚══════════════════════════════════════════════════════════╝")

mc_df = mc_load_all_files()

# Check for NaN/missing in feature columns
println("\nNaN check on feature columns:")
for col in AUGMENTED_FEATURES
    if hasproperty(mc_df, col)
        v = mc_df[!, col]
        n_nan = eltype(v) <: AbstractFloat ? count(isnan, v) : 0
        n_nan > 0 && println("  $col: $n_nan NaN values")
    else
        println("  $col: MISSING from DataFrame!")
    end
end
println("  (no output above = all clean)")

# Replace NaN with 0.0 for model features to avoid LightGBM/probit issues
for col in AUGMENTED_FEATURES
    if hasproperty(mc_df, col) && eltype(mc_df[!, col]) <: AbstractFloat
        mc_df[!, col] = replace(mc_df[!, col], NaN => 0.0)
    end
end

results = run_model_comparison(mc_df)

println("\nDone! Plot displayed.")
println("Access the combined DataFrame as `mc_df`.")
println("Re-run: results = run_model_comparison(mc_df)")
