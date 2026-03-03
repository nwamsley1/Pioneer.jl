#=
FRAGCORR Focused Analysis — Narrow RT Window Fragment Correlation
=================================================================
Self-contained script that:
1. Loads focused FRAGCORR traces (±15s around best PSM RT) and computes unsupervised scores
2. Joins with PSM scoring features
3. Trains LightGBM/probit models with and without FRAGCORR features
4. Validates by checking DIA-NN precursor recovery

Usage:
    include("scripts/fragcorr_focused_analysis.jl")
=#

# ============================================================
# Step 1: Include + Setup
# ============================================================

include(joinpath(@__DIR__, "fragcorr_unsupervised_scores.jl"))

using LightGBM
using SpecialFunctions: erf
using CSV

const FC_DATA_DIR = "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/OlsenEclipse_fragcorr_focused/temp_data"
const FC_PSM_DIR  = joinpath(FC_DATA_DIR, "first_pass_psms")
const FC_CORR_DIR = joinpath(FC_DATA_DIR, "first_pass_corr")

# ============================================================
# Step 2: Load & Join Data
# ============================================================

function fc_load_and_join(filename::String)
    psm_path  = joinpath(FC_PSM_DIR, filename)
    corr_path = joinpath(FC_CORR_DIR, filename)

    psm_df = DataFrame(Tables.columntable(Arrow.Table(psm_path)))
    gdf = load_fragcorr_data(corr_path)
    fragcorr_df = score_all_precursors(gdf)

    combined = innerjoin(psm_df, fragcorr_df;
                         on=[:precursor_idx, :target],
                         makeunique=true)
    println("  Joined: $(nrow(psm_df)) PSMs x $(nrow(fragcorr_df)) FRAGCORR -> $(nrow(combined)) combined")
    return combined
end

function fc_load_all_files()
    files = filter(f -> endswith(f, ".arrow"), readdir(FC_PSM_DIR))
    corr_files = Set(readdir(FC_CORR_DIR))
    files = filter(f -> f in corr_files, files)

    println("\nLoading $(length(files)) file pairs...")
    dfs = DataFrame[]
    for f in files
        println("\n-- $f --")
        push!(dfs, fc_load_and_join(f))
    end

    combined = vcat(dfs...; cols=:union)
    n_target = sum(combined.target)
    n_decoy  = sum(.!combined.target)
    println("\n========================================")
    println("Combined: $(nrow(combined)) rows ($n_target targets, $n_decoy decoys)")
    println("========================================")
    return combined
end

# ============================================================
# Step 3: Probit IRLS (standalone, no external dependencies)
# ============================================================

function probit_train(X::Matrix{Float64}, y::Vector{Bool};
                      max_iter::Int=30, tol::Float64=1e-2,
                      z_bounds::Tuple{Float64,Float64}=(-8.0, 8.0))
    n, p = size(X)
    beta = zeros(p)
    eta = zeros(n)
    Z = zeros(n)
    W = zeros(n)
    beta_old = zeros(p)
    old_loss = -Inf

    for _ in 1:max_iter
        copyto!(beta_old, beta)

        mul!(eta, X, beta)
        @inbounds for i in 1:n
            eta[i] = clamp(eta[i], z_bounds[1], z_bounds[2])
        end

        loss = 0.0
        @inbounds for i in 1:n
            phi = exp(-eta[i]^2 / 2) / sqrt(2pi)
            mu = (1 + erf(eta[i] / sqrt(2))) / 2
            mu = clamp(mu, 1e-12, 1 - 1e-12)
            if y[i]
                Z[i] = eta[i] + (1 - mu) / phi
                loss += log(mu)
            else
                Z[i] = eta[i] - mu / phi
                loss += log1p(-mu)
            end
            W[i] = phi^2 / (mu * (1 - mu))
        end

        XW = X' * Diagonal(W)
        XWX = XW * X
        XWZ = XW * Z
        beta_new = XWX \ XWZ

        delta_beta = beta_new .- beta
        step = 1.0
        accepted = false
        while step > 1e-4
            beta_trial = beta .+ step .* delta_beta
            eta_trial = X * beta_trial
            @inbounds for i in 1:n
                eta_trial[i] = clamp(eta_trial[i], z_bounds[1], z_bounds[2])
            end
            trial_loss = 0.0
            @inbounds for i in 1:n
                mu_t = clamp((1 + erf(eta_trial[i] / sqrt(2))) / 2, 1e-12, 1 - 1e-12)
                trial_loss += y[i] ? log(mu_t) : log1p(-mu_t)
            end
            if trial_loss >= loss
                beta .= beta_trial
                loss = trial_loss
                accepted = true
                break
            end
            step *= 0.5
        end
        if !accepted
            beta .+= step .* delta_beta
        end

        if norm(beta .- beta_old) < tol || abs(loss - old_loss) < tol
            break
        end
        old_loss = loss
    end

    return beta
end

function probit_predict(X::Matrix{Float64}, beta::Vector{Float64})
    eta = X * beta
    probs = Vector{Float64}(undef, length(eta))
    @inbounds for i in eachindex(eta)
        probs[i] = (1 + erf(eta[i] / sqrt(2))) / 2
    end
    return probs
end

# ============================================================
# Step 4: Cross-Validated Model Training
# ============================================================

const FC_BASELINE_FEATURES = [:scribe, :spectral_contrast, :city_block, :entropy_score, :y_count, :err_norm]
const FC_FRAGCORR_FEATURES = [:eigengap_raw, :lambda1_frac_raw, :ev1_raw, :median_corr_raw, :best_frag_corr_mean]
const FC_AUGMENTED_FEATURES = vcat(FC_BASELINE_FEATURES, FC_FRAGCORR_FEATURES)

fc_cv_fold_assignment(precursor_idx::AbstractVector) = Int.(precursor_idx .% 3 .+ 1)

function fc_run_cv_lightgbm(df::DataFrame, features::Vector{Symbol}; n_folds::Int=3)
    folds = fc_cv_fold_assignment(df.precursor_idx)
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
            num_class = 1,
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
        y_hat = ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
        preds[test_mask] .= Float64.(y_hat)
        fold_models[k] = estimator
    end

    return preds, fold_models
end

function fc_run_cv_probit(df::DataFrame, features::Vector{Symbol}; n_folds::Int=3)
    folds = fc_cv_fold_assignment(df.precursor_idx)
    preds = fill(NaN, nrow(df))

    for k in 1:n_folds
        train_mask = folds .!= k
        test_mask  = folds .== k

        X_train = Matrix{Float64}(df[train_mask, features])
        y_train = Vector{Bool}(df[train_mask, :target])
        X_test  = Matrix{Float64}(df[test_mask, features])

        beta = probit_train(X_train, y_train)
        preds[test_mask] .= probit_predict(X_test, beta)
    end

    return preds
end

# ============================================================
# Step 5: Model Comparison
# ============================================================

function fc_run_model_comparison(df::DataFrame;
                                 thresholds=[0.001, 0.005, 0.01, 0.02, 0.05, 0.10])
    is_target = df.target
    n_targets = sum(is_target)
    n_decoys  = sum(.!is_target)

    conditions = Dict{String, Vector{Float64}}()

    println("\n-- Training LightGBM (baseline) --")
    preds_lgbm_base, _ = fc_run_cv_lightgbm(df, FC_BASELINE_FEATURES)
    conditions["LightGBM baseline"] = preds_lgbm_base

    println("\n-- Training LightGBM (augmented) --")
    preds_lgbm_aug, lgbm_aug_models = fc_run_cv_lightgbm(df, FC_AUGMENTED_FEATURES)
    conditions["LightGBM augmented"] = preds_lgbm_aug

    println("\n-- Training Probit (baseline) --")
    preds_probit_base = fc_run_cv_probit(df, FC_BASELINE_FEATURES)
    conditions["Probit baseline"] = preds_probit_base

    println("\n-- Training Probit (augmented) --")
    preds_probit_aug = fc_run_cv_probit(df, FC_AUGMENTED_FEATURES)
    conditions["Probit augmented"] = preds_probit_aug

    conditions["Pipeline prob (ref)"] = Float64.(df.prob)
    conditions["eigengap_raw (unsup)"] = Float64.(df.eigengap_raw)

    # Compute q-values
    all_qvals = Dict{String, Vector{Float64}}()
    for (name, scores) in conditions
        all_qvals[name] = compute_qvalues(scores, is_target; higher_is_better=true)
    end

    # Print summary table
    ordered_names = [
        "LightGBM baseline", "LightGBM augmented",
        "Probit baseline", "Probit augmented",
        "Pipeline prob (ref)", "eigengap_raw (unsup)",
    ]
    thresh_strs = [string(round(t * 100, digits=1)) * "%" for t in thresholds]
    w = 32 + 10 * length(thresholds)
    println("\n\nTargets: $n_targets  |  Decoys: $n_decoys")
    println("="^w)
    println(rpad("condition", 32), join(rpad.(thresh_strs, 10)))
    println("-"^w)

    for name in ordered_names
        qvals = all_qvals[name]
        counts = [targets_at_qvalue(qvals, is_target, t) for t in thresholds]
        println(rpad(name, 32), join(rpad.(string.(counts), 10)))
    end
    println("="^w)

    # LightGBM feature importances
    println("\n-- LightGBM Feature Importances (augmented, fold 3) --")
    try
        gains = LightGBM.gain_importance(lgbm_aug_models[end])
        feat_imp = collect(zip(FC_AUGMENTED_FEATURES, gains))
        sort!(feat_imp; by=x -> x[2], rev=true)
        for (feat, gain) in feat_imp
            println("  ", rpad(string(feat), 28), round(gain, digits=2))
        end
    catch e
        println("  Could not extract importances: $e")
    end

    return (conditions=conditions, qvals=all_qvals, lgbm_aug_models=lgbm_aug_models,
            ordered_names=ordered_names)
end

# ============================================================
# Step 6: DIA-NN Recovery
# ============================================================

const FC_ModTuple = Tuple{Int, Int}
const FC_CanonicalKey = Tuple{String, Vector{FC_ModTuple}, Int}

function fc_parse_diann_modified_seq(mod_seq::AbstractString, charge::Int)::FC_CanonicalKey
    stripped = IOBuffer()
    mods = FC_ModTuple[]
    pos = 0
    i = 1
    while i <= ncodeunits(mod_seq)
        c = mod_seq[i]
        if c == '('
            close_idx = findnext(')', mod_seq, i)
            mod_str = mod_seq[i+1:close_idx-1]
            colon_idx = findlast(':', mod_str)
            unimod_id = parse(Int, mod_str[colon_idx+1:end])
            push!(mods, (pos, unimod_id))
            i = close_idx + 1
        else
            pos += 1
            write(stripped, c)
            i += 1
        end
    end
    sort!(mods)
    return (String(take!(stripped)), mods, charge)
end

function fc_parse_pioneer_structural_mods(mods_str::AbstractString)::Vector{FC_ModTuple}
    isempty(mods_str) && return FC_ModTuple[]
    mods = FC_ModTuple[]
    for m in eachmatch(r"\((\d+),\w,(\w+:\d+)\)", mods_str)
        pos = parse(Int, m.captures[1])
        mod_id_str = m.captures[2]
        colon_idx = findlast(':', mod_id_str)
        unimod_id = parse(Int, mod_id_str[colon_idx+1:end])
        push!(mods, (pos, unimod_id))
    end
    sort!(mods)
    return mods
end

function fc_pioneer_canonical_key(seq::AbstractString, structural_mods::AbstractString, charge::Integer)::FC_CanonicalKey
    mods = fc_parse_pioneer_structural_mods(structural_mods)
    return (String(seq), mods, Int(charge))
end

const FC_LIBRARY_PATH = "/Users/nathanwamsley/Data/SPEC_LIBS/3P_rank_based.poin"
const FC_DIANN_PR_MATRIX = joinpath("/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse",
    "DIA-NN_results", "OlsenExplorisThreeProteome500ng-11-24-2025-report.pr_matrix.tsv")
const FC_PIONEER_FILES = ["E45H50Y5_2", "E45H50Y5_3", "E45H50Y5_4",
                          "E5H50Y45_1", "E5H50Y45_2", "E5H50Y45_4"]

function fc_diann_col_to_pioneer_file(col::AbstractString)::Union{String, Nothing}
    m = match(r"(E\d+H\d+Y\d+)_\d+SPD_DIA_(\d+)", col)
    m === nothing && return nothing
    return "$(m.captures[1])_$(m.captures[2])"
end

function fc_run_diann_recovery(fc_df::DataFrame, results)
    is_target = fc_df.target
    ordered_names = results.ordered_names

    # --- Pioneer library: precursor_idx -> canonical key ---
    println("\nLoading Pioneer library...")
    lib_df = DataFrame(Arrow.Table(joinpath(FC_LIBRARY_PATH, "precursors_table.arrow")))

    mc_pidx_set = Set{UInt32}(fc_df[fc_df.target, :precursor_idx])
    println("  Unique target precursor_idx in fc_df: $(length(mc_pidx_set))")

    pidx_to_key = Dict{UInt32, FC_CanonicalKey}()
    for i in 1:nrow(lib_df)
        pidx = UInt32(i)
        pidx in mc_pidx_set || continue
        lib_df.is_decoy[i] && continue
        seq = lib_df.sequence[i]
        smods_raw = lib_df.structural_mods[i]
        smods = ismissing(smods_raw) ? "" : smods_raw
        charge = Int(lib_df.prec_charge[i])
        pidx_to_key[pidx] = fc_pioneer_canonical_key(seq, smods, charge)
    end
    println("  Built $(length(pidx_to_key)) precursor_idx -> key mappings")

    # --- DIA-NN fair key set ---
    println("\nLoading DIA-NN pr_matrix...")
    diann_df = CSV.read(FC_DIANN_PR_MATRIX, DataFrame; delim='\t')
    n_diann = nrow(diann_df)
    println("  Loaded $n_diann precursors")

    metadata_cols = Set(["Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                         "First.Protein.Description", "Proteotypic",
                         "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge", "Precursor.Id"])
    run_cols = [c for c in names(diann_df) if c ∉ metadata_cols]

    pioneer_run_cols = String[]
    for col in run_cols
        pf = fc_diann_col_to_pioneer_file(col)
        pf !== nothing && pf in Set(FC_PIONEER_FILES) && push!(pioneer_run_cols, col)
    end
    println("  Pioneer run columns matched: $(length(pioneer_run_cols))")

    diann_fair_keys = Set{FC_CanonicalKey}()
    diann_all_keys = Set{FC_CanonicalKey}()
    for i in 1:n_diann
        key = fc_parse_diann_modified_seq(
            diann_df[i, "Modified.Sequence"],
            Int(diann_df[i, "Precursor.Charge"])
        )
        push!(diann_all_keys, key)
        for col in pioneer_run_cols
            val = diann_df[i, col]
            if !ismissing(val) && val isa Number && val > 0
                push!(diann_fair_keys, key)
                break
            end
        end
    end
    println("  DIA-NN total unique keys: $(length(diann_all_keys))")
    println("  DIA-NN fair keys (detected in Pioneer runs): $(length(diann_fair_keys))")

    # --- Annotate fc_df with is_diann ---
    println("\nAnnotating fc_df...")
    is_diann = falses(nrow(fc_df))
    let n_mapped = 0, n_matched = 0
        for i in 1:nrow(fc_df)
            fc_df.target[i] || continue
            pidx = UInt32(fc_df.precursor_idx[i])
            key = get(pidx_to_key, pidx, nothing)
            key === nothing && continue
            n_mapped += 1
            if key in diann_fair_keys
                is_diann[i] = true
                n_matched += 1
            end
        end
        println("  Targets mapped to keys: $n_mapped")
        println("  Targets matching DIA-NN fair set: $n_matched")
    end
    fc_df.is_diann = is_diann

    # --- Recovery tables ---
    thresholds = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10]
    is_diann_col = fc_df.is_diann

    println("\n" * "="^70)
    println("DIA-NN OVERLAP AT MULTIPLE Q-VALUE THRESHOLDS")
    println("="^70)

    thresh_strs = [string(round(t * 100, digits=1)) * "%" for t in thresholds]
    w = 32 + 18 * length(thresholds)
    println("\n" * rpad("Condition", 32) * join([rpad("$t (tot/diann/%)", 18) for t in thresh_strs]))
    println("-"^w)

    recovery_data = Dict{String, Dict{Float64, NamedTuple}}()

    for name in ordered_names
        qvals = results.qvals[name]
        recovery_data[name] = Dict{Float64, NamedTuple}()

        row_parts = String[rpad(name, 32)]
        for t in thresholds
            passing = qvals .<= t
            total_targets = count(passing .& is_target)
            diann_targets = count(passing .& is_target .& is_diann_col)
            pct = total_targets > 0 ? 100.0 * diann_targets / total_targets : 0.0
            push!(row_parts, rpad("$total_targets/$diann_targets/$(round(pct, digits=1))%", 18))
            recovery_data[name][t] = (total=total_targets, diann=diann_targets, pct=pct)
        end
        println(join(row_parts))
    end
    println("-"^w)

    # Clean summary at 1% FDR
    println("\n\n" * "="^60)
    println("TABLE 1: DIA-NN Overlap at 1% FDR")
    println("="^60)
    println(rpad("Condition", 32) * rpad("Total", 10) * rpad("DIA-NN", 10) * "DIA-NN%")
    println("-"^60)
    for name in ordered_names
        r = recovery_data[name][0.01]
        println(rpad(name, 32) * rpad(string(r.total), 10) * rpad(string(r.diann), 10) *
                "$(round(r.pct, digits=1))%")
    end
    println("-"^60)

    # --- Venn decomposition ---
    println("\n\n" * "="^60)
    println("TABLE 2: Venn Decomposition (LightGBM baseline vs augmented, 1% FDR)")
    println("="^60)

    qv_base = results.qvals["LightGBM baseline"]
    qv_aug  = results.qvals["LightGBM augmented"]

    pass_base = (qv_base .<= 0.01) .& is_target
    pass_aug  = (qv_aug  .<= 0.01) .& is_target

    both_models = pass_base .& pass_aug
    aug_only    = pass_aug .& .!pass_base
    base_only   = pass_base .& .!pass_aug

    n_both      = count(both_models)
    n_aug_only  = count(aug_only)
    n_base_only = count(base_only)

    diann_both      = count(both_models .& is_diann_col)
    diann_aug_only  = count(aug_only .& is_diann_col)
    diann_base_only = count(base_only .& is_diann_col)

    pct_aug_only_diann  = n_aug_only > 0 ? 100.0 * diann_aug_only / n_aug_only : 0.0
    pct_base_only_diann = n_base_only > 0 ? 100.0 * diann_base_only / n_base_only : 0.0
    pct_both_diann      = n_both > 0 ? 100.0 * diann_both / n_both : 0.0

    println("Both models:          $(rpad(n_both, 8))  ($diann_both = $(round(pct_both_diann, digits=1))% are DIA-NN validated)")
    println("Aug only (FRAGCORR):  $(rpad(n_aug_only, 8))  ($diann_aug_only = $(round(pct_aug_only_diann, digits=1))% are DIA-NN validated)")
    println("Base only:            $(rpad(n_base_only, 8))  ($diann_base_only = $(round(pct_base_only_diann, digits=1))% are DIA-NN validated)")
    println("-"^60)

    return recovery_data
end

# ============================================================
# Step 7: Comparison with Broad-Window Results
# ============================================================

function fc_print_broad_comparison()
    println("\n\n" * "="^70)
    println("COMPARISON: Focused (+-15s) vs Broad (full fragment index window)")
    println("="^70)
    println("""
    To compare, run both scripts and note targets at 1% FDR:

    Broad window:  include("scripts/fragcorr_model_comparison.jl")
    Focused:       include("scripts/fragcorr_focused_analysis.jl")  (this script)

    Key questions:
    1. Does focused FRAGCORR improve LightGBM augmented targets at 1% FDR?
    2. Is the DIA-NN overlap % higher for focused FRAGCORR?
    3. Are aug-only precursors more likely to be DIA-NN validated?
    """)
end

# ============================================================
# Step 8: Save Plots
# ============================================================

function fc_save_plots(fc_df::DataFrame, results)
    is_target = fc_df.target
    is_diann_col = hasproperty(fc_df, :is_diann) ? fc_df.is_diann : falses(nrow(fc_df))
    ordered_names = results.ordered_names

    results_dir = joinpath(FC_DATA_DIR, "..", "results")
    mkpath(results_dir)

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

    # Q-value curve: all conditions
    p_total = plot(xlabel="q-value", ylabel="Targets passing",
                   title="Focused FRAGCORR: targets at q-value threshold",
                   legend=:topleft, size=(900, 550))
    for name in ordered_names
        qvals = results.qvals[name]
        counts = [targets_at_qvalue(qvals, is_target, t) for t in qv_range]
        plot!(p_total, collect(qv_range), counts;
              label=name, lw=2.5, color=colors_map[name], ls=styles_map[name])
    end
    savefig(p_total, joinpath(results_dir, "focused_qvalue_curves.png"))
    println("  Saved focused_qvalue_curves.png")

    # LightGBM baseline vs augmented with DIA-NN overlay
    p_combined = plot(xlabel="q-value", ylabel="Targets passing",
                      title="Focused FRAGCORR: LightGBM total vs DIA-NN validated",
                      legend=:topleft, size=(900, 550))
    for name in ["LightGBM baseline", "LightGBM augmented"]
        qvals = results.qvals[name]
        col = colors_map[name]
        total_counts = [count((qvals .<= t) .& is_target) for t in qv_range]
        diann_counts = [count((qvals .<= t) .& is_target .& is_diann_col) for t in qv_range]
        plot!(p_combined, collect(qv_range), total_counts;
              label="$name (total)", lw=2.5, color=col, ls=:solid)
        plot!(p_combined, collect(qv_range), diann_counts;
              label="$name (DIA-NN)", lw=2.0, color=col, ls=:dash)
    end
    savefig(p_combined, joinpath(results_dir, "focused_lgbm_diann_overlay.png"))
    println("  Saved focused_lgbm_diann_overlay.png")

    display(p_combined)
    return (p_total=p_total, p_combined=p_combined)
end

# ============================================================
# Main
# ============================================================

println("\n+========================================================+")
println("|  FRAGCORR Focused Analysis: +-15s RT Window            |")
println("+========================================================+")

fc_df = fc_load_all_files()

# NaN check and cleanup
println("\nNaN check on feature columns:")
for col in FC_AUGMENTED_FEATURES
    if hasproperty(fc_df, col)
        v = fc_df[!, col]
        n_nan = eltype(v) <: AbstractFloat ? count(isnan, v) : 0
        n_nan > 0 && println("  $col: $n_nan NaN values")
    else
        println("  $col: MISSING from DataFrame!")
    end
end
println("  (no output above = all clean)")

for col in FC_AUGMENTED_FEATURES
    if hasproperty(fc_df, col) && eltype(fc_df[!, col]) <: AbstractFloat
        fc_df[!, col] = replace(fc_df[!, col], NaN => 0.0)
    end
end

# Run model comparison
fc_results = fc_run_model_comparison(fc_df)

# Run DIA-NN recovery analysis
println("\n+========================================================+")
println("|  DIA-NN Recovery Analysis                              |")
println("+========================================================+")
fc_recovery = fc_run_diann_recovery(fc_df, fc_results)

# Print broad comparison notes
fc_print_broad_comparison()

# Save plots
println("\nSaving plots...")
fc_plots = fc_save_plots(fc_df, fc_results)

println("\nDone! Access results:")
println("  fc_df       -- Combined DataFrame (with :is_diann column)")
println("  fc_results  -- Model comparison results (conditions, qvals, models)")
println("  fc_recovery -- DIA-NN recovery data dict")
println("  fc_plots    -- Saved plot objects")
