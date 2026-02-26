"""
LightGBM-based cross-validated rescoring of first-pass PSMs.

After per-file probit scoring, this pools PSMs across files and trains a
LightGBM model using 2-fold cross-validation (train fold 0 → predict fold 1,
train fold 1 → predict fold 0). The improved probabilities replace the probit
scores and feed into `get_best_precursors_accross_runs`.
"""

const FIRST_PASS_LGBM_FEATURES = [
    :spectral_contrast, :city_block, :entropy_score, :scribe,
    :charge2, :poisson, :irt_error, :missed_cleavage, :Mox,
    :b_count, :TIC, :y_count, :err_norm, :spectrum_peak_count,
    #:sp_raw
]

const FIRST_PASS_LGBM_DERIVED_FEATURES = Symbol[#=:gof_chi2, :gof_log_pval=#]

const FIRST_PASS_LGBM_SAMPLE_SIZE = 1_000_000
const FIRST_PASS_LGBM_TRAIN_QVAL = 0.01f0
const FIRST_PASS_LGBM_ITER_SCHEME = [150, 300, 300]
const FIRST_PASS_LGBM_MIN_PEP_THRESHOLD = 0.90f0

"""
    _log_psm_diagnostics(targets, q_values, label)

Log target counts at multiple q-value thresholds for pooled PSMs.
"""
function _log_psm_diagnostics(targets::AbstractVector{Bool},
                              q_values::AbstractVector,
                              label::String)
    thresholds = (0.01f0, 0.05f0, 0.10f0, 0.15f0)
    parts = String[]
    for t in thresholds
        nt = count(i -> targets[i] && Float32(q_values[i]) <= t, eachindex(targets))
        push!(parts, "q≤$(t): $nt")
    end
    @user_info "  PSM targets ($label): " * join(parts, " | ")
end

"""
    _compute_global_precursor_diagnostics(all_psms, all_psms_paths, valid_indices,
                                          prec_is_decoy, fdr_scale_factor, label)

Replicate the log-odds global probability calculation from
`get_best_precursors_accross_runs` and log diagnostics at multiple FDR thresholds.
"""
function _compute_global_precursor_diagnostics(
    all_psms::DataFrame,
    all_psms_paths::Vector{String},
    valid_indices::Vector{Int},
    prec_is_decoy::AbstractVector{Bool},
    fdr_scale_factor::Float32,
    label::String
)
    # Collect best prob per precursor per file
    prec_probs_by_run = Dictionary{UInt32, Vector{Float32}}()
    n_valid_files = 0

    for ms_file_idx in valid_indices
        path = all_psms_paths[ms_file_idx]
        if isempty(path)
            continue
        end

        file_mask = all_psms.ms_file_idx .== UInt32(ms_file_idx)
        if !any(file_mask)
            continue
        end
        n_valid_files += 1

        prec_col = all_psms.precursor_idx[file_mask]
        prob_col = all_psms.prob[file_mask]

        file_best = Dictionary{UInt32, Float32}()
        for i in eachindex(prec_col)
            pid = prec_col[i]
            p = prob_col[i]
            if haskey(file_best, pid)
                file_best[pid] = max(file_best[pid], p)
            else
                insert!(file_best, pid, p)
            end
        end
        for (pid, p) in pairs(file_best)
            if haskey(prec_probs_by_run, pid)
                push!(prec_probs_by_run[pid], p)
            else
                insert!(prec_probs_by_run, pid, Float32[p])
            end
        end
    end

    if n_valid_files == 0
        @user_info "  Global precursors ($label): no valid files"
        return
    end

    sqrt_n_files = max(1, floor(Int, sqrt(n_valid_files)))
    n_unique = length(prec_probs_by_run)

    global_probs = Vector{Float32}(undef, n_unique)
    global_targets = Vector{Bool}(undef, n_unique)
    for (i, (pid, probs)) in enumerate(pairs(prec_probs_by_run))
        global_probs[i] = _logodds_combine(probs, sqrt_n_files)
        global_targets[i] = !prec_is_decoy[pid]
    end

    # Global q-values
    global_qvals = Vector{Float32}(undef, n_unique)
    get_qvalues!(global_probs, global_targets, global_qvals;
                 doSort=true, fdr_scale_factor=fdr_scale_factor)

    # Global PEP
    global_peps = Vector{Float32}(undef, n_unique)
    get_PEP!(global_probs, global_targets, global_peps;
             doSort=true, fdr_scale_factor=fdr_scale_factor)

    n_targets = count(global_targets)
    n_decoys = n_unique - n_targets

    # FDR thresholds
    fdr_parts = String[]
    for t in (0.01f0, 0.05f0, 0.10f0, 0.15f0)
        n = count(<=(t), global_qvals)
        push!(fdr_parts, "$(round(Int, t*100))%: $n")
    end

    # Global PEP thresholds
    pep_parts = String[]
    for t in (0.1f0, 0.5f0, 0.75f0, 0.9f0)
        nt = count(i -> global_peps[i] <= t && global_targets[i], eachindex(global_peps))
        nd = count(i -> global_peps[i] <= t && !global_targets[i], eachindex(global_peps))
        push!(pep_parts, "≤$t: T=$nt D=$nd")
    end

    @user_info "  Global precursors ($label): $n_unique unique ($n_targets T, $n_decoys D) " *
        "from $n_valid_files files (top_n=$sqrt_n_files)"
    @user_info "    FDR: " * join(fdr_parts, " | ")
    @user_info "    PEP: " * join(pep_parts, " | ")
end

"""
    _log_feature_importances(models, features)

Log gain-based feature importances averaged across CV fold models.
"""
function _log_feature_importances(models::Vector{LightGBMModel},
                                  features::Vector{Symbol})
    # Collect gain vectors from each fold
    all_gains = Vector{Vector{Float64}}()
    for m in models
        imp = importance(m)
        if imp !== nothing
            push!(all_gains, Float64[last(x) for x in imp])
        end
    end
    if isempty(all_gains)
        return
    end

    n_features = length(features)
    avg_gains = zeros(Float64, n_features)
    for gains in all_gains
        for j in 1:min(length(gains), n_features)
            avg_gains[j] += gains[j]
        end
    end
    avg_gains ./= length(all_gains)

    # Sort by importance descending
    order = sortperm(avg_gains; rev=true)
    total = sum(avg_gains)
    parts = String[]
    for idx in order
        g = avg_gains[idx]
        pct = total > 0 ? round(g / total * 100; digits=1) : 0.0
        push!(parts, "$(features[idx])=$(round(g; digits=1)) ($(pct)%)")
    end
    @user_info "  Feature importance (gain): " * join(parts, ", ")
end

"""
    _select_first_pass_training_data(sample_psms, iter_prob, iter_qval,
                                      test_fold, iteration, min_pep_threshold,
                                      max_q_value) -> (BitVector, Vector{Bool})

Select training data for one iteration of first-pass LightGBM rescoring.
Returns (train_mask, train_labels) where train_labels may relabel low-confidence
targets as negatives (negative mining). Inline equivalent of `QValueNegativeMining`.
"""
function _select_first_pass_training_data(
    sample_psms::DataFrame,
    iter_prob::Vector{Float32},
    iter_qval::Vector{Float32},
    test_fold::UInt8,
    iteration::Int,
    min_pep_threshold::Float32,
    max_q_value::Float32
)
    other_fold = sample_psms.cv_fold .!= test_fold

    if iteration == 1
        # Iteration 1: use all data from the other fold, no filtering
        train_mask = other_fold
        train_labels = Vector{Bool}(sample_psms.target[train_mask])
        return train_mask, train_labels
    end

    # Iteration 2+: negative mining on the other fold
    fold_indices = findall(other_fold)
    fold_probs = iter_prob[fold_indices]
    fold_targets = Vector{Bool}(sample_psms.target[fold_indices])

    # Compute PEP to relabel worst-scoring targets as negatives
    order = sortperm(fold_probs; rev=true)
    sorted_scores = fold_probs[order]
    sorted_targets = fold_targets[order]
    peps = Vector{Float32}(undef, length(order))
    get_PEP!(sorted_scores, sorted_targets, peps; doSort=false)

    idx_cutoff = findfirst(>=(min_pep_threshold), peps)
    if idx_cutoff !== nothing
        worst_idxs = order[idx_cutoff:end]
        fold_targets[worst_idxs] .= false
    end

    # Filter: keep all decoys + targets with q-value <= threshold
    fold_qvals = iter_qval[fold_indices]
    keep = BitVector([(!t) || (t && q <= max_q_value)
                      for (t, q) in zip(fold_targets, fold_qvals)])

    # Build full-size mask
    train_mask = falses(nrow(sample_psms))
    train_mask[fold_indices[keep]] .= true

    train_labels = fold_targets[keep]
    return train_mask, train_labels
end

"""
    _compute_gof_chi2_features!(all_psms)

Derive χ² goodness-of-fit features from the stored `sp_raw` column.

1. Recover S_p from log2 storage: S_p = 2^sp_raw - 1
2. Compute degrees of freedom: df = max(b_count + y_count - 1, 1)
3. Estimate proportionality constant ĉ via robust median on high-confidence targets:
   W_p = S_p / χ²_median(df), where χ²_median(k) ≈ k*(1-2/(9k))³
   ĉ = median(W_p)
4. For all PSMs: gof_chi2 = S_p / ĉ
5. gof_log_pval = log10(ccdf(Chisq(df), gof_chi2)), clamped at -30
"""
function _compute_gof_chi2_features!(all_psms::DataFrame)
    n = nrow(all_psms)

    # Recover S_p from log2(1 + S_p) storage
    sp_raw_col = Float64.(all_psms.sp_raw)
    s_p = @. exp2(sp_raw_col) - 1.0

    # Degrees of freedom per PSM
    df_vec = Int[max(Int(all_psms.b_count[i]) + Int(all_psms.y_count[i]) - 1, 1) for i in 1:n]

    # χ² median approximation: median(χ²_k) ≈ k*(1 - 2/(9k))³
    function chi2_median_approx(k::Int)
        kf = Float64(k)
        kf * (1.0 - 2.0 / (9.0 * kf))^3
    end

    # Estimate ĉ from high-confidence targets (q ≤ 0.01)
    good_mask = all_psms.target .& (Float32.(all_psms.q_value) .<= FIRST_PASS_LGBM_TRAIN_QVAL)
    good_idx = findall(good_mask)

    if length(good_idx) < 50
        @user_warn "GOF χ² features: only $(length(good_idx)) targets at q≤0.01, skipping"
        all_psms[!, :gof_chi2] = zeros(Float32, n)
        all_psms[!, :gof_log_pval] = zeros(Float32, n)
        return
    end

    # W_p = S_p / χ²_median(df) for good targets
    w_p = Float64[s_p[i] / max(chi2_median_approx(df_vec[i]), 1e-10) for i in good_idx]
    c_hat = median(w_p)
    c_hat = max(c_hat, 1e-10)  # safety floor

    @user_info "  GOF χ² features: ĉ = $(round(c_hat; digits=4)) from $(length(good_idx)) targets"

    # Compute features for all PSMs
    gof_chi2 = Vector{Float32}(undef, n)
    gof_log_pval = Vector{Float32}(undef, n)

    for i in 1:n
        chi2_val = Float64(s_p[i] / c_hat)
        gof_chi2[i] = Float32(chi2_val)

        # log10(ccdf(Chisq(df), chi2_val)), clamped at -30
        df_i = df_vec[i]
        if chi2_val <= 0.0 || df_i < 1
            gof_log_pval[i] = 0.0f0
        else
            pval = ccdf(Chisq(df_i), chi2_val)
            gof_log_pval[i] = Float32(max(log10(max(pval, 1e-30)), -30.0))
        end
    end

    all_psms[!, :gof_chi2] = gof_chi2
    all_psms[!, :gof_log_pval] = gof_log_pval
end

"""
    rescore_first_pass_with_lightgbm!(search_context, fdr_scale_factor)

Read all first-pass Arrow files, pool PSMs, train 2-fold CV LightGBM models,
rescore all PSMs, and write updated Arrow files back.
"""
function rescore_first_pass_with_lightgbm!(
    search_context::SearchContext,
    fdr_scale_factor::Float32
)
    # Collect paths from all valid files
    all_psms_paths = getFirstPassPsms(getMSData(search_context))
    valid_indices = get_valid_file_indices(search_context)
    valid_paths = [all_psms_paths[i] for i in valid_indices]

    if isempty(valid_paths)
        @user_warn "LightGBM rescore: no valid first-pass files, skipping"
        return
    end

    # 1. Read & pool all per-file Arrow files into a single DataFrame
    dfs = DataFrame[]
    for path in valid_paths
        tbl = Arrow.Table(path)
        if isempty(tbl[:ms_file_idx])
            continue
        end
        push!(dfs, DataFrame(tbl))
    end

    if isempty(dfs)
        @user_warn "LightGBM rescore: all files empty, skipping"
        return
    end

    all_psms = vcat(dfs...; cols=:intersect)
    dfs = nothing # free memory

    # 2. Add cv_fold column from library
    precursors = getPrecursors(getSpecLib(search_context))
    all_psms[!, :cv_fold] = UInt8[getCvFold(precursors, pid) for pid in all_psms.precursor_idx]

    # Derive target column from library if not present
    if !hasproperty(all_psms, :target)
        is_decoy = getIsDecoy(precursors)
        all_psms[!, :target] = Bool[!is_decoy[pid] for pid in all_psms.precursor_idx]
    end

    # 3. Compute derived GOF χ² features from sp_raw (needs target & q_value)
    #if hasproperty(all_psms, :sp_raw)
    #    _compute_gof_chi2_features!(all_psms)
    #end

    # 4. Determine which features are actually present (base + derived)
    available_features = Symbol[f for f in vcat(FIRST_PASS_LGBM_FEATURES, FIRST_PASS_LGBM_DERIVED_FEATURES)
                                if hasproperty(all_psms, f)]

    if length(available_features) < 3
        @user_warn "LightGBM rescore: only $(length(available_features)) features available, skipping"
        return
    end

    n_total = nrow(all_psms)
    n_targets_pre = count(all_psms.target)
    prec_is_decoy = getIsDecoy(precursors)

    # ── Before-rescore diagnostics ──
    @user_info "LightGBM rescore: $n_total pooled PSMs ($n_targets_pre targets)"
    _log_psm_diagnostics(all_psms.target, all_psms.q_value, "probit")
    _compute_global_precursor_diagnostics(
        all_psms, all_psms_paths, valid_indices,
        prec_is_decoy, fdr_scale_factor, "probit")

    # 5. Sample for training (keep all for prediction)
    n_sample = min(FIRST_PASS_LGBM_SAMPLE_SIZE, n_total)
    if n_sample < n_total
        sample_idx = sort(randperm(n_total)[1:n_sample])
        sample_psms = all_psms[sample_idx, :]
    else
        sample_psms = all_psms
    end

    # 6. Cross-validated training & prediction (2-fold CV, iterative with negative mining)
    new_probs = Vector{Float32}(undef, n_total)
    fill!(new_probs, 0.5f0)

    folds = sort(unique(all_psms.cv_fold))
    if length(folds) < 2
        @user_warn "LightGBM rescore: fewer than 2 CV folds found, skipping"
        return
    end

    trained_models = LightGBMModel[]

    for test_fold in folds
        # Initialize per-fold iteration vectors from probit scores
        # Use sample_psms indices for training-side prob/qval tracking
        iter_prob = Vector{Float32}(sample_psms.prob)
        iter_qval = Vector{Float32}(sample_psms.q_value)

        local final_model::LightGBMModel

        for (itr, num_round) in enumerate(FIRST_PASS_LGBM_ITER_SCHEME)
            train_mask, train_labels = _select_first_pass_training_data(
                sample_psms, iter_prob, iter_qval,
                test_fold, itr,
                FIRST_PASS_LGBM_MIN_PEP_THRESHOLD,
                FIRST_PASS_LGBM_TRAIN_QVAL
            )

            n_train_targets = count(train_labels)
            n_train_decoys = length(train_labels) - n_train_targets

            if n_train_targets < 100 || n_train_decoys < 100
                @user_warn "LightGBM rescore fold $test_fold iter $itr: insufficient training data " *
                    "(targets=$n_train_targets, decoys=$n_train_decoys), skipping all rescoring"
                return
            end

            @user_info "  Fold $test_fold iter $itr/$( length(FIRST_PASS_LGBM_ITER_SCHEME)): " *
                "$n_train_targets targets, $n_train_decoys decoys, $num_round rounds"

            # Build and train fresh model (matches SimpleLightGBM second-pass config)
            classifier = build_lightgbm_classifier(
                num_iterations = num_round,
                max_depth = 4,
                num_leaves = 15,
                learning_rate = 0.1,
                feature_fraction = 0.8,
                bagging_fraction = 0.8,
                bagging_freq = 1,
                min_data_in_leaf = 20,
                min_gain_to_split = 0.1
            )

            train_data = sample_psms[train_mask, available_features]
            model = fit_lightgbm_model(classifier, train_data, train_labels)
            final_model = model

            # If not the last iteration, predict on the training fold's sample data
            # to update iter_prob/iter_qval for the next iteration's negative mining
            if itr < length(FIRST_PASS_LGBM_ITER_SCHEME)
                train_fold_mask = sample_psms.cv_fold .!= test_fold
                train_fold_data = sample_psms[train_fold_mask, available_features]
                train_fold_preds = predict(model, train_fold_data)

                iter_prob[train_fold_mask] .= train_fold_preds

                # Recompute q-values on the training fold for negative mining
                tf_indices = findall(train_fold_mask)
                tf_probs = iter_prob[tf_indices]
                tf_targets = Vector{Bool}(sample_psms.target[tf_indices])
                tf_qvals = Vector{Float32}(undef, length(tf_indices))
                get_qvalues!(tf_probs, tf_targets, tf_qvals; doSort=true)
                iter_qval[tf_indices] .= tf_qvals
            end
        end

        push!(trained_models, final_model)

        # Predict on ALL PSMs in the test fold using the final model
        test_mask = all_psms.cv_fold .== test_fold
        test_data = all_psms[test_mask, available_features]
        preds = predict(final_model, test_data)

        new_probs[test_mask] .= preds
    end

    # Log feature importances (averaged across folds)
    _log_feature_importances(trained_models, available_features)

    # 7. Update scores: replace prob, recompute q_value and PEP
    all_psms[!, :prob] = new_probs

    # Recompute q-values (must stay Float16 for downstream readPSMs! compatibility)
    q_values = Vector{Float32}(undef, n_total)
    get_qvalues!(new_probs, all_psms.target, q_values;
                 doSort=true, fdr_scale_factor=fdr_scale_factor)
    all_psms[!, :q_value] = Float16.(q_values)

    # Recompute PEP
    pep_values = Vector{Float16}(undef, n_total)
    get_PEP!(new_probs, all_psms.target, pep_values;
             doSort=true, fdr_scale_factor=fdr_scale_factor)
    all_psms[!, :PEP] = pep_values

    # ── After-rescore diagnostics ──
    _log_psm_diagnostics(all_psms.target, all_psms.q_value, "LightGBM")
    _compute_global_precursor_diagnostics(
        all_psms, all_psms_paths, valid_indices,
        prec_is_decoy, fdr_scale_factor, "LightGBM")

    # 8. Write back to per-file Arrow files — drop cv_fold and derived GOF columns
    drop_cols = Symbol[:cv_fold]
    for f in FIRST_PASS_LGBM_DERIVED_FEATURES
        hasproperty(all_psms, f) && push!(drop_cols, f)
    end
    select!(all_psms, Not(drop_cols))

    grouped = groupby(all_psms, :ms_file_idx)
    for (key, file_psms) in pairs(grouped)
        ms_file_idx = Int(key.ms_file_idx)
        if ms_file_idx < 1 || ms_file_idx > length(all_psms_paths)
            continue
        end
        path = all_psms_paths[ms_file_idx]
        if isempty(path)
            continue
        end
        Arrow.write(path, file_psms)
    end

    return
end
