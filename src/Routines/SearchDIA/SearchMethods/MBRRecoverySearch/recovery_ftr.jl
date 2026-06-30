# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# Recovery FTR — paired-counterfactual false-transfer-rate model over the
# MBRRecoverySearch recovery seeds (V1.5).
#
# Each recovery seed already carries the MBR_*_true features (the receiver
# extraction compared to the genuine cross-file donor of the SAME precursor).
# Here we compute the MBR_*_false counterfactual contrast (the receiver compared
# to the nearest opposite-class donor of a DIFFERENT precursor, reusing
# PrecursorScoringSearch's partner-pool machinery), train a 2-fold CV LightGBM on
# the doubled (true/false) frame with is_real labels — mirroring
# apply_mbr_filter_paired! in mbr_ftr.jl — and report the funnel:
#
#     candidates  →  passing FTR (PEP ≤ α)  →  target/decoy q ≤ fdr
#
# plus the model's gain feature importances.
#
# Recovery rows have no main-search score (trace_prob = 0), so unlike the main
# MBR path — which leans on the downstream :global_qval for target/decoy FDR —
# the recovery set's only FDR handle is its own target/decoy composition. The
# "<1% FDR" count is therefore a target/decoy q-value computed directly on the
# FTR is_real score.

const RECOVERY_FTR_FEATURES_TRUE = Symbol[
    :MBR_max_pair_prob_true,
    :MBR_log2_weight_ratio_true,
    :MBR_log2_explained_ratio_true,
    :MBR_best_irt_diff_true,
    :MBR_best_rt_diff_true,
    :MBR_log_by_diff_true,
    :MBR_smoothed_frag_hellinger_true,
]
const RECOVERY_FTR_FEATURES_FALSE = Symbol[
    :MBR_max_pair_prob_false,
    :MBR_log2_weight_ratio_false,
    :MBR_log2_explained_ratio_false,
    :MBR_best_irt_diff_false,
    :MBR_best_rt_diff_false,
    :MBR_log_by_diff_false,
    :MBR_smoothed_frag_hellinger_false,
]

"""
    _collect_recovery_ftr_rows(psm_paths, donor_dict, partner_pools, fragment_keys) -> DataFrame

Reads each per-file recovery seed sidecar (in `psm_paths` order so the file
index matches the donor pool / `my_file` exclusion), reconstructs the
`ReceiverExtraction`, resolves a counterfactual donor, and computes the
MBR_*_false features. Returns one row per seed with the 7 _true (from sidecar)
and 7 _false (computed) features, the library target flag, cv_fold, and a
`MBR_is_missing_false` flag for seeds with no counterfactual partner.
"""
function _collect_recovery_ftr_rows(
    psm_paths::Vector{String},
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    fragment_keys::_MBRFragmentAnnotationKeys,
)
    pid_c = UInt32[]; tgt_c = Bool[]; cv_c = UInt8[]
    t1 = Float32[]; t2 = Float32[]; t3 = Float32[]; t4 = Float32[]; t5 = Float32[]; t6 = Float32[]; t7 = Float32[]
    f1 = Float32[]; f2 = Float32[]; f3 = Float32[]; f4 = Float32[]; f5 = Float32[]; f6 = Float32[]; f7 = Float32[]
    miss_f = Bool[]

    false_cache = Dict{UInt32, Union{Nothing, _MBRDonorEntry}}()
    for (file_idx, psm_path) in enumerate(psm_paths)
        sc = read_recovery_seed_sidecar(recovery_seed_sidecar_path(psm_path))
        nrow(sc) == 0 && continue
        my_file = UInt32(file_idx)
        @inbounds for i in 1:nrow(sc)
            P = UInt32(sc.precursor_idx[i])
            receiver = ReceiverExtraction(
                UInt32(sc.scan_idx[i]),
                Float32(sc.weight[i]),
                Float32(sc.log2_intensity_explained[i]),
                Float32(sc.irt_obs[i]),
                Float32(sc.rt_obs[i]),
                Float32(sc.log_by_ratio[i]),
                (Float32(sc.frag1_smoothed_intensity[i]), Float32(sc.frag2_smoothed_intensity[i]),
                 Float32(sc.frag3_smoothed_intensity[i]), Float32(sc.frag4_smoothed_intensity[i]),
                 Float32(sc.frag5_smoothed_intensity[i]), Float32(sc.frag6_smoothed_intensity[i]),
                 Float32(sc.frag7_smoothed_intensity[i]), Float32(sc.frag8_smoothed_intensity[i])),
            )
            push!(pid_c, P); push!(tgt_c, Bool(sc.target[i])); push!(cv_c, UInt8(sc.cv_fold[i]))
            push!(t1, Float32(sc.MBR_max_pair_prob_true[i]))
            push!(t2, Float32(sc.MBR_log2_weight_ratio_true[i]))
            push!(t3, Float32(sc.MBR_log2_explained_ratio_true[i]))
            push!(t4, Float32(sc.MBR_best_irt_diff_true[i]))
            push!(t5, Float32(sc.MBR_best_rt_diff_true[i]))
            push!(t6, Float32(sc.MBR_log_by_diff_true[i]))
            push!(t7, Float32(sc.MBR_smoothed_frag_hellinger_true[i]))

            fdonor = _false_donor_for_pid(false_cache, donor_dict, partner_pools, P, my_file)
            if fdonor === nothing
                push!(miss_f, true)
                # Sentinels mirror the MBR feature convention (-1 generic,
                # 1 for Hellinger). These rows are excluded as non-candidates.
                push!(f1, -1f0); push!(f2, -1f0); push!(f3, -1f0); push!(f4, -1f0)
                push!(f5, -1f0); push!(f6, -1f0); push!(f7, 1f0)
            else
                ff = compute_seed_mbr_features(receiver, fdonor, Float32(sc.irt_pred[i]), fragment_keys, P)
                push!(miss_f, false)
                push!(f1, ff.MBR_max_pair_prob_true)
                push!(f2, ff.MBR_log2_weight_ratio_true)
                push!(f3, ff.MBR_log2_explained_ratio_true)
                push!(f4, ff.MBR_best_irt_diff_true)
                push!(f5, ff.MBR_best_rt_diff_true)
                push!(f6, ff.MBR_log_by_diff_true)
                push!(f7, ff.MBR_smoothed_frag_hellinger_true)
            end
        end
    end

    return DataFrame(
        precursor_idx = pid_c, target = tgt_c, cv_fold = cv_c,
        MBR_max_pair_prob_true = t1, MBR_log2_weight_ratio_true = t2,
        MBR_log2_explained_ratio_true = t3, MBR_best_irt_diff_true = t4,
        MBR_best_rt_diff_true = t5, MBR_log_by_diff_true = t6,
        MBR_smoothed_frag_hellinger_true = t7,
        MBR_max_pair_prob_false = f1, MBR_log2_weight_ratio_false = f2,
        MBR_log2_explained_ratio_false = f3, MBR_best_irt_diff_false = f4,
        MBR_best_rt_diff_false = f5, MBR_log_by_diff_false = f6,
        MBR_smoothed_frag_hellinger_false = f7,
        MBR_is_missing_false = miss_f,
    )
end

"""
    run_recovery_ftr!(search_context, results; alpha=0.01, fdr=0.01) -> NamedTuple | nothing

Trains the recovery FTR over all emitted recovery seeds and logs the funnel +
feature importances. Returns the funnel counts (or `nothing` if there are no
seeds / counterfactuals). Requires `results` to be initialized (donor pool +
fragment keys built during the per-file extraction).
"""
function run_recovery_ftr!(
    search_context::SearchContext,
    results::MBRRecoverySearchResults;
    alpha::Float32 = 0.01f0,
    fdr::Float32 = 0.01f0,
)
    isempty(results.donor_pool) && return nothing
    psm_paths = list_post_scoring_psm_files(search_context)
    isempty(psm_paths) && return nothing
    precursors = getPrecursors(getSpecLib(search_context))
    return recovery_ftr_from_inputs(
        psm_paths, results.donor_pool, results.fragment_keys, precursors;
        alpha = alpha, fdr = fdr,
    )
end

"""
    recovery_ftr_from_inputs(psm_paths, donor_dict, fragment_keys, precursors; alpha, fdr)

Core recovery-FTR pass over pre-loaded inputs. Separated from
`run_recovery_ftr!` so it can be driven directly from preserved on-disk sidecars
(post-scoring PSM files + recovery seed sidecars) without re-running a search.
"""
function recovery_ftr_from_inputs(
    psm_paths::Vector{String},
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    fragment_keys::_MBRFragmentAnnotationKeys,
    precursors;
    alpha::Float32 = 0.01f0,
    fdr::Float32 = 0.01f0,
)
    partner_pools = build_counterfactual_partner_pools(psm_paths, precursors)
    df = _collect_recovery_ftr_rows(psm_paths, donor_dict, partner_pools, fragment_keys)
    n = nrow(df)
    if n == 0
        @user_info "Recovery FTR: no recovery seeds — skipping"
        return nothing
    end

    cand = .!df.MBR_is_missing_false
    n_cand = count(cand)
    if n_cand == 0
        @user_info "Recovery FTR: $n seeds but 0 have a counterfactual partner — skipping"
        return nothing
    end
    sub = df[cand, :]

    # Doubled (true/false) frame; label = is_real (real donor vs counterfactual).
    X = vcat(feature_matrix(sub, RECOVERY_FTR_FEATURES_TRUE),
             feature_matrix(sub, RECOVERY_FTR_FEATURES_FALSE))
    y = vcat(trues(n_cand), falses(n_cand))
    cv_double = vcat(sub.cv_fold, sub.cv_fold)

    pos0 = findall(cv_double .== UInt8(0))
    pos1 = findall(cv_double .== UInt8(1))
    ftr_double = fill(NaN32, 2 * n_cand)
    last_cls = nothing
    for (train_pos, test_pos) in ((pos1, pos0), (pos0, pos1))
        (isempty(train_pos) || isempty(test_pos)) && continue
        y_tr = y[train_pos]
        if length(unique(y_tr)) == 1
            ftr_double[test_pos] .= y_tr[1] ? 1f0 : 0f0
            continue
        end
        cls = build_lightgbm_classifier(; SHARED_LGBM_HP...)
        LightGBM.fit!(cls, X[train_pos, :], _prepare_labels(y_tr); verbosity = -1)
        raw = LightGBM.predict(cls, X[test_pos, :])
        s = ndims(raw) == 2 ? dropdims(raw; dims = 2) : raw
        ftr_double[test_pos] .= Float32.(s)
        last_cls = cls
    end

    # FTR pass: top-half (real) rows with PEP ≤ α (is_real labels on doubled frame).
    pep_double = Vector{Float32}(undef, 2 * n_cand)
    get_PEP!(ftr_double, y, pep_double)
    passing = pep_double[1:n_cand] .<= alpha
    n_passing = count(passing)
    tgt = Vector{Bool}(sub.target)
    n_pass_t = count(passing .& tgt)
    n_pass_d = n_passing - n_pass_t

    # <1% FDR: target/decoy q-value on the FTR is_real score over the recovery
    # rows (top half). Recovery rows carry no main score, so this is their only
    # FDR control.
    score_top = ftr_double[1:n_cand]
    qv = Vector{Float32}(undef, n_cand)
    get_qvalues!(score_top, tgt, qv; doSort = true, fdr_scale_factor = 1.0f0)
    fdr_pass = qv .<= fdr
    n_fdr_t = count(fdr_pass .& tgt)
    n_fdr_d = count(fdr_pass) - n_fdr_t

    importances = Tuple{Symbol, Float64}[]
    if last_cls !== nothing
        imp = importance(LightGBMModel(last_cls, RECOVERY_FTR_FEATURES_TRUE, nothing))
        imp !== nothing && (importances = sort([(Symbol(f), Float64(g)) for (f, g) in imp], by = x -> -x[2]))
    end

    emp_fdr = n_pass_t > 0 ? round(100 * n_pass_d / n_pass_t, digits=2) : 0.0
    @user_info "Recovery FTR — paired-counterfactual model:"
    if !isempty(importances)
        @user_info "  feature importances (gain):"
        for (fname, gain) in importances
            @user_info "    $(rpad(string(fname), 36)) $(round(gain, digits=2))"
        end
    end
    @user_info "  recovery candidates (seed + counterfactual): $n_cand / $n seeds"
    @user_info "  passing FTR (PEP ≤ $alpha): $n_passing  (targets=$n_pass_t decoys=$n_pass_d; " *
               "empirical target/decoy FDR=$emp_fdr%)"
    @user_info "  becoming PSMs at q ≤ $fdr (target/decoy on FTR score): $n_fdr_t targets (decoys=$n_fdr_d)"

    return (
        n_seeds = n, n_candidates = n_cand,
        n_passing = n_passing, n_pass_t = n_pass_t, n_pass_d = n_pass_d,
        n_fdr_t = n_fdr_t, n_fdr_d = n_fdr_d,
        importances = importances,
    )
end
