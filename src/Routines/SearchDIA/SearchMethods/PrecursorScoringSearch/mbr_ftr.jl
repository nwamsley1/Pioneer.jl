# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# MBR false-transfer controller as a recovery mechanism scoped to the
# transfer-candidate cohort.
#
# The counterfactual frame contains each candidate once for the real transfer
# and once per counterfactual:
#   * top block: candidate features using the same-precursor donor
#   * counterfactual blocks: the same receiver with alternate donor features
# The FTR label is true for top-half same-precursor transfers and false for
# counterfactual transfers. Semi-supervised iterations use only the
# top-half rows passing the 3% counterfactual FTR threshold as positives for the
# next training round and for the next stopping evaluation; unselected top-half
# rows are omitted from the semi-supervised objective.
# The resulting OOF score drives FTR recovery.
#
# Architecture unchanged otherwise:
#   * The qval pipeline downstream uses :trace_prob (= :trace_prob_prepass,
#     non-MBR pass-1). Non-candidate rows can never be score-inflated by MBR.
#   * Recovered candidates get :mbr_recovered = true inside the post-qvalue
#     filtered set; downstream PEP and run-level q-values are recomputed
#     without a recovery bypass.

# Donor q-value threshold for candidate eligibility. Defines the prepass-score
# floor a donor must clear (= score-at-q≤MBR_DONOR_Q_THRESHOLD among targets).
# Tighter (0.001) → only very confident donors → fewer candidates, safer FTR.
# Looser (1.0) → any cross-run PSM is a valid donor → more candidates, broader recovery.
const MBR_DONOR_Q_THRESHOLD = Float32(0.01)
const MBR_SEMISUPERVISED_FTR_THRESHOLD = Float32(0.03)
const MBR_RECOVERY_ALPHA_DEFAULT = Float32(0.01)

function _mbr_env_flag(name::AbstractString, env = ENV)::Bool
    raw = lowercase(strip(get(env, name, "")))
    return raw in ("1", "true", "yes", "y", "on")
end

function _mbr_debug_dir_from_env(env = ENV)::Union{Nothing, String}
    raw = strip(get(env, "PIONEER_MBR_DEBUG_DIR", ""))
    isempty(raw) && return nothing
    return raw
end

function _mbr_keep_temporary_sidecars_from_env(env = ENV)::Bool
    return _mbr_env_flag("PIONEER_MBR_KEEP_SIDECARS", env)
end

function _mbr_use_top_counterfactual_qvalues_from_env(env = ENV)::Bool
    raw = lowercase(strip(get(env, "PIONEER_MBR_USE_TOP_COUNTERFACTUAL_QVALUES", "")))
    isempty(raw) && return true
    raw in ("1", "true", "yes", "y", "on") && return true
    raw in ("0", "false", "no", "n", "off") && return false
    error("PIONEER_MBR_USE_TOP_COUNTERFACTUAL_QVALUES must be a boolean flag, got \"$raw\"")
end

function _mbr_recovery_alpha_from_env(env = ENV)::Float32
    raw = get(env, "PIONEER_MBR_FTR_ALPHA", "")
    isempty(raw) && return MBR_RECOVERY_ALPHA_DEFAULT
    value = tryparse(Float32, raw)
    (value === nothing || !isfinite(value) || value <= 0.0f0 || value > 1.0f0) &&
        error("PIONEER_MBR_FTR_ALPHA must be a finite number in (0, 1], got \"$raw\"")
    return value
end

# ─────────────────────────────────────────────────────────────────────────────
# Legacy Phase 6 (single-frame) FTR controller removed 2026-05-18.
# The Batch F paired controller (apply_mbr_filter_paired! below) supplanted it;
# the legacy FTR_FEATURES list + _train_ftr_lgbm + apply_mbr_filter! were dead
# code (no callers). Helpers (_transfer_class_counts, _pct, _log_transfer_counts)
# were only used by the legacy path and were removed with it.
# ─────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
# MBR Batch F — paired-counterfactual FTR controller with q-value threshold
# ─────────────────────────────────────────────────────────────────────────────

# FTR features for the Batch F controller. Same as Phase 8f's FTR_FEATURES but
# - drops :trace_prob_mbr (Pass 2 no longer exists in Batch F)
# - uses the _true variants for the in-place MBR features (the _false copies
#   are swapped into the counterfactual blocks of the training frame).
# rtv3 (2026-05-13) + rtv4 (2026-05-14) — slim FTR feature set.
# The pre-rtv3 list had ~50 features. The ~40 non-MBR features were
# row-level (identical between top and counterfactual blocks of the
# training frame), so they couldn't directly discriminate real-vs-
# counterfactual. Empirically the FTR LGBM was using them for
# interaction modeling. Dropping them gained +1,015 IDs (rtv3 doc)
# and dropped Dennis FTR by 0.31 pp.
#
# rtv4 dropped `trace_prob_prepass` (OOF) too — the in-fold score
# absorbs ~88% of the combined OOF+infold LGBM gain. The trace_prob_*
# columns are still computed and used downstream (qval pipeline,
# MBR top-K aggregation) — they just aren't exposed to the FTR LGBM
# as features.
#
# Dead features dropped (gain ≈ 0 across all rtv3/rtv4 runs):
#   - MBR_best_is_missing_true (the -1 sentinel on every other MBR feature
#     already encodes missingness)
#   - MBR_frag_top1_match_true (was redundant with MBR_frag_rank_corr_true;
#     both subsequently dropped along with the other M0 frag-6-vector features)
const MBR_CROSS_RUN_FTR_FEATURE_DROPS = Set([
    :total_ions,
    :missed_cleavage,
    :y_count,
    :Mox,
    :n_frags_detected_union,
    :n_frags_detected_intersection,
    :frag1_smoothed_intensity,
    :frag2_smoothed_intensity,
    :frag3_smoothed_intensity,
    :frag4_smoothed_intensity,
    :frag5_smoothed_intensity,
    :frag6_smoothed_intensity,
    :frag7_smoothed_intensity,
    :frag8_smoothed_intensity,
    :n_scans_other_windows,
    :other_window_weight_corr,
    :other_window_apex_delta_irt,
    :smoothness,
    :n_contiguous_scans,
    :scan_prec_mz_n_precursors,
    :delta_frame_peak_center,
    :ms1_m1_intensity,
    :ms1_m1_to_m0_ratio,
    :ms1_m1_to_m0_pred,
    :ms1_m0_mass_err_ppm,
    :ms1_m0_peak_n_precursors,
    :ms1_m0_peak_frag_intensity_fraction,
    :flanking_ms1_m0_candidate_fraction,
    :flanking_frag_candidate_fraction,
    :flanking_ms1_frag_sum_corr,
    :flanking_frag_corr_mean,
    :flanking_frag_corr_strength,
    :flanking_frag_corr_effective_n,
    :flanking_frag_corr_best_m0,
    :flanking_signal_support,
    :spectrum_peak_count,
    :prec_mz,
    :longest_y,
    :top3_ms2_mass_error_mean,
    :ms1_m0_intensity,
    :ms1_isotope_dotp_m0_m1_m2,
    :ms1_m0_m1_m2_window_fraction,
    :ms1_m0_m1_m2_window_fraction_pc,
    :ms1_ms2_explained_delta,
    :ms1_ms2_explained_delta_pc,
    :n_scans,
    :irt_fwhm,
    :frag_apex_gt2x_flank_bitvec_rank,
    :frag_apex_dispersion_irt,
    :frag_corr_strength,
    :n_correlated_fragments,
    :n_correlated_fragments_bitvec_rank,
    :n_frags_detected_union_bitvec_rank,
    :n_frags_detected_intersection_bitvec_rank,
    :frag_corr_effective_n,
    :frag_corr_best_m0,
])

const MBR_CROSS_RUN_FTR_FEATURES = Symbol[
    feature for feature in ADVANCED_FEATURE_SET
    if !(feature in MBR_CROSS_RUN_FTR_FEATURE_DROPS)
]

const FTR_FEATURES_F_TRUE = Symbol[
    MBR_CROSS_RUN_FTR_FEATURES...,
    :MBR_best_pair_prob_true,
    :MBR_worst_pair_prob_true,
    :MBR_best_run_similarity_true,
    :MBR_log2_weight_lod_ratio,
    :MBR_best_log2_weight_ratio_true,
    :MBR_worst_log2_weight_ratio_true,
    :MBR_best_log2_explained_ratio_true,
    :MBR_worst_log2_explained_ratio_true,
    :MBR_best_abs_n_scans_diff_true,
    :MBR_worst_abs_n_scans_diff_true,
    :MBR_best_log2_n_scans_ratio_true,
    :MBR_worst_log2_n_scans_ratio_true,
    :MBR_best_irt_diff_true,
    :MBR_worst_irt_diff_true,
    :MBR_best_observed_irt_diff_true,
    :MBR_worst_observed_irt_diff_true,
    :MBR_single_donor_true,
    :MBR_best_hellinger_source_prob_true,
    :MBR_best_smoothed_frag_hellinger_true,
    :MBR_best_smoothed_frag_hellinger_rank_true,
    :MBR_best_smoothed_frag_hellinger_margin_true,
    :MBR_best_corr_frag_hellinger_true,
    :MBR_best_corr_frag_hellinger_rank_true,
    :MBR_best_corr_frag_hellinger_margin_true,
    :MBR_best_donor_frag_corr_bitvec_rank_true,
    :MBR_best_receiver_corr_frag_hellinger_true,
    :MBR_best_receiver_corr_frag_hellinger_rank_true,
    :MBR_best_receiver_corr_frag_hellinger_margin_true,
    :MBR_receiver_frag_corr_bitvec_rank,
    :MBR_best_shared_corr_frag_hellinger_true,
    :MBR_best_shared_corr_frag_hellinger_rank_true,
    :MBR_best_shared_corr_frag_hellinger_margin_true,
    :MBR_best_shared_corr_frag_bitvec_rank_true,
]

const MBR_HELLINGER_CONTRAST_FEATURES = Symbol[
    :MBR_best_smoothed_frag_hellinger_rank_true,
    :MBR_best_smoothed_frag_hellinger_margin_true,
    :MBR_best_corr_frag_hellinger_rank_true,
    :MBR_best_corr_frag_hellinger_margin_true,
    :MBR_best_receiver_corr_frag_hellinger_rank_true,
    :MBR_best_receiver_corr_frag_hellinger_margin_true,
    :MBR_best_shared_corr_frag_hellinger_rank_true,
    :MBR_best_shared_corr_frag_hellinger_margin_true,
]
# Dropped 2026-05-16:
#   :MBR_top_n_median_score_true / :MBR_top_n_irt_diff_true
#   A/B on Olsen 23-file + MTAC 6-file (top-K vs top-1 donor): no material ID,
#   PG, or EFDR-MAE change. K trimmed to 2 for same-file fallback only.
# Dropped 2026-05-17:
#   :MBR_frag_pattern_cosine_true / :_scribe_true / :MBR_frag_rank_corr_true
#   A/B on YeastMBR 40-file + HelaOnly 20-file: dropping the three M0 fragment
#   6-vector features lowered Dennis FTR 2.59 % → 2.40 % and raised PGs +0.16 %,
#   at a cost of −0.06 % precursor IDs / −2.1 % MBR recoveries (YeastMBR);
#   HelaOnly pure-FTR 0.252 % → 0.242 %, PGs +0.33 %, IDs −0.65 %. Net quality
#   win; donor entries no longer need frag1_int..frag8_int.
#
# Ablated 2026-07-03:
#   :MBR_best_rt_diff_true, :MBR_best_log_by_diff_true, worst Hellinger inputs,
#   and donor-library Hellinger inputs. Non-Hellinger best/worst precursor
#   summaries are kept, with :MBR_single_donor_* marking rows whose worst donor
#   columns are sentinel-valued because only one donor was available.

# Same features but with the MBR columns swapped to a counterfactual suffix.
# Each entry must mirror
# FTR_FEATURES_F_TRUE position-for-position so the LGBM sees the same
# semantic feature columns in the same order. Row-level features pass through
# the default `f` branch unchanged.
function _mbr_counterfactual_feature_name(f::Symbol, counterfactual_idx::Int)
    s = String(f)
    if startswith(s, "MBR_") && endswith(s, "_true")
        stem = s[1:(end - length("_true"))]
        suffix = counterfactual_idx == 1 ? "_false" : "_false$(counterfactual_idx)"
        return Symbol(stem * suffix)
    end
    return f
end

const FTR_FEATURES_F_FALSE_BY_COUNTERFACTUAL = ntuple(
    counterfactual_idx -> Symbol[
        _mbr_counterfactual_feature_name(f, counterfactual_idx)
        for f in FTR_FEATURES_F_TRUE
    ],
    MBR_MAX_COUNTERFACTUALS,
)

const FTR_FEATURES_F_FALSE = FTR_FEATURES_F_FALSE_BY_COUNTERFACTUAL[1]
const FTR_FEATURES_F_FALSE2 = FTR_FEATURES_F_FALSE_BY_COUNTERFACTUAL[2]
const FTR_FEATURES_F_FALSE3 = FTR_FEATURES_F_FALSE_BY_COUNTERFACTUAL[3]
const FTR_FEATURES_F_FALSE4 = FTR_FEATURES_F_FALSE_BY_COUNTERFACTUAL[4]

function _mbr_ftr_features_true_from_env(env = ENV)
    features = copy(FTR_FEATURES_F_TRUE)
    if _mbr_env_flag("PIONEER_MBR_DISABLE_HELLINGER_CONTRAST", env)
        filter!(f -> !(f in MBR_HELLINGER_CONTRAST_FEATURES), features)
    end
    return features
end

function _mbr_best_hellinger_col(
    counterfactual_idx::Int,
    base::AbstractString = "MBR_best_smoothed_frag_hellinger",
)
    true_col = Symbol(base * "_true")
    return counterfactual_idx == 0 ?
        true_col :
        _mbr_counterfactual_feature_name(true_col, counterfactual_idx)
end

function _mbr_best_hellinger_rank_col(
    counterfactual_idx::Int,
    base::AbstractString = "MBR_best_smoothed_frag_hellinger",
)
    true_col = Symbol(base * "_rank_true")
    return counterfactual_idx == 0 ?
        true_col :
        _mbr_counterfactual_feature_name(true_col, counterfactual_idx)
end

function _mbr_best_hellinger_margin_col(
    counterfactual_idx::Int,
    base::AbstractString = "MBR_best_smoothed_frag_hellinger",
)
    true_col = Symbol(base * "_margin_true")
    return counterfactual_idx == 0 ?
        true_col :
        _mbr_counterfactual_feature_name(true_col, counterfactual_idx)
end

function _mbr_add_best_hellinger_rank_features!(
    psms::DataFrame;
    n_counterfactuals::Int,
    base::AbstractString = "MBR_best_smoothed_frag_hellinger",
)
    (1 <= n_counterfactuals <= MBR_MAX_COUNTERFACTUALS) ||
        error("n_counterfactuals must be in 1:$(MBR_MAX_COUNTERFACTUALS), got $n_counterfactuals")

    n = nrow(psms)
    best_h_cols = Symbol[_mbr_best_hellinger_col(0, base)]
    for counterfactual_idx in 1:n_counterfactuals
        push!(best_h_cols, _mbr_best_hellinger_col(counterfactual_idx, base))
    end
    for col in best_h_cols
        hasproperty(psms, col) || error("Missing MBR Hellinger feature column $col")
    end

    @inbounds for counterfactual_idx in 0:n_counterfactuals
        source_col = _mbr_best_hellinger_col(counterfactual_idx, base)
        rank_col = _mbr_best_hellinger_rank_col(counterfactual_idx, base)
        source_values = psms[!, source_col]
        ranks = Vector{Float32}(undef, n)
        for i in 1:n
            source_value = Float32(source_values[i])
            if isfinite(source_value) && source_value >= 0.0f0
                rank = 1
                for col in best_h_cols
                    value = Float32(psms[i, col])
                    if isfinite(value) && value >= 0.0f0 && value < source_value
                        rank += 1
                    end
                end
                ranks[i] = Float32(rank)
            else
                ranks[i] = -1.0f0
            end
        end
        psms[!, rank_col] = ranks
    end

    return psms
end

function _mbr_add_best_hellinger_margin_features!(
    psms::DataFrame;
    n_counterfactuals::Int,
    base::AbstractString = "MBR_best_smoothed_frag_hellinger",
)
    (1 <= n_counterfactuals <= MBR_MAX_COUNTERFACTUALS) ||
        error("n_counterfactuals must be in 1:$(MBR_MAX_COUNTERFACTUALS), got $n_counterfactuals")

    n = nrow(psms)
    best_h_cols = Symbol[_mbr_best_hellinger_col(0, base)]
    for counterfactual_idx in 1:n_counterfactuals
        push!(best_h_cols, _mbr_best_hellinger_col(counterfactual_idx, base))
    end
    for col in best_h_cols
        hasproperty(psms, col) || error("Missing MBR Hellinger feature column $col")
    end

    @inbounds for counterfactual_idx in 0:n_counterfactuals
        source_col = _mbr_best_hellinger_col(counterfactual_idx, base)
        margin_col = _mbr_best_hellinger_margin_col(counterfactual_idx, base)
        source_values = psms[!, source_col]
        margins = Vector{Float32}(undef, n)
        for i in 1:n
            source_value = Float32(source_values[i])
            if isfinite(source_value) && source_value >= 0.0f0
                best_other = Inf32
                for (other_idx, col) in enumerate(best_h_cols)
                    other_idx == counterfactual_idx + 1 && continue
                    value = Float32(psms[i, col])
                    if isfinite(value) && value >= 0.0f0 && value < best_other
                        best_other = value
                    end
                end
                margins[i] = isfinite(best_other) ? best_other - source_value : -1.0f0
            else
                margins[i] = -1.0f0
            end
        end
        psms[!, margin_col] = margins
    end

    return psms
end

function _mbr_repeat_blocks(values, n_blocks::Int)
    return vcat((copy(values) for _ in 1:n_blocks)...)
end

function _mbr_write_ftr_debug_tables!(
    debug_dir::AbstractString,
    sub::DataFrame,
    available_true::Vector{Symbol},
    x_blocks::Vector{Matrix{Float32}},
    ftr_score_double::AbstractVector{<:Real},
    qvals_double::AbstractVector{<:Real},
    pep_double::AbstractVector{<:Real},
    eval_labels::AbstractVector{Bool},
    eval_mask::AbstractVector{Bool},
    recovered_in_cand::AbstractVector{Bool};
    n_counterfactuals::Int,
    alpha::Float32,
    q_thresh::Float32,
    prob_thresh::Float32,
    best_iter::Int,
    n_positive::Int,
    combined_error_qvals_top = nothing,
    combined_error_rates_top = nothing,
    combined_diagnostics = nothing,
)
    mkpath(debug_dir)
    n_cand = nrow(sub)
    n_blocks = 1 + n_counterfactuals
    n_frame = n_blocks * n_cand
    @assert length(x_blocks) == n_blocks
    @assert length(ftr_score_double) == n_frame
    @assert length(qvals_double) == n_frame
    @assert length(pep_double) == n_frame
    @assert length(eval_labels) == n_frame
    @assert length(eval_mask) == n_frame

    candidate_idx = vcat((UInt32.(1:n_cand) for _ in 1:n_blocks)...)
    block_idx = vcat((
        fill(UInt8(block), n_cand)
        for block in 0:n_counterfactuals
    )...)
    frame = DataFrame(
        mbr_candidate_idx = candidate_idx,
        mbr_counterfactual_idx = block_idx,
        mbr_is_true_block = block_idx .== UInt8(0),
    )

    id_cols = Symbol[
        :precursor_idx,
        :ms_file_idx,
        :scan_idx,
        :cv_fold,
        :target,
        :qval,
        :global_qval,
        :trace_prob_prepass,
        :trace_prob_infold,
    ]
    for col in id_cols
        hasproperty(sub, col) || continue
        frame[!, col] = _mbr_repeat_blocks(collect(sub[!, col]), n_blocks)
    end

    X = vcat(x_blocks...)
    @assert size(X, 1) == n_frame
    @assert size(X, 2) == length(available_true)
    for (feature_idx, feature_name) in enumerate(available_true)
        frame[!, feature_name] = X[:, feature_idx]
    end

    recovered_frame = vcat(Bool.(recovered_in_cand), falses(n_counterfactuals * n_cand))
    frame[!, :mbr_ftr_score] = Float32.(ftr_score_double)
    frame[!, :mbr_ftr_qval] = Float32.(qvals_double)
    frame[!, :mbr_ftr_pep] = Float32.(pep_double)
    frame[!, :mbr_eval_label] = Bool.(eval_labels)
    frame[!, :mbr_eval_mask] = Bool.(eval_mask)
    frame[!, :mbr_recovered_top] = recovered_frame
    writeArrow(joinpath(debug_dir, "mbr_ftr_frame.arrow"), frame)

    candidates = copy(sub)
    candidates[!, :mbr_candidate_idx] = UInt32.(1:n_cand)
    candidates[!, :mbr_ftr_score_true] = Float32.(ftr_score_double[1:n_cand])
    candidates[!, :mbr_ftr_qval_true_debug] = Float32.(qvals_double[1:n_cand])
    candidates[!, :mbr_ftr_pep_true_debug] = Float32.(pep_double[1:n_cand])
    if combined_error_qvals_top !== nothing
        @assert length(combined_error_qvals_top) == n_cand
        candidates[!, :mbr_total_error_qval_true_debug] = Float32.(combined_error_qvals_top)
    end
    if combined_error_rates_top !== nothing
        @assert length(combined_error_rates_top) == n_cand
        candidates[!, :mbr_total_error_rate_true_debug] = Float32.(combined_error_rates_top)
    end
    candidates[!, :mbr_recovered_debug] = Bool.(recovered_in_cand)
    writeArrow(joinpath(debug_dir, "mbr_ftr_candidates.arrow"), candidates)

    summary = DataFrame(
        n_candidates = Int[n_cand],
        n_counterfactuals = Int[n_counterfactuals],
        n_frame_rows = Int[n_frame],
        alpha = Float32[alpha],
        q_thresh = Float32[q_thresh],
        prob_thresh = Float32[prob_thresh],
        best_iter = Int[best_iter],
        n_positive = Int[n_positive],
        n_recovered = Int[count(recovered_in_cand)],
        features = [join(string.(available_true), ",")],
    )
    if combined_diagnostics !== nothing
        summary[!, :base_targets] = Int[combined_diagnostics.base_targets]
        summary[!, :base_decoys] = Int[combined_diagnostics.base_decoys]
        summary[!, :baseline_error_rate] = Float32[combined_diagnostics.baseline_error_rate]
        summary[!, :mbr_targets] = Int[combined_diagnostics.mbr_targets]
        summary[!, :mbr_decoys] = Int[combined_diagnostics.mbr_decoys]
        summary[!, :mbr_false_transfers] = Int[combined_diagnostics.mbr_false_transfers]
        summary[!, :internal_ftr_targets] = Int[combined_diagnostics.internal_ftr_targets]
        summary[!, :internal_ftr_errors] = Int[combined_diagnostics.internal_ftr_errors]
        summary[!, :internal_ftr_estimate] = Float32[combined_diagnostics.internal_ftr_estimate]
        summary[!, :total_errors] = Int[combined_diagnostics.total_errors]
        summary[!, :total_targets] = Int[combined_diagnostics.total_targets]
        summary[!, :combined_error_rate] = Float32[combined_diagnostics.combined_error_rate]
    end
    writeArrow(joinpath(debug_dir, "mbr_ftr_summary.arrow"), summary)

    open(joinpath(debug_dir, "mbr_ftr_features.txt"), "w") do io
        for feature_name in available_true
            println(io, feature_name)
        end
    end
    return nothing
end

function _mbr_n_counterfactuals_from_env(env = ENV)::Int
    raw = get(env, "PIONEER_MBR_N_COUNTERFACTUALS", "")
    isempty(raw) && return MBR_DEFAULT_N_COUNTERFACTUALS
    value = tryparse(Int, raw)
    (value === nothing || value < 1 || value > MBR_MAX_COUNTERFACTUALS) &&
        error("PIONEER_MBR_N_COUNTERFACTUALS must be an integer in 1:$(MBR_MAX_COUNTERFACTUALS), got \"$raw\"")
    return value
end

function _mbr_transfer_training_labels(
    positive_top::AbstractVector{Bool};
    n_counterfactuals::Int = 1,
)
    n = length(positive_top)
    labels = falses((1 + n_counterfactuals) * n)
    @inbounds for i in 1:n
        labels[i] = positive_top[i]
    end
    return labels
end

function _mbr_transfer_training_mask(
    positive_top::AbstractVector{Bool};
    negative_top::Union{Nothing, AbstractVector{Bool}} = nothing,
    n_counterfactuals::Int = 1,
    frame_mask::Union{Nothing, AbstractVector{Bool}} = nothing,
)
    n = length(positive_top)
    negative_top !== nothing && @assert length(negative_top) == n
    mask = trues((1 + n_counterfactuals) * n)
    @inbounds for i in 1:n
        mask[i] = positive_top[i] || (negative_top !== nothing && negative_top[i])
    end
    if frame_mask !== nothing
        @assert length(frame_mask) == length(mask)
        mask .&= Bool.(frame_mask)
    end
    return mask
end

function _mbr_target_positive_top(
    positive_top::AbstractVector{Bool},
    target_top::AbstractVector{Bool},
)
    @assert length(positive_top) == length(target_top)
    return Bool.(positive_top) .& Bool.(target_top)
end

function _mbr_counterfactual_missing_col(counterfactual_idx::Int)::Symbol
    return counterfactual_idx == 1 ?
        :MBR_best_is_missing_false :
        Symbol("MBR_best_is_missing_false$(counterfactual_idx)")
end

function _mbr_counterfactual_present_matrix(
    psms::DataFrame,
    row_idx::AbstractVector{<:Integer};
    n_counterfactuals::Int,
)
    (1 <= n_counterfactuals <= MBR_MAX_COUNTERFACTUALS) ||
        error("n_counterfactuals must be in 1:$(MBR_MAX_COUNTERFACTUALS), got $n_counterfactuals")
    n = length(row_idx)
    present = falses(n, n_counterfactuals)
    @inbounds for counterfactual_idx in 1:n_counterfactuals
        miss_col = _mbr_counterfactual_missing_col(counterfactual_idx)
        hasproperty(psms, miss_col) ||
            error("Missing MBR counterfactual sidecar column $miss_col")
        miss = Bool.(psms[!, miss_col])
        for (j, row) in enumerate(row_idx)
            present[j, counterfactual_idx] = !miss[row]
        end
    end
    return present
end

function _mbr_transfer_frame_mask(counterfactual_present::AbstractMatrix{Bool})
    n_cand, n_counterfactuals = size(counterfactual_present)
    mask = trues((1 + n_counterfactuals) * n_cand)
    @inbounds for counterfactual_idx in 1:n_counterfactuals
        offset = counterfactual_idx * n_cand
        for i in 1:n_cand
            mask[offset + i] = counterfactual_present[i, counterfactual_idx]
        end
    end
    return mask
end

function _mbr_transfer_positive_top(
    ftr_qvals_top::AbstractVector{<:Real};
    ftr_threshold::Float32 = MBR_SEMISUPERVISED_FTR_THRESHOLD,
)
    n = length(ftr_qvals_top)
    positive_top = falses(n)
    @inbounds for i in 1:n
        positive_top[i] = Float32(ftr_qvals_top[i]) <= ftr_threshold
    end
    return positive_top
end

function _mbr_initial_transfer_positive_top(
    receiver_scores::AbstractVector{<:Real},
    target_top::AbstractVector{Bool},
)
    n = length(receiver_scores)
    @assert length(target_top) == n
    return Bool.(target_top)
end

function _mbr_transfer_counterfactual_labels(
    n_cand::Int;
    n_counterfactuals::Int = 1,
)
    labels = falses((1 + n_counterfactuals) * n_cand)
    @inbounds for i in 1:n_cand
        labels[i] = true
    end
    return labels
end

function _mbr_top_counterfactual_eval_mask(
    scores_double::AbstractVector{<:Real},
    n_cand::Int;
    n_counterfactuals::Int = 1,
    eval_mask::Union{Nothing, AbstractVector{Bool}} = nothing,
)
    n_double = (1 + n_counterfactuals) * n_cand
    @assert length(scores_double) == n_double
    if eval_mask !== nothing
        @assert length(eval_mask) == n_double
    end

    source_mask = eval_mask === nothing ? trues(n_double) : Bool.(eval_mask)
    top_mask = falses(n_double)
    @inbounds for i in 1:n_cand
        top_mask[i] = source_mask[i]
        best_idx = 0
        best_score = -Inf32
        for counterfactual_idx in 1:n_counterfactuals
            idx = counterfactual_idx * n_cand + i
            source_mask[idx] || continue
            score = Float32(scores_double[idx])
            rank_score = isfinite(score) ? score : -Inf32
            if best_idx == 0 || rank_score > best_score
                best_idx = idx
                best_score = rank_score
            end
        end
        best_idx != 0 && (top_mask[best_idx] = true)
    end
    return top_mask
end

function _mbr_transfer_iteration_metrics(
    scores_double::AbstractVector{<:Real},
    n_cand::Int;
    ftr_threshold::Float32 = MBR_SEMISUPERVISED_FTR_THRESHOLD,
    n_counterfactuals::Int = 1,
    eval_labels::AbstractVector{Bool} = _mbr_transfer_counterfactual_labels(
        n_cand;
        n_counterfactuals = n_counterfactuals,
    ),
    eval_mask::Union{Nothing, AbstractVector{Bool}} = nothing,
    use_top_counterfactual_qvalues::Bool = _mbr_use_top_counterfactual_qvalues_from_env(),
    target_top::Union{Nothing, AbstractVector{Bool}} = nothing,
)
    n_double = (1 + n_counterfactuals) * n_cand
    @assert length(scores_double) == n_double
    @assert length(eval_labels) == n_double
    target_top !== nothing && @assert length(target_top) == n_cand
    if eval_mask !== nothing
        @assert length(eval_mask) == n_double
    end

    qvals_double = fill(Inf32, n_double)
    pep_double = fill(Inf32, n_double)
    qvalue_eval_mask = if use_top_counterfactual_qvalues
        _mbr_top_counterfactual_eval_mask(
            scores_double,
            n_cand;
            n_counterfactuals = n_counterfactuals,
            eval_mask = eval_mask,
        )
    else
        eval_mask === nothing ? trues(n_double) : Bool.(eval_mask)
    end
    eval_idx = findall(qvalue_eval_mask)
    if !isempty(eval_idx)
        scores_eval = Float32.(scores_double[eval_idx])
        labels_eval = eval_labels[eval_idx]
        qvals_eval = Vector{Float32}(undef, length(eval_idx))
        pep_eval = Vector{Float32}(undef, length(eval_idx))
        get_qvalues!(scores_eval, labels_eval, qvals_eval)
        get_PEP!(scores_eval, labels_eval, pep_eval)
        @inbounds for (j, i) in enumerate(eval_idx)
            qvals_double[i] = qvals_eval[j]
            pep_double[i] = pep_eval[j]
        end
    end

    qvals_top = qvals_double[1:n_cand]
    pep_top = pep_double[1:n_cand]
    positive_top = _mbr_transfer_positive_top(
        qvals_top,
        ftr_threshold = ftr_threshold,
    )
    if target_top !== nothing
        positive_top = _mbr_target_positive_top(positive_top, target_top)
    end

    return (
        qvals_double = qvals_double,
        pep_double = pep_double,
        qvals_top = qvals_top,
        pep_top = pep_top,
        qvalue_eval_mask = qvalue_eval_mask,
        positive_top = positive_top,
        n_positive = count(positive_top),
    )
end

"""
    apply_mbr_filter_paired!(psms; alpha=0.01f0, q_thresh=0.01f0) -> NamedTuple

MBR Batch F — paired-counterfactual FTR controller.

Algorithm (see BATCH_F_PLAN.md):
1. Candidate gate = (pre_qvals > q_thresh) & (MBR_best_pair_prob_true ≥
   prob_thresh) & true donor present. Missing counterfactual donor blocks are
   masked out of FTR training/evaluation rather than dropping the candidate.
2. Build a training frame with one true block and one block per counterfactual:
   - top half: candidate rows with FTR_FEATURES_F_TRUE values.
   - counterfactual blocks: same receivers with the nearest
     FTR_FEATURES_F_FALSE* values when available.
3. Train iterative 2-fold CV LightGBM. Each iteration uses the previous round's
   3% counterfactual-FTR top-half rows as positives, plus available
   counterfactual rows as negatives.
4. Compute q-values and PEPs on the counterfactual frame. Recover candidate rows whose
   top-half row has q-value ≤ alpha.

Sets `:mbr_recovered`, `:MBR_transfer_candidate`, `:mbr_target_decoy_prob`,
`:ftr_qval_true` on `psms` in-place. Does NOT mutate `:trace_prob`.
"""
function _mbr_recovery_mask(qvals_top, pep_top, alpha::Float32)
    return Float32.(qvals_top) .<= alpha
end

@inline _mbr_rank_score(x) = isfinite(Float32(x)) ? Float32(x) : -Inf32
@inline _mbr_internal_ftr_estimate(errors::Integer, targets::Integer) =
    targets > 0 ? Float32(errors / targets) : NaN32

function _mbr_top_counterfactual_scores(
    scores_double::AbstractVector{<:Real},
    n_cand::Int;
    n_counterfactuals::Int,
    eval_mask::Union{Nothing, AbstractVector{Bool}} = nothing,
)
    n_double = (1 + n_counterfactuals) * n_cand
    @assert length(scores_double) == n_double
    if eval_mask !== nothing
        @assert length(eval_mask) == n_double
    end

    source_mask = eval_mask === nothing ? trues(n_double) : Bool.(eval_mask)
    scores = fill(-Inf32, n_cand)
    @inbounds for i in 1:n_cand
        best_score = -Inf32
        for counterfactual_idx in 1:n_counterfactuals
            idx = counterfactual_idx * n_cand + i
            source_mask[idx] || continue
            score = _mbr_rank_score(scores_double[idx])
            score > best_score && (best_score = score)
        end
        scores[i] = best_score
    end
    return scores
end

function _mbr_combined_error_recovery(
    scores_top::AbstractVector{<:Real},
    target_top::AbstractVector{Bool},
    top_counterfactual_scores::AbstractVector{<:Real};
    base_targets::Integer,
    base_decoys::Integer,
    alpha::Float32,
)
    n = length(scores_top)
    @assert length(target_top) == n
    @assert length(top_counterfactual_scores) == n

    combined_error_qvals = fill(Inf32, n)
    combined_error_rates = fill(Inf32, n)
    recovered = falses(n)
    base_targets_i = Int(base_targets)
    base_decoys_i = Int(base_decoys)
    baseline_error_rate = base_targets_i > 0 ? Float32(base_decoys_i / base_targets_i) : Inf32
    if n == 0 || base_targets_i <= 0 || baseline_error_rate >= alpha
        return (
            recovered = recovered,
            combined_error_qvals = combined_error_qvals,
            combined_error_rates = combined_error_rates,
            threshold = Float32(Inf),
            n_recovered = 0,
            base_targets = base_targets_i,
            base_decoys = base_decoys_i,
            baseline_error_rate = baseline_error_rate,
            mbr_targets = 0,
            mbr_decoys = 0,
            mbr_false_transfers = 0,
            internal_ftr_targets = 0,
            internal_ftr_errors = 0,
            internal_ftr_estimate = NaN32,
            total_errors = base_decoys_i,
            total_targets = base_targets_i,
            combined_error_rate = baseline_error_rate,
        )
    end

    valid_idx = Int[]
    sizehint!(valid_idx, n)
    @inbounds for i in 1:n
        isfinite(Float32(scores_top[i])) && push!(valid_idx, i)
    end
    if isempty(valid_idx)
        return (
            recovered = recovered,
            combined_error_qvals = combined_error_qvals,
            combined_error_rates = combined_error_rates,
            threshold = Float32(Inf),
            n_recovered = 0,
            base_targets = base_targets_i,
            base_decoys = base_decoys_i,
            baseline_error_rate = baseline_error_rate,
            mbr_targets = 0,
            mbr_decoys = 0,
            mbr_false_transfers = 0,
            internal_ftr_targets = 0,
            internal_ftr_errors = 0,
            internal_ftr_estimate = NaN32,
            total_errors = base_decoys_i,
            total_targets = base_targets_i,
            combined_error_rate = baseline_error_rate,
        )
    end

    order = sort(valid_idx, by = i -> Float32(scores_top[i]), rev = true)
    false_event_scores = Float32[]
    sizehint!(false_event_scores, count(target_top))
    @inbounds for i in 1:n
        target_top[i] || continue
        score = Float32(scores_top[i])
        cf_score = Float32(top_counterfactual_scores[i])
        if isfinite(score) && isfinite(cf_score)
            push!(false_event_scores, min(score, cf_score))
        end
    end
    sort!(false_event_scores, rev = true)

    group_starts = Int[]
    group_ends = Int[]
    group_rates = Float32[]
    n_targets = 0
    n_decoys = 0
    false_ptr = 0
    pos = 1
    @inbounds while pos <= length(order)
        group_start = pos
        threshold = Float32(scores_top[order[pos]])
        while pos <= length(order) && Float32(scores_top[order[pos]]) == threshold
            idx = order[pos]
            if target_top[idx]
                n_targets += 1
            else
                n_decoys += 1
            end
            pos += 1
        end
        while false_ptr < length(false_event_scores) && false_event_scores[false_ptr + 1] >= threshold
            false_ptr += 1
        end
        total_targets = base_targets_i + n_targets
        total_errors = base_decoys_i + n_decoys + false_ptr
        rate = total_targets > 0 ? Float32(total_errors / total_targets) : Inf32
        push!(group_starts, group_start)
        push!(group_ends, pos - 1)
        push!(group_rates, rate)
    end

    running_min = Inf32
    @inbounds for group_idx in length(group_rates):-1:1
        running_min = min(running_min, group_rates[group_idx])
        for pos_idx in group_starts[group_idx]:group_ends[group_idx]
            idx = order[pos_idx]
            combined_error_qvals[idx] = running_min
            combined_error_rates[idx] = group_rates[group_idx]
        end
    end

    recovered .= combined_error_qvals .<= alpha
    n_recovered = count(recovered)
    if n_recovered == 0
        return (
            recovered = recovered,
            combined_error_qvals = combined_error_qvals,
            combined_error_rates = combined_error_rates,
            threshold = Float32(Inf),
            n_recovered = 0,
            base_targets = base_targets_i,
            base_decoys = base_decoys_i,
            baseline_error_rate = baseline_error_rate,
            mbr_targets = 0,
            mbr_decoys = 0,
            mbr_false_transfers = 0,
            internal_ftr_targets = 0,
            internal_ftr_errors = 0,
            internal_ftr_estimate = NaN32,
            total_errors = base_decoys_i,
            total_targets = base_targets_i,
            combined_error_rate = baseline_error_rate,
        )
    end

    threshold = minimum(Float32.(scores_top[recovered]))
    mbr_targets = 0
    mbr_decoys = 0
    mbr_false_transfers = 0
    @inbounds for i in 1:n
        if recovered[i]
            if target_top[i]
                mbr_targets += 1
            else
                mbr_decoys += 1
            end
        end
        if target_top[i]
            score = Float32(scores_top[i])
            cf_score = Float32(top_counterfactual_scores[i])
            if isfinite(score) && isfinite(cf_score) && score >= threshold && cf_score >= threshold
                mbr_false_transfers += 1
            end
        end
    end
    total_targets = base_targets_i + mbr_targets
    total_errors = base_decoys_i + mbr_decoys + mbr_false_transfers
    combined_error_rate = total_targets > 0 ? Float32(total_errors / total_targets) : Inf32
    internal_ftr_targets = mbr_targets
    internal_ftr_errors = mbr_false_transfers
    internal_ftr_estimate = _mbr_internal_ftr_estimate(internal_ftr_errors, internal_ftr_targets)

    return (
        recovered = recovered,
        combined_error_qvals = combined_error_qvals,
        combined_error_rates = combined_error_rates,
        threshold = threshold,
        n_recovered = n_recovered,
        base_targets = base_targets_i,
        base_decoys = base_decoys_i,
        baseline_error_rate = baseline_error_rate,
        mbr_targets = mbr_targets,
        mbr_decoys = mbr_decoys,
        mbr_false_transfers = mbr_false_transfers,
        internal_ftr_targets = internal_ftr_targets,
        internal_ftr_errors = internal_ftr_errors,
        internal_ftr_estimate = internal_ftr_estimate,
        total_errors = total_errors,
        total_targets = total_targets,
        combined_error_rate = combined_error_rate,
    )
end

function apply_mbr_filter_paired!(
    psms::DataFrame;
    alpha::Float32    = 0.01f0,
    q_thresh::Float32 = 0.01f0,
    prob_thresh_override::Union{Nothing, Float32} = nothing,
)
    t0 = time()
    n_counterfactuals = _mbr_n_counterfactuals_from_env()

    n = nrow(psms)
    if n == 0
        @debug_l1 "MBR Batch F — apply_mbr_filter_paired!: no PSMs to filter"
        psms[!, :mbr_recovered]          = falses(0)
        psms[!, :MBR_transfer_candidate] = falses(0)
        psms[!, :mbr_target_decoy_prob]  = Float32[]
        psms[!, :ftr_qval_true]          = Float32[]
        psms[!, :ftr_pep_true]           = Float32[]
        psms[!, :mbr_total_error_qval_true] = Float32[]
        psms[!, :mbr_total_error_rate_true] = Float32[]
        return (
            n_candidates=0,
            threshold=Float32(Inf),
            n_recovered=0,
            prob_thresh=Float32(Inf),
            base_targets=0,
            base_decoys=0,
            baseline_error_rate=NaN32,
            mbr_targets=0,
            mbr_decoys=0,
            mbr_false_transfers=0,
            internal_ftr_targets=0,
            internal_ftr_errors=0,
            internal_ftr_estimate=NaN32,
            total_errors=0,
            total_targets=0,
            combined_error_rate=NaN32,
            elapsed_s=0.0,
        )
    end

    # ── 1. Pre-MBR q-values from the non-MBR (pass-1) score ──
    pre_score  = Float32.(psms[!, :trace_prob_prepass])
    target_col = Bool.(psms[!, :target])
    pre_qvals = if hasproperty(psms, :qval)
        Float32.(psms[!, :qval])
    else
        qvals = Vector{Float32}(undef, n)
        get_qvalues!(pre_score, target_col, qvals)
        qvals
    end
    global_pass = if hasproperty(psms, :global_qval)
        Float32.(psms[!, :global_qval]) .<= q_thresh
    else
        trues(n)
    end

    donor_target_pass = (pre_qvals .<= MBR_DONOR_Q_THRESHOLD) .& target_col
    prob_thresh = if prob_thresh_override !== nothing
        prob_thresh_override
    elseif any(donor_target_pass)
        minimum(pre_score[donor_target_pass])
    else
        Float32(Inf)
    end

    # ── 2. Candidates: failed q-cut, have a true donor. Counterfactual rows
    # can be partial; missing counterfactual blocks are masked below.
    mbr_pp_t       = Float32.(psms[!, :MBR_best_pair_prob_true])
    mbr_miss_t     = Bool.(psms[!, :MBR_best_is_missing_true])
    candidate_mask = global_pass .&
                     (pre_qvals .> q_thresh) .&
                     (mbr_pp_t  .>= prob_thresh) .&
                     (.!mbr_miss_t)
    n_cand = count(candidate_mask)
    cand_idx = findall(candidate_mask)

    @debug_l1 "MBR Batch F — paired FTR recovery:"
    @debug_l1 "  MBR_DONOR_Q_THRESHOLD = $MBR_DONOR_Q_THRESHOLD; prob_thresh = $(round(prob_thresh, digits=4))"
    @debug_l1 "  candidates (pre-MBR q>$q_thresh + true donor; $(n_counterfactuals) counterfactuals): $n_cand / $n ($(round(100*n_cand/n, digits=2))%)"
    @debug_l1 "  α (q-value FTR budget): $alpha"

    if n_cand == 0
        psms[!, :mbr_recovered]          = falses(n)
        psms[!, :MBR_transfer_candidate] = candidate_mask
        psms[!, :mbr_target_decoy_prob]  = fill(NaN32, n)
        psms[!, :ftr_qval_true]          = fill(NaN32, n)
        psms[!, :ftr_pep_true]           = fill(NaN32, n)
        psms[!, :mbr_total_error_qval_true] = fill(NaN32, n)
        psms[!, :mbr_total_error_rate_true] = fill(NaN32, n)
        base_pass = global_pass .& (pre_qvals .<= q_thresh)
        base_targets = count(base_pass .& target_col)
        base_decoys = count(base_pass .& .!target_col)
        baseline_error_rate = base_targets > 0 ? Float32(base_decoys / base_targets) : Inf32
        return (
            n_candidates=0,
            threshold=Float32(Inf),
            n_recovered=0,
            prob_thresh=prob_thresh,
            base_targets=base_targets,
            base_decoys=base_decoys,
            baseline_error_rate=baseline_error_rate,
            mbr_targets=0,
            mbr_decoys=0,
            mbr_false_transfers=0,
            internal_ftr_targets=0,
            internal_ftr_errors=0,
            internal_ftr_estimate=NaN32,
            total_errors=base_decoys,
            total_targets=base_targets,
            combined_error_rate=baseline_error_rate,
            elapsed_s=time()-t0,
        )
    end

    # ── 3. Build the true + selected-counterfactual training frame ──
    # Top half: true-donor rows, positive under the counterfactual FTR label.
    # Counterfactual blocks: same receivers with counterfactual MBR values,
    # always negative.
    counterfactual_present = _mbr_counterfactual_present_matrix(
        psms,
        cand_idx;
        n_counterfactuals = n_counterfactuals,
    )
    available_counterfactual_counts = vec(sum(counterfactual_present; dims = 2))
    n_full_counterfactuals = count(==(n_counterfactuals), available_counterfactual_counts)
    n_zero_counterfactuals = count(==(0), available_counterfactual_counts)
    n_available_counterfactual_rows = sum(available_counterfactual_counts)
    @debug_l1 "  counterfactual coverage: available rows=$n_available_counterfactual_rows / $(n_cand * n_counterfactuals); " *
              "full=$n_full_counterfactuals; zero=$n_zero_counterfactuals"

    sub = psms[cand_idx, :]
    true_features_all = _mbr_ftr_features_true_from_env()
    if :MBR_best_smoothed_frag_hellinger_rank_true in true_features_all
        _mbr_add_best_hellinger_rank_features!(
            sub;
            n_counterfactuals = n_counterfactuals,
            base = "MBR_best_smoothed_frag_hellinger",
        )
    end
    if :MBR_best_smoothed_frag_hellinger_margin_true in true_features_all
        _mbr_add_best_hellinger_margin_features!(
            sub;
            n_counterfactuals = n_counterfactuals,
            base = "MBR_best_smoothed_frag_hellinger",
        )
    end
    if :MBR_best_corr_frag_hellinger_rank_true in true_features_all
        _mbr_add_best_hellinger_rank_features!(
            sub;
            n_counterfactuals = n_counterfactuals,
            base = "MBR_best_corr_frag_hellinger",
        )
    end
    if :MBR_best_corr_frag_hellinger_margin_true in true_features_all
        _mbr_add_best_hellinger_margin_features!(
            sub;
            n_counterfactuals = n_counterfactuals,
            base = "MBR_best_corr_frag_hellinger",
        )
    end
    if :MBR_best_receiver_corr_frag_hellinger_rank_true in true_features_all
        _mbr_add_best_hellinger_rank_features!(
            sub;
            n_counterfactuals = n_counterfactuals,
            base = "MBR_best_receiver_corr_frag_hellinger",
        )
    end
    if :MBR_best_receiver_corr_frag_hellinger_margin_true in true_features_all
        _mbr_add_best_hellinger_margin_features!(
            sub;
            n_counterfactuals = n_counterfactuals,
            base = "MBR_best_receiver_corr_frag_hellinger",
        )
    end
    if :MBR_best_shared_corr_frag_hellinger_rank_true in true_features_all
        _mbr_add_best_hellinger_rank_features!(
            sub;
            n_counterfactuals = n_counterfactuals,
            base = "MBR_best_shared_corr_frag_hellinger",
        )
    end
    if :MBR_best_shared_corr_frag_hellinger_margin_true in true_features_all
        _mbr_add_best_hellinger_margin_features!(
            sub;
            n_counterfactuals = n_counterfactuals,
            base = "MBR_best_shared_corr_frag_hellinger",
        )
    end
    available_true  = filter(f -> hasproperty(sub, f), true_features_all)
    # Position-aligned _false sets (only matters that the MBR cols flip; non-MBR
    # cols are identical between TRUE and FALSE lists, so available_true tells
    # us which non-MBR cols actually exist).
    feat_to_idx = Dict{Symbol, Int}()
    for (i, f) in enumerate(true_features_all); feat_to_idx[f] = i; end
    x_blocks = Matrix{Float32}[feature_matrix(sub, available_true)]
    for counterfactual_idx in 1:n_counterfactuals
        available_false = Symbol[]
        false_features = Symbol[
            _mbr_counterfactual_feature_name(f, counterfactual_idx)
            for f in true_features_all
        ]
        for f in available_true
            push!(available_false, false_features[feat_to_idx[f]])
        end
        push!(x_blocks, feature_matrix(sub, available_false))
    end
    X = vcat(x_blocks...)
    eval_labels = _mbr_transfer_counterfactual_labels(
        n_cand;
        n_counterfactuals = n_counterfactuals,
    )
    frame_eval_mask = _mbr_transfer_frame_mask(counterfactual_present)
    eval_mask = copy(frame_eval_mask)
    cv_double = vcat((sub[!, :cv_fold] for _ in 1:(1 + n_counterfactuals))...)
    target_top = Bool.(sub[!, :target])
    initial_positive_top = _mbr_initial_transfer_positive_top(
        Float32.(sub[!, :trace_prob_prepass]),
        target_top,
    )
    receiver_decoy_top = .!target_top

    # ── 4. Semi-supervised 2-fold CV LightGBM on the counterfactual frame ──
    pos0 = findall(cv_double .== 0)
    pos1 = findall(cv_double .== 1)
    training_labels = _mbr_transfer_training_labels(
        initial_positive_top;
        n_counterfactuals = n_counterfactuals,
    )
    training_mask = _mbr_transfer_training_mask(
        initial_positive_top;
        negative_top = receiver_decoy_top,
        n_counterfactuals = n_counterfactuals,
        frame_mask = frame_eval_mask,
    )
    best_state = nothing
    previous_positive_count = -1
    @debug_l1 "  initial MBR training positives: $(count(initial_positive_top)) / $n_cand " *
              "(target candidates; target-decoys trained as negatives)"

    for iter_idx in 1:SCORING_SEMISUPERVISED_MAX_ITERATIONS
        ftr_score_double = fill(NaN32, (1 + n_counterfactuals) * n_cand)
        last_cls_iter = nothing
        n_train_targets = 0
        n_train_decoys = 0

        for (train_pos_all, test_pos) in ((pos1, pos0), (pos0, pos1))
            train_pos = train_pos_all[training_mask[train_pos_all]]
            (isempty(train_pos) || isempty(test_pos)) && continue
            y_tr = training_labels[train_pos]
            n_targets_fold = count(y_tr)
            n_train_targets += n_targets_fold
            n_train_decoys += length(y_tr) - n_targets_fold
            if length(unique(y_tr)) == 1
                ftr_score_double[test_pos] .= y_tr[1] ? 1f0 : 0f0
                continue
            end
            cls = build_lightgbm_classifier(; SHARED_LGBM_HP...)
            LightGBM.fit!(cls, X[train_pos, :], _prepare_labels(y_tr); verbosity = -1)
            raw = LightGBM.predict(cls, X[test_pos, :])
            s = ndims(raw) == 2 ? dropdims(raw; dims = 2) : raw
            ftr_score_double[test_pos] .= Float32.(s)
            last_cls_iter = cls
        end

        metrics = _mbr_transfer_iteration_metrics(
            ftr_score_double,
            n_cand,
            eval_labels = eval_labels,
            eval_mask = eval_mask,
            n_counterfactuals = n_counterfactuals,
            target_top = target_top,
        )
        state = (
            iter = iter_idx,
            scores = ftr_score_double,
            metrics = metrics,
            last_cls = last_cls_iter,
            n_positive = metrics.n_positive,
            n_train_targets = n_train_targets,
            n_train_decoys = n_train_decoys,
            eval_labels = copy(eval_labels),
            eval_mask = copy(eval_mask),
        )
        if best_state === nothing || state.n_positive >= best_state.n_positive
            best_state = state
        end

        @debug_l1 "  MBR transfer semi-supervised iter $iter_idx: " *
                   "train positives=$n_train_targets negatives=$n_train_decoys; " *
                   "FTR≤$(round(100 * MBR_SEMISUPERVISED_FTR_THRESHOLD, digits=2))% " *
                   "positives=$(state.n_positive)"

        if iter_idx > 1 && !_scoring_target_gain_sufficient(
            previous_positive_count,
            state.n_positive,
        )
            @debug_l1 "  MBR transfer semi-supervised stopping: " *
                       "iter $iter_idx positives=$(state.n_positive) did not improve by " *
                       "$(round(100 * SCORING_SEMISUPERVISED_MIN_TARGET_GAIN, digits=2))% " *
                       "over $previous_positive_count; using iter $(best_state.iter) " *
                       "with positives=$(best_state.n_positive)"
            break
        elseif iter_idx == SCORING_SEMISUPERVISED_MAX_ITERATIONS
            @debug_l1 "  MBR transfer semi-supervised stopping: hit max iterations " *
                       "$SCORING_SEMISUPERVISED_MAX_ITERATIONS; using iter $(best_state.iter) " *
                       "with positives=$(best_state.n_positive)"
            break
        end

        previous_positive_count = state.n_positive
        valid_transfer_labels = _mbr_transfer_training_labels(
            metrics.positive_top;
            n_counterfactuals = n_counterfactuals,
        )
        valid_transfer_mask = _mbr_transfer_training_mask(
            metrics.positive_top;
            negative_top = receiver_decoy_top,
            n_counterfactuals = n_counterfactuals,
            frame_mask = frame_eval_mask,
        )
        training_labels = valid_transfer_labels
        training_mask = valid_transfer_mask
        eval_labels = valid_transfer_labels
        eval_mask = valid_transfer_mask
    end

    # ── 5. Q-value AND PEP on the counterfactual frame ──
    ftr_score_double = best_state.scores
    qvals_double = best_state.metrics.qvals_double
    pep_double = best_state.metrics.pep_double
    last_cls = best_state.last_cls
    # qvals_double[i] = (# negatives ranked ≥ row i) / (# positives ranked ≥ row i),
    # monotonized to be non-increasing as score increases.
    # pep_double[i]   = per-row P(negative | score), isotonic-regression-derived.

    # ── 6. Recovery: top-half rows (real-MBR) with q-value ≤ alpha.
    # PEP is still retained as a diagnostic/reporting column, but the recovery
    # decision uses the paired FTR q-value.
    qvals_top = qvals_double[1:n_cand]
    pep_top   = pep_double[1:n_cand]
    ftr_score_top = ftr_score_double[1:n_cand]
    base_pass = global_pass .& (pre_qvals .<= q_thresh)
    base_targets = count(base_pass .& target_col)
    base_decoys = count(base_pass .& .!target_col)
    combined_score_top = copy(ftr_score_top)
    eval_top = Bool.(best_state.metrics.qvalue_eval_mask[1:n_cand])
    @inbounds for i in 1:n_cand
        eval_top[i] || (combined_score_top[i] = NaN32)
    end
    top_counterfactual_scores = _mbr_top_counterfactual_scores(
        ftr_score_double,
        n_cand;
        n_counterfactuals = n_counterfactuals,
        eval_mask = best_state.metrics.qvalue_eval_mask,
    )
    combined = _mbr_combined_error_recovery(
        combined_score_top,
        Bool.(sub[!, :target]),
        top_counterfactual_scores;
        base_targets = base_targets,
        base_decoys = base_decoys,
        alpha = alpha,
    )
    recovered_in_cand = combined.recovered
    n_recovered = count(recovered_in_cand)

    debug_dir = _mbr_debug_dir_from_env()
    if debug_dir !== nothing
        _mbr_write_ftr_debug_tables!(
            debug_dir,
            sub,
            available_true,
            x_blocks,
            ftr_score_double,
            qvals_double,
            pep_double,
            best_state.eval_labels,
            best_state.metrics.qvalue_eval_mask,
            recovered_in_cand;
            n_counterfactuals = n_counterfactuals,
            alpha = alpha,
            q_thresh = q_thresh,
            prob_thresh = prob_thresh,
            best_iter = best_state.iter,
            n_positive = best_state.n_positive,
            combined_error_qvals_top = combined.combined_error_qvals,
            combined_error_rates_top = combined.combined_error_rates,
            combined_diagnostics = combined,
        )
        @debug_l1 "  Wrote MBR FTR debug tables to $debug_dir"
    end

    mbr_recovered_full = falses(n)
    target_decoy_prob_full = fill(NaN32, n)
    ftr_qval_full      = fill(NaN32, n)
    ftr_pep_full       = fill(NaN32, n)
    combined_qval_full = fill(NaN32, n)
    combined_rate_full = fill(NaN32, n)
    @inbounds for (k, i) in enumerate(cand_idx)
        ftr_qval_full[i] = qvals_top[k]
        ftr_pep_full[i]  = pep_top[k]
        combined_qval_full[i] = combined.combined_error_qvals[k]
        combined_rate_full[i] = combined.combined_error_rates[k]
        if recovered_in_cand[k]
            mbr_recovered_full[i] = true
            target_decoy_prob_full[i] = ftr_score_top[k]
        end
    end
    psms[!, :mbr_recovered]          = mbr_recovered_full
    psms[!, :MBR_transfer_candidate] = candidate_mask
    psms[!, :mbr_target_decoy_prob]  = target_decoy_prob_full
    psms[!, :ftr_qval_true]          = ftr_qval_full
    psms[!, :ftr_pep_true]           = ftr_pep_full
    psms[!, :mbr_total_error_qval_true] = combined_qval_full
    psms[!, :mbr_total_error_rate_true] = combined_rate_full

    # ── 7. Threshold for log: highest ftr_score on a top-half row within the
    # combined target/decoy + counterfactual error budget.
    τ = combined.threshold

    # ── 8. Recovered transfer composition ──
    # Use library label + MBR_is_best_decoy-equivalent — donor polarity from
    # the true-donor's target/decoy class. The dual-features file does not
    # currently log donor target/decoy; infer from the precursor itself (T←T
    # for targets, D←D for decoys, since the donor is the SAME precursor).
    # T←D / D←T can only happen via the counterfactual, which is never used
    # for recovery; here counts are simpler: T or D recoveries.
    n_t_rec = count(i -> mbr_recovered_full[i] && target_col[i], 1:n)
    n_d_rec = n_recovered - n_t_rec

    @debug_l1 "  FTR frame rows: $((1 + n_counterfactuals) * n_cand)  (true donors=1, counterfactuals=$(n_counterfactuals))"
    @debug_l1 "  τ (combined-error score, total error ≤ α): $(round(τ, digits=4))"
    @debug_l1 "  baseline errors: decoys=$(combined.base_decoys) / targets=$(combined.base_targets) " *
              "($(round(100 * combined.baseline_error_rate, digits=4))%)"
    @debug_l1 "  combined errors: total=$(combined.total_errors) / targets=$(combined.total_targets) " *
              "($(round(100 * combined.combined_error_rate, digits=4))%); " *
              "MBR target-decoys=$(combined.mbr_decoys), MBR false-transfers=$(combined.mbr_false_transfers)"
    @debug_l1 "  internal FTR estimate for final accepted MBR target candidates: " *
              "$(combined.internal_ftr_errors)/$(combined.internal_ftr_targets) " *
              "($(round(100 * combined.internal_ftr_estimate, digits=4))%)"
    @debug_l1 "  internal FTR-q≤$(alpha) would recover $(count(_mbr_recovery_mask(qvals_top, pep_top, alpha))) candidates"
    @debug_l1 "  RECOVERED (combined total error ≤ α): $n_recovered ($(round(100*n_recovered/max(n_cand,1), digits=2))% of candidates)"
    @debug_l1 "  recovered targets: $n_t_rec   recovered decoys: $n_d_rec"

    if last_cls !== nothing
        lgbm_model = LightGBMModel(last_cls, available_true, nothing)
        imp = importance(lgbm_model)
        if imp !== nothing
            sorted_imp = sort(imp, by = x -> -x[2])
            @debug_l1 "  Batch F FTR-model feature importances (gain):"
            for (fname, gain) in sorted_imp
                @debug_l1 "    $(rpad(string(fname), 36)) $(round(gain, digits=2))"
            end
        end
    end

    elapsed = time() - t0
    @debug_l1 "  apply_mbr_filter_paired! elapsed: $(round(elapsed, digits=2))s"

    return (
        n_candidates = n_cand,
        threshold    = τ,
        n_recovered  = n_recovered,
        prob_thresh  = prob_thresh,
        base_targets = combined.base_targets,
        base_decoys = combined.base_decoys,
        baseline_error_rate = combined.baseline_error_rate,
        mbr_targets = combined.mbr_targets,
        mbr_decoys = combined.mbr_decoys,
        mbr_false_transfers = combined.mbr_false_transfers,
        internal_ftr_targets = combined.internal_ftr_targets,
        internal_ftr_errors = combined.internal_ftr_errors,
        internal_ftr_estimate = combined.internal_ftr_estimate,
        total_errors = combined.total_errors,
        total_targets = combined.total_targets,
        combined_error_rate = combined.combined_error_rate,
        elapsed_s    = elapsed,
    )
end
