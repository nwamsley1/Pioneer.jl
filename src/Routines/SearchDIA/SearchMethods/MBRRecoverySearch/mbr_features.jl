# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# Inline MBR_*_true feature compute for recovered seed rows.
#
# Formulas mirror PrecursorScoringSearch/mbr_streaming.jl's _compute_mbr_inner!
# 1:1 (see lines 530-585 there). The only difference is that our receiver +
# donor share the same precursor_idx — we're "recovering" the donor's pid in
# a different file — which lets us skip the donor-dict lookup and the
# counterfactual partner pairing at compute time.
#
# Sentinel = -1.0f0 (matches the existing MBR feature sentinel convention).
#   Used when the receiver's PMM weight is non-positive, or when fragment
#   intensities are all zero (Hellinger uses 1.0f0 in that case to match the
#   existing _mbr_smoothed_spectrum_hellinger_from_sqrt early-return).

"""
    ReceiverExtraction

Per-cell PMM extraction output that mbr_features.jl consumes. Mirrors the
schema of `_MBRDonorEntry`'s payload (without trace_prob, ms_file_idx, is_decoy)
so the inline MBR feature compute is symmetric.
"""
struct ReceiverExtraction
    scan_idx::UInt32
    weight::Float32
    log2_intensity_explained::Float32
    irt_obs::Float32                          # mapped through rt_to_irt at chosen scan
    rt_obs::Float32                           # raw RT (min) at chosen scan
    log_by_ratio::Float32
    smoothed_frag_sqrt::NTuple{8, Float32}
end

"""
    compute_seed_mbr_features(receiver, donor, irt_pred, fragment_keys) -> NamedTuple

Computes the seven MBR_*_true features + MBR_is_missing_true (= false) for a
recovered seed. Receiver and donor share the same precursor_idx, so the
Hellinger fragment-key lookup gets the same pid both sides.

`irt_pred` is the library iRT for the shared pid. It cancels in
`MBR_best_irt_diff_true` because (irt_pred − receiver.irt_obs) −
(irt_pred − donor.irt_obs) = donor.irt_obs − receiver.irt_obs, but we
keep the same formula shape to mirror the existing semantics exactly.
"""
function compute_seed_mbr_features(
    receiver::ReceiverExtraction,
    donor::_MBRDonorEntry,
    irt_pred::Float32,
    fragment_keys::_MBRFragmentAnnotationKeys,
)
    pid = donor.precursor_idx

    # MBR_max_pair_prob_true = donor's pass-1 LightGBM trace_prob
    max_pair_prob = donor.trace_prob

    # MBR_log2_weight_ratio_true = log2(receiver_weight / donor_weight),
    # sentinel if either is non-positive
    log2_weight_ratio = if donor.weight > 0.0f0 && receiver.weight > 0.0f0
        log2(receiver.weight / donor.weight)
    else
        -1.0f0
    end

    # MBR_log2_explained_ratio_true = receiver_log2ie - donor_log2ie
    # (already in log2 space, so it's a difference not a ratio)
    log2_explained_ratio = receiver.log2_intensity_explained - donor.log2_intensity_explained

    # MBR_best_irt_diff_true: receiver iRT residual minus donor iRT residual
    receiver_irt_residual = irt_pred - receiver.irt_obs
    best_irt_diff = abs(receiver_irt_residual - donor.irt_residual)

    # MBR_log_by_diff_true: difference in log(b)/log(y) ratio
    log_by_diff = receiver.log_by_ratio - donor.log_by_ratio

    # MBR_best_rt_diff_true: absolute RT-minutes difference
    best_rt_diff = abs(receiver.rt_obs - donor.rt_obs)

    # MBR_smoothed_frag_hellinger_true: Hellinger on sqrt-L1-normalized
    # top-8 fragment intensities, matching ranks via library fragment keys
    smoothed_frag_hellinger = _mbr_smoothed_spectrum_hellinger_from_sqrt(
        receiver.smoothed_frag_sqrt,
        donor.smoothed_frag_sqrt,
        fragment_keys,
        pid,         # recipient_pid
        pid,         # donor_pid (same)
    )

    return (
        MBR_max_pair_prob_true            = max_pair_prob,
        MBR_log2_weight_ratio_true        = log2_weight_ratio,
        MBR_log2_explained_ratio_true     = log2_explained_ratio,
        MBR_best_irt_diff_true            = best_irt_diff,
        MBR_best_rt_diff_true             = best_rt_diff,
        MBR_log_by_diff_true              = log_by_diff,
        MBR_smoothed_frag_hellinger_true  = smoothed_frag_hellinger,
        MBR_is_missing_true               = false,
    )
end
