# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

"""
Internal schemas and constants for post-integration match-between-runs
rescoring.

MBR is deliberately implemented outside either PrecursorScoringSearch or
IntegrateChromatogramsSearch: precursor scoring stages candidates, while
chromatogram integration supplies the receiver evidence and triggers final
rescoring.
"""

const MBR_DONOR_Q_THRESHOLD = 0.01f0
const MBR_N_COUNTERFACTUALS = 3
const MBR_SEMISUPERVISED_FTR_THRESHOLD = 0.03f0
const MBR_SEMISUPERVISED_MAX_ITERATIONS = 5
# Per OOF model. The natural training ratio is at most three
# counterfactual negatives per accepted positive transfer, and the combined
# 2.5M-row ceiling matches experiment-wide precursor scoring.
const MBR_MAX_POSITIVE_TRAIN_PER_FOLD = 625_000
const MBR_MAX_NEGATIVE_TRAIN_PER_FOLD = 1_875_000

const PASS1_SIDECAR_SUFFIX = ".pass1_sidecar.arrow"
const MBR_SIDECAR_SUFFIX = ".mbr_sidecar.arrow"
const RECOVERY_SIDECAR_SUFFIX = ".recovery_sidecar.arrow"

const MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT = ntuple(_ -> 0.0f0, 8)
const MBR_INTEGRATED_FRAGMENT_SQRT_COLUMNS = ntuple(
    rank -> Symbol("MBR_integrated_frag$(rank)_sqrt"),
    8,
)
const MBR_INTEGRATED_TEMPORAL_MEAN_SQRT_COLUMNS = ntuple(
    rank -> Symbol("MBR_integrated_temporal_mean_frag$(rank)_sqrt"),
    8,
)
const MBR_INTEGRATED_TEMPORAL_TRACE_COLUMN =
    :MBR_integrated_temporal_shadow_trace
const MBR_TEMPORAL_TRACE_STRIDE = 9
const MBR_INTEGRATED_APEX_IRT_COLUMN = :MBR_integrated_apex_irt_obs
const MBR_INTEGRATED_WEIGHT_COLUMN = :MBR_integrated_weight
const MBR_INTEGRATED_LOG2_INTENSITY_EXPLAINED_COLUMN =
    :MBR_integrated_log2_intensity_explained
const MBR_INTEGRATED_FITTED_MANHATTAN_DISTANCE_COLUMN =
    :MBR_integrated_fitted_manhattan_distance
const MBR_INTEGRATED_FITTED_HELLINGER_COLUMN =
    :MBR_integrated_fitted_hellinger
const MBR_INTEGRATED_SMOOTHED_2D_SHADOW_HELLINGER_COLUMN =
    :MBR_integrated_smoothed_2d_shadow_hellinger
const MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_COLUMN =
    :MBR_integrated_n_correlated_fragments
const MBR_INTEGRATED_FRAG_CORR_BITVEC_COLUMN =
    :MBR_integrated_frag_corr_bitvec
const MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_BITVEC_RANK_COLUMN =
    :MBR_integrated_n_correlated_fragments_bitvec_rank
const MBR_INTEGRATED_FRAG_CORR_STRENGTH_COLUMN =
    :MBR_integrated_frag_corr_strength
const MBR_INTEGRATED_FRAG_CORR_EFFECTIVE_N_COLUMN =
    :MBR_integrated_frag_corr_effective_n
const MBR_INTEGRATED_FRAG_CORR_BEST_WEIGHT_COLUMN =
    :MBR_integrated_frag_corr_best_weight
const MBR_INTEGRATED_N_SCANS_COLUMN = :points_integrated
const MBR_INTERNAL_INTEGRATED_COLUMNS = Symbol[
    MBR_INTEGRATED_FRAGMENT_SQRT_COLUMNS...,
    MBR_INTEGRATED_TEMPORAL_MEAN_SQRT_COLUMNS...,
    MBR_INTEGRATED_TEMPORAL_TRACE_COLUMN,
    MBR_INTEGRATED_APEX_IRT_COLUMN,
    MBR_INTEGRATED_WEIGHT_COLUMN,
    MBR_INTEGRATED_LOG2_INTENSITY_EXPLAINED_COLUMN,
    MBR_INTEGRATED_FITTED_MANHATTAN_DISTANCE_COLUMN,
    MBR_INTEGRATED_FITTED_HELLINGER_COLUMN,
    MBR_INTEGRATED_SMOOTHED_2D_SHADOW_HELLINGER_COLUMN,
    MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_COLUMN,
    MBR_INTEGRATED_FRAG_CORR_BITVEC_COLUMN,
    MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_BITVEC_RANK_COLUMN,
    MBR_INTEGRATED_FRAG_CORR_STRENGTH_COLUMN,
    MBR_INTEGRATED_FRAG_CORR_EFFECTIVE_N_COLUMN,
    MBR_INTEGRATED_FRAG_CORR_BEST_WEIGHT_COLUMN,
    :integration_start_scan,
    :integration_stop_scan,
]

const _MBRIrtPool = NamedTuple{
    (:pids, :irts),
    Tuple{Vector{UInt32}, Vector{Float32}},
}

_empty_mbr_irt_pool() = (pids = UInt32[], irts = Float32[])

struct _MBRCounterfactualEligibility
    global_passed::BitSet
    run_passed_by_file::Dict{UInt32, BitSet}
end

struct _MBRPartnerPools
    charge_by_pid::Vector{UInt8}
    length_by_pid::Vector{UInt8}
    irt_by_pid::Vector{Float32}
    charge_length_pool::Dict{Tuple{Int, Int}, _MBRIrtPool}
    file_charge_length_pool::Dict{Tuple{UInt32, Int, Int}, _MBRIrtPool}
end

"""
Compact donor evidence retained after chromatogram integration.

Only baseline rows that passed both the global and run-level q-value gates can
become donors. Counterfactuals reuse the same structure, which guarantees that
real and false transfers are compared through identical feature code.
"""
struct _MBRDonorEntry
    trace_prob::Float32
    precursor_idx::UInt32
    weight::Float32
    log2_intensity_explained::Float32
    irt_obs::Float32
    n_scans::Float32
    integrated_frag_sqrt::NTuple{8, Float32}
    frag_corr_bitvec::UInt8
    frag_corr_bitvec_rank::UInt16
    ms_file_idx::UInt32
end

const MBR_RECEIVER_FEATURES = Symbol[
    :trace_prob_infold,
    MBR_INTEGRATED_FITTED_MANHATTAN_DISTANCE_COLUMN,
    MBR_INTEGRATED_FITTED_HELLINGER_COLUMN,
    MBR_INTEGRATED_SMOOTHED_2D_SHADOW_HELLINGER_COLUMN,
    MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_COLUMN,
    MBR_INTEGRATED_FRAG_CORR_STRENGTH_COLUMN,
    MBR_INTEGRATED_FRAG_CORR_EFFECTIVE_N_COLUMN,
    MBR_INTEGRATED_FRAG_CORR_BEST_WEIGHT_COLUMN,
    MBR_INTEGRATED_N_SCANS_COLUMN,
]

const MBR_PAIRED_FEATURE_STEMS = String[
    "MBR_best_pair_prob",
    "MBR_best_run_similarity",
    "MBR_best_log2_weight_ratio",
    "MBR_best_log2_explained_ratio",
    "MBR_best_observed_irt_diff",
    "MBR_best_abs_n_scans_diff",
    "MBR_best_log2_n_scans_ratio",
    "MBR_best_temporal_frag_hellinger",
    "MBR_best_temporal_corr_frag_hellinger",
    "MBR_best_temporal_receiver_corr_frag_hellinger",
    "MBR_best_temporal_shared_corr_frag_hellinger",
    "MBR_best_donor_frag_corr_bitvec_rank",
    "MBR_best_shared_corr_frag_bitvec_rank",
]

@inline _mbr_true_feature(stem::AbstractString) = Symbol(stem * "_true")

@inline function _mbr_false_feature(
    stem::AbstractString,
    counterfactual_idx::Int,
)
    suffix = counterfactual_idx == 1 ? "_false" : "_false$(counterfactual_idx)"
    return Symbol(stem * suffix)
end

@inline function _mbr_missing_feature(counterfactual_idx::Int)
    return counterfactual_idx == 0 ?
        :MBR_best_is_missing_true :
        _mbr_false_feature("MBR_best_is_missing", counterfactual_idx)
end

const MBR_FTR_FEATURES_TRUE = Symbol[
    MBR_RECEIVER_FEATURES...,
    (_mbr_true_feature(stem) for stem in MBR_PAIRED_FEATURE_STEMS)...,
]

function _mbr_ftr_features_false(counterfactual_idx::Int)
    return Symbol[
        MBR_RECEIVER_FEATURES...,
        (_mbr_false_feature(stem, counterfactual_idx)
         for stem in MBR_PAIRED_FEATURE_STEMS)...,
    ]
end
