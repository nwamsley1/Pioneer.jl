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
const MBR_SEMISUPERVISED_MAX_ITERATIONS =
    SCORING_SEMISUPERVISED_MAX_ITERATIONS
const MBR_RUN_CLUSTER_MIN_CHILD_SIZE = 2
const MBR_RUN_CLUSTER_SPLIT_ITERATIONS = 4
const MBR_RUN_CLUSTER_MIN_GAIN = 0.01f0
const MBR_LOD_WEIGHT_QUANTILE = 0.05f0
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
    fold_by_pid::Vector{UInt8}
    mz_bin_by_pid::Vector{UInt8}
    charge_by_pid::Vector{UInt8}
    length_by_pid::Vector{UInt8}
    irt_by_pid::Vector{Float32}
    fold_mz_charge_length_pool::Dict{
        Tuple{Int, Int, Int, Int},
        _MBRIrtPool,
    }
    fold_charge_length_pool::Dict{Tuple{Int, Int, Int}, _MBRIrtPool}
    charge_length_pool::Dict{Tuple{Int, Int}, _MBRIrtPool}
    file_charge_length_pool::Dict{Tuple{UInt32, Int, Int}, _MBRIrtPool}
end

struct _MBRIntensityPosting
    file_idx::UInt32
    log2_weight::Float32
end

struct _MBRRunIntensityProfile
    file_idx::UInt32
    precursor_columns::Vector{Int32}
    values::Vector{Float32}
end

struct _MBRSparseCentroid
    precursor_columns::Vector{Int32}
    values::Vector{Float32}
end

struct _MBRSphericalRunCluster
    members::Vector{Int}
    centroid::_MBRSparseCentroid
    similarity_sum::Float64
end

struct _MBRReceiverRunClusters
    cluster_by_file::Dict{UInt32, UInt32}
    cluster_sizes::Vector{UInt32}
    support_by_precursor::Dict{UInt32, Dict{UInt32, UInt32}}
    passed_by_file::Dict{UInt32, BitSet}
end

function _MBRReceiverRunClusters()
    return _MBRReceiverRunClusters(
        Dict{UInt32, UInt32}(),
        UInt32[],
        Dict{UInt32, Dict{UInt32, UInt32}}(),
        Dict{UInt32, BitSet}(),
    )
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
    irt_residual::Float32
    irt_obs::Float32
    n_scans::Float32
    integrated_frag_sqrt::NTuple{8, Float32}
    frag_corr_bitvec::UInt8
    frag_corr_bitvec_rank::UInt16
    ms_file_idx::UInt32
end

const MBR_RECEIVER_FEATURES = Symbol[
    :trace_prob_infold,
    :fitted_manhattan_distance,
    :irt_error,
    :poisson,
    :err_norm,
    :weight,
    :gof,
    :sequence_length,
    :fitted_hellinger,
    :weight_ratio_at_scan,
    :weight_rank_at_scan,
    :ms1_weight_apex_to_m0_apex_irt,
    :smoothed_2d_shadow_hellinger,
    :precursor_fraction_transmitted,
    :log_by_ratio_m0,
    :log2_intensity_explained,
]

const MBR_SHARED_FEATURES = Symbol[
    :MBR_log2_weight_lod_ratio,
    :MBR_receiver_frag_corr_bitvec_rank,
]

const MBR_PAIRED_FEATURE_STEMS = String[
    "MBR_best_pair_prob",
    "MBR_worst_pair_prob",
    "MBR_best_run_similarity",
    "MBR_receiver_cluster_support_count",
    "MBR_receiver_cluster_peer_count",
    "MBR_receiver_cluster_support_fraction",
    "MBR_best_log2_weight_ratio",
    "MBR_worst_log2_weight_ratio",
    "MBR_best_log2_explained_ratio",
    "MBR_worst_log2_explained_ratio",
    "MBR_best_abs_n_scans_diff",
    "MBR_worst_abs_n_scans_diff",
    "MBR_best_log2_n_scans_ratio",
    "MBR_worst_log2_n_scans_ratio",
    "MBR_best_irt_diff",
    "MBR_worst_irt_diff",
    "MBR_best_observed_irt_diff",
    "MBR_worst_observed_irt_diff",
    "MBR_single_donor",
    "MBR_best_hellinger_source_prob",
    "MBR_best_temporal_frag_hellinger",
    "MBR_best_temporal_corr_frag_hellinger",
    "MBR_best_temporal_receiver_corr_frag_hellinger",
    "MBR_best_temporal_shared_corr_frag_hellinger",
    "MBR_best_donor_frag_corr_bitvec_rank",
    "MBR_best_shared_corr_frag_bitvec_rank",
]

const MBR_HELLINGER_CONTRAST_BASE_STEMS = String[
    "MBR_best_temporal_frag_hellinger",
    "MBR_best_temporal_corr_frag_hellinger",
    "MBR_best_temporal_receiver_corr_frag_hellinger",
    "MBR_best_temporal_shared_corr_frag_hellinger",
]

const MBR_CONTRAST_FEATURE_STEMS = String[
    feature
    for base in MBR_HELLINGER_CONTRAST_BASE_STEMS
    for feature in (base * "_rank", base * "_margin")
]

const MBR_MODEL_PAIRED_FEATURE_STEMS = vcat(
    MBR_PAIRED_FEATURE_STEMS,
    MBR_CONTRAST_FEATURE_STEMS,
)

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
    MBR_SHARED_FEATURES...,
    (_mbr_true_feature(stem) for stem in MBR_MODEL_PAIRED_FEATURE_STEMS)...,
]

function _mbr_ftr_features_false(counterfactual_idx::Int)
    return Symbol[
        MBR_RECEIVER_FEATURES...,
        MBR_SHARED_FEATURES...,
        (_mbr_false_feature(stem, counterfactual_idx)
         for stem in MBR_MODEL_PAIRED_FEATURE_STEMS)...,
    ]
end
