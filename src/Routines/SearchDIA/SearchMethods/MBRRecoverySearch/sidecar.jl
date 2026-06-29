# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# Recovery seed sidecar I/O.
#
# Per-file Arrow file written next to main_search_psms/{file}.arrow with
# suffix RECOVERY_SEED_SIDECAR_SUFFIX (".recovery_seed.arrow"). Format is
# vertical-concat-friendly: one row per RecoveredSeedRow, columns flatten
# the NTuple{8} into eight scalar Float32 columns matching the existing
# `frag{1..8}_smoothed_intensity` schema in the main file.
#
# The sidecar is the source-of-truth for what V0 emits. V1.5 will read it
# back, pair each row with a counterfactual donor partner, train a paired
# FTR LightGBM, and gate downstream consumption (IntegrateChromatogramSearch
# + MaxLFQSearch) on the resulting `mbr_recovered_seed_qval`.

"""
    recovery_seed_sidecar_path(main_psm_path) -> String

Sidecar path for a given post-PrecursorScoringSearch main PSM file path.
"""
function recovery_seed_sidecar_path(main_psm_path::String)
    dir = dirname(main_psm_path)
    base = replace(basename(main_psm_path), r"\.arrow$" => "")
    return joinpath(dir, base * RECOVERY_SEED_SIDECAR_SUFFIX)
end

"""
    write_recovery_seed_sidecar(path, seeds::Vector{RecoveredSeedRow})

Materializes the seed vector as a DataFrame and writes an Arrow file.
No-op when `seeds` is empty (writes an empty file so PrecursorScoringSearch's
downstream merge can rely on the sidecar's existence).
"""
function write_recovery_seed_sidecar(path::String, seeds::Vector{RecoveredSeedRow})
    n = length(seeds)
    df = DataFrame(
        # Identity
        precursor_idx                     = Vector{UInt32}(undef, n),
        target                            = Vector{Bool}(undef, n),
        scan_idx                          = Vector{UInt32}(undef, n),
        cv_fold                           = Vector{UInt8}(undef, n),
        # Receiver raw features
        weight                            = Vector{Float32}(undef, n),
        log2_intensity_explained          = Vector{Float32}(undef, n),
        irt_obs                           = Vector{Float32}(undef, n),
        irt_pred                          = Vector{Float32}(undef, n),
        rt_obs                            = Vector{Float32}(undef, n),
        log_by_ratio                      = Vector{Float32}(undef, n),
        prec_mz                           = Vector{Float32}(undef, n),
        # Top-8 smoothed fragments (flattened from NTuple{8})
        frag1_smoothed_intensity          = Vector{Float32}(undef, n),
        frag2_smoothed_intensity          = Vector{Float32}(undef, n),
        frag3_smoothed_intensity          = Vector{Float32}(undef, n),
        frag4_smoothed_intensity          = Vector{Float32}(undef, n),
        frag5_smoothed_intensity          = Vector{Float32}(undef, n),
        frag6_smoothed_intensity          = Vector{Float32}(undef, n),
        frag7_smoothed_intensity          = Vector{Float32}(undef, n),
        frag8_smoothed_intensity          = Vector{Float32}(undef, n),
        # Precomputed MBR_*_true features
        MBR_max_pair_prob_true            = Vector{Float32}(undef, n),
        MBR_log2_weight_ratio_true        = Vector{Float32}(undef, n),
        MBR_log2_explained_ratio_true     = Vector{Float32}(undef, n),
        MBR_best_irt_diff_true            = Vector{Float32}(undef, n),
        MBR_best_rt_diff_true             = Vector{Float32}(undef, n),
        MBR_log_by_diff_true              = Vector{Float32}(undef, n),
        MBR_smoothed_frag_hellinger_true  = Vector{Float32}(undef, n),
        MBR_is_missing_true               = Vector{Bool}(undef, n),
        # Markers
        mbr_recovered_seed                = trues(n),
        trace_prob                        = zeros(Float32, n),
    )

    @inbounds for i in 1:n
        s = seeds[i]
        df.precursor_idx[i]                    = s.precursor_idx
        df.target[i]                           = s.target
        df.scan_idx[i]                         = s.scan_idx
        df.cv_fold[i]                          = s.cv_fold
        df.weight[i]                           = s.weight
        df.log2_intensity_explained[i]         = s.log2_intensity_explained
        df.irt_obs[i]                          = s.irt_obs
        df.irt_pred[i]                         = s.irt_pred
        df.rt_obs[i]                           = s.rt_obs
        df.log_by_ratio[i]                     = s.log_by_ratio
        df.prec_mz[i]                          = s.prec_mz
        df.frag1_smoothed_intensity[i]         = s.smoothed_frag_sqrt[1]
        df.frag2_smoothed_intensity[i]         = s.smoothed_frag_sqrt[2]
        df.frag3_smoothed_intensity[i]         = s.smoothed_frag_sqrt[3]
        df.frag4_smoothed_intensity[i]         = s.smoothed_frag_sqrt[4]
        df.frag5_smoothed_intensity[i]         = s.smoothed_frag_sqrt[5]
        df.frag6_smoothed_intensity[i]         = s.smoothed_frag_sqrt[6]
        df.frag7_smoothed_intensity[i]         = s.smoothed_frag_sqrt[7]
        df.frag8_smoothed_intensity[i]         = s.smoothed_frag_sqrt[8]
        df.MBR_max_pair_prob_true[i]           = s.MBR_max_pair_prob_true
        df.MBR_log2_weight_ratio_true[i]       = s.MBR_log2_weight_ratio_true
        df.MBR_log2_explained_ratio_true[i]    = s.MBR_log2_explained_ratio_true
        df.MBR_best_irt_diff_true[i]           = s.MBR_best_irt_diff_true
        df.MBR_best_rt_diff_true[i]            = s.MBR_best_rt_diff_true
        df.MBR_log_by_diff_true[i]             = s.MBR_log_by_diff_true
        df.MBR_smoothed_frag_hellinger_true[i] = s.MBR_smoothed_frag_hellinger_true
        df.MBR_is_missing_true[i]              = s.MBR_is_missing_true
    end

    Arrow.write(path, df)
    return path
end

"""
    read_recovery_seed_sidecar(path) -> DataFrame

Returns an empty DataFrame if the sidecar doesn't exist.
"""
function read_recovery_seed_sidecar(path::String)
    isfile(path) || return DataFrame()
    return DataFrame(Arrow.Table(path))
end
