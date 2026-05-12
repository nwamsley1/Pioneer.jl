# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

abstract type UnscoredPSM{T<:AbstractFloat} <: PSM end

"""
    MainUnscoredPSM{T<:AbstractFloat} <: UnscoredPSM{T}

Per-(scan, precursor) accumulator updated inline during the fused scan
pass via `apply_main_scoring!`. After the scan completes, `Score!`
reads these into a `MainSearchScoredPSM` row.
"""
struct MainUnscoredPSM{T<:AbstractFloat} <: UnscoredPSM{T}
    best_rank::UInt8 #Highest ranking predicted framgent that was observed
    best_rank_iso::UInt8
    topn::UInt8 #How many of the topN predicted fragments were observed.
    topn_iso::UInt8
    longest_y::UInt8
    longest_b::UInt8
    longest_y_iso::UInt8
    isotope_count::UInt8
    b_count::UInt8
    b_int::T
    y_count::UInt8
    y_int::T
    y_count_iso::UInt8
    p_count::UInt8
    non_cannonical_count::UInt8
    error::T
    matched_rank_mask::UInt8   # bit i = matched M0 fragment of rank (i+1); bits 0..7 → ranks 1..8
    # Per-rank M0 fragment intensities (rank 1-6); for fragment-correlation features
    frag1_int::T
    frag2_int::T
    frag3_int::T
    frag4_int::T
    frag5_int::T
    frag6_int::T
    # Per-rank M+1 (first isotopologue) fragment intensities; captures only iso_idx==1.
    # Added 2026-05-12 to enable per-fragment M0/M+1 correlation features.
    frag1_int_m1::T
    frag2_int_m1::T
    frag3_int_m1::T
    frag4_int_m1::T
    frag5_int_m1::T
    frag6_int_m1::T
    # E7 (Batch E, 2026-05-12): top-3 fragment ppm-error capture. Accumulates
    # |ppm_err| for M0 matches at ranks 1-3 only; mean computed in Score!.
    top3_abs_ppm_err_sum::T
    top3_ppm_err_count::UInt8
    # E1/E2 (Batch E, 2026-05-12): per-rank post-NCE/iso-spline predicted
    # intensity captures (M0 only). Captures the max across this precursor's
    # scans — predicted intensity varies with quad transmission per scan; max
    # gives a stable per-fragment value for spectral_contrast / scribe.
    frag1_pred::T
    frag2_pred::T
    frag3_pred::T
    frag4_pred::T
    frag5_pred::T
    frag6_pred::T
    precursor_idx::UInt32
    ms_file_idx::UInt32
end

MainUnscoredPSM{Float32}() = MainUnscoredPSM(UInt8(255), UInt8(255), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), Float32(0), zero(UInt8), Float32(0), zero(UInt8), zero(UInt8), zero(UInt8), Float32(0), zero(UInt8), Float32(0), Float32(0), Float32(0), Float32(0), Float32(0), Float32(0), Float32(0), Float32(0), Float32(0), Float32(0), Float32(0), Float32(0), Float32(0), zero(UInt8), Float32(0), Float32(0), Float32(0), Float32(0), Float32(0), Float32(0), UInt32(0), UInt32(0))

"""
    TuningUnscoredPSM{T<:AbstractFloat} <: UnscoredPSM{T}

Slim per-(scan, precursor) accumulator for tuning paths
(ParameterTuningSearch, QuadTuningSearch, IntegrateChromatogramsSearch).
Same structural shape as `MainUnscoredPSM` but omits the MainSearch-only
fragment-chromatogram captures (`matched_rank_mask`, `frag1_int..frag6_int`)
since tuning code never reads them. Written by `apply_tuning_scoring!`.
"""
struct TuningUnscoredPSM{T<:AbstractFloat} <: UnscoredPSM{T}
    best_rank::UInt8
    best_rank_iso::UInt8
    topn::UInt8
    topn_iso::UInt8
    longest_y::UInt8
    longest_b::UInt8
    longest_y_iso::UInt8
    isotope_count::UInt8
    b_count::UInt8
    b_int::T
    y_count::UInt8
    y_int::T
    y_count_iso::UInt8
    p_count::UInt8
    non_cannonical_count::UInt8
    error::T
    precursor_idx::UInt32
    ms_file_idx::UInt32
end

TuningUnscoredPSM{Float32}() = TuningUnscoredPSM(UInt8(255), UInt8(255), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), Float32(0), zero(UInt8), Float32(0), zero(UInt8), zero(UInt8), zero(UInt8), Float32(0), UInt32(0), UInt32(0))
