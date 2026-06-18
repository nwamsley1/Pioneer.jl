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
    # Per-rank fragment trace intensities (rank 1-8); sums matched isotope peaks
    # predicted at >=25% of the fragment's most abundant isotope.
    frag1_int::T
    frag2_int::T
    frag3_int::T
    frag4_int::T
    frag5_int::T
    frag6_int::T
    frag7_int::T
    frag8_int::T
    # E7 (Batch E, 2026-05-12): top-3 fragment ppm-error capture. Accumulates
    # |ppm_err| for M0 matches at ranks 1-3 only; mean computed in Score!.
    top3_abs_ppm_err_sum::T
    top3_ppm_err_count::UInt8
    # E3 (Batch E, 2026-05-12): sum of post-spline predicted intensity over
    # ALL matched M0 fragments (any rank). Paired with b_int + y_int to give
    # matched_ratio = log2((obs matched) / (pred matched)) in Score!.
    pred_int_sum_m0::T
    precursor_idx::UInt32
    ms_file_idx::UInt32
end

MainUnscoredPSM{Float32}() = MainUnscoredPSM{Float32}(
    UInt8(255), zero(UInt8), zero(UInt8), zero(UInt8),
    zero(UInt8), zero(UInt8), Float32(0),
    zero(UInt8), Float32(0), zero(UInt8),
    zero(UInt8), zero(UInt8), Float32(0),
    Float32(0), Float32(0), Float32(0), Float32(0),
    Float32(0), Float32(0), Float32(0), Float32(0),
    Float32(0), zero(UInt8), Float32(0),
    UInt32(0), UInt32(0),
)

"""
    TuningUnscoredPSM{T<:AbstractFloat} <: UnscoredPSM{T}

Slim per-(scan, precursor) accumulator for tuning paths
(ParameterTuningSearch, QuadTuningSearch, IntegrateChromatogramsSearch).
Same structural shape as `MainUnscoredPSM` but omits the MainSearch-only
fragment-chromatogram captures (`matched_rank_mask`, `frag1_int..frag8_int`)
since tuning code never reads them. Written by `apply_tuning_scoring!`.
"""
struct TuningUnscoredPSM{T<:AbstractFloat} <: UnscoredPSM{T}
    best_rank::UInt8
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

TuningUnscoredPSM{Float32}() = TuningUnscoredPSM(UInt8(255), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), zero(UInt8), Float32(0), zero(UInt8), Float32(0), zero(UInt8), zero(UInt8), zero(UInt8), Float32(0), UInt32(0), UInt32(0))
