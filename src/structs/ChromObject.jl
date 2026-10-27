# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

abstract type ChromObject end
struct MS2ChromObject <: ChromObject
    rt::Float32
    intensity::Float32
    scan_idx::UInt32
    precursor_idx::UInt32
end

struct MS2MBRChromObject <: ChromObject
    rt::Float32
    intensity::Float32
    scan_idx::UInt32
    precursor_idx::UInt32
    shadow_frag1_int::Float32
    shadow_frag2_int::Float32
    shadow_frag3_int::Float32
    shadow_frag4_int::Float32
    shadow_frag5_int::Float32
    shadow_frag6_int::Float32
    shadow_frag7_int::Float32
    shadow_frag8_int::Float32
    fitted_frag1_int::Float32
    fitted_frag2_int::Float32
    fitted_frag3_int::Float32
    fitted_frag4_int::Float32
    fitted_frag5_int::Float32
    fitted_frag6_int::Float32
    fitted_frag7_int::Float32
    fitted_frag8_int::Float32
end
struct MS1ChromObject <: ChromObject
    rt::Float32
    intensity::Float32
    m0::Bool
    n_iso::UInt8
    scan_idx::UInt32
    precursor_idx::UInt32
end

function growChromObjects!(chromatograms::Vector{ChromObject}, block_size::Int64)
    chromatograms = append!(chromatograms, Vector{ChromObject}(undef, block_size))
end
