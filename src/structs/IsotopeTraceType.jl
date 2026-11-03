# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

abstract type IsotopeTraceType end

struct CombineTraces <: IsotopeTraceType
    min_ft::Float32
end

seperateTraces(itt::CombineTraces) = false

function getPsmGroupbyCols(itt::CombineTraces)
    return [:precursor_idx]
end

struct SeperateTraces <: IsotopeTraceType
end
 
seperateTraces(itt::SeperateTraces) = true

function getPsmGroupbyCols(itt::SeperateTraces)
    return [:precursor_idx,:isotopes_captured]
end




