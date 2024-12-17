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




