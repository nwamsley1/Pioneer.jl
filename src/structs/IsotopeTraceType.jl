abstract type IsotopeTraceType end

struct CombineTraces <: IsotopeTraceType
    min_ft::Float32
end

separateTraces(itt::CombineTraces) = false

function getPsmGroupbyCols(itt::CombineTraces)
    return [:precursor_idx]
end

struct SeparateTraces <: IsotopeTraceType
end
 
separateTraces(itt::SeparateTraces) = true

function getPsmGroupbyCols(itt::SeparateTraces)
    return [:precursor_idx,:isotopes_captured]
end




