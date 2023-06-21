function selectTransitions(fragment_list::Vector{Vector{LibraryFragment{T}}}, pep_ids::Base.Iterators.Take{Indices{UInt32}}, ppm::AbstractFloat = 20.0) where {T<:AbstractFloat}
    transitions = Vector{LibraryFragment{T}}()
    for pep_id in pep_ids
        append!(transitions, fragment_list[pep_id])
    end
    return sort!(transitions, by = x->getFragMZ(x))
end