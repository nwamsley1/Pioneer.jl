function selectTransitions(fragment_list::Vector{Vector{LibraryFragment{T}}}, pep_ids::Vector{UInt32}, ppm::AbstractFloat = 20.0) where {T<:AbstractFloat}
    transitions = Vector{LibraryFragment{T}}()
    for pep_id in pep_ids
        append!(transitions, fragment_list[pep_id])
    end
    return sort!(transitions, by = x->getFragMZ(x))
end

#=a = [1, 2, 3, 4, 5]
b = [10, 10, 10, 10, 1]
filter!(x->b[x]>1, @view(a[1:5]))
=#