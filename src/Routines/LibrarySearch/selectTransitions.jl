function selectTransitions(fragment_list::Vector{Vector{LibraryFragment{T}}}, pep_ids::Vector{UInt32}, ppm::AbstractFloat = 20.0) where {T<:AbstractFloat}
    transitions = Vector{LibraryFragment{T}}()
    for pep_id in pep_ids
        append!(transitions, fragment_list[pep_id])
    end
    return sort!(transitions, by = x->getFragMZ(x))
end

function selectTransitions(fragment_list::Vector{Vector{LibraryFragment{T}}}, pep_ids::Base.Generator, ppm::AbstractFloat = 20.0) where {T<:AbstractFloat}
    transitions = Vector{LibraryFragment{T}}()
    i = 1
    for pep_id::UInt32 in pep_ids
        #println(pep_id)
        append!(transitions, fragment_list[pep_id])
        i += 1
        if i > 20
            break
        end
    end
    return sort!(transitions, by = x->getFragMZ(x))
end

#=a = [1, 2, 3, 4, 5]
b = [10, 10, 10, 10, 1]
filter!(x->b[x]>1, @view(a[1:5]))
=#