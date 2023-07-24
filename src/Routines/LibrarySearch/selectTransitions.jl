function selectTransitions(fragment_list::Vector{Vector{LibraryFragment{T}}}, counter::Counter{I,C}, topN::Int, ppm::AbstractFloat = 20.0) where {T<:AbstractFloat, I,C<:Unsigned}
    transitions = Vector{LibraryFragment{T}}()
    i = 1
    while i <= min(topN, counter.matches)
        append!(transitions, fragment_list[getID(counter, i)])
        i += 1
    end
    return sort!(transitions, by = x->getFragMZ(x))
end


#Get relevant framgents given a retention time and precursor mass using a retentionTimeIndex object
function selectTransitions(fragment_list::Vector{Vector{LibraryFragment{V}}}, rt_index::retentionTimeIndex{T, U}, rt::T, rt_tol::T, prec_mz::U, prec_tol::U) where {T,U,V<:AbstractFloat}
    
    #Get matching precursors within an RT bin 
    function addTransitions!(transitions::Vector{LibraryFragment{V}}, fragment_list::Vector{Vector{LibraryFragment{V}}}, precs::Vector{Tuple{UInt32, U}}, prec_mz::U, prec_tol::U)
        start = searchsortedfirst(precs, by = x->last(x), prec_mz - prec_tol) #First precursor in the isolation window
        stop = searchsortedlast(precs, by = x->last(x), prec_mz + prec_tol) #Last precursor in the isolation window
        for i in start:stop #Get transitions for each precursor
            append!(transitions, fragment_list[first(precs[i])])
        end
    end

    transitions = Vector{LibraryFragment{V}}() #Initialize list of transitions/fragment Ions

    i = 1
    rt_start = max(searchsortedfirst(rt_index.rt_bins, rt - rt_tol, lt=(r,x)->r.lb<x) - 1, 1) #First RT bin to search
    rt_stop = min(searchsortedlast(rt_index.rt_bins, rt + rt_tol, lt=(x, r)->r.ub>x) + 1, length(rt_index.rt_bins)) #Last RT bin to search 
    for i in rt_start:rt_stop #Add transitions
        addTransitions!(transitions, fragment_list, rt_index.rt_bins[i].prec, prec_mz, prec_tol)
    end

    return sort!(transitions, by = x->getFragMZ(x)) #Sort transitions by their fragment m/z. 
end
#=function selectTransitions(fragment_list::Vector{Vector{LibraryFragment{T}}}, pep_ids::Base.Generator, ppm::AbstractFloat = 20.0) where {T<:AbstractFloat}
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

function selectTransitions(fragment_list::Vector{Vector{LibraryFragment{T}}}, pep_ids::Vector{UInt32}, ppm::AbstractFloat = 20.0) where {T<:AbstractFloat}
    transitions = Vector{LibraryFragment{T}}()4
    for pep_id in pep_ids
        append!(transitions, fragment_list[pep_id])
    end
    return sort!(transitions, by = x->getFragMZ(x))
end=#

#=a = [1, 2, 3, 4, 5]
b = [10, 10, 10, 10, 1]
filter!(x->b[x]>1, @view(a[1:5]))
=#