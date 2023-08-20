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
function selectTransitions!(transitions::Vector{LibraryFragment{V}}, prec_ids::Vector{UInt32}, fragment_list::Vector{Vector{LibraryFragment{V}}}, rt_index::retentionTimeIndex{T, U}, rt::T, rt_tol::T, prec_mz::U, prec_tol::U) where {T,U,V<:AbstractFloat}
    
    #Get matching precursors within an RT bin 
    function addTransitions!(transitions::Vector{LibraryFragment{V}}, prec_ids::Vector{UInt32}, transition_idx::Int64, prec_idx::Int64, fragment_list::Vector{Vector{LibraryFragment{V}}}, precs::Vector{Tuple{UInt32, U}}, prec_mz::U, prec_tol::U)
        start = searchsortedfirst(precs, by = x->last(x), prec_mz - prec_tol) #First precursor in the isolation window
        stop = searchsortedlast(precs, by = x->last(x), prec_mz + prec_tol) #Last precursor in the isolation window
        for i in start:stop #Get transitions for each precursor
            for frag in fragment_list[first(precs[i])]
                transitions[transition_idx] = frag
                transition_idx += 1 
            end
            prec_ids[prec_idx] = first(precs[i])
        end
        return transition_idx, prec_idx
    end

    #transitions = Vector{LibraryFragment{V}}() #Initialize list of transitions/fragment Ions

    i = 1
    transition_idx = 1
    prec_idx = 1
    rt_start = max(searchsortedfirst(rt_index.rt_bins, rt - rt_tol, lt=(r,x)->r.lb<x) - 1, 1) #First RT bin to search
    rt_stop = min(searchsortedlast(rt_index.rt_bins, rt + rt_tol, lt=(x, r)->r.ub>x) + 1, length(rt_index.rt_bins)) #Last RT bin to search 
    #prec_ids = UInt32[]
    for i in rt_start:rt_stop #Add transitions
        transition_idx, prec_idx = addTransitions!(transitions, prec_ids, transition_idx, prec_idx, fragment_list, rt_index.rt_bins[i].prec, prec_mz, prec_tol)
    end

    return sort!(transitions, by = x->getFragMZ(x)), prec_ids, transition_idx, prec_idx #Sort transitions by their fragment m/z. 
end
