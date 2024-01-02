




function selectTransitions!(transitions::Vector{LibraryFragment{V}},
                            precursors::Vector{LibraryPrecursor{Float32}},
                            fragment_list::Vector{Vector{LibraryFragment{V}}}, 
                            iso_splines::IsotopeSplineModel{Float64},
                            isotopes::Vector{Float64},
                            #precursors::Vector{LibraryPrecursor{Float32}}
                            counter::Counter{I,Float32}, 
                            topN::Int, 
                            iRT::U, 
                            iRT_tol::U, 
                            mz_bounds::Tuple{Float32, Float32};
                            block_size = 10000) where {T,V,U<:AbstractFloat, I<:Unsigned}
    #transitions = Vector{LibraryFragment{T}}()
    i = 1
    transition_idx = 0
    prec_idx = false
    while i <= min(topN, counter.matches)

        prec = precursors[getID(counter, i)]
        #Enforce iRT tolerance on precursors
        if abs(getIRT(prec) - iRT) > iRT_tol
            i += 1
            continue
        end

        #Manage isotope errors
        mz_low = - first(mz_bounds) - 3*NEUTRON/getCharge(prec)
        mz_high = last(mz_bounds) + NEUTRON/getCharge(prec)

        if (getMz(prec) < mz_low) | (getMz(prec) > mz_high)
            i += 1
            continue
        end

        prec_isotope_set = getPrecursorIsotopeSet(getMz(prec), getCharge(prec), mz_bounds)
        #println("prec.sequence ", prec.sequence)
        for frag in fragment_list[getID(counter, i)]
            #if abs((getiRT(precursors[getPrecID(frag)]) - rt)) > 5.0
            #    continue
            #end
            getFragAbundance!(
                        isotopes,
                        iso_splines,
                        isotope(Float64(frag.frag_mz*frag.frag_charge), Int64(frag.sulfur_count), 0),
                        isotope(Float64(prec.mz*prec.charge), Int64(prec.sulfur_count), 0),
                        prec_isotope_set)
            #println(isotopes)
            if length(isotopes) == 0
                continue
            end
            #println("isotopes $isotopes")
            #isotope_count = min(length(isotopes), 2)
            norm_fac = isotopes[1] + isotopes[2]
            isotopes[1] = isotopes[1]/(norm_fac)# = isotopes[1:isotope_count]/sum(isotopes[1:isotope_count])
            isotopes[2] = isotopes[2]/(norm_fac)
            for iso_idx in (0, )#(0, 1)

                transition_idx += 1

                transitions[transition_idx] = LibraryFragment(
                    Float32(frag.frag_mz + iso_idx*NEUTRON/frag.frag_charge),
                    frag.frag_charge,
                    frag.is_y_ion,
                    Bool(iso_idx), #Is the fragment an isotope?
                    frag.ion_position,
                    frag.ion_index,
                    max(Float32(frag.intensity*isotopes[iso_idx + 1]),zero(Float32)),
                    frag.prec_charge,
                    frag.prec_id,
                    frag.rank,
                    frag.sulfur_count
                )

                #Grow array if exceeds length
                if transition_idx >= length(transitions)
                    append!(transitions, [LibraryFragment{V}() for _ in range(1, block_size)])
                end
            end

            for i in range(1, length(isotopes))
                isotopes[i] = 0.0
            end

        end

        i += 1
    end

    sort!(@view(transitions[1:transition_idx]), 
            by = x->getFragMZ(x),
            alg=PartialQuickSort(1:transition_idx)
            #alg = TimSort)
    )

    reset!(counter)

    #return transition_idx, 0#sort!(transitions, by = x->getFragMZ(x))
    return transition_idx, false
end

#Get relevant framgents given a retention time and precursor mass using a retentionTimeIndex object
function selectRTIndexedTransitions!(transitions::Vector{LibraryFragment{V}}, 
                            precursors::Vector{LibraryPrecursor{Float32}},
                            fragment_list::Vector{Vector{LibraryFragment{V}}}, 
                            iso_splines::IsotopeSplineModel{Float64},
                            isotopes::Vector{Float64},
                            prec_ids::Vector{UInt32}, 
                            rt_index::Union{retentionTimeIndex{Float32, Float32}, Missing}, 
                            rt::U, 
                            rt_tol::U, 
                            prec_mz::U, 
                            prec_tol::U,
                            mz_bounds::Tuple{Float32, Float32};
                            block_size = 10000) where {T,U,V<:AbstractFloat}
    
    #Get matching precursors within an RT bin 
    function addTransitions!(transitions::Vector{LibraryFragment{V}}, prec_ids::Vector{UInt32}, transition_idx::Int64, prec_idx::Int64, fragment_list::Vector{Vector{LibraryFragment{V}}}, precs::Vector{Tuple{UInt32, U}}, prec_mz::U, prec_tol::U)
        start = searchsortedfirst(precs, by = x->last(x), prec_mz - prec_tol - 3*NEUTRON/2) #First precursor in the isolation window
        stop = searchsortedlast(precs, by = x->last(x), prec_mz + prec_tol + NEUTRON/2) #Last precursor in the isolation window
        for i in start:stop #Get transitions for each precursor
            prec_idx += 1
            prec = precursors[first(precs[i])]
            mz_low = - first(mz_bounds) - 3*NEUTRON/getCharge(prec)
            mz_high = last(mz_bounds) + NEUTRON/getCharge(prec)
            #$mz_low = - first(mz_bounds) #- 3*NEUTRON/getCharge(prec)
            #mz_high = last(mz_bounds) #+ NEUTRON/getCharge(prec)
            if (getMz(prec) < mz_low) | (getMz(prec) > mz_high)
                continue
            end
            prec_isotope_set = getPrecursorIsotopeSet(getMz(prec), getCharge(prec), mz_bounds)

            for frag in fragment_list[first(precs[i])]
                getFragAbundance!(
                        isotopes,
                        iso_splines,
                        isotope(Float64(frag.frag_mz*frag.frag_charge), Int64(frag.sulfur_count), 0),
                        isotope(Float64(prec.mz*prec.charge), Int64(prec.sulfur_count), 0),
                        prec_isotope_set)
                #println(isotopes)
                if length(isotopes) == 0
                    continue
                end
                #println("isotopes $isotopes")
                #isotope_count = min(length(isotopes), 2)
                #norm_fac = isotopes[1] + isotopes[2]
                #isotopes[1] = isotopes[1]/(norm_fac)# = isotopes[1:isotope_count]/sum(isotopes[1:isotope_count])
                #isotopes[2] = isotopes[2]/(norm_fac)
                for iso_idx in (0, )#(0, 1)

                    transition_idx += 1

                    transitions[transition_idx] = LibraryFragment(
                        Float32(frag.frag_mz + iso_idx*NEUTRON/frag.frag_charge),
                        frag.frag_charge,
                        frag.is_y_ion,
                        Bool(iso_idx), #Is the fragment an isotope?
                        frag.ion_position,
                        frag.ion_index,
            
                        max(Float32(frag.intensity*isotopes[iso_idx + 1]),zero(Float32)),
                        frag.prec_charge,
                        frag.prec_id,
                        frag.rank,
                        frag.sulfur_count
                    )

                    #Grow array if exceeds length
                    if transition_idx >= length(transitions)
                        append!(transitions, [LibraryFragment{V}() for _ in range(1, block_size)])
                    end
                end

                for i in range(1, length(isotopes)) #Critical step
                    isotopes[i] = 0.0
                end

            end
            #println("precursor ", prec.sequence)
            #println("prec_isotope_set $prec_isotope_set")
            #println("transitions ", transitions[1:transition_idx])
            #Grow array if exceeds length
            (prec_idx > length(prec_ids)) ? append!(prec_ids, zeros(UInt32, block_size)) : nothing

            prec_ids[prec_idx] = first(precs[i])
        end
        return transition_idx, prec_idx
    end

    #transitions = Vector{LibraryFragment{V}}() #Initialize list of transitions/fragment Ions

    i = 1
    transition_idx = 0
    prec_idx = 0
    rt_start = max(searchsortedfirst(rt_index.rt_bins, rt - rt_tol, lt=(r,x)->r.lb<x) - 1, 1) #First RT bin to search
    rt_stop = min(searchsortedlast(rt_index.rt_bins, rt + rt_tol, lt=(x, r)->r.ub>x) + 1, length(rt_index.rt_bins)) #Last RT bin to search 
    #prec_ids = UInt32[]
    for i in rt_start:rt_stop #Add transitions
        transition_idx, prec_idx = addTransitions!(transitions, prec_ids, transition_idx, prec_idx, fragment_list, rt_index.rt_bins[i].prec, prec_mz, prec_tol)
    end

    sort!(@view(transitions[1:transition_idx]), 
          by = x->getFragMZ(x),
          alg=PartialQuickSort(1:transition_idx)) #Optimize?
    #return sort!(transitions, by = x->getFragMZ(x)), prec_ids, transition_idx, prec_idx #Sort transitions by their fragment m/z. 
    return transition_idx, prec_idx
end

function selectIsotopes!(isotopes::Vector{Isotope{T}},
                        prec_list::Vector{Tuple{Union{V, Missing}, UInt32}}, 
                        isotope_dict::UnorderedDictionary{UInt32, Vector{Isotope{T}}}, 
                        prec_ids::Vector{UInt32}, 
                        rt::U, 
                        rt_tol::U) where {T,U,V<:AbstractFloat}
    i = 0
    ion_idx = 0
    prec_idx = 0
    rt_start = searchsortedfirst(prec_list, rt - rt_tol, lt=(r,x)->first(r)<x) #First RT bin to search
    rt_stop = searchsortedlast(prec_list, rt + rt_tol, lt=(x, r)->first(r)>x) #Last RT bin to search 
    #return rt_start, rt_stop
    for i in range(rt_start, rt_stop)
        prec_idx += 1
        for iso in isotope_dict[last(prec_list[i])]
            ion_idx += 1
            (ion_idx > length(isotopes)) ?  append!(isotopes, [Isotope{T}() for _ in range(1, block_size)]) : nothing

            isotopes[ion_idx] = iso
            #append!(isotopes, isotope_dict[last(prec_list[i])])
        end
        (prec_idx > length(prec_ids)) ? append!(prec_ids, zeros(UInt32, block_size)) : nothing
        prec_ids[prec_idx] = last(prec_list[i])
    end
    sort!(@view(isotopes[1:ion_idx]), 
        by = x->getMZ(x),
        alg=PartialQuickSort(1:ion_idx))
    return ion_idx, prec_idx#sort(isotopes, by = x->getMZ(x))
end
