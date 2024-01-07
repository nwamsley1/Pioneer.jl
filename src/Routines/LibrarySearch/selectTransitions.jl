



function selectTransitions!(transitions::Vector{LibraryFragment{T}},
                            precursors::Vector{LibraryPrecursor{T}},
                            fragment_list::Vector{Vector{LibraryFragment{T}}}, 
                            iso_splines::IsotopeSplineModel{U},
                            isotopes::Vector{U},
                            #precursors::Vector{LibraryPrecursor{Float32}}
                            counter::Counter{I,T}, 
                            topN::Int, 
                            iRT::V, 
                            iRT_tol::V, 
                            mz_bounds::Tuple{W, W};
                            isotope_err_bounds::Tuple{Int64,Int64} = (3, 1),
                            block_size = 10000) where {T,V,U,W<:AbstractFloat, I<:Unsigned}

    i = 1
    transition_idx = 0

    while i <= min(topN, counter.matches)

       
        prec_idx = getID(counter, i)
        prec = precursors[prec_idx]
        i += 1
        #Enforce iRT tolerance on precursors
        abs(getIRT(prec) - iRT) > iRT_tol ? continue : nothing

        #Manage isotope errors
        mz_low = - first(mz_bounds) - first(isotope_err_bounds)*NEUTRON/getCharge(prec)
        mz_high = last(mz_bounds) + last(isotope_err_bounds)*NEUTRON/getCharge(prec)

        (getMz(prec) < mz_low) | (getMz(prec) > mz_high) ? continue : nothing

        #println("prec.sequence ", prec.sequence)
        transition_idx = fillTransitionList!(transitions, 
                                            prec_idx,
                                            transition_idx,
                                            fragment_list, 
                                            isotopes, iso_splines, prec, mz_bounds, block_size)

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
function selectRTIndexedTransitions!(transitions::Vector{LibraryFragment{T}}, 
                            precursors::Vector{LibraryPrecursor{T}},
                            fragment_list::Vector{Vector{LibraryFragment{T}}}, 
                            iso_splines::IsotopeSplineModel{U},
                            isotopes::Vector{U},
                            prec_ids::Vector{UInt32}, 
                            rt_index::Union{retentionTimeIndex{T, T}, Missing}, 
                            rt::T, 
                            rt_tol::T, 
                            mz_bounds::Tuple{T, T};
                            isotope_err_bounds::Tuple{I,I} = (3, 1),
                            block_size = 10000) where {T,U<:AbstractFloat,I<:Integer}
    i = 1
    transition_idx = 0
    prec_idx = 0
    rt_start = max(searchsortedfirst(rt_index.rt_bins, rt - rt_tol, lt=(r,x)->r.lb<x) - 1, 1) #First RT bin to search
    rt_stop = min(searchsortedlast(rt_index.rt_bins, rt + rt_tol, lt=(x, r)->r.ub>x) + 1, length(rt_index.rt_bins)) #Last RT bin to search 

    for i in rt_start:rt_stop #Add transitions
        precs = rt_index.rt_bins[i].prec
        start = searchsortedfirst(precs, by = x->last(x), first(mz_bounds) - first(isotope_err_bounds)*NEUTRON/2) #First precursor in the isolation window
        stop = searchsortedlast(precs, by = x->last(x), last(mz_bounds) + last(isotope_err_bounds)*NEUTRON/2) #Last precursor in the isolation window
        for i in start:stop #Get transitions for each precursor
            prec_idx += 1 #Keep track of number of precursors 
            prec = precursors[first(precs[i])]
            mz_low = - first(mz_bounds) - first(isotope_err_bounds)*NEUTRON/getCharge(prec)
            mz_high = last(mz_bounds) + last(isotope_err_bounds)*NEUTRON/getCharge(prec)

            #If precursor m/z (with isotope error) out of qaudrupole isolation bounds 
            (getMz(prec) < mz_low) | (getMz(prec) > mz_high) ? continue : nothing

            #Which precursor isotopes where captured in the quadrupole isolation window? 
            #For example, return (0, 3) if M+0 through M+3 isotopes were captured 
            transition_idx = fillTransitionList!(transitions, 
                                                first(precs[i]), #prec_idx
                                                transition_idx,
                                                fragment_list, 
                                                isotopes, iso_splines, prec, mz_bounds, block_size)

            #Grow array if exceeds length
            (prec_idx > length(prec_ids)) ? append!(prec_ids, zeros(UInt32, block_size)) : nothing

            prec_ids[prec_idx] = first(precs[i])
        end
    end

    sort!(@view(transitions[1:transition_idx]), 
          by = x->getFragMZ(x),
          alg=PartialQuickSort(1:transition_idx))

    #return sort!(transitions, by = x->getFragMZ(x)), prec_ids, transition_idx, prec_idx #Sort transitions by their fragment m/z. 
    return transition_idx, prec_idx
end

function fillTransitionList!(transitions::Vector{LibraryFragment{T}}, prec_idx::UInt32, transition_idx::Int64, fragment_list::Vector{Vector{LibraryFragment{T}}}, isotopes::Vector{U}, iso_splines::IsotopeSplineModel{U}, prec::LibraryPrecursor{V},  mz_bounds::Tuple{W,W}, block_size::Int64)::Int64 where {T,U,V,W<:AbstractFloat,I<:Integer}
    
    prec_isotope_set = getPrecursorIsotopeSet(getMz(prec), getCharge(prec), mz_bounds)
    #for frag in fragment_list[getID(counter, i)]
    for frag in fragment_list[prec_idx]

        #Estimate isotope abundances 
        getFragIsotopes!(isotopes, iso_splines, prec, frag, prec_isotope_set)

        for iso_idx in range(0, length(isotopes) - 1)

            #Skip if missing
            isnan(isotopes[iso_idx + 1]) ? continue : nothing

            #if iso_idx > 0
            #    continue
            #end
            if isotopes[iso_idx + 1]/first(isotopes) < 0.5
                continue
            end

            transition_idx += 1

            transitions[transition_idx] = LibraryFragment(
                Float32(frag.frag_mz + iso_idx*NEUTRON/frag.frag_charge), #Estimated isotopic m/z
                frag.frag_charge,
                frag.is_y_ion,
                iso_idx>0, #Is the fragment an isotope?
                frag.ion_position,
                frag.ion_index,
                Float32(isotopes[iso_idx + 1]), #Estimated relative abundance 
                frag.prec_charge,
                frag.prec_id,
                frag.rank,
                frag.sulfur_count
            )::LibraryFragment{T}
            

            #Grow array if exceeds length
            
            if transition_idx >= length(transitions)
                append!(transitions, [LibraryFragment{T}() for _ in range(1, block_size)])
            end
        end
    end
    return transition_idx
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
