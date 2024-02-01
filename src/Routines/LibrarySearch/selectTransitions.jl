function selectTransitions!(transitions::Vector{DetailedFrag{Float32}},
                            precursors::Vector{LibraryPrecursorIon{Float32}},
                            #fragment_list::Vector{Vector{DetailedFrag{Float32}}}, 
                            library_fragment_lookup::LibraryFragmentLookup{Float32}, 
                            iso_splines::IsotopeSplineModel{Float64},
                            isotopes::Vector{Float64},
                            #precursors::Vector{LibraryPrecursor{Float32}}
                            counter::Counter{UInt32,UInt8}, 
                            topN::Int64, 
                            iRT::Float32, 
                            iRT_tol::Float32, 
                            mz_bounds::Tuple{Float32, Float32};
                            isotope_err_bounds::Tuple{Int64,Int64} = (3, 1),
                            block_size::Int64 = 10000)# where {V,W<:AbstractFloat}

    i = 1
    transition_idx = 0
    #NEUTRON = Float64(1.00335)
    while i <= min(topN, counter.matches)

       
        prec_idx = getID(counter, i)
        prec = precursors[prec_idx]
        #Enforce iRT tolerance on precursors
        if abs(getIRT(prec) - iRT) > iRT_tol
            i += 1
             continue
        end

        #Manage isotope errors
        mz_low = - first(mz_bounds) - first(isotope_err_bounds)*NEUTRON/getPrecCharge(prec)
        mz_high = last(mz_bounds) + last(isotope_err_bounds)*NEUTRON/getPrecCharge(prec)

        if (getMZ(prec) < mz_low) | (getMZ(prec) > mz_high)
            i += 1
            continue
        end

        for frag_idx in getPrecFragRange(library_fragment_lookup, prec_idx)#fragment_list[getID(counter, i)]
            transition_idx += 1
            transitions[transition_idx] = getFrag(library_fragment_lookup, frag_idx) #fragment_list.frags[frag_idx]
            
            if transition_idx > length(transitions)
                append!(transitions, [DetailedFrag{Float32}() for _ in range(1, block_size)])
            end
            #Grow array if exceeds length
        end
        i += 1
    end

    sort!(@view(transitions[1:transition_idx]), 
            by = x->getMZ(x),
            alg=PartialQuickSort(1:transition_idx)
            #alg = TimSort)
    )

    reset!(counter)

    #return transition_idx, 0#sort!(transitions, by = x->getFragMZ(x))
    return transition_idx, false
end

#Get relevant framgents given a retention time and precursor mass using a retentionTimeIndex object
function selectRTIndexedTransitions!(transitions::Vector{DetailedFrag{Float32}}, 
                            precursors::Vector{LibraryPrecursorIon{Float32}},
                            library_fragment_lookup::LibraryFragmentLookup{Float32}, 
                            iso_splines::IsotopeSplineModel{Float64},
                            isotopes::Vector{Float64},
                            prec_ids::Vector{UInt32}, 
                            rt_index::Union{retentionTimeIndex{Float32, Float32}, Missing}, 
                            rt::Float32, 
                            rt_tol::Float32, 
                            mz_bounds::Tuple{Float32, Float32};
                            isotope_err_bounds::Tuple{Int64, Int64} = (3, 1),
                            block_size = 10000)
    transition_idx = 0
    prec_idx = 0
    rt_start = max(searchsortedfirst(rt_index.rt_bins, rt - rt_tol, lt=(r,x)->r.lb<x) - 1, 1) #First RT bin to search
    rt_stop = min(searchsortedlast(rt_index.rt_bins, rt + rt_tol, lt=(x, r)->r.ub>x) + 1, length(rt_index.rt_bins)) #Last RT bin to search 
    #rt_start = max(searchsortedfirst(rt_index.rt_bins, rt, lt=(r,x)->r.lb<x) - 1, 1)
    #rt_stop = rt_start
    for rt_bin_idx in rt_start:rt_stop #Add transitions
        precs = rt_index.rt_bins[rt_bin_idx].prec
        start = searchsortedfirst(precs, by = x->last(x), first(mz_bounds) - first(isotope_err_bounds)*NEUTRON/2) #First precursor in the isolation window
        stop = searchsortedlast(precs, by = x->last(x), last(mz_bounds) + last(isotope_err_bounds)*NEUTRON/2) #Last precursor in the isolation window
        for i in start:stop #Get transitions for each precursor
            prec_idx += 1 #Keep track of number of precursors 
            prec = precursors[first(precs[i])]
            mz_low = - first(mz_bounds) - first(isotope_err_bounds)*NEUTRON/getPrecCharge(prec)
            mz_high = last(mz_bounds) + last(isotope_err_bounds)*NEUTRON/getPrecCharge(prec)

            #If precursor m/z (with isotope error) out of qaudrupole isolation bounds 
            (getMZ(prec) < mz_low) | (getMZ(prec) > mz_high) ? continue : nothing

            #Which precursor isotopes where captured in the quadrupole isolation window? 
            #For example, return (0, 3) if M+0 through M+3 isotopes were captured 
            transition_idx = @inline fillTransitionList!(transitions, 
                                                first(precs[i]), #prec_idx
                                                transition_idx,
                                                library_fragment_lookup, 
                                                isotopes, iso_splines, prec, mz_bounds, block_size)

            #Grow array if exceeds length
            (prec_idx > length(prec_ids)) ? append!(prec_ids, zeros(UInt32, block_size)) : nothing

            prec_ids[prec_idx] = first(precs[i])
        end
    end

    sort!(@view(transitions[1:transition_idx]), 
          by = x->getMZ(x),
          alg=PartialQuickSort(1:transition_idx))

    #return sort!(transitions, by = x->getFragMZ(x)), prec_ids, transition_idx, prec_idx #Sort transitions by their fragment m/z. 
    return transition_idx, prec_idx
end

function fillTransitionList!(transitions::Vector{DetailedFrag{Float32}}, 
                            prec_idx::UInt32, 
                            transition_idx::Int64, 
                            library_fragment_lookup::LibraryFragmentLookup{Float32}, 
                            isotopes::Vector{Float64}, iso_splines::IsotopeSplineModel{Float64}, 
                            prec::LibraryPrecursorIon{Float32},  
                            mz_bounds::Tuple{Float32,Float32}, block_size::Int64)::Int64 #where {T,U,V,W<:AbstractFloat,I<:Integer}
    NEUTRON = Float64(1.00335)
    prec_isotope_set = getPrecursorIsotopeSet(getMZ(prec), getPrecCharge(prec), mz_bounds)
    #for frag in fragment_list[getID(counter, i)]
    for frag_idx in getPrecFragRange(library_fragment_lookup, prec_idx)

        frag = getFrag(library_fragment_lookup, frag_idx)
        #Estimate isotope abundances 
        getFragIsotopes!(isotopes, iso_splines, prec, frag, prec_isotope_set)

        for iso_idx in range(0, length(isotopes) - 1)

            #Skip if missing
            isnan(isotopes[iso_idx + 1]) ? continue : nothing
            #iso_idx > 0 ? continue : nothing
            #isotopes[iso_idx + 1]/first(isotopes) < 0.5 ? continue : nothing

            
            
            transition_idx += 1
            transitions[transition_idx] = DetailedFrag(
                frag.prec_id,

                Float32(frag.mz + iso_idx*NEUTRON/frag.frag_charge), #Estimated isotopic m/z
                Float16(isotopes[iso_idx + 1]), #Estimated relative abundance 

                frag.is_y_ion,
                iso_idx>0, #Is the fragment an isotope?

                frag.frag_charge,
                frag.ion_position,
                frag.prec_charge,
                frag.rank,
                frag.sulfur_count
            )#::LibraryFragment{T}
            

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