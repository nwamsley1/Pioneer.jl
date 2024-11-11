function selectTransitions!(transitions::Vector{DetailedFrag{Float32}},
                            scan_to_prec_idx::UnitRange{Int64},
                            precursors_passed_scoring::Vector{UInt32},
                            prec_mzs::Arrow.Primitive{Float32, Vector{Float32}},
                            prec_charges::Arrow.Primitive{UInt8, Vector{UInt8}},
                            prec_irts::Arrow.Primitive{Float32, Vector{Float32}},
                            prec_sulfur_counts::Arrow.Primitive{UInt8, Vector{UInt8}},
                            iso_splines::IsotopeSplineModel,
                            quad_transmission_func::QuadTransmissionFunction,
                            precursor_transmission::Vector{Float32},
                            isotopes::Vector{Float32},
                            n_frag_isotopes::Int64,
                            library_fragment_lookup::LibraryFragmentLookup{Float32}, 
                            iRT::Float32, 
                            iRT_tol::Float32, 
                            frag_mz_bounds::Tuple{Float32, Float32};
                            isotope_err_bounds::Tuple{Int64,Int64} = (3, 1),
                            block_size::Int64 = 10000)# where {V,W<:AbstractFloat}
    #isotopes = Float32[0.0f0]
    transition_idx = 0
    for i in scan_to_prec_idx

       
        prec_idx =  precursors_passed_scoring[i]
        prec_charge = prec_charges[prec_idx]
        prec_mz = prec_mzs[prec_idx]
        #Enforce iRT tolerance on precursors
        if abs(prec_irts[prec_idx] - iRT) > iRT_tol
             continue
        end

        #Manage isotope errors. NEUTRON is global constant. 
        mz_low = getPrecMinBound(quad_transmission_func) - NEUTRON*first(isotope_err_bounds)/prec_charge
        mz_high = getPrecMaxBound(quad_transmission_func) + NEUTRON*last(isotope_err_bounds)/prec_charge
        if (prec_mz < mz_low) | (prec_mz > mz_high)
            continue
        end
        prec_sulfur_count, prec_charge, prec_mz = prec_sulfur_counts[prec_idx], prec_charges[prec_idx], prec_mzs[prec_idx]
        transition_idx = @inline fillTransitionList!(
                                                    transitions, 
                                                    getPrecFragRange(library_fragment_lookup, prec_idx),
                                                    getFragments(library_fragment_lookup),
                                                    prec_mz,
                                                    prec_charge,
                                                    prec_sulfur_count,
                                                    transition_idx,
                                                    quad_transmission_func,
                                                    precursor_transmission,
                                                    isotopes, 
                                                    n_frag_isotopes,
                                                    iso_splines, 
                                                    frag_mz_bounds,
                                                    block_size
                                                    )
        
        #=
        for frag_idx in getPrecFragRange(library_fragment_lookup, prec_idx)#fragment_list[getID(counter, i)]
            frag = getFrag(library_fragment_lookup, frag_idx) 
            if (getMZ(frag) < first(frag_mz_bounds)) |  (getMZ(frag) > last(frag_mz_bounds))
               continue
            end
            transition_idx += 1
            transitions[transition_idx] = getFrag(library_fragment_lookup, frag_idx) #fragment_list.frags[frag_idx]
            if transition_idx + 1 > length(transitions)
                append!(transitions, [DetailedFrag{Float32}() for _ in range(1, block_size)])
            end
            #Grow array if exceeds length
        end
        =#
    end

    sort!(@view(transitions[1:transition_idx]), 
            by = x->getMZ(x),
            alg=PartialQuickSort(1:transition_idx)
            #alg = TimSort)
    )

    #reset!(counter)

    #return transition_idx, 0#sort!(transitions, by = x->getFragMZ(x))
    return transition_idx, false
end

function selectTransitions!(
                            transitions::Vector{DetailedFrag{Float32}},
                            library_fragment_lookup::LibraryFragmentLookup{Float32},
                            scan_to_prec_idx::UnitRange{Int64},
                            precursors_passed_scoring::Vector{UInt32};
                            max_rank::Int64 = 5,
                            block_size::Int64 = 10000)# where {V,W<:AbstractFloat}

    transition_idx = 0
    for i in scan_to_prec_idx
        prec_idx =  precursors_passed_scoring[i]
        for frag_idx in getPrecFragRange(library_fragment_lookup, prec_idx)
            frag = getFrag(library_fragment_lookup, frag_idx) 
            if getRank(frag) <= max_rank
                transition_idx += 1
                transitions[transition_idx] = frag
            end
            if transition_idx + 1 > length(transitions)
                append!(transitions, [DetailedFrag{Float32}() for _ in range(1, block_size)])
            end
            #Grow array if exceeds length
        end
    end

    sort!(@view(transitions[1:transition_idx]), 
            by = x->getMZ(x),
            alg=PartialQuickSort(1:transition_idx)
            #alg = TimSort)
    )

    #reset!(counter)

    #return transition_idx, 0#sort!(transitions, by = x->getFragMZ(x))
    return transition_idx, false
end
#Get relevant framgents given a retention time and precursor mass using a retentionTimeIndex object
function selectRTIndexedTransitions!(
                            transitions::Vector{DetailedFrag{Float32}}, 
                            precs_temp::Vector{UInt32},
                            precs_temp_size::Int64,
                            library_fragment_lookup::LibraryFragmentLookup{Float32}, 
                            prec_mzs::AbstractArray{Float32},
                            prec_charges::AbstractArray{UInt8},
                            prec_sulfur_counts::AbstractArray{UInt8},
                            iso_splines::IsotopeSplineModel,
                            quad_transmission_func::QuadTransmissionFunction,
                            precursor_transmission::Vector{Float32},
                            isotopes::Vector{Float32},
                            n_frag_isotopes::Int64,
                            rt_index::Union{retentionTimeIndex{Float32, Float32}, Missing}, 
                            rt_start_idx::Int64, 
                            rt_stop_idx::Int64,
                            frag_mz_bounds::Tuple{Float32, Float32},
                            isotope_err_bounds::Tuple{Int64, Int64},
                            block_size::Int64)
    transition_idx = 0
    n = 0
    min_prec_mz = getPrecMinBound(quad_transmission_func)
    max_prec_mz  = getPrecMaxBound(quad_transmission_func)
    for rt_bin_idx in range(rt_start_idx, rt_stop_idx) #Add transitions
        precs = rt_index.rt_bins[rt_bin_idx].prec
        start = searchsortedfirst(precs, by = x->last(x), min_prec_mz - first(isotope_err_bounds)*NEUTRON/2) #First precursor in the isolation window
        stop = searchsortedlast(precs, by = x->last(x), max_prec_mz + last(isotope_err_bounds)*NEUTRON/2) #Last precursor in the isolation window
        for i in start:stop #Get transitions for each precursor
            #println("entered loop")
            precs_temp_size += 1
            n += 1 #Keep track of number of precursors 
            prec_idx = first(precs[i])
            precs_temp[precs_temp_size] = prec_idx
            prec_sulfur_count, prec_charge, prec_mz = prec_sulfur_counts[prec_idx], prec_charges[prec_idx], prec_mzs[prec_idx]
            mz_low = min_prec_mz - first(isotope_err_bounds)*NEUTRON/prec_charge
            mz_high = max_prec_mz + last(isotope_err_bounds)*NEUTRON/prec_charge
            #If precursor m/z (with isotope error) out of qaudrupole isolation bounds 
            (prec_mz < mz_low) | (prec_mz > mz_high) ? continue : nothing
            #println("passed ternary")
            #Which precursor isotopes where captured in the quadrupole isolation window? 
            #For example, return (0, 3) if M+0 through M+3 isotopes were captured 
            transition_idx = @inline fillTransitionList!(transitions, 
                                            getPrecFragRange(library_fragment_lookup, prec_idx),
                                            getFragments(library_fragment_lookup),
                                            prec_mz,
                                            prec_charge,
                                            prec_sulfur_count,
                                            transition_idx,
                                            quad_transmission_func,
                                            precursor_transmission,
                                            isotopes, 
                                            n_frag_isotopes,
                                            iso_splines, 
                                            frag_mz_bounds,
                                            block_size
                                            )
        end
    end

    sort!(@view(transitions[1:transition_idx]), 
          by = x->getMZ(x),
          alg=PartialQuickSort(1:transition_idx))
          #println("transition_idx $transition_idx, n $n precs_temp_size $precs_temp_size")
    #return sort!(transitions, by = x->getFragMZ(x)), prec_ids, transition_idx, prec_idx #Sort transitions by their fragment m/z. 
    return transition_idx, n, precs_temp_size
end
function selectRTIndexedTransitions!(
                            transitions::Vector{DetailedFrag{Float32}}, 
                            precs_temp::Vector{UInt32},
                            precs_temp_size::Int64,
                            precursors_passing::Set{UInt32},
                            library_fragment_lookup::LibraryFragmentLookup{Float32}, 
                            prec_mzs::AbstractArray{Float32},
                            prec_charges::AbstractArray{UInt8},
                            prec_sulfur_counts::AbstractArray{UInt8},
                            iso_splines::IsotopeSplineModel,
                            quad_transmission_func::QuadTransmissionFunction,
                            precursor_transmission::Vector{Float32},
                            isotopes::Vector{Float32},
                            n_frag_isotopes::Int64,
                            rt_index::Union{retentionTimeIndex{Float32, Float32}, Missing}, 
                            rt_start_idx::Int64, 
                            rt_stop_idx::Int64,
                            frag_mz_bounds::Tuple{Float32, Float32},
                            isotope_err_bounds::Tuple{Int64, Int64},
                            block_size::Int64)
    transition_idx = 0
    n = 0
    min_prec_mz = getPrecMinBound(quad_transmission_func)
    max_prec_mz  = getPrecMaxBound(quad_transmission_func)
    for rt_bin_idx in range(rt_start_idx, rt_stop_idx) #Add transitions
        precs = rt_index.rt_bins[rt_bin_idx].prec
        start = searchsortedfirst(precs, by = x->last(x), min_prec_mz - first(isotope_err_bounds)*NEUTRON/2) #First precursor in the isolation window
        stop = searchsortedlast(precs, by = x->last(x), max_prec_mz + last(isotope_err_bounds)*NEUTRON/2) #Last precursor in the isolation window
        for i in start:stop #Get transitions for each precursor
            prec_idx = first(precs[i])
            if prec_idx ∉ precursors_passing
                continue
            end
            precs_temp_size += 1
            n += 1 #Keep track of number of precursors 
            precs_temp[precs_temp_size] = prec_idx
            prec_sulfur_count, prec_charge, prec_mz = prec_sulfur_counts[prec_idx], prec_charges[prec_idx], prec_mzs[prec_idx]
            mz_low = min_prec_mz - first(isotope_err_bounds)*NEUTRON/prec_charge
            mz_high = max_prec_mz + last(isotope_err_bounds)*NEUTRON/prec_charge

            #If precursor m/z (with isotope error) out of qaudrupole isolation bounds 
            (prec_mz < mz_low) | (prec_mz > mz_high) ? continue : nothing

            #Which precursor isotopes where captured in the quadrupole isolation window? 
            #For example, return (0, 3) if M+0 through M+3 isotopes were captured 
            transition_idx = @inline fillTransitionList!(transitions, 
                                            getPrecFragRange(library_fragment_lookup, prec_idx),
                                            getFragments(library_fragment_lookup),
                                            prec_mz,
                                            prec_charge,
                                            prec_sulfur_count,
                                            transition_idx,
                                            quad_transmission_func,
                                            precursor_transmission,
                                            isotopes, 
                                            n_frag_isotopes,
                                            iso_splines, 
                                            frag_mz_bounds,
                                            block_size
                                            )
        end
    end

    sort!(@view(transitions[1:transition_idx]), 
          by = x->getMZ(x),
          alg=PartialQuickSort(1:transition_idx))
    #return sort!(transitions, by = x->getFragMZ(x)), prec_ids, transition_idx, prec_idx #Sort transitions by their fragment m/z. 
    return transition_idx, n, precs_temp_size
end

function selectTransitionsForQuadEstimation!(
                            prec_idxs::AbstractVector{UInt32},
                            transitions::Vector{DetailedFrag{Float32}}, 
                            library_fragment_lookup::LibraryFragmentLookup{Float32}, 
                            prec_mzs::AbstractArray{Float32},
                            prec_charges::AbstractArray{UInt8},
                            prec_sulfur_counts::AbstractArray{UInt8},
                            iso_splines::IsotopeSplineModel,
                            precursor_transmission::Vector{Float32},
                            isotopes::Vector{Float32},
                            frag_mz_bounds::Tuple{Float32, Float32},
                            block_size::Int64)
    transition_idx = 0
    n = 0
    for prec_idx in prec_idxs
        for prec_iso_idx in range(0, 2)
            n += 1 #Keep track of number of precursors 
            prec_sulfur_count, prec_charge, prec_mz = prec_sulfur_counts[prec_idx], prec_charges[prec_idx], prec_mzs[prec_idx]
            #Which precursor isotopes where captured in the quadrupole isolation window? 
            #For example, return (0, 3) if M+0 through M+3 isotopes were captured 
            transition_idx = @inline fillTransitionListForQuadEstimation!(transitions, 
                                            getPrecFragRange(library_fragment_lookup, prec_idx),
                                            getFragments(library_fragment_lookup),
                                            prec_mz,
                                            prec_charge,
                                            prec_sulfur_count,
                                            transition_idx,
                                            precursor_transmission,
                                            isotopes, 
                                            prec_iso_idx,
                                            iso_splines, 
                                            frag_mz_bounds,
                                            block_size
                                            )
        end
    end

    sort!(@view(transitions[1:transition_idx]), 
          by = x->getMZ(x),
          alg=PartialQuickSort(1:transition_idx))
          #println("transition_idx $transition_idx, n $n precs_temp_size $precs_temp_size")
    #return sort!(transitions, by = x->getFragMZ(x)), prec_ids, transition_idx, prec_idx #Sort transitions by their fragment m/z. 
    return transition_idx, n
end


function fillTransitionList!(transitions::Vector{DetailedFrag{Float32}}, 
                            precursor_fragment_range::UnitRange{UInt64},
                            fragment_ions::Vector{DetailedFrag{Float32}},
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            transition_idx::Int64, 
                            quad_transmission_func::QuadTransmissionFunction,
                            precursor_transmission::Vector{Float32},
                            isotopes::Vector{Float32}, 
                            n_frag_isotopes::Int64,
                            iso_splines::IsotopeSplineModel, 
                            frag_mz_bounds::Tuple{Float32, Float32},
                            block_size::Int64)::Int64 #where {T,U,V,W<:AbstractFloat,I<:Integer}

    #NEUTRON = Float64(1.00335)

    #precursor_transmission = zeros(Float32, 5)
    getPrecursorIsotopeTransmission!(
        precursor_transmission,
        prec_mz,
        prec_charge,
        quad_transmission_func
    )
    prec_isotope_set = getPrecursorIsotopeSet(
        prec_mz,
        prec_charge,
        quad_transmission_func
    )
    if (first(prec_isotope_set) > 0) | (last(prec_isotope_set) < 2)
        transition_idx = fillPrecursorFragments!(
            prec_isotope_set,
            transitions,
            precursor_fragment_range,
            fragment_ions,
            prec_mz,
            prec_charge,
            prec_sulfur_count,
            transition_idx,
            precursor_transmission,
            isotopes,
            n_frag_isotopes,
            iso_splines,
            frag_mz_bounds,
            block_size)
    else
        transition_idx = fillPrecursorFragments!(
            prec_isotope_set,
            transitions,
            precursor_fragment_range,
            fragment_ions,
            transition_idx,
            n_frag_isotopes,
            iso_splines,
            frag_mz_bounds,
            block_size)
    end
    return transition_idx
end

function fillPrecursorFragments!(
    prec_isotope_set::Tuple{Int64, Int64},
    transitions::Vector{DetailedFrag{Float32}}, 
    precursor_fragment_range::UnitRange{UInt64},
    fragment_ions::Vector{DetailedFrag{Float32}},
    prec_mz::Float32,
    prec_charge::UInt8,
    prec_sulfur_count::UInt8,
    transition_idx::Int64, 
    precursor_transmission::Vector{Float32},
    isotopes::Vector{Float32}, 
    n_frag_isotopes::Int64,
    iso_splines::IsotopeSplineModel, 
    frag_mz_bounds::Tuple{Float32, Float32},
    block_size::Int64)
    for frag_idx in precursor_fragment_range
        frag = fragment_ions[frag_idx]
        if frag.rank > 25
            continue
        end
        #Estimate isotope abundances 
        
        getFragIsotopes!(isotopes, 
                        precursor_transmission,
                        iso_splines, 
                        prec_mz,
                        prec_charge, 
                        prec_sulfur_count,
                        frag,
                        )
        
        for iso_idx in range(0, min(n_frag_isotopes - 1, last(prec_isotope_set)))

            frag_mz = Float32(frag.mz + iso_idx*NEUTRON/frag.frag_charge)
            if (frag_mz < first(frag_mz_bounds)) |  (frag_mz > last(frag_mz_bounds))
                continue
            end

            transition_idx += 1
            transitions[transition_idx] = DetailedFrag(
                frag.prec_id,

                frag_mz, #Estimated isotopic m/z
                Float16(isotopes[iso_idx + 1]), #Estimated relative abundance 

                
                frag.ion_type,
                frag.is_y,
                frag.is_b,
                frag.is_p,
                iso_idx>0, #Is the fragment an isotope?

                frag.frag_charge,
                frag.ion_position,
                frag.prec_charge,
                frag.rank,
                frag.sulfur_count
            )#::LibraryFragment{T}
            

            #Grow array if exceeds length
            
            if transition_idx >= length(transitions)
                append!(transitions, [DetailedFrag{Float32}() for _ in range(1, block_size)])
            end
        end
    end
    return transition_idx
end

function fillPrecursorFragments!(
    prec_isotope_set::Tuple{Int64, Int64},
    transitions::Vector{DetailedFrag{Float32}}, 
    precursor_fragment_range::UnitRange{UInt64},
    fragment_ions::Vector{DetailedFrag{Float32}},
    transition_idx::Int64, 
    n_frag_isotopes::Int64,
    iso_splines::IsotopeSplineModel, 
    frag_mz_bounds::Tuple{Float32, Float32},
    block_size::Int64)
    for frag_idx in precursor_fragment_range
        frag = fragment_ions[frag_idx]
        if frag.rank > 25
            continue
        end
        #Estimate isotope abundances 
        for iso_idx in range(0, min(n_frag_isotopes - 1, last(prec_isotope_set)))

            frag_mz = Float32(frag.mz + iso_idx*NEUTRON/frag.frag_charge)
            if (frag_mz < first(frag_mz_bounds)) |  (frag_mz > last(frag_mz_bounds))
                continue
            end
            intensity = Float16(frag.intensity*iso_splines(min(Int64(frag.sulfur_count), 5), iso_idx, frag_mz*frag.frag_charge))
            transition_idx += 1
            transitions[transition_idx] = DetailedFrag(
                frag.prec_id,

                frag_mz, #Estimated isotopic m/z
                intensity, #Estimated relative abundance 

                
                frag.ion_type,
                frag.is_y,
                frag.is_b,
                frag.is_p,
                iso_idx>0, #Is the fragment an isotope?

                frag.frag_charge,
                frag.ion_position,
                frag.prec_charge,
                frag.rank,
                frag.sulfur_count
            )#::LibraryFragment{T}
            

            #Grow array if exceeds length
            
            if transition_idx >= length(transitions)
                append!(transitions, [DetailedFrag{Float32}() for _ in range(1, block_size)])
            end
        end
    end
    return transition_idx
end

function fillTransitionListForQuadEstimation!(transitions::Vector{DetailedFrag{Float32}}, 
                            precursor_fragment_range::UnitRange{UInt64},
                            fragment_ions::Vector{DetailedFrag{Float32}},
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            transition_idx::Int64, 
                            precursor_transmission::Vector{Float32},
                            isotopes::Vector{Float32}, 
                            prec_iso_idx::Int64,
                            iso_splines::IsotopeSplineModel, 
                            frag_mz_bounds::Tuple{Float32, Float32},
                            block_size::Int64)::Int64 #where {T,U,V,W<:AbstractFloat,I<:Integer}

    NEUTRON = Float64(1.00335)
    fill!(precursor_transmission, zero(Float32))
    precursor_transmission[prec_iso_idx+1] = one(Float32)#
    iso_fac = one(Float32)/(iso_splines(min(Int64(prec_sulfur_count), 5), prec_iso_idx, prec_mz*prec_charge))#one(Float32)
    for frag_idx in precursor_fragment_range
        frag = fragment_ions[frag_idx]
        if frag.rank > 5
            continue
        end
        #Estimate isotope abundances 
        getFragIsotopes!(isotopes, 
                        precursor_transmission,
                        iso_splines, 
                        prec_mz,
                        prec_charge, 
                        prec_sulfur_count,
                        frag,
                        )
        #Mono and mono+1 fragment isotopes 
        for iso_idx in range(0, 2)#prec_iso_idx)
            frag_mz = Float32(frag.mz + iso_idx*NEUTRON/frag.frag_charge)
            if (frag_mz < first(frag_mz_bounds)) |  (frag_mz > last(frag_mz_bounds))
                continue
            end
            transition_idx += 1
            transitions[transition_idx] = DetailedFrag(
                UInt32((frag.prec_id -1)*3 + (iso_idx + 1)), #will turn out odd for M+1 precursor and even for M+0 precursor
                frag_mz, #Estimated isotopic m/z
                Float16(isotopes[iso_idx + 1]*iso_fac), #Estimated relative abundance 
                frag.ion_type,
                frag.is_y,
                frag.is_b,
                frag.is_p,
                iso_idx>0, #Is the fragment an isotope?
                frag.frag_charge,
                frag.ion_position,
                frag.prec_charge,
                frag.rank,
                frag.sulfur_count
            )
            #Grow array if exceeds length
            if transition_idx >= length(transitions)
                append!(transitions, [DetailedFrag{Float32}() for _ in range(1, block_size)])
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