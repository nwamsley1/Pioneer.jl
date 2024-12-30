
"""
    buildPionLib(spec_lib_path::String,
                y_start_index::UInt8,
                y_start::UInt8,
                b_start_index::UInt8,
                b_start::UInt8,
                include_p_index::Bool,
                include_p::Bool,
                include_isotope::Bool,
                include_immonium::Bool,
                include_internal::Bool,
                include_neutral_diff::Bool,
                max_frag_charge::UInt8,
                max_frag_rank::UInt8,
                min_frag_intensity::Float32,
                rank_to_score::Vector{UInt8},
                frag_bounds::FragBoundModel,
                frag_bin_tol_ppm::Float32,
                rt_bin_tol_ppm::Float32,
                model_type::KoinaModelType = InstrumentSpecificModel("default"))

Build Pioneer spectral library including fragment indices.

Parameters:
- spec_lib_path: Path to library directory
- y/b_start_index: Minimum index for index-building
- y/b_start: Minimum index for detailed fragments 
- include_*: Flags for fragment type inclusion
- max_frag_charge: Maximum fragment charge to include
- max_frag_rank: Maximum number of fragments per precursor
- min_frag_intensity: Minimum relative intensity threshold
- rank_to_score: Mapping from rank to scoring value
- frag_bounds: Model for fragment m/z bounds
- frag/rt_bin_tol_ppm: Binning tolerances
- model_type: Type of prediction model used
"""
function buildPionLib(spec_lib_path::String,
                      y_start_index::UInt8,
                      y_start::UInt8,
                      b_start_index::UInt8,
                      b_start::UInt8,
                      include_p_index::Bool,
                      include_p::Bool,
                      include_isotope::Bool,
                      include_immonium::Bool,
                      include_internal::Bool,
                      include_neutral_diff::Bool,
                      max_frag_charge::UInt8,
                      max_frag_rank::UInt8,
                      min_frag_intensity::Float32,
                      rank_to_score::Vector{UInt8},
                      frag_bounds::FragBoundModel,
                      frag_bin_tol_ppm::Float32,
                      rt_bin_tol_ppm::Float32,
                      model_type::KoinaModelType,
                      )
    fragments_table, prec_to_frag, precursors_table = nothing, nothing, nothing
    try
        fragments_table = Arrow.Table(joinpath(spec_lib_path,"fragments_table.arrow"));
        prec_to_frag = Arrow.Table(joinpath(spec_lib_path,"prec_to_frag.arrow"));
        precursors_table = Arrow.Table(joinpath(spec_lib_path,"precursors_table.arrow"));
    catch e 
        @error "could not find library..."
        return nothing
    end

    #Simple fragments that go into the fragment index 
    println("Get index fragments...")
    simple_frags = getSimpleFrags(
        fragments_table[:mz],
        fragments_table[:is_y],
        fragments_table[:is_b],
        fragments_table[:is_p],
        fragments_table[:fragment_index],
        fragments_table[:charge],
        fragments_table[:isotope],
        fragments_table[:is_internal],
        fragments_table[:is_immonium],
        fragments_table[:has_neutral_diff],
        precursors_table[:mz],
        precursors_table[:irt],
        precursors_table[:prec_charge],#precursor_charge],
        prec_to_frag[:start_idx],
        y_start_index,
        b_start_index,
        include_p_index,
        include_isotope,
        include_immonium,
        include_internal,
        include_neutral_diff,
        max_frag_charge,
        frag_bounds,
        rank_to_score
    );

    println("Build fragment index...")
    ##########
    #Builds fragment indexes and saves them for the spec_lib_path
    sort!(simple_frags, by = x->x.prec_irt)
    buildFragmentIndex!(
        spec_lib_path,
        simple_frags,
        frag_bin_tol_ppm, #frag_tol_ppm
        rt_bin_tol_ppm,  #irt_bin_width
        index_name = ""
    );
    println("Build presearch fragment index...")
    sort!(simple_frags, by = x->x.prec_irt)
    buildFragmentIndex!(
        spec_lib_path,
        simple_frags,
        frag_bin_tol_ppm, #frag_tol_ppm
        typemax(Float32),  #irt_tol
        index_name = "presearch_"
    );
    simple_frags = nothing
    GC.gc()

    println("Get full fragments list...")
    detailed_frags, pid_to_fid = getDetailedFrags(
    fragments_table[:mz],
    fragments_table[:intensity],
    fragments_table[:is_y],
    fragments_table[:is_b],
    fragments_table[:is_p],
    fragments_table[:fragment_index],
    fragments_table[:charge],
    fragments_table[:sulfur_count],
    fragments_table[:ion_type],
    fragments_table[:isotope],
    fragments_table[:is_internal],
    fragments_table[:is_immonium],
    fragments_table[:has_neutral_diff],
    precursors_table[:mz],
    precursors_table[:prec_charge],#:precursor_charge],
    prec_to_frag[:start_idx],
    y_start,
    b_start,
    include_p,
    include_isotope,
    include_immonium,
    include_internal,
    include_neutral_diff,
    max_frag_charge,
    frag_bounds,
    max_frag_rank,
    min_frag_intensity,
    model_type
    );
    
    save_detailed_frags(
        joinpath(spec_lib_path, "detailed_fragments.jld2"),
        detailed_frags
    )

    jldsave(
        joinpath(spec_lib_path, "precursor_to_fragment_indices.jld2");
        pid_to_fid
    )

    rm(joinpath(spec_lib_path,"fragments_table.arrow"));
    rm(joinpath(spec_lib_path,"prec_to_frag.arrow"));
    rm(joinpath(spec_lib_path,"precursors.arrow"));

    return nothing
end


"""
    buildPionLib(spec_lib_path::String,
                y_start_index::UInt8,
                y_start::UInt8,
                b_start_index::UInt8,
                b_start::UInt8,
                include_p_index::Bool,
                include_p::Bool,
                include_isotope::Bool,
                include_immonium::Bool,
                include_internal::Bool,
                include_neutral_diff::Bool,
                max_frag_charge::UInt8,
                max_frag_rank::UInt8,
                min_frag_intensity::Float32,
                rank_to_score::Vector{UInt8},
                frag_bounds::FragBoundModel,
                frag_bin_tol_ppm::Float32,
                rt_bin_tol_ppm::Float32,
                model_type::KoinaModelType = InstrumentSpecificModel("default"))

Build Pioneer spectral library including fragment indices.

Parameters:
- spec_lib_path: Path to library directory
- y/b_start_index: Minimum index for index-building
- y/b_start: Minimum index for detailed fragments 
- include_*: Flags for fragment type inclusion
- max_frag_charge: Maximum fragment charge to include
- max_frag_rank: Maximum number of fragments per precursor
- min_frag_intensity: Minimum relative intensity threshold
- rank_to_score: Mapping from rank to scoring value
- frag_bounds: Model for fragment m/z bounds
- frag/rt_bin_tol_ppm: Binning tolerances
- model_type: Type of prediction model used
"""
function buildPionLib(spec_lib_path::String,
                      y_start_index::UInt8,
                      y_start::UInt8,
                      b_start_index::UInt8,
                      b_start::UInt8,
                      include_p_index::Bool,
                      include_p::Bool,
                      include_isotope::Bool,
                      include_immonium::Bool,
                      include_internal::Bool,
                      include_neutral_diff::Bool,
                      max_frag_charge::UInt8,
                      max_frag_rank::UInt8,
                      min_frag_intensity::Float32,
                      rank_to_score::Vector{UInt8},
                      frag_bounds::FragBoundModel,
                      frag_bin_tol_ppm::Float32,
                      rt_bin_tol_ppm::Float32,
                      model_type::SplineCoefficientModel,
                      )
    fragments_table, prec_to_frag, precursors_table = nothing, nothing, nothing
    try
        fragments_table = Arrow.Table(joinpath(spec_lib_path,"fragments_table.arrow"));
        prec_to_frag = Arrow.Table(joinpath(spec_lib_path,"prec_to_frag.arrow"));
        precursors_table = Arrow.Table(joinpath(spec_lib_path,"precursors_table.arrow"));
    catch e 
        @error "could not find library..."
        return nothing
    end

    #Simple fragments that go into the fragment index 
    println("Get index fragments...")
    simple_frags = getSimpleFrags(
        fragments_table[:mz],
        fragments_table[:is_y],
        fragments_table[:is_b],
        fragments_table[:is_p],
        fragments_table[:fragment_index],
        fragments_table[:charge],
        fragments_table[:isotope],
        fragments_table[:is_internal],
        fragments_table[:is_immonium],
        fragments_table[:has_neutral_diff],
        precursors_table[:mz],
        precursors_table[:irt],
        precursors_table[:prec_charge],#precursor_charge],
        prec_to_frag[:start_idx],
        y_start_index,
        b_start_index,
        include_p_index,
        include_isotope,
        include_immonium,
        include_internal,
        include_neutral_diff,
        max_frag_charge,
        frag_bounds,
        rank_to_score
    );

    println("Build fragment index...")
    ##########
    #Builds fragment indexes and saves them for the spec_lib_path
    sort!(simple_frags, by = x->x.prec_irt)
    buildFragmentIndex!(
        spec_lib_path,
        simple_frags,
        frag_bin_tol_ppm, #frag_tol_ppm
        rt_bin_tol_ppm,  #irt_bin_width
        index_name = ""
    );
    println("Build presearch fragment index...")
    sort!(simple_frags, by = x->x.prec_irt)
    buildFragmentIndex!(
        spec_lib_path,
        simple_frags,
        frag_bin_tol_ppm, #frag_tol_ppm
        typemax(Float32),  #irt_tol
        index_name = "presearch_"
    );
    simple_frags = nothing
    GC.gc()

    println("Get full fragments list...")
    detailed_frags, pid_to_fid = getDetailedFrags(
    fragments_table[:mz],
    fragments_table[:coefficients],
    fragments_table[:intensity],
    fragments_table[:is_y],
    fragments_table[:is_b],
    fragments_table[:is_p],
    fragments_table[:fragment_index],
    fragments_table[:charge],
    fragments_table[:sulfur_count],
    fragments_table[:ion_type],
    fragments_table[:isotope],
    fragments_table[:is_internal],
    fragments_table[:is_immonium],
    fragments_table[:has_neutral_diff],
    precursors_table[:mz],
    precursors_table[:prec_charge],#:precursor_charge],
    prec_to_frag[:start_idx],
    y_start,
    b_start,
    include_p,
    include_isotope,
    include_immonium,
    include_internal,
    include_neutral_diff,
    max_frag_charge,
    frag_bounds,
    max_frag_rank,
    min_frag_intensity,
    model_type
    );
    
    save_detailed_frags(
        joinpath(spec_lib_path, "detailed_fragments.jld2"),
        detailed_frags
    )

    jldsave(
        joinpath(spec_lib_path, "precursor_to_fragment_indices.jld2");
        pid_to_fid
    )

    #rm(joinpath(spec_lib_path,"fragments_table.arrow"));
    #rm(joinpath(spec_lib_path,"prec_to_frag.arrow"));
    #rm(joinpath(spec_lib_path,"precursors.arrow"));
    

    return nothing
end

function fragFilter(
    frag_is_y::Bool,
    frag_is_b::Bool,
    frag_is_p::Bool,
    frag_index::UInt8,
    frag_charge::UInt8,
    frag_isotope::UInt8,
    frag_internal::Bool,
    frag_immonium::Bool,
    frag_neutral_diff::Bool,
    frag_mz::Float32,
    frag_bounds::FragBoundModel,
    prec_mz::Float32,
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8)
    
    #println("frag_neutral_diff $frag_neutral_diff")
    min_frag_mz, max_frag_mz = frag_bounds(prec_mz)
    if (frag_mz < min_frag_mz) | (frag_mz > max_frag_mz)
        return false
    end
    if frag_is_y
        if frag_index < y_start
            return false
        end
    end
    if frag_is_b
        if frag_index < b_start
            return false
        end
    end
    if frag_is_p
        if !include_p
            return false
        end
    end
    if frag_immonium
        if !include_immonium
            return false
        end
    end
    if frag_internal
        if !include_internal
            return false
        end
    end
    if frag_neutral_diff
        if !include_neutral_diff
            return false
        end
    end
    if include_isotope
        if !iszero(frag_isotope)
            return false
        end
    end
    if frag_index < b_start
        return false
    end
    if frag_charge > max_frag_charge
        return false
    end
    return true
end

function getSimpleFrags(
    frag_mz::AbstractVector{Float32},
    frag_is_y::AbstractVector{Bool},
    frag_is_b::AbstractVector{Bool},
    frag_is_p::AbstractVector{Bool},
    frag_index::AbstractVector{UInt8},
    frag_charge::AbstractVector{UInt8},
    frag_isotope::AbstractVector{UInt8},
    frag_internal::AbstractVector{Bool},
    frag_immonium::AbstractVector{Bool},
    frag_neutral_diff::AbstractVector{Bool},
    precursor_mz::AbstractVector{Float32},
    precursor_irt::AbstractVector{Float32},
    precursor_charge::AbstractVector{UInt8},
    prec_to_frag_idx::AbstractVector{UInt64},
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8,
    frag_bounds::FragBoundModel,
    rank_to_score::Vector{UInt8},
    )
    if (length(prec_to_frag_idx) - 1) != (length(precursor_mz))
        println("mistake")
    end
    #Maximum ranked fragment that can be included in the fragment index
    max_rank_index = length(rank_to_score)
    #Number of precursors 
    n_precursors = UInt32(length(precursor_mz))
    simple_frags = Vector{SimpleFrag{Float32}}(undef, n_precursors*max_rank_index)
    simple_frag_idx = 0
    for pid in range(one(UInt32), n_precursors)
        prec_mz = precursor_mz[pid]
        frag_start_idx, frag_stop_idx = prec_to_frag_idx[pid], prec_to_frag_idx[pid+1] - 1
        rank = 1
        for frag_idx in range(frag_start_idx, frag_stop_idx)
            if fragFilter(
                    frag_is_y[frag_idx],
                    frag_is_b[frag_idx],
                    frag_is_p[frag_idx],
                    frag_index[frag_idx],
                    frag_charge[frag_idx],
                    frag_isotope[frag_idx],
                    frag_internal[frag_idx],
                    frag_immonium[frag_idx],
                    frag_neutral_diff[frag_idx],
                    frag_mz[frag_idx],
                    frag_bounds,
                    prec_mz,
                    y_start,
                    b_start,
                    include_p,
                    include_isotope,
                    include_immonium,
                    include_internal,
                    include_neutral_diff,
                    max_frag_charge)==false
                continue
            end
            simple_frag_idx += 1
            simple_frags[simple_frag_idx] = SimpleFrag(
                frag_mz[frag_idx],
                pid,
                precursor_mz[pid],
                precursor_irt[pid],
                precursor_charge[pid],
                rank_to_score[rank]
            )
            rank += 1
            if rank > max_rank_index
                break
            end
        end

    end
    return simple_frags[1:simple_frag_idx]
end

function buildFragmentIndex!(
                            folder_out::String,
                            frag_ions::Vector{SimpleFrag{T}}, 
                            frag_bin_tol_ppm::AbstractFloat, 
                            rt_bin_tol::AbstractFloat;
                            index_name::String = ""
                            ) where {T<:AbstractFloat}

    function buildFragIndex!(
        index_fragments::Vector{IndexFragment{T}}, 
        rt_bins::Vector{FragIndexBin{T}},
        frag_bins::Vector{FragIndexBin{T}},
        frag_ions::Vector{SimpleFrag{T}}, 
        frag_bin_tol_ppm::AbstractFloat,
        rt_bin_tol::AbstractFloat) where {T<:AbstractFloat}

        start_irt, stop_irt = getIRT(first(frag_ions)), getIRT(first(frag_ions))
        start_idx, stop_idx = 1, 1
        rt_bin_idx = 1
        frag_bin_idx = 1
        #Within each iRT bin (defined by `rt_size`) sort the fragments by precursor_mz

        for (i, frag_ion) in enumerate(frag_ions)
            stop_irt = getIRT(frag_ion)
            stop_idx = i
            #diff_mz = stop_fragmz - start_fragmz
            #mean_mz = (stop_fragmz + start_fragmz)/2
            #Is the difference between the first and last fragment greater than the bin_ppm?
            if ((stop_irt - start_irt) > rt_bin_tol) & (stop_idx != start_idx)#(diff_mz/(mean_mz/1e6) > frag_bin_tol_ppm) & (stop_idx != start_idx)
                stop_idx = i - 1 #i - 1 fragment is the last that should be incluced in the bin
                stop_irt = getIRT(frag_ions[stop_idx])
                #Within the fragment bin Sort by IRT
                sort!(@view(frag_ions[start_idx:stop_idx]), by = x->getMZ(x), alg=QuickSort)
                #Add new fragbin
                #Build fragment bins for the retention time bin 
                first_frag_bin_idx = frag_bin_idx
                frag_bin_idx = buildFragBins!(index_fragments,
                            frag_bins,
                            frag_bin_idx,
                            frag_ions,
                            start_idx,stop_idx,
                            frag_bin_tol_ppm)
                rt_bins[rt_bin_idx] = FragIndexBin(start_irt, 
                                                    stop_irt, #important that stop_idx is i - 1 and not i
                                                        UInt32(first_frag_bin_idx),
                                                        UInt32(frag_bin_idx-1)
                                                        ) #-1 is critical
                rt_bin_idx += 1

                start_idx, stop_idx = i, i
                start_irt = getIRT(frag_ion)
            end
        end


        #Last bin is special case 
        if start_idx != length(frag_ions)
            stop_irt =  getIRT(frag_ions[stop_idx])
            sort!(@view(frag_ions[start_idx:stop_idx]), by = x->getMZ(x), alg=QuickSort)
            #Add new fragbin
            #Build RT bins for the frag bin
            first_frag_bin_idx = frag_bin_idx
            frag_bin_idx = buildFragBins!(index_fragments,
                        frag_bins,
                        frag_bin_idx,
                        frag_ions,
                        start_idx,stop_idx,frag_bin_tol_ppm)

            rt_bins[rt_bin_idx] = FragIndexBin(start_irt, 
                                                    stop_irt, #important that stop_idx is i - 1 and not i
                                                    UInt32(first_frag_bin_idx),
                                                    UInt32(frag_bin_idx-1)) #-1 is critical
            rt_bin_idx += 1
        else
            first_frag_bin_idx = frag_bin_idx

            frag_bin_idx = buildFragBins!(index_fragments,
                        frag_bins,
                        frag_bin_idx,
                        frag_ions,
                        start_idx,stop_idx,frag_bin_tol_ppm)

            #Add new fragbin
            rt_bins[rt_bin_idx] = FragIndexBin(start_irt, 
                                                    getIRT(frag_ions[stop_idx]), #important that stop_idx is i - 1 and not i
                                                    UInt32(first_frag_bin_idx),
                                                    UInt32(frag_bin_idx-1)) #-1 is critical
            rt_bin_idx += 1
        end

        return frag_bin_idx, rt_bin_idx
    end

    function buildFragBins!(index_fragments::Vector{IndexFragment{T}}, 
                            frag_bins::Vector{FragIndexBin{T}},
                            frag_bin_idx::Int64,
                            frag_ions::Vector{SimpleFrag{T}}, 
                            start::Int64, 
                            stop::Int64, 
                            frag_bin_tol_ppm::AbstractFloat) where {T<:AbstractFloat}
        start_idx, stop_idx = start, start
        start_fragmz, stop_fragmz = getMZ(frag_ions[start]), getMZ(frag_ions[start])
        for i in range(start, stop)
            frag_ion = frag_ions[i]
            stop_fragmz = getMZ(frag_ion)
            stop_idx = i

            diff_mz = stop_fragmz - start_fragmz
            mean_mz = (stop_fragmz + start_fragmz)/2
            if (diff_mz/(mean_mz/1e6) > frag_bin_tol_ppm) & (stop_idx != start_idx)
                stop_idx = i - 1
                stop_fragmz =  getMZ(frag_ions[stop_idx]) #Need to set before sorting 
                sort!(@view(frag_ions[start_idx:stop_idx]), by = x->getPrecMZ(x), alg=QuickSort)
                #Add new rt bin
                frag_bins[frag_bin_idx] = FragIndexBin(start_fragmz, 
                                                        stop_fragmz, #important that stop_idx is i - 1 and not i
                                                    UInt32(start_idx),
                                                    UInt32(stop_idx)
                                                )
                frag_bin_idx += 1
                for idx in range(start_idx, stop_idx)
                    index_fragments[idx] = IndexFragment(
                                                    getPrecID(frag_ions[idx]),
                                                    getPrecMZ(frag_ions[idx]),
                                                    getScore(frag_ions[idx]),
                                                    getPrecCharge(frag_ions[idx])
                                                    )
                end
                start_idx, stop_idx = i, i
                start_fragmz = getMZ(frag_ions[stop_idx])
            end
        end

        #Last bin is special case 
        if start_idx != stop
            stop_fragmz = getMZ(frag_ions[stop_idx])
            sort!(@view(frag_ions[start_idx:stop_idx]), by = x->getPrecMZ(x), alg=QuickSort)
            #Add new fragbin
            frag_bins[frag_bin_idx] = FragIndexBin(start_fragmz, 
                        stop_fragmz, #important that stop_idx is i - 1 and not i
                        UInt32(start_idx),
                        UInt32(stop_idx)
                    )
            frag_bin_idx += 1
            for idx in range(start_idx, stop_idx)
                index_fragments[idx] = IndexFragment(
                                                getPrecID(frag_ions[idx]),
                                                getPrecMZ(frag_ions[idx]),
                                                getScore(frag_ions[idx]),
                                                getPrecCharge(frag_ions[idx])
                                                )
            end
        else
            frag_bins[frag_bin_idx] = FragIndexBin(start_fragmz, 
                        getMZ(frag_ions[stop]), #important that stop_idx is i - 1 and not i
                        UInt32(start_idx),
                        UInt32(stop_idx)
                    )
            frag_bin_idx += 1
            index_fragments[stop_idx] = IndexFragment(
                getPrecID(frag_ions[stop_idx]),
                getPrecMZ(frag_ions[stop_idx]),
                getScore(frag_ions[stop_idx]),
                getPrecCharge(frag_ions[stop_idx])
                )
        end
        return frag_bin_idx
    end

    sort!(frag_ions, by = x->getIRT(x), alg=QuickSort)
    #diff = getPPM(getMZ(frag_ions[start]), bin_ppm) #ppm tolerance of the current fragment bin
    #Get smallest iRT in the library
    index_fragments = Vector{IndexFragment{T}}(undef, length(frag_ions))
    rt_bins = Vector{FragIndexBin{T}}(undef, length(frag_ions))
    frag_bins = Vector{FragIndexBin{T}}(undef, length(frag_ions))
    println("building fragment index...")
    frag_bin_idx, rt_bin_idx = buildFragIndex!(index_fragments,
                    rt_bins,
                    frag_bins,
                    frag_ions,
                    frag_bin_tol_ppm,
                    rt_bin_tol)
                    
    fragments = (IndexFragment = index_fragments,)
    rt_bins  = (FragIndexBin = rt_bins[1:rt_bin_idx-1],)
    frag_bins = (FragIndexBin = frag_bins[1:frag_bin_idx-1],)
    println("writing tables...")
    if !isdir(folder_out)
        mkdir(folder_out)
    end
    Arrow.write(joinpath(folder_out, index_name*"f_index_fragments.arrow"), fragments)
    Arrow.write(joinpath(folder_out, index_name*"f_index_rt_bins.arrow"), rt_bins)
    Arrow.write(joinpath(folder_out, index_name*"f_index_fragment_bins.arrow"), frag_bins)

    return 
end

function getDetailedFrags(
    frag_mz::AbstractVector{Float32},
    frag_intensity::AbstractVector{Float16},
    frag_is_y::AbstractVector{Bool},
    frag_is_b::AbstractVector{Bool},
    frag_is_p::AbstractVector{Bool},
    frag_index::AbstractVector{UInt8},
    frag_charge::AbstractVector{UInt8},
    frag_sulfur_count::AbstractVector{UInt8},
    frag_ion_type::AbstractVector{UInt16},
    frag_isotope::AbstractVector{UInt8},
    frag_internal::AbstractVector{Bool},
    frag_immonium::AbstractVector{Bool},
    frag_neutral_diff::AbstractVector{Bool},
    precursor_mz::AbstractVector{Float32},
    precursor_charge::AbstractVector{UInt8},
    prec_to_frag_idx::AbstractVector{UInt64},
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8,
    frag_bounds::FragBoundModel,
    max_frag_rank::UInt8,
    min_frag_intensity::AbstractFloat,
    koina_model::KoinaModelType)
    println("in this function yo")
    if (length(prec_to_frag_idx) - 1) != (length(precursor_mz))
        println("mistake")
    end

    #Maximum ranked fragment that can be included in the fragment index
    max_frag_rank = UInt8(min(255, max_frag_rank))
    #Number of precursors 
    n_precursors = UInt32(length(precursor_mz))
    #Keep track of number of fragments to allocate 
    n_frags = zero(UInt64)
    println("counting fragments...")
    for pid in ProgressBar(range(one(UInt32), n_precursors))
        prec_mz = precursor_mz[pid] #Filter on precursor mass
        frag_start_idx, frag_stop_idx = prec_to_frag_idx[pid], prec_to_frag_idx[pid+1] - 1
        rank = 1
        #count fragments for the current precursor
        for frag_idx in range(frag_start_idx, frag_stop_idx)
            #Filter on fragment properties
            if min_frag_intensity > frag_intensity[frag_idx]
                continue
            end
            if !fragFilter(
                    frag_is_y[frag_idx],
                    frag_is_b[frag_idx],
                    frag_is_p[frag_idx],
                    frag_index[frag_idx],
                    frag_charge[frag_idx],
                    frag_isotope[frag_idx],
                    frag_internal[frag_idx],
                    frag_immonium[frag_idx],
                    frag_neutral_diff[frag_idx],
                    frag_mz[frag_idx],
                    frag_bounds,
                    prec_mz,
                    y_start,
                    b_start,
                    include_p,
                    include_isotope,
                    include_immonium,
                    include_internal,
                    include_neutral_diff,
                    max_frag_charge)
                continue
            end
            #update counters 
            n_frags += one(UInt64)
            rank += 1
            if rank > max_frag_rank
                break
            end
        end
    end
    println("n_frags $n_frags")
    detailed_frags = Vector{DetailedFrag{Float32}}(
                                    undef, 
                                    n_frags)   
    prec_to_frag_idx_new = Vector{UInt64}(undef, n_precursors + 1)
    detailed_frag_idx = 1
    println("writing fragments...")
    for pid in ProgressBar(range(one(UInt32), n_precursors))
        prec_mz = precursor_mz[pid]
        #Index of the first fragment for the precursor 
        prec_to_frag_idx_new[pid] = UInt64(detailed_frag_idx)
        frag_start_idx, frag_stop_idx = prec_to_frag_idx[pid], prec_to_frag_idx[pid+1] - 1
        rank = 1
        for frag_idx in range(frag_start_idx, frag_stop_idx)
            if !fragFilter(
                    frag_is_y[frag_idx],
                    frag_is_b[frag_idx],
                    frag_is_p[frag_idx],
                    frag_index[frag_idx],
                    frag_charge[frag_idx],
                    frag_isotope[frag_idx],
                    frag_internal[frag_idx],
                    frag_immonium[frag_idx],
                    frag_neutral_diff[frag_idx],
                    frag_mz[frag_idx],
                    frag_bounds,
                    prec_mz,
                    y_start,
                    b_start,
                    include_p,
                    include_isotope,
                    include_immonium,
                    include_internal,
                    include_neutral_diff,
                    max_frag_charge)
                continue
            end
            if min_frag_intensity > frag_intensity[frag_idx]
                continue
            end
            is_y, is_internal, is_immonium = frag_is_y[frag_idx], frag_internal[frag_idx], frag_immonium[frag_idx]
            is_b, is_p = frag_is_b[frag_idx], frag_is_p[frag_idx]
            detailed_frags[detailed_frag_idx] = DetailedFrag(
                UInt32(pid),

                frag_mz[frag_idx],
                frag_intensity[frag_idx],

                frag_ion_type[frag_idx],
                is_y,
                is_b,
                is_p,
                frag_isotope[frag_idx]>0,

                frag_charge[frag_idx],
                frag_index[frag_idx],
                precursor_charge[pid],
                UInt8(rank),
                frag_sulfur_count[frag_idx]
            )
            detailed_frag_idx += 1
            rank += 1
            if rank > max_frag_rank
                break
            end
        end

    end
    prec_to_frag_idx_new[end] = UInt64(detailed_frag_idx)
    return detailed_frags[1:detailed_frag_idx-1], prec_to_frag_idx_new
end

function getDetailedFrags(
    frag_mz::AbstractVector{Float32},
    frag_coef::AbstractVector{NTuple{N, Float32}},
    frag_intensity::AbstractVector{Float16},
    frag_is_y::AbstractVector{Bool},
    frag_is_b::AbstractVector{Bool},
    frag_is_p::AbstractVector{Bool},
    frag_index::AbstractVector{UInt8},
    frag_charge::AbstractVector{UInt8},
    frag_sulfur_count::AbstractVector{UInt8},
    frag_ion_type::AbstractVector{UInt16},
    frag_isotope::AbstractVector{UInt8},
    frag_internal::AbstractVector{Bool},
    frag_immonium::AbstractVector{Bool},
    frag_neutral_diff::AbstractVector{Bool},
    precursor_mz::AbstractVector{Float32},
    precursor_charge::AbstractVector{UInt8},
    prec_to_frag_idx::AbstractVector{UInt64},
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8,
    frag_bounds::FragBoundModel,
    max_frag_rank::UInt8,
    min_frag_intensity::AbstractFloat,
    koina_model::SplineCoefficientModel) where {N}

    if (length(prec_to_frag_idx) - 1) != (length(precursor_mz))
        println("mistake")
    end

    #Maximum ranked fragment that can be included in the fragment index
    max_frag_rank = UInt8(min(255, max_frag_rank))
    #Number of precursors 
    n_precursors = UInt32(length(precursor_mz))
    #Keep track of number of fragments to allocate 
    n_frags = zero(UInt64)
    println("counting fragments...")
    for pid in ProgressBar(range(one(UInt32), n_precursors))
        prec_mz = precursor_mz[pid] #Filter on precursor mass
        frag_start_idx, frag_stop_idx = prec_to_frag_idx[pid], prec_to_frag_idx[pid+1] - 1
        rank = 1
        #count fragments for the current precursor
        for frag_idx in range(frag_start_idx, frag_stop_idx)
            #Filter on fragment properties
            if !fragFilter(
                    frag_is_y[frag_idx],
                    frag_is_b[frag_idx],
                    frag_is_p[frag_idx],
                    frag_index[frag_idx],
                    frag_charge[frag_idx],
                    frag_isotope[frag_idx],
                    frag_internal[frag_idx],
                    frag_immonium[frag_idx],
                    frag_neutral_diff[frag_idx],
                    frag_mz[frag_idx],
                    frag_bounds,
                    prec_mz,
                    y_start,
                    b_start,
                    include_p,
                    include_isotope,
                    include_immonium,
                    include_internal,
                    include_neutral_diff,
                    max_frag_charge)
                continue
            end
            #update counters 
            n_frags += one(UInt64)
            rank += 1
            if rank > max_frag_rank
                break
            end
        end
    end
    println("n_frags $n_frags")
    n_tuple_p = typeof(first(frag_coef))
    n_tuple_size = length(n_tuple_p.parameters)
    n_tuple_type = eltype(n_tuple_p)
    detailed_frags = Vector{SplineDetailedFrag{n_tuple_size, n_tuple_type}}(
                                    undef, 
                                    n_frags)   
    prec_to_frag_idx_new = Vector{UInt64}(undef, n_precursors + 1)
    detailed_frag_idx = 1
    println("writing fragments...")
    for pid in ProgressBar(range(one(UInt32), n_precursors))
        prec_mz = precursor_mz[pid]
        #Index of the first fragment for the precursor 
        prec_to_frag_idx_new[pid] = UInt64(detailed_frag_idx)
        frag_start_idx, frag_stop_idx = prec_to_frag_idx[pid], prec_to_frag_idx[pid+1] - 1
        rank = 1
        for frag_idx in range(frag_start_idx, frag_stop_idx)
            if !fragFilter(
                    frag_is_y[frag_idx],
                    frag_is_b[frag_idx],
                    frag_is_p[frag_idx],
                    frag_index[frag_idx],
                    frag_charge[frag_idx],
                    frag_isotope[frag_idx],
                    frag_internal[frag_idx],
                    frag_immonium[frag_idx],
                    frag_neutral_diff[frag_idx],
                    frag_mz[frag_idx],
                    frag_bounds,
                    prec_mz,
                    y_start,
                    b_start,
                    include_p,
                    include_isotope,
                    include_immonium,
                    include_internal,
                    include_neutral_diff,
                    max_frag_charge)
                continue
            end
            is_y, is_internal, is_immonium = frag_is_y[frag_idx], frag_internal[frag_idx], frag_immonium[frag_idx]
            is_b, is_p = frag_is_b[frag_idx], frag_is_p[frag_idx]
            detailed_frags[detailed_frag_idx] = SplineDetailedFrag(
                UInt32(pid),

                frag_mz[frag_idx],
                frag_coef[frag_idx],

                frag_ion_type[frag_idx],
                is_y,
                is_b,
                is_p,
                frag_isotope[frag_idx]>0,

                frag_charge[frag_idx],
                frag_index[frag_idx],
                precursor_charge[pid],
                UInt8(rank),
                frag_sulfur_count[frag_idx]
            )
            detailed_frag_idx += 1
            rank += 1
            if rank > max_frag_rank
                break
            end
        end

    end
    prec_to_frag_idx_new[end] = UInt64(detailed_frag_idx)
    return detailed_frags[1:detailed_frag_idx-1], prec_to_frag_idx_new
end
