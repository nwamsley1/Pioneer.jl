"""
    sortPrositLib(prec_mz::Vector{T},prec_irt::Vector{T},rt_bin_tol::AbstractFloat) where {T<:AbstractFloat}

Given a list of precursor m/z ratios, corresponding precursor indexed retention times, and a retention time bin tolerance
does a multi level sorts and retruns the indices of the sorted permutation. 
    1) Sort precursors by retention time and divide into retention time bins. 
    2) Within each retention time bin sort by precursor mz
    3) Filter precursors by the min/max precursor mz and remove indices from 
    the sort that correspond to these precursors
    4) Return list of inices correponding to the correct order 

- `prec_mz::Vector{T}` -- N length vector of precursor mz's
- `prec_irt::Vector{T} -- N length vector of precursor retention times (arbitrary scale)
- `rt_bin_tol::AbstractFloat` -- width of the retention time bins

"""
function sortPrositLib(
                        prec_mz::AbstractArray{T},
                        prec_irt::AbstractArray{T},
                        rt_bin_tol::AbstractFloat
                        ) where {T<:AbstractFloat}

    #Precursor indices sorted by retention time 
    sorted_indices = sortperm(prec_irt)
    #First and last indices of the current irt bin 
    start_irt_idx, stop_irt_idx = 1, 1
    start_irt, stop_irt = prec_irt[start_irt_idx], prec_irt[stop_irt_idx]
    for i in range(1, length(prec_irt))
        stop_irt_idx = i
        stop_irt = prec_irt[stop_irt_idx]
        #If the bin width has been exceeded, step back one index and sort the bin
        #by precusor mz
        if ((stop_irt - start_irt) > rt_bin_tol) & (start_irt_idx != stop_irt_idx)
            sort!(
                    @view(sorted_indices[start_irt_idx:(stop_irt_idx - 1)]),
                    by = x->prec_mz[x],
                    alg=QuickSort
                )
            #Get first index and rt of the next bin
            start_irt_idx = i
            start_irt = prec_irt[start_irt_idx]
        end
    end
    #Sort the final irt bin by precursor m/z 
    sort!(@view(sorted_indices[start_irt_idx:(stop_irt_idx)]),
            by = x->prec_mz[x],
            alg=QuickSort)
    return sorted_indices
end
#=
 prec_frags::Arrow.List{SubArray{PrositFrag, 1, 
                                    Arrow.Struct{PrositFrag, Tuple{Arrow.Primitive{Float32, Vector{Float32}}, 
                                    Arrow.Primitive{Float32, Vector{Float32}}, 
                                    Arrow.Primitive{Char, Vector{UInt32}}, 
                                    Arrow.Primitive{UInt8, Vector{UInt8}}, 
                                    Arrow.Primitive{UInt8, Vector{UInt8}}, 
                                    Arrow.Primitive{UInt8, Vector{UInt8}}}, (:intensity, :mass, :type, :index, :charge, :sulfur_count)}, Tuple{UnitRange{Int64}}, true}, Int32, 
                                        Arrow.Struct{PrositFrag, Tuple{Arrow.Primitive{Float32, Vector{Float32}}, 
                                        Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{Char, Vector{UInt32}}, 
                                        Arrow.Primitive{UInt8, Vector{UInt8}},
                                         Arrow.Primitive{UInt8, Vector{UInt8}}, 
                                         Arrow.Primitive{UInt8, Vector{UInt8}}}, (:intensity, :mass, :type, :index, :charge, :sulfur_count)}}
=#
"""
    parsePioneerLib(
        folder_out::String,
        modified_sequence::AbstractArray{String},
        accession_numbers::AbstractArray{String},
        decoy::AbstractArray{Bool},
        charge::AbstractArray{UInt8},
        irt::AbstractArray{Float32},
        prec_mz::AbstractArray{Float32},
        prec_frags::Arrow.List{SubArray{PrositFrag, 1, 
                    Arrow.Struct{PrositFrag, Tuple{Arrow.Primitive{Float32, Vector{Float32}}, 
                    Arrow.Primitive{Float32, Vector{Float32}}, 
                    Arrow.Primitive{Char, Vector{UInt32}}, 
                    Arrow.Primitive{UInt8, Vector{UInt8}}, 
                    Arrow.Primitive{UInt8, Vector{UInt8}}, 
                    Arrow.Primitive{UInt8, Vector{UInt8}}}, (:intensity, :mass, :type, :index, :charge, :sulfur_count)}, Tuple{UnitRange{Int64}}, true}, Int32, 
                        Arrow.Struct{PrositFrag, Tuple{Arrow.Primitive{Float32, Vector{Float32}}, 
                        Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{Char, Vector{UInt32}}, 
                        Arrow.Primitive{UInt8, Vector{UInt8}},
                        Arrow.Primitive{UInt8, Vector{UInt8}}, 
                        Arrow.Primitive{UInt8, Vector{UInt8}}}, (:intensity, :mass, :type, :index, :charge, :sulfur_count)}},
        frag_mz_min::T,
        frag_mz_max::T,
        prec_mz_min::T,
        prec_mz_max::T;
        rank_to_score::Vector{UInt8}=UInt8[4, 2, 2, 1, 1],
        rt_bin_tol::T = 1.0f0,
        y_start_index::Int64 = 4, 
        b_start_index::Int64 = 3,
        y_start::Int64 = 3, 
        b_start::Int64 = 2,
        missed_cleavage_regex::Regex = r"[KR][^P|\$]") where {T<:AbstractFloat}

Converts a pioneer formated spectral library to a fragment index according to specified parameters.

Outputs: Writes the following to `folder_out/pioneer_lib`
- `precursor_table.arrow` -- Arrow table with one row for each precursor 
    Arrow.Table with 12412074 rows, 10 columns, and schema:
    :irt                Float32
    :mz                 Float32
    :is_decoy           Bool
    :accession_numbers  String
    :sequence           String
    :prec_charge        UInt8
    :missed_cleavages   UInt8
    :variable_mods      UInt8
    :length             UInt8
    :sulfur_count       UInt8
- `simple_fragments.jld2` -- JLD2 formated Vector{SimpleFrag}
- `precursor_to_fragment_indices.jld2` -- JLD2 formated Vector{UnitRange{UInt32}}
    Each range is a range of indices of `detailed_fragments` corresponding to the complete framgent
    list for a precursor
- `detailed_fragments.jld2` -- JLD2 formated Vector{DetailedFrag} 

returns the following 
- `f_simp::Vector{SimpleFrag}` -- used to generate the fragment binary search index in a subsequent step 
- `folder_out::String` -- path to the folder to which the outputs are written 

Arguments: 
- `folder_out::String` -- N length vector of precursor mz's
- `modified_sequence::AbstractArray{String}` -- Sequence annotated with UniMod IDs: "PEPTIDEM(Unimod:4)K"
- `accession_numbers::AbstractArray{String}` -- Uniprot accession numbers for all proteins the precursor maps to. semi-colon delimited: "Q7Z572;Q7Z2Y5"
- `decoy::AbstractArray{Bool}` -- true if the precursor is a decoy and false otherwise 
- `charge::AbstractArray{UInt8}` -- precursor charge state
- `irt::AbstractArray{Float32}` -- library retention time 
- `irt::AbstractArray{Float32}` -- library retention time 
- `prec_mz::AbstractArray{Float32}` -- Sequence annotated with UniMod IDs: "PEPTIDEM(Unimod:4)K"
- `prec_frags::Arrow.List{}` -- Essentially Vector{Vector{PrositFrag}}. One subarray containing the fragments for each precursor. 
        prec_frags[1]
        17-element view(::Arrow.Struct{PrositFrag, Tuple{Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{Float32, Vector{Float32}}, Arrow.Primitive{Char, Vector{UInt32}}, Arrow.Primitive{UInt8, Vector{UInt8}}, Arrow.Primitive{UInt8, Vector{UInt8}}, Arrow.Primitive{UInt8, Vector{UInt8}}}, (:intensity, :mass, :type, :index, :charge, :sulfur_count)}, 1:17) with eltype PrositFrag:
        PrositFrag(1.0f0, 147.11281f0, 'y', 0x01, 0x01, 0x00)
        PrositFrag(0.32908705f0, 204.13426f0, 'y', 0x02, 0x01, 0x00)
        PrositFrag(0.21963413f0, 261.15573f0, 'y', 0x03, 0x01, 0x00)
        ⋮
        PrositFrag(0.00079682295f0, 96.04305f0, 'b', 0x05, 0x03, 0x00)
        PrositFrag(0.0005110648f0, 115.0502f0, 'b', 0x06, 0x03, 0x00)
        PrositFrag(0.0013098607f0, 134.05736f0, 'b', 0x07, 0x03, 0x00)
- `frag_mz_min::T` -- Exclude fragments with a lower m/z 
- `frag_mz_max::T` -- Exclude fragments with a higher m/z
- `prec_mz_min::T` -- Exclude precursors with a lower m/z 
- `prec_mz_max::T` -- Exclude precursors with a higher m/z
- `rank_to_score::Vector{UInt8}=UInt8[4, 2, 2, 1, 1]` -- fragments are ranked in descending order by intensity. 
    this vector maps fragment rank to fragment score in the fragment index. Fragments with rank higher than the 
    length of this vector are not included in the fragment index. 
- `rt_bin_tol::T=1.0f0` -- Size of the library retention time bins. 
- `y_start_index::Int64 = 4` -- excludes y ions from the fragment index before this index (y1, y2, y3 excluded in this example) 
- `b_start_index::Int64 = 3` -- see above 
- `y_start::Int64 = 3` -- excludes y ions from the complete fragment list before this index (y1, y2 excluded in this example)
- `b_start::Int64 = 2` -- excludes b ions from the complete fragment list before this index (y1, y2 excluded in this example)
- `missed_cleavage_regex::Regex = r"[KR][^P|\$]"` -- set to recognize missed tryptic cleavages. Failed cleavage before a proline doesn't count (but really it should?)
"""
function parsePioneerLib(
                        folder_out::String,
                        id_to_annotation::Vector{PioneerFragAnnotation},
                        irt::AbstractArray{Float32},
                        prec_mz::AbstractArray{Float32},
                        decoy::AbstractArray{Bool},

                        frag_mz::AbstractArray{Float32},
                        frag_intensity::AbstractArray{Float16},
                        ion_types::AbstractArray{UInt16},
                        frag_is_y::AbstractArray{Bool},
                        frag_index::AbstractArray{UInt8},
                        frag_charge::AbstractArray{UInt8},
                        frag_isotope::AbstractArray{UInt8},
                        frag_internal::AbstractArray{Bool},
                        frag_immonium::AbstractArray{Bool},
                        frag_internal_ind::AbstractArray{Tuple{UInt8, UInt8}},
                        frag_sulfur_count::AbstractArray{UInt8},

                        proteome_names::AbstractArray{Union{Missing, String}},#Arrow.List{Union{Missing, String}, Int32, Vector{UInt8}},
                        accession_numbers::AbstractArray{Union{Missing, String}},
                        sequence::AbstractArray{Union{Missing, String}},
                        structural_mods::AbstractArray{Union{Missing, String}},
                        isotopic_mods::AbstractArray{Union{Missing, String}},

                        charge::AbstractArray{UInt8},
                        sulfur_count::AbstractArray{UInt8},
                        sequence_length::AbstractArray{UInt8},
                        prec_frag_ranges::Vector{UnitRange{UInt32}},
                        #prec_frags::Vector{PioneerFrag},
                        frag_bounds::FragBoundModel,
                        prec_mz_min::T,
                        prec_mz_max::T;
                        rank_to_score::Vector{UInt8}=UInt8[4, 2, 2, 1, 1],
                        rt_bin_tol::T = 1.0f0,
                        exclude_from_index::Set{Char} = Set{Char}(),
                        y_start_index::Int64 = 4, 
                        b_start_index::Int64 = 3,
                        y_start::Int64 = 3, 
                        b_start::Int64 = 2,
                        missed_cleavage_regex::Regex = r"[KR][^U|\$]",
                        ) where {T<:AbstractFloat}

    max_rank_index = length(rank_to_score)

    println("sorting precursors...")
    @time sorted_indices = sortPrositLib(
        prec_mz,
        irt,
        rt_bin_tol
    )

    detailed_fragments_path = joinpath(folder_out, "detailed_fragments.jld2")
    precursor_to_fragment_indices_path = joinpath(folder_out, "precursor_to_fragment_indices.jld2")
    simple_fragments_path = joinpath(folder_out, "simple_fragments.arrow")
    precursor_table_path = joinpath(folder_out, "precursor_table.arrow")
    
    println("allocating memory for output...")            
    N_precs = length(sorted_indices)
    #How many fragments? Determines size of pre-allocated arrays
    frag_count = length(frag_mz)
    ###########
    #Initialize Containers
    println("frag_count $frag_count")
    frags_detailed = Vector{DetailedFrag{Float32}}(undef, frag_count)
    #Keeps track of with fragments in "frags_detailed" correspond to which precursor
    precursor_indices = Vector{UnitRange{UInt32}}(undef, N_precs)
    #Detailed specification for each precursor. 
    precursors = Vector{LibraryPrecursorIon{Float32}}(undef, N_precs)
    #Used to fuild the framgment index for MSFragger-like searching. 
    frags_simple = Vector{SimpleFrag{Float32}}(undef, N_precs*max_rank_index)
    #maps precursor ids to the original row of the data framgent
    prec_id_to_library_row = Vector{UInt32}(undef, N_precs)
    #Is the fragment within the constraints? m/z scan range and exceeding minimum b,y ion index?
    function inScanRange(frag::PioneerFrag, ion_type::Char, low_frag_mz::T, high_frag_mz::T, y_start_ind::Int64, b_start_ind::Int64) where {T<:AbstractFloat}
        mz, index = getMZ(frag), getIndex(frag)
        if ion_type == 'p'
            return true
        elseif ion_type == 'y'
            return (mz>low_frag_mz)&(mz<high_frag_mz)&(index>=y_start_ind)
        else
            return (mz>low_frag_mz)&(mz<high_frag_mz)&(index>=b_start_ind)
        end
    end

    #Given a vector of fragments for a precursor `Vector{PrositFrag}}`
    #Write the top n (determined by max_rank_index)
    #to the simple frags vector `Vector{SimpleFrag{T}}`
    #This is used to build the framgent index. 
    function addSimpleFrags!(frags_simple::Vector{SimpleFrag{T}}, 
                                id_to_annotation::Vector{PioneerFragAnnotation},
                                exclude_from_index::Set{Char},
                                frags_simple_idx::Int64, 
                                rank_to_score::Vector{UInt8},
                                prec_idx::Int64,
                                prec::LibraryPrecursorIon{T},
                                prosit_frags::Vector{PioneerFrag}, #Critically, must be sorted in descending order by predicted intensity
                                max_frag_idx::Int64,
                                max_rank_index::Int64,
                                frag_mz_min::T,
                                frag_mz_max::T,
                                y_start::Int64,
                                b_start::Int64) where {T<:AbstractFloat}

        #Is assumed that frags are sorted in ascending rank order 
        simple_frags_added = 0
        for i in range(1, max_frag_idx)
            frag = prosit_frags[i]
            #If the fragment meets the constraints add to the simple frags list 
            if getBaseType(id_to_annotation[getType(frag)]) ∈ exclude_from_index
                continue
            end
            if inScanRange(frag, getBaseType(id_to_annotation[getType(frag)]),frag_mz_min, frag_mz_max, y_start, b_start)
                score = rank_to_score[simple_frags_added+1]#getScore(simple_frags_added+1) #score added to precursor for matching this fragment
                frags_simple[frags_simple_idx + simple_frags_added] = SimpleFrag(
                    getMZ(frag),
                    UInt32(prec_idx), #prec_id
                    getMZ(prec),
                    getIRT(prec),
                    getCharge(prec),
                    score
                )
                simple_frags_added += 1
            end
            if simple_frags_added == max_rank_index #Finished adding the top fragments 
                return frags_simple_idx + simple_frags_added
            end
        end
        return frags_simple_idx + simple_frags_added
    end

    #Add all fragments for a given precursor to the 
    #detailed frags list. Under m/z and minimum y/b index constraints 
    function addDetailedFrags!(frags_detailed::Vector{DetailedFrag{T}},
                                frags_detailed_idx::Int64,
                                prec_idx::Int64,
                                prec::LibraryPrecursorIon{T},
                                prosit_frags::Vector{PioneerFrag},
                                id_to_annotation::Vector{PioneerFragAnnotation},
                                max_frag_idx::Int64,
                                frag_mz_min::T,
                                frag_mz_max::T,
                                y_start::Int64, b_start::Int64) where {T<:AbstractFloat}
        detailed_frags_added = 0
        for i in range(1, max_frag_idx)
            frag = prosit_frags[i]
            annotation = id_to_annotation[getType(frag)]

            ion_type = zero(UInt8)
            if  annotation.base_type=='b'
                ion_type = one(UInt8)
            elseif annotation.base_type=='y'
                ion_type = UInt8(2)
            elseif annotation.base_type=='p'
                ion_type = UInt8(3)
            end
    
             #If the fragment meets the constraints add to the detailed frags list 
            if inScanRange(frag,annotation.base_type,frag_mz_min,frag_mz_max,y_start,b_start)
                f_charge = getCharge(frag)
                if annotation.base_type=='p'
                    f_charge = getCharge(prec)
                end
                frags_detailed[frags_detailed_idx + detailed_frags_added] = DetailedFrag(
                    UInt32(prec_idx), #prec_id

                    getMZ(frag), #mz
                    Float16(getIntensity(frag)), #intensity

                    ion_type,
                    false, #not isotope

                    f_charge, #frag_charge
                    getIndex(frag), #ion_position
                    getCharge(prec), #prec_charge
                    UInt8(detailed_frags_added+1), #rank
                    getSulfurCount(frag), #sulfur count
                )
                detailed_frags_added += 1
            end
        end
        return frags_detailed_idx, frags_detailed_idx + detailed_frags_added - 1
    end

    #Add a precursor ion to the precursors list 
    function addPrecursor!(precursors::Vector{LibraryPrecursorIon{Float32}}, 
                            precursor_idx::Int64,
                            irt::AbstractFloat,
                            prec_mz::AbstractFloat,

                            decoy::Bool,

                            proteome_name::String,
                            accession_numbers::String,
                            sequence::String,
                            structural_mods::Union{Missing, String},
                            isotopic_mods::Union{Missing, String},

                            charge::UInt8,
                            sulfur_count::UInt8,
                            sequence_length::UInt8)

        precursor_idx += 1

        function countMissedCleavages(seq::String, pattern::Regex)
            #pattern = r"[KR][^P|$]"
            # Count the number of matches of the pattern in the input string
            return length(collect((eachmatch(pattern, seq))))
        end

        missed_cleavage = countMissedCleavages(sequence, missed_cleavage_regex)
        precursors[precursor_idx] = LibraryPrecursorIon(
            Float32(irt),
            Float32(prec_mz), #need to fix chronologer output to include precursor mz from the prosit table
            decoy,

            proteome_name,
            accession_numbers,
            sequence,
            structural_mods,
            isotopic_mods,

            UInt8(charge),
            UInt8(missed_cleavage),
            sequence_length,
            sulfur_count
        )

        return precursor_idx
    end
    #thread_tasks = map(tasks) do task

    #Threads.@spawn begin 
    frag_det_idx_start, frag_det_idx_stop = 1, 1
    frag_simple_idx = 1
    frags_simple_idx, precursor_idx = zero(Int64), zero(Int64)
    max_N_frags = 1000
    frags = Vector{PioneerFrag}(undef, max_N_frags)
    #Loop through rows of prosit library (1-1 correspondence between rows and precursors)
    for row_idx in ProgressBar(sorted_indices)
        if (prec_mz[row_idx] < prec_mz_min) | (prec_mz[row_idx] > prec_mz_max)
            continue
        end
        #Get frag_mz bounds 
        frag_mz_min, frag_mz_max = frag_bounds(prec_mz[row_idx])
        #Grow placeholder fragment array if necessary 
        n_frags = length(prec_frag_ranges[row_idx]) 
        if n_frags > max_N_frags
            max_N_frags = n_frags
            frags = Vector{PioneerFrag}(undef, max_N_frags)
        end
        #Write precursor frags into the temporary array 
        for (i, frag_idx) in enumerate(prec_frag_ranges[row_idx])
            try
            frags[i] = PioneerFrag(
                frag_mz[frag_idx],
                frag_intensity[frag_idx],
                ion_types[frag_idx],
                frag_is_y[frag_idx],
                frag_index[frag_idx],
                frag_charge[frag_idx],
                frag_isotope[frag_idx],
                frag_internal[frag_idx],
                frag_immonium[frag_idx],
                frag_internal_ind[frag_idx],
                frag_sulfur_count[frag_idx],
            )#prec_frags[frag_idx]    
            catch
                println(typeof(frag_mz))
                println(typeof(frag_intensity))
                println(typeof(ion_types))
                println(typeof(frag_is_y))
                println(typeof(frag_index))
                println(typeof(frag_charge))
                println(typeof(frag_isotope))
                println(typeof(frag_internal))
                println(typeof(frag_internal_ind))
                println(typeof(frag_sulfur_count))
                error("urmom")
            end      
        end 
        #Sort sublist of fragments by descending order of intensity. 
        sort!(@view(frags[1:n_frags]), by = x -> getIntensity(x), rev = true, alg = QuickSort)
       
        precursor_idx = addPrecursor!(precursors, 
                        precursor_idx,
                        irt[row_idx],
                        prec_mz[row_idx],
                        decoy[row_idx],
                        proteome_names[row_idx],
                        accession_numbers[row_idx],
                        sequence[row_idx],
                        structural_mods[row_idx],
                        isotopic_mods[row_idx],
                        charge[row_idx],
                        sulfur_count[row_idx],
                        sequence_length[row_idx]
                        )
        #Map the new precursor id to the data table row 
        prec_id_to_library_row[precursor_idx] = UInt32(row_idx)

        frag_det_idx_start, frag_det_idx_stop = addDetailedFrags!(frags_detailed,
                                                                    frag_det_idx_start,
                                                                    precursor_idx,
                                                                    precursors[precursor_idx],
                                                                    frags,
                                                                    id_to_annotation,
                                                                    n_frags,
                                                                    frag_mz_min,
                                                                    frag_mz_max,
                                                                    y_start, b_start)
        #println("frag_det_idx_start $frag_det_idx_start, frag_det_idx_stop $frag_det_idx_stop")
        precursor_indices[precursor_idx] = frag_det_idx_start:frag_det_idx_stop
        frag_det_idx_start = frag_det_idx_stop + 1
        frag_simple_idx = addSimpleFrags!(frags_simple,
                                            id_to_annotation,
                                            exclude_from_index,
                                            frag_simple_idx,
                                            rank_to_score,
                                            precursor_idx,
                                            precursors[precursor_idx],
                                            frags,
                                            n_frags,
                                            max_rank_index,
                                            frag_mz_min,
                                            frag_mz_max,
                                            y_start_index, b_start_index)
    end
    
    println("writing detailed_frags...")
    @time jldopen(detailed_fragments_path, "w") do file
        file["detailed_fragments"] = frags_detailed[1:frag_det_idx_stop-1]
    end

    println("writing prec to fragment indices ...")
    @time jldopen(precursor_to_fragment_indices_path, "w") do file
        file["precursor_to_fragment_indices"] = precursor_indices[1:precursor_idx]
    end

    println("writing simple fragments... ")
    @time jldopen(simple_fragments_path, "w") do file
        file["simple_fragments"] = frags_simple[1:frag_simple_idx-1]
    end

    println("writing precursor tables...")
    @time Arrow.write(precursor_table_path,DataFrame(precursors[1:precursor_idx]))

    return frags_simple[1:frag_simple_idx-1], folder_out
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
    Arrow.write(joinpath(folder_out, index_name*"f_index_fragments.arrow"), fragments)
    Arrow.write(joinpath(folder_out, index_name*"f_index_rt_bins.arrow"), rt_bins)
    Arrow.write(joinpath(folder_out, index_name*"f_index_fragment_bins.arrow"), frag_bins)

    return 
end