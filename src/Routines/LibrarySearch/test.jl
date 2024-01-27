abstract type IonType end
abstract type FragmentIndexType <: IonType end
getMZ(f::FragmentIndexType) = f.frag_mz
getFragMZ(f::FragmentIndexType) = f.frag_mz
getPrecID(f::FragmentIndexType) = f.prec_id
getPrecCharge(f::FragmentIndexType) = f.prec_charge

import Base.<
import Base.>

<(y::FragmentIndexType, x::T) where {T<:Real} = getFragMZ(y) < x
<(x::T, y::FragmentIndexType) where {T<:Real} = <(y, x)
>(y::FragmentIndexType, x::T) where {T<:Real} = getFragMZ(y) > x
>(x::T, y::FragmentIndexType) where {T<:Real} = >(y, x)

struct FragmentIon{T<:AbstractFloat} <: FragmentIndexType
    frag_mz::T
    prec_id::UInt32
    prec_mz::T
    prec_intensity::T #Needs to be updated
    prec_rt::T
    prec_charge::UInt8
end

getPrecMZ(f::FragmentIon) = f.prec_mz
getIntensity(f::FragmentIon) = f.prec_intensity[]
getRT(f::FragmentIon) = f.prec_rt

struct LibraryFragment{T<:AbstractFloat} <: FragmentIndexType
    frag_mz::T
    frag_charge::UInt8
    is_y_ion::Bool
    is_isotope::Bool
    ion_position::UInt8
    ion_index::UInt8
    intensity::Float32
    prec_charge::UInt8
    prec_id::UInt32
    rank::UInt8
    sulfur_count::UInt8
end

getIntensity(f::LibraryFragment) = f.intensity
isyIon(f::LibraryFragment) = f.is_y_ion
getIonIndex(f::LibraryFragment) = f.ion_index
getIonPosition(f::LibraryFragment) = f.ion_position
getFragCharge(f::LibraryFragment) = f.frag_charge
getRank(f::LibraryFragment) = f.rank
sulfurCount(f::LibraryFragment) = f.sulfur
#LibraryFragment{T}() where {T<:AbstractFloat} = LibraryFragment(zero(T), zero(UInt8), false, zero(UInt8), zero(UInt8), zero(Float32), zero(UInt8), zero(UInt32), zero(UInt8), zero(UInt8))
LibraryFragment{T}() where {T<:AbstractFloat} = LibraryFragment(zero(T), 
zero(UInt8), 
false, 
false,
zero(UInt8), 
zero(UInt8), 
zero(Float32), 
zero(UInt8),
zero(UInt32), 
zero(UInt8),
zero(UInt8)
 )

struct LibraryPrecursor{T<:AbstractFloat}
    iRT::T
    struct LibraryPrecursor{T<:AbstractFloat}
    iRT::T
    mz::T
    total_intensity::T
    base_peak_intensity::T
    isDecoy::Bool
    charge::UInt8
    pep_id::UInt32
    prot_ids::Vector{UInt32}
    accession_numbers::String
    sequence::String
    missed_cleavages::UInt8
    variable_mods::UInt8
    length::UInt8
    sulfur_count::UInt8
end

isDecoy(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.isDecoy
getIRT(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.iRT
getCharge(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.charge
getMz(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.mz
getTotalIntensity(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.total_intensity
getPepID(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.pep_id
getBasePeakInt(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.base_peak_intensity
sulfurCount(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.sulfur_count
    mz::T
    total_intensity::T
    base_peak_intensity::T
    isDecoy::Bool
    charge::UInt8
    pep_id::UInt32
    prot_ids::Vector{UInt32}
    accession_numbers::String
    sequence::String
    missed_cleavages::UInt8
    variable_mods::UInt8
    length::UInt8
    sulfur_count::UInt8
end

isDecoy(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.isDecoy
getIRT(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.iRT
getCharge(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.charge
getMz(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.mz
getTotalIntensity(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.total_intensity
getPepID(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.pep_id
getBasePeakInt(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.base_peak_intensity
sulfurCount(p::LibraryPrecursor{T}) where {T<:AbstractFloat} = p.sulfur_count
#addIntensity!(p::LibraryPrecursor{T}, intensity::T) where {T<:AbstractFloat} = p.total_intensity[] += intensity
#ArrowTypes.arrowname(::Type{LibraryFragment{Float32}}) = :LibraryFragment
#ArrowTypes.JuliaType(::Val{:LibraryFragment}) = LibraryFragment
#fixed_mods = [(p=r"C", r="C[Carb]")]
#mods_dict = Dict("Carb" => Float64(57.021464),
#                 "Ox" => Float64(15.994915)
#                 )

function parseSequence(seq::String, charge::UInt8, fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, mods_dict::Dict{String, T}) where {T<:AbstractFloat}

    seq = fixedMods(replace(seq, "M(ox)"=>"M[Ox]"), fixed_mods)

    mz = getMZ(Precursor(
                            seq, 
                            charge = charge,
                            mods_dict = mods_dict
                        )
                )

    function getPlainSeq(seq::String)
        pattern = r"\[.*?\]"
        # Remove all occurrences of the pattern from the input string
        return replace(seq, pattern => "")
    end

    function countMissedCleavages(seq::String)
        pattern = r"[KR][^P|$]"
        # Count the number of matches of the pattern in the input string
        return length(collect((eachmatch(pattern, seq))))
    end
    seq = getPlainSeq(seq)

    return Float32(mz), UInt8(length(seq)), UInt8(countMissedCleavages(seq))
end

function split_array_into_chunks(arr::Vector, num_chunks::Int)
    len = length(arr)
    chunk_size = div(len, num_chunks)
    
    # Create a list of tuples for each chunk
    chunks = [[start_idx, min(start_idx + chunk_size - 1, len)] for start_idx in 1:chunk_size:len]
    
    # Adjust the last chunk's stopping index to the end of the array
    last_chunk = last(chunks)
    if last_chunk[2] < len
        last_chunk = [last_chunk[1], len]
        chunks[end] = last_chunk
    end
    
    return chunks
end
#CSV.write("outputfile.csv",test_prosit[1000000:1000036,:],delim=',')

function getMZBounds(mz::T) where {T<:AbstractFloat}
    #low_frag_mz = Float32(((mz - 404.4337)÷8.0037 + 2)*1.667 + 83.667)
    #high_frag_mz = Float32(((mz - 404.4337)÷8.0037 - 2)*25 + 1265.0)
    return 250.0f0, 1500.0f0
end

#test = CSV.File("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/Prosit_HumanYeastEcoli_NCE33_corrected_092823.csv")
#test_df = DataFrame(test)


function parsePrositLib(prosit_csv_path::String, fixed_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}}, 
                        mods_dict::Dict{String, T}, getMZBounds::Function; 
                        max_rank_index::Int64 = 5, y_start::Int64 = 2, b_start::Int64 = 1) where {T<:AbstractFloat}
    println("START")
    #"/Users/n.t.wamsley/Projects/PROSIT/prosit1/my_prosit/examples/peptidelist.msms"
    prosit_library = CSV.File(prosit_csv_path)
    charge_pattern = r"\((\w+)\+\)"
    index_pattern = r"[yb]([0-9]{1,2})[\(]*"
    println("LOADED")
    #How many fragments? Determines size of pre-allocated arrays
    frag_count = 0
    for i in ProgressBar(range(1, length(prosit_library)))
        frag_count += count(";", prosit_library[i].intensities) + 1
    end

    ###########
    #Initialize Containers
    #List of lists. The n'th element of the outer list is indexed by the precursor id.
    #Each inner list contains a list of "LibraryFragment"s. These are a detailed specificaion of a framgent ion
    frags_detailed = Vector{LibraryFragment{Float32}}(undef, frag_count)
    #Keeps track of with fragments in "frags_detailed" correspond to which precursor
    precursor_indices = [(zero(UInt32), zero(UInt32)) for _ in range(1, length(prosit_library))]
    #Detailed specification for each precursor. 
    precursors = Vector{LibraryPrecursor{Float32}}(undef, length(prosit_library))

    #Used to fuild the framgment index for MSFragger-like searching. 
    frags_simple = Vector{FragmentIon{Float32}}(undef, length(prosit_library)*3)
    #Loop through rows of prosit library (1-1 correspondence between rows and precursors)
    #lk = ReentrantLock()
    #N = split_array_into_chunks(frags_detailed, Threads.nthreads())

    function inScanRange(mass::T, ion_index::UInt8, low_frag_mz::T, high_frag_mz::T, is_y::Bool, y_start_ind::Int64, b_start_ind::Int64) where {T<:AbstractFloat}
        if is_y
            return (mass>low_frag_mz)&(mass<high_frag_mz)&(ion_index>=y_start_ind)
        else
            return (mass>low_frag_mz)&(mass<high_frag_mz)&(ion_index>=b_start_ind)
        end
    end

    #thread_tasks = map(tasks) do task

    #Threads.@spawn begin 
    frag_det_idx_start, frag_det_idx_stop = 1, 1
    frag_simple_idx = 1
    for (prec_idx, precursor) in ProgressBar(enumerate(prosit_library))
        #############
        #Parse Column/Precursor
        #############
        matches =  split(precursor[:matched_ions],';')::Vector{SubString{String}}
        masses = parse.(Float32, split(precursor[:masses],';')::Vector{SubString{String}})::Vector{Float32}
        intensities = parse.(Float32, split(precursor[:intensities],';')::Vector{SubString{String}})::Vector{Float32}
        base_peak_intensity = maximum(intensities) #Base peak intensity
        #intensities = intensities./sum(intensities)
        charge = UInt8(precursor[:Charge]::Int64)
        pep_id =  UInt32(precursor[:pep_id]::Int64)
        decoy =  precursor[:decoy]::Bool
        iRT = Float32(precursor[:iRT]::Float64)
        sequence = precursor[:modified_sequence]::String31

        accession_numbers = precursor[:accession_numbers]::String
        mz, len_AA, missed_cleavages = parseSequence(String(sequence), charge, fixed_mods, mods_dict)

        #pre-allocate library fragments for the n'th precursor
        #nth_precursor_frags = Vector{LibraryFragment{Float32}}()#undef, length(matches))
        total_intensity = zero(Float32)

        #ms2 scan range 
        low_frag_mz, high_frag_mz = getMZBounds(mz)

        n_frags = 0
        for (i, _match) in enumerate(matches)

            ion_index = parse(UInt8, match(index_pattern, _match).captures[1]) 

            if !inScanRange(masses[i], ion_index, low_frag_mz, high_frag_mz, (_match[1] == 'y'), y_start, b_start)
                intensities[i] = zero(Float32)
            else
                n_frags += 1
                total_intensity += intensities[i]
            end
        end
        
        #nth_precursor_frags = Vector{LibraryFragment{Float32}}(undef, n_frags)
        ranks = ordinalrank(intensities, rev = true)
        #########
        #Parse each fragment ion 
        for i in eachindex(matches)
            fragment_name = matches[i]
            #ion_type = match(pattern, fragment_name)
            ion_index = parse(UInt8, match(index_pattern, fragment_name).captures[1]) #for y1 or b12 this is "1" and "12 respectively
            ion_charge = match(charge_pattern, fragment_name)
            frag_charge = ion_charge === nothing ? one(UInt8) : parse(UInt8, ion_charge.captures[1])

            if iszero(intensities[i])#!inScanRange(masses[i], ion_index, low_frag_mz, high_frag_mz, (matches[i][1] == 'y'), y_start, b_start)
                #frag_det_idx += 1
                continue
            else
                frag_det_idx_stop += 1
            end

            #Add the n'th fragment 
            frags_detailed[frag_det_idx_stop] =  LibraryFragment(
                                                masses[i], #frag_mz
                                                frag_charge, #frag_charge
                                                fragment_name[1] == 'y' ? true : false, #is_y_ion
                                                false, #is_isotope
                                                ion_index, #ion_position
                                                UInt8(i), #ion_index
                                                intensities[i], #intensity
                                                charge, #prec_charge
                                                UInt32(prec_idx), #prec_id
                                                UInt8(ranks[i]),
                                                one(UInt8) #sulfurr_count
                                                )
            ###########
            #Only retain most important fragments for the index!
            if (ranks[i]<=max_rank_index)


                intensity = zero(Float32)
                if ranks[i] == 1
                    intensity = one(Float32)
                elseif ranks[i] < 4
                    intensity = Float32(0.7)
                end

               frags_simple[frag_simple_idx] = FragmentIon(
                            masses[i], #frag_mz
                            UInt32(prec_idx), #prec_id
                            mz, #prec_mz 
                            intensity, #prec_intensity
                            iRT, #prec_rt
                            charge, #prec_charge
                        )
                frag_simple_idx += 1
            end

        end

        precursor_indices[prec_idx] = (frag_det_idx_start, frag_det_idx_stop)
        frag_det_idx_start = frag_det_idx_stop + 1

        precursors[prec_idx] = LibraryPrecursor(
                                            iRT, #iRT
                                            mz, #mz
                                            total_intensity, #total_intensity
                                            base_peak_intensity,
                                            decoy, #isDecoy
                                            charge, #charge
                                            pep_id, #pep_id
                                            UInt32[1], #prot_ids
                                            String(accession_numbers), #accession_numbers
                                            String(sequence),
                                            missed_cleavages, #missed_cleavages
                                            UInt8(1), #variable_mods
                                            len_AA, #length
                                            one(UInt8) #sulfur count
        )
    end

    #end

    #return vcat(frags_simple...), frags_detailed, precursors
    return frags_simple, frags_detailed, precursors
end

"""
    PrecursorBinItem{T<:AbstractFloat}

Item in a precursor bin. Minimal information required to know if the fragment has a correct precursor mass, and if so, 
what the precursor ID is.  

### Fields

- prec_id::UInt32 -- Unique identifier for the precursor
- prec_mz::T -- m/z of the precursor

### Examples

### GetterMethods

- getPrecID(pbi::PrecursorBinItem) = pbi.prec_id
- getPrecMZ(pbi::PrecursorBinItem) = pbi.prec_mz
"""
struct PrecursorBinItem{T<:AbstractFloat}
    prec_id::UInt32
    prec_mz::T
    intensity::T
    charge::UInt8
end

getPrecID(pbi::PrecursorBinItem) = pbi.prec_id
getPrecMZ(pbi::PrecursorBinItem) = pbi.prec_mz
getPrecCharge(pbi::PrecursorBinItem) = pbi.charge
getIntensity(pbi::PrecursorBinItem) = pbi.intensity


"""
    PrecursorBin{T<:AbstractFloat}

Container for precursor Items. 

### Fields

- precs::Vector{PrecursorBinItem{T}} -- Precursor items

### Examples

- PrecursorBin(T::DataType, N::Int) = PrecursorBin(Vector{PrecursorBinItem{T}}(undef, N)) -- Constructor

### GetterMethods

- getPrecursors(pb::PrecursorBin) = pb.precs

### Methods

- setPrecursor!(pb::PrecursorBin, index::Int, pbi::PrecursorBinItem)
"""
struct PrecursorBin{T<:AbstractFloat}
    precs::Vector{PrecursorBinFragment{T}}
end

getPrecursors(pb::PrecursorBin)  = pb.precs
getPrecursor(pb::PrecursorBin, i::Int64) = getPrecursors(pb)[i] 
getLength(pb::PrecursorBin)= length(pb.precs)

function setPrecursor!(pb::PrecursorBin, index::Int, pbi::PrecursorBinFragment)
    pb.precs[index] = pbi
end

PrecursorBin(T::DataType, N::Int) = PrecursorBin(Vector{PrecursorBinFragment{T}}(undef, N))

"""
    FragBin{T<:AbstractFloat}

Represents a bin of sorted fragment ions. Gives the lowest and highest m/z ratio of ions in the bin. 
Also contains a unique identifier for the corresponding precursor bin. 

### Fields

- lb::T -- Lowest m/z of a fragment ion in the `FragBin`
- ub::T -- Highest m/z of a fragment ion in the `FragBin`
- prec_bin::UInt32 -- Identifier of the corresponding `PrecursorBin`

### Examples

- FragBin() = FragBin(0.0, 0.0, UInt32(0)) -- Constructor

### GetterMethods

- getLowMZ(fb::FragBin) = fb.lb
- getHighMZ(fb::FragBin) = fb.ub
- getPrecBinID(fb::FragBin) = fb.prec_bin

### Methods

- setPrecursor!(pb::PrecursorBin, index::Int, pbi::PrecursorBinItem)
"""
struct FragBin{T<:AbstractFloat}
    lb::T
    ub::T
    sub_bin::UInt32
end

getLowMZ(fb::FragBin) = fb.lb
getHighMZ(fb::FragBin) = fb.ub
getPrecBinID(fb::FragBin) = fb.sub_bin
FragBin() = FragBin(0.0, 0.0, UInt32(0))

struct RTBin{T<:AbstractFloat}
    lb::T
    ub::T
    sub_bin::UInt32
end
RTBin() = RTBin(0.0, 0.0, UInt32(0))
getLow(rb::RTBin) = rb.lb
getHigh(rb::RTBin) = rb.ub
getPrecBinID(rb::RTBin) = rb.sub_bin


"""
    FragmentIndex{T<:AbstractFloat}

A fragment index for an MSFragger-style/fragment-centric search. Contains fragment bins 
that indicate the highest and lowest m/z of a fragment ion in the bin. Each `FragBin` links
to a `PrecursorBin`. A precursor bin has the precursor m/z's and peptide ID's for each fragment ion. 

### Fields

- fragment_bins::Vector{FragBin{T}} -- Lowest m/z of a fragment ion in the `FragBin`
- precursor_bins::Vector{PrecursorBin{T}} -- Highest m/z of a fragment ion in the `FragBin`

### Examples

- FragmentIndex(T::DataType, M::Int, N::Int) = FragmentIndex(fill(FragBin(), N), fill(PrecursorBin(T, M), N))

### GetterMethods

- getFragmentBin(fi::FragmentIndex, bin::Int) = fi.fragment_bins[bin]
- getPrecursorBin(fi::FragmentIndex, bin::Int64) = fi.precursor_bins[bin]

### Methods

- setFragmentBin!(fi::FragmentIndex, bin::Int64, frag_bin::FragBin)
- setPrecursorBinItem!(fi::FragmentIndex{T}, bin::Int64, index::Int64, prec_bin_item::PrecursorBinItem{T}) where {T<:AbstractFloat}
"""
struct FragmentIndex{T<:AbstractFloat}
    fragment_bins::Vector{FragBin{T}}
    rt_bins::Vector{Vector{RTBin{T}}}
    precursor_bins::Vector{PrecursorBin{T}}
end

getFragBins(fi::FragmentIndex) = fi.fragment_bins
getRTBins(fi::FragmentIndex) = fi.rt_bins
getPrecBins(fi::FragmentIndex) = fi.precursor_bins


#FragmentIndex(T::DataType, M::Int, N::Int) =  FragmentIndex(fill(FragBin(), N), fill(PrecursorBin(T, M), N))

function FragmentIndex(T::DataType) 
    return FragmentIndex(Vector{FragBin{T}}(), Vector{Vector{RTBin{T}}}(), Vector{PrecursorBin{T}}())
end

getFragmentBin(fi::FragmentIndex, bin::Int) = fi.fragment_bins[bin]

function setFragmentBin!(fi::FragmentIndex{T}, bin::Int64, frag_bin::FragBin{T}) where {T<:AbstractFloat}
    fi.fragment_bins[bin] = frag_bin
end

function addFragmentBin!(fi::FragmentIndex{T}, frag_bin::FragBin{T}) where {T<:AbstractFloat}
    push!(getFragBins(fi), frag_bin)
end

#function setRTBin!(fi::FragmentIndex{T}, bin::Int64, rt_bin::RTBin{T}) where {T<:AbstractFloat}
#    fi.rt_bins[bin] = rt_bin
#end

function addRTBin!(fi::FragmentIndex{T}) where {T<:AbstractFloat}
    push!(getRTBins(fi), Vector{RTBin{T}}())
end


getPrecursorBin(fi::FragmentIndex, bin::UInt32) = fi.precursor_bins[bin]
getPrecursorBinLength(fi::FragmentIndex, bin::Int64) = getLength(fi.precursor_bins[bin])

function setPrecursorBinItem!(fi::FragmentIndex{T}, bin::UInt32, index::Int, prec_bin_item::PrecursorBinItem{T}) where {T<:AbstractFloat}
    setPrecursor!(fi.precursor_bins[bin], index, prec_bin_item)
end

function addPrecursorBinItem!(fi::FragmentIndex{T}, bin::UInt32, prec_bin_item::PrecursorBinItem{T}) where {T<:AbstractFloat}
    push!(getPrecursors(fi.precursor_bins[bin]), prec_bin_item)
end

function addPrecursorBin!(fi::FragmentIndex{T}, prec_bin::PrecursorBin{T}) where {T<:AbstractFloat}
    push!(getPrecBins(fi), prec_bin)
end

getPPM(frag_mz::T, ppm::T) where {T<:AbstractFloat} = ppm*frag_mz/1e6


function buildFragmentIndex!(frag_ions::Vector{FragmentIon{T}}, bin_ppm::AbstractFloat, rt_size::AbstractFloat; low_frag_mz::AbstractFloat = 80.0, high_frag_mz::AbstractFloat = 3000.0, low_prec_mz::AbstractFloat = 300.0, high_prec_mz::AbstractFloat = 1100.0) where {T<:AbstractFloat}
    println("sorting...")
    sort!(frag_ions, by = x->x.frag_mz)
    println("sorted")
    #The fragment ions are divided into bins of roughtly equal m/z width.
    #That should correspond to roughly half the fragment mass accuracy of the detector?
    frag_index = FragmentIndex(T) 

    function fillPrecursorBin!(frag_index::FragmentIndex{<:AbstractFloat}, frag_ions::Vector{FragmentIon{T}}, frag_bin_idx::UInt32, start::Int, stop::Int, low_prec_mz::AbstractFloat, high_prec_mz::AbstractFloat)
        for ion_index in range(start, stop)
            pep_id = getPrecID(frag_ions[ion_index])
                prec_mz = getPrecMZ(frag_ions[ion_index])

                if (prec_mz < low_prec_mz) | (prec_mz > high_prec_mz) #Precursor m/z outside the bounds
                    continue
                end
                #Add precursor corresponding to the charge state
                addPrecursorBinItem!(frag_index,
                                     frag_bin_idx,
                                    PrecursorBinItem(pep_id, prec_mz, getIntensity(frag_ions[ion_index]), getPrecCharge(frag_ions[ion_index]))
                                    )
        end
    end

    rt_bin_idx = UInt32(1) #Current fragment bin index
    prec_bin_idx = one(UInt32)
    start = 1 #Fragment index of the first fragment in the current bin
    stop = 2

    diff = getPPM(getFragMZ(frag_ions[start]), bin_ppm) #ppm tolerance of the current fragment bin
    #Build bins 
    while stop < (length(frag_ions)) + 1

        if (stop % 1_000_000) == 0
            println(stop/1_000_000)
        end
        #Haven't reached minimum fragment m/z yet
        if getFragMZ(frag_ions[start]) < low_frag_mz
            start += 1
            stop += 1
            continue
        end

        #Does the difference between the highest and lowest m/z in the bin 
        #enought exceed 10 ppm of the lowest m/z?
        if ((getFragMZ(frag_ions[min(stop + 1, length(frag_ions))]) - getFragMZ(frag_ions[start])) > diff) | ((stop + 1) >= length(frag_ions))
            
            #Add a new fragment bin
            addFragmentBin!(frag_index, 
                            FragBin(getFragMZ(frag_ions[start]),
                                    getFragMZ(frag_ions[stop]),
                                    rt_bin_idx
                                    )
                            )

            #Holds the RT bins for these fragments. 
            addRTBin!(frag_index)
            
            #Sort fragments in the frag_bin by retention time
            frag_ions[start:stop] = sort(frag_ions[start:stop], by = frag -> frag.prec_rt)
            start_rt_idx = start
            #Make new precursor and retention time bins. 
            for i in start:stop
                #The first and last fragment differ in retention time by a threshold
                if ((frag_ions[i].prec_rt + (-1)*frag_ions[start_rt_idx].prec_rt) > rt_size) | (i == stop)

                    #Create a retention time bin to hold the start_idx:i-1 fragments 
                    push!(frag_index.rt_bins[rt_bin_idx] , RTBin(frag_ions[start_rt_idx].prec_rt,
                                                            frag_ions[max(i-1, start_rt_idx)].prec_rt,
                                                            prec_bin_idx)
                            )

                    #Create a precursor bin to hold the start_idx:i-1fragments. 
                    addPrecursorBin!(frag_index, 
                                        #Preallocate an empty precursor bin of the correct length 
                                        PrecursorBin(Vector{PrecursorBinItem{T}}())#undef, (last_frag_in_bin - start + 1)*length(charges)))
                                        )
                    if i == stop
                        fillPrecursorBin!(frag_index, frag_ions, prec_bin_idx, start_rt_idx, i, low_prec_mz, high_prec_mz)
                    else
                        fillPrecursorBin!(frag_index, frag_ions, prec_bin_idx, start_rt_idx, i-1, low_prec_mz, high_prec_mz)
                    end

                    start_rt_idx = i
                    prec_bin_idx += one(UInt32)
                end
            end
            #Update counters and ppm tolerance 
            rt_bin_idx += UInt32(1)
            start = stop + 1
            diff = getPPM(getFragMZ(frag_ions[min(start, length(frag_ions))]), bin_ppm)
            stop = start
            #Maximum fragment m/z reached. Stop adding bins.
            if getFragMZ(frag_ions[min(start, length(frag_ions))]) > high_frag_mz
                break
            end
            
        else
            stop += 1
        end
    end

    function sortPrecursorBins!(frag_index::FragmentIndex{<:AbstractFloat})
        for i in 1:length(frag_index.precursor_bins)
            sort!(frag_index.precursor_bins[i].precs, by = x->getPrecMZ(x));
        end
        return nothing
    end
    sortPrecursorBins!(frag_index)
    return frag_index
end