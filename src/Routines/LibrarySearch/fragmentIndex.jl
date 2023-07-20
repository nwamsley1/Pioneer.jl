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
#=struct PrecursorBinItem{T<:AbstractFloat}
    prec_id::UInt32
    prec_mz::T
    charge::UInt8
end

getPrecID(pbi::PrecursorBinItem) = pbi.prec_id
getPrecMZ(pbi::PrecursorBinItem) = pbi.prec_mz
getPrecCharge(pbi::PrecursorBinItem) = pbi.charge
getIntensity(pbi::PrecursorBinItem) = pbi.intensity=#

struct PrecursorBinItem{T<:AbstractFloat}
    prec_id::UInt32
    prec_mz::T
    intensity::Float32
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
    precs::Vector{PrecursorBinItem{T}}
end

getPrecursors(pb::PrecursorBin)  = pb.precs
getPrecursor(pb::PrecursorBin, i::Int64) = getPrecursors(pb)[i] 
getLength(pb::PrecursorBin)= length(pb.precs)

function setPrecursor!(pb::PrecursorBin, index::Int, pbi::PrecursorBinItem)
    pb.precs[index] = pbi
end

PrecursorBin(T::DataType, N::Int) = PrecursorBin(Vector{PrecursorBinItem{T}}(undef, N))

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

function buildFragmentIndex!(frag_ions::Vector{FragmentIon{T}}, bin_ppm::AbstractFloat, rt_size::AbstractFloat; low_frag_mz::AbstractFloat = 150.0, high_frag_mz::AbstractFloat = 1700.0, low_prec_mz::AbstractFloat = 300.0, high_prec_mz::AbstractFloat = 1100.0) where {T<:AbstractFloat}
    println("sorting...")
    sort!(frag_ions, by = x->x.frag_mz)
    println("sorted")
    #The fragment ions are divided into bins of roughtly equal m/z width.
    #That should correspond to roughly half the fragment mass accuracy of the detector?
    frag_index = FragmentIndex(T) 

    function fillPrecursorBin!(frag_index::FragmentIndex{<:AbstractFloat}, frag_ions::Vector{FragmentIon{T}}, frag_bin_idx::UInt32, start::Int, stop::Int, low_prec_mz::AbstractFloat, high_prec_mz::AbstractFloat)
        for ion_index in range(start, stop)
            pep_id = getPrecID(frag_ions[ion_index])
                prec_mz = getPrecMZ(frag_ions[ion_index])#(getPrecMZ(frag_ions[ion_index]) + PROTON*(charge-1))/charge #m/z of the precursor

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

    bin = UInt32(1) #Current fragment bin index
    prec_bin_idx = one(UInt32)
    start = 1 #Fragment index of the first fragment in the current bin

    getPPM(frag_mz::T, ppm::T) = ppm*frag_mz/1e6

    diff = getPPM(getFragMZ(frag_ions[start]), bin_ppm) #ppm tolerance of the current fragment bin

    #Build bins 
    for stop in 2:length(frag_ions)

        if (stop % 1_000_000) == 0
            println(stop/1_000_000)
        end
        #Haven't reached minimum fragment m/z yet
        if getFragMZ(frag_ions[stop]) < low_frag_mz
            start += 1
            continue
        end

        #Doex the difference between the highest and lowest m/z in the bin 
        #enought exceed 10 ppm of the lowest m/z?
        if (getFragMZ(frag_ions[stop]) - getFragMZ(frag_ions[start])) > diff

            #Nedds to be stop - 1 to gaurantee the smallest and largest fragment
            #in the bin differ by less than diff 
            last_frag_in_bin = stop - 1

            #Add a new fragment bin
            addFragmentBin!(frag_index, 
                            FragBin(getFragMZ(frag_ions[start]),
                                    getFragMZ(frag_ions[last_frag_in_bin]),
                                    bin
                                    )
                            )
            addRTBin!(frag_index)
            #Sort fragments in the frag_bin by retention time
            frag_ions[start:last_frag_in_bin] = sort(frag_ions[start:last_frag_in_bin], by = frag -> frag.prec_rt)
            start_rt = frag_ions[start].prec_rt
            start_rt_idx = start
            for i in start:last_frag_in_bin
                if (frag_ions[i].prec_rt - frag_ions[start_rt_idx].prec_rt) > rt_size

                    push!(frag_index.rt_bins[bin] , RTBin(frag_ions[start_rt_idx].prec_rt,
                                                            frag_ions[i].prec_rt,
                                                            prec_bin_idx)
                            )

                    addPrecursorBin!(frag_index, 
                                        #Preallocate an empty precursor bin of the correct length 
                                        PrecursorBin(Vector{PrecursorBinItem{T}}())#undef, (last_frag_in_bin - start + 1)*length(charges)))
                                        )
                
                    fillPrecursorBin!(frag_index, frag_ions, prec_bin_idx, start_rt_idx, i, low_prec_mz, high_prec_mz)

                    #Sort the precursor bin by precursor m/z
                    #sort!(getPrecursors(getPrecursorBin(frag_index, prec_bin_idx)), by = x->getPrecMZ(x));
                    #sort!(frag_index.precursor_bins[prec_bin_idx].precs, by = x->getPrecMZ(x));
                    start_rt_idx = i+1
                    prec_bin_idx += one(UInt32)
                end
            end
            #Update counters and ppm tolerance 
            bin += UInt32(1)
            start = stop
            diff = getPPM(getFragMZ(frag_ions[start]), bin_ppm)

            #Maximum fragment m/z reached. Stop adding bins. 
            if getFragMZ(frag_ions[stop]) > high_frag_mz
                break
            end
        end
    end
    return frag_index
end

