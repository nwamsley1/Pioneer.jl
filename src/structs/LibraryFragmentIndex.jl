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
getFragmentBin(fi::FragmentIndex, bin::Int) = fi.fragment_bins[bin]
getPrecursorBin(fi::FragmentIndex, bin::UInt32) = fi.precursor_bins[bin]
getPrecursorBinLength(fi::FragmentIndex, bin::Int64) = getLength(fi.precursor_bins[bin])

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

function setPrecursorBinItem!(fi::FragmentIndex{T}, bin::UInt32, index::Int, prec_bin_item::PrecursorBinFragment{T}) where {T<:AbstractFloat}
    setPrecursor!(fi.precursor_bins[bin], index, prec_bin_item)
end

function addPrecursorBinFragment!(fi::FragmentIndex{T}, bin::UInt32, prec_bin_item::PrecursorBinFragment{T}) where {T<:AbstractFloat}
    push!(getPrecursors(fi.precursor_bins[bin]), prec_bin_item)
end

function addPrecursorBin!(fi::FragmentIndex{T}, prec_bin::PrecursorBin{T}) where {T<:AbstractFloat}
    push!(getPrecBins(fi), prec_bin)
end

getPPM(frag_mz::T, ppm::T) where {T<:AbstractFloat} = ppm*frag_mz/1e6
