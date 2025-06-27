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

abstract type FragmentIndexBin{T<:AbstractFloat} end 

getLow(fb::FragmentIndexBin{T}) where {T<:AbstractFloat} = fb.lb
getHigh(fb::FragmentIndexBin{T}) where {T<:AbstractFloat} = fb.ub
getSubBinRange(fb::FragmentIndexBin{T}) where {T<:AbstractFloat} = fb.first_bin:fb.last_bin

struct FragIndexBin{T<:AbstractFloat} <: FragmentIndexBin{T}
    lb::T
    ub::T
    first_bin::UInt32
    last_bin::UInt32
end
ArrowTypes.arrowname(::Type{FragIndexBin{Float32}}) = :FragIndexBin
ArrowTypes.JuliaType(::Val{:FragIndexBin}) = FragIndexBin

struct IndexFragment{T<:AbstractFloat} <: LibraryFragmentIon{T}
    prec_id::UInt32
    prec_mz::T #Only need to tell if the peptide is in the quad isolation window
    score::UInt8 
    charge::UInt8
end
ArrowTypes.arrowname(::Type{IndexFragment{Float32}}) = :IndexFragment
ArrowTypes.JuliaType(::Val{:IndexFragment}) = IndexFragment

getScore(ind_frag::IndexFragment{T}) where {T<:AbstractFloat} = ind_frag.score

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
    fragment_bins::Arrow.Struct{FragIndexBin, Tuple{Arrow.Primitive{T, Vector{T}}, Arrow.Primitive{T, Vector{T}}, Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{UInt32, Vector{UInt32}}}, (:lb, :ub, :first_bin, :last_bin)}
    rt_bins::Arrow.Struct{FragIndexBin, Tuple{Arrow.Primitive{T, Vector{T}}, Arrow.Primitive{T, Vector{T}}, Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{UInt32, Vector{UInt32}}}, (:lb, :ub, :first_bin, :last_bin)}
    fragments::Arrow.Struct{IndexFragment, Tuple{Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{T, Vector{T}}, Arrow.Primitive{UInt8, Vector{UInt8}}, Arrow.Primitive{UInt8, Vector{UInt8}}}, (:prec_id, :prec_mz, :score, :charge)}
end
#struct FragmentIndex{T<:AbstractFloat}
#    fragment_bins::Vector{FragIndexBin{T}}
#    rt_bins::Vector{FragIndexBin{T}}
#    fragments::Vector{IndexFragment{T}}
#end

getFragBins(fi::FragmentIndex{T}) where {T<:AbstractFloat} = fi.fragment_bins
getRTBins(fi::FragmentIndex{T}) where {T<:AbstractFloat} = fi.rt_bins
getFragmentBin(fi::FragmentIndex{T}, frag_bin_idx::I) where {T<:AbstractFloat,I<:Integer} = getFragBins(fi)[frag_bin_idx]
getRTBin(fi::FragmentIndex{T}, rt_bin_idx::I) where {T<:AbstractFloat,I<:Integer} = getRTBins(fi)[rt_bin_idx]
getFragments(fi::FragmentIndex{T}) where {T<:AbstractFloat} = fi.fragments


abstract type SpectralLibrary end

struct FragmentIndexLibrary <: SpectralLibrary
    presearch_fragment_index::FragmentIndex{Float32}
    fragment_index::FragmentIndex{Float32}
    precursors::LibraryPrecursors
    proteins::LibraryProteins
    fragment_lookup_table::StandardFragmentLookup
end

struct SplineFragmentIndexLibrary <: SpectralLibrary
    presearch_fragment_index::FragmentIndex{Float32}
    fragment_index::FragmentIndex{Float32}
    precursors::LibraryPrecursors
    proteins::LibraryProteins
    fragment_lookup_table::SplineFragmentLookup
end


getPresearchFragmentIndex(sl::SpectralLibrary) = sl.presearch_fragment_index
getFragmentIndex(sl::SpectralLibrary) = sl.fragment_index
getPrecursors(sl::SpectralLibrary) = sl.precursors
getFragmentLookupTable(sl::SpectralLibrary) = sl.fragment_lookup_table 
getProteins(sl::SpectralLibrary) = sl.proteins
