"""
    FragmentIon{T<:AbstractFloat}

Simple representation of a fragment ion for use in building fragment indices. 

### Fields

- frag_mz::T -- m/z ratio of the fragment ion
- pep_id::UInt32 -- Identifier of the peptide from which the fragment came
- prec_mz::T  -- m/z ratio of the precursor from which the fragment came

### Examples

### GetterMethods

- getFragMZ(f::FragmentIon) = f.frag_mz
- getPepID(f::FragmentIon) = f.pep_id
- getPrecMZ(f::FragmentIon) = f.prec_mz
- <(y::FragmentIon, x::T) where {T<:Real} = getFragMZ(y) < x
- >(y::FragmentIon, x::T) where {T<:Real} = getFragMZ(y) > x
"""
struct FragmentIon{T<:AbstractFloat}
    frag_mz::T
    pep_id::UInt32
    prec_mz::T
end

getFragMZ(f::FragmentIon) = f.frag_mz
getPepID(f::FragmentIon) = f.pep_id
getPrecMZ(f::FragmentIon) = f.prec_mz

#Make sortable by fragment mz. 
import Base.<
import Base.>
<(y::FragmentIon, x::T) where {T<:Real} = getFragMZ(y) < x
>(y::FragmentIon, x::T) where {T<:Real} = getFragMZ(y) > x

function getSortedFragmentList(peptides::UnorderedDictionary{UInt32, Peptide}, mods_dict::Dict{String, T}; 
                                frag_charges::Vector{UInt8} = UInt8[1, 2], frag_isotopes::Vector{UInt8} = UInt8[0],
                                y_start::Int = 3, b_start::Int = 3, low_mz::Float64 = 300.0, high_mz::Float64 = 1700.0)::Vector{FragmentIon{T}} where {T<:AbstractFloat}
    fragment_list = Vector{FragmentIon{T}}()

    for (id, peptide) in pairs(peptides) #Each iteration of the loop gets adds the fragments for one peptide
        residues = getResidues(getSeq(peptide), mods_dict)::Vector{Residue{T}} #AA Residues for the given peptide. 

        #IMPORTANT: Should only consider +1 precursor charge here. Can calculate m/z for other precursor charge states later. 
        prec_mz = getIonMZ(residues, UInt8(1)) #Charge state for the peptide with +1 charge

        #Get fragments for each precursor that derives from the given peptide 
        #Charge and isotopic state distinguish precursors deriving from the same peptide. 
        for (frag_charge, frag_isotope) in zip(frag_charges, frag_isotopes)
            #Get the fragment ion mz's
            for frag_mz in getFragIons(getResidues(getSeq(peptide), mods_dict), charge = frag_charge, isotope = frag_isotope, b_start = b_start, y_start = y_start)
                if (frag_mz > low_mz) & (frag_mz < high_mz) #Filter fragments with high or low m/z's 
                    push!(fragment_list, FragmentIon(frag_mz, id, prec_mz)) #Add a fragment. 
                end
            end
        end
    end
    sort!(fragment_list, by = x -> getFragMZ(x)) #Sort fragment List. 
    return fragment_list
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
end

getPrecID(pbi::PrecursorBinItem) = pbi.prec_id
getPrecMZ(pbi::PrecursorBinItem) = pbi.prec_mz

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
getPrecursor(pb::PrecursorBin, i::Int64) = pb.precs[i] 
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
    prec_bin::UInt32
end

getLowMZ(fb::FragBin) = fb.lb
getHighMZ(fb::FragBin) = fb.ub
getPrecBinID(fb::FragBin) = fb.prec_bin
FragBin() = FragBin(0.0, 0.0, UInt32(0))


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
    precursor_bins::Vector{PrecursorBin{T}}
end

getFragBins(fi::FragmentIndex) = fi.fragment_bins
getPrecBins(fi::FragmentIndex) = fi.precursor_bins


#FragmentIndex(T::DataType, M::Int, N::Int) =  FragmentIndex(fill(FragBin(), N), fill(PrecursorBin(T, M), N))

function FragmentIndex(T::DataType) 
    return FragmentIndex(Vector{FragBin{T}}(), Vector{PrecursorBin{T}}())
end

getFragmentBin(fi::FragmentIndex, bin::Int) = fi.fragment_bins[bin]

function setFragmentBin!(fi::FragmentIndex{T}, bin::Int64, frag_bin::FragBin{T}) where {T<:AbstractFloat}
    fi.fragment_bins[bin] = frag_bin
end

getPrecursorBin(fi::FragmentIndex, bin::Int64) = fi.precursor_bins[bin]
getPrecursorBinLength(fi::FragmentIndex, bin::Int64) = getLength(fi.precursor_bins[bin])

function setPrecursorBinItem!(fi::FragmentIndex{T}, bin::Int64, index::Int64, prec_bin_item::PrecursorBinItem{T}) where {T<:AbstractFloat}
    setPrecursor!(fi.precursor_bins[bin], index, prec_bin_item)
end

function makeFragmentIndex!(frag_ions::Vector{FragmentIon{T}}, charges::Vector{UInt8} = UInt8[2, 3, 4]) where {T<:AbstractFloat}
   
    #The fragment ions are divided into bins of size N*length(charges). 
    #bin_count = lengthfrag_ions)÷N + 1
    #bin_size = N*length(charges)

    #Pre-allocate the FragmentIndex to the correct size
    #println("0")
    frag_index = FragmentIndex(T) 

    function fillPrecursorBin!(frag_index::FragmentIndex, frag_ions::Vector{FragmentIon{T}}, bin::Int, start::Int, stop::Int)
        i = 1
        for ion_index in range(start, stop)
            
            for charge in charges
                ion = frag_ions[ion_index]
                #Add precursor corresponding to the charge state
                #=setPrecursorBinItem!(frag_index, 
                                    bin, 
                                    i, 
                                    PrecursorBinItem(getPepID(ion), (getPrecMZ(ion)+ PROTON*charge)/charge)
                                    )=#
                frag_index.precursor_bins[bin].precs[i] = PrecursorBinItem(getPepID(ion), (getPrecMZ(ion)+ PROTON*charge)/charge)
                i += 1
            end
        end
    end
    bin = 1
    start = 1
    for stop in 2:length(frag_ions)#bin in 1:(bin_count - 1)
        #Ready to make another precursor bin. 
        if (getFragMZ(frag_ions[stop]) - getFragMZ(frag_ions[start])) > 0.001
            
            #Add a new fragment bin
            #println("A")
            push!(frag_index.fragment_bins, FragBin(getFragMZ(frag_ions[start]),getFragMZ(frag_ions[stop]),UInt32(bin)))
            #setFragmentBin!(frag_index, bin, FragBin(getFragMZ(frag_ions[start]),getFragMZ(frag_ions[stop]),UInt32(bin)));

            #Add a new precursor bin
            #println("B")
            push!(frag_index.precursor_bins, PrecursorBin(Vector{PrecursorBinItem{T}}(undef, (stop - start + 1)*length(charges)))) #PrecursorBin(T, start - stop))
            #println("C")
            fillPrecursorBin!(frag_index, frag_ions, bin, start, stop)

            
            sort!(getPrecursors(getPrecursorBin(frag_index, bin)), by = x->getPrecMZ(x));
            bin += 1
            start = stop 
        end
    end
    return frag_index
end


using FASTX
using CodecZlib
using Dictionaries
file_path = "/Users/n.t.wamsley/RIS_temp/HAMAD_MAY23/mouse_SIL_List/UP000000589_10090.fasta.gz"
@time begin
    @time begin 
        peptides_fasta = digestFasta(parseFasta(file_path))
        test_table = PrecursorTable()
        fixed_mods = [(p=r"C", r="C[Carb]")]
        var_mods = [(p=r"(K$)", r="[Hlys]"), (p=r"(R$)", r="[Harg]")]
        buildPrecursorTable!(test_table, peptides_fasta, fixed_mods, var_mods, 2)
    end

    const mods_dict = Dict("Carb" => Float64(57.021464),
    "Harg" => Float64(10.008269),
    "Hlys" => Float64(8.014199),
    "Hglu" => Float64(6))

    @time f_list = getSortedFragmentList(test_table.id_to_pep, mods_dict);
    @time f_index = makeFragmentIndex!(test_table, f_list, 256);
end

