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
<(x::T, y::FragmentIon) where {T<:Real} = <(y, x)
>(y::FragmentIon, x::T) where {T<:Real} = getFragMZ(y) > x
>(x::T, y::FragmentIon) where {T<:Real} = >(y, x)

function getSortedFragmentList(peptides::UnorderedDictionary{UInt32, Peptide}, mods_dict::Dict{String, T}; 
                                frag_charges::Vector{UInt8} = UInt8[1], frag_isotopes::Vector{UInt8} = UInt8[0],
                                y_start::Int = 3, b_start::Int = 3, low_mz::Float64 = 150.0, high_mz::Float64 = 1700.0)::Vector{FragmentIon{T}} where {T<:AbstractFloat}
    
    #List of all fragment ions from the `peptides`
    fragment_list = Vector{FragmentIon{T}}()

    for (id, peptide) in pairs(peptides) #Each iteration of the loop adds the fragments for one peptide
        residues = getResidues(getSeq(peptide), mods_dict)::Vector{Residue{T}} #AA Residues for the given peptide. 

        #IMPORTANT: Should only consider +1 precursor charge here. Can calculate m/z for other precursor charge states later. 
        prec_mz = getIonMZ(residues, UInt8(1)) #Charge state for the Precursor with +1 charge (monoisotopic). See "src/precursor.jl"

        #Get fragments for each precursor that derives from the given peptide 
        #Charge and isotopic state distinguish precursors deriving from the same peptide. 
        for (frag_charge, frag_isotope) in zip(frag_charges, frag_isotopes)
            #Get the fragment ion mz's
            for frag_mz in getFragIons(residues, charge = frag_charge, isotope = frag_isotope, b_start = b_start, y_start = y_start)
                if (frag_mz > low_mz) & (frag_mz < high_mz) #Filter fragments with high or low m/z's 
                    push!(fragment_list, FragmentIon(frag_mz, id, prec_mz)) #Add a fragment. 
                end
            end
        end
    end

    sort!(fragment_list, by = x -> getFragMZ(x)) #Sort fragment list by m/z

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
    charge::UInt8
end

getPrecID(pbi::PrecursorBinItem) = pbi.prec_id
getPrecMZ(pbi::PrecursorBinItem) = pbi.prec_mz
getPrecCharge(pbi::PrecursorBinItem) = pbi.charge
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


function addFragmentBin!(fi::FragmentIndex{T}, frag_bin::FragBin{T}) where {T<:AbstractFloat}
    push!(getFragBins(fi), frag_bin)
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


function buildFragmentIndex!(frag_ions::Vector{FragmentIon{T}}, bin_ppm::T; max_charge::UInt8 = UInt8(4), min_charge::UInt8 = UInt8(2), low_prec_mz::Float64 = 300.0, high_prec_mz::Float64 = 1100.0) where {T<:AbstractFloat}
   
    #The fragment ions are divided into bins of roughtly equal m/z width.
    #That should correspond to roughly half the fragment mass accuracy of the detector?
    frag_index = FragmentIndex(T) 

    function fillPrecursorBin!(frag_index::FragmentIndex, frag_ions::Vector{FragmentIon{T}}, max_charge::UInt8, min_charge::UInt8, bin::UInt32, start::Int, stop::Int, low_prec_mz::Float64, high_prec_mz::Float64)
        i = 1 #Index of current fragme nt
        for ion_index in range(start, stop)
            pep_id = getPepID(frag_ions[ion_index])
            for charge in UInt8[1]#charges #For each precursor charge 
                prec_mz = getPrecMZ(frag_ions[ion_index])#(getPrecMZ(frag_ions[ion_index]) + PROTON*(charge-1))/charge #m/z of the precursor
                if (prec_mz < low_prec_mz/max_charge) | (prec_mz > high_prec_mz*min_charge) #Precursor m/z outside the bounds
                    continue
                end
                #Add precursor corresponding to the charge state
                addPrecursorBinItem!(frag_index,
                                    bin,
                                    #i,
                                    PrecursorBinItem(pep_id, prec_mz, charge)
                                    )
                #frag_index.precursor_bins[bin].precs[i] = PrecursorBinItem(getPepID(ion), (getPrecMZ(ion) + PROTON*(charge-1))/charge)
                i += 1 #Move to the next fragment 
            end
        end
    end

    bin = UInt32(1) #Current fragment bin index
    start = 1 #Fragment index of the first fragment in the current bin

    getPPM(frag_mz::T, ppm::T) = ppm*frag_mz/1e6

    diff = getPPM(getFragMZ(frag_ions[start]), bin_ppm) #ppm tolerance of the current fragment bin

    #Build bins 
    for stop in 2:length(frag_ions)

        #Ready to make another precursor bin. 
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
            addPrecursorBin!(frag_index, 
                                #Preallocate an empty precursor bin of the correct length 
                                PrecursorBin(Vector{PrecursorBinItem{T}}())#undef, (last_frag_in_bin - start + 1)*length(charges)))
                                )
        
            fillPrecursorBin!(frag_index, frag_ions, max_charge, min_charge, bin, start, last_frag_in_bin, low_prec_mz, high_prec_mz)

            #Sort the precursor bin by precursor m/z
            sort!(getPrecursors(getPrecursorBin(frag_index, bin)), by = x->getPrecMZ(x));

            #Update counters and ppm tolerance 
            bin += UInt32(1)
            start = stop
            diff = getPPM(getFragMZ(frag_ions[start]), bin_ppm)
        end
    end
    return frag_index
end


using FASTX
using CodecZlib
using Dictionaries
#Example gzipped mouse proteome. 
file_path = "/Users/n.t.wamsley/RIS_temp/HAMAD_MAY23/mouse_SIL_List/UP000000589_10090.fasta.gz"
@time begin
    @time begin 
        peptides_fasta = digestFasta(parseFasta(file_path))
        test_table = PrecursorTable()
        fixed_mods = [(p=r"C", r="C[Carb]")]
        var_mods = [(p=r"(K$)", r="[Hlys]"), (p=r"(R$)", r="[Harg]"), (p=r"(M)", r="[MOx]")]
        buildPrecursorTable!(test_table, peptides_fasta, fixed_mods, var_mods, 2)
    end

    const mods_dict = Dict("Carb" => Float64(57.021464),
    "Harg" => Float64(10.008269),
    "Hlys" => Float64(8.014199),
    "MOx" => Float64(15.9949))


    @time f_list = getSortedFragmentList(test_table.id_to_pep, mods_dict);
    @time f_index = buildFragmentIndex!(f_list, 10.0);
end

println("Size of precursor bins in GB: ", (sum([sum([sizeof(prec) for prec in prec_bin.precs]) for prec_bin in f_index.precursor_bins]))/1e9)
println("Size of fragment list in GB: ", sum([sizeof(x) for x in f_list])/1e9)
println("Size of fragment bins in GB: ", sum([sizeof(x) for x in f_index.fragment_bins])/1e9)

using Plots
bin_sizes = [length(prec_bin.precs) for prec_bin in f_index.precursor_bins]
n = length(bin_sizes)
p = plot((1:n), cumsum(sort(bin_sizes))/sum(bin_sizes),
    xlabel = "sample", ylabel = "Probability", 
    title = "Empirical Cumluative Distribution", label = "")
histogram(bin_sizes[(bin_sizes.<4e4) .& (bin_sizes.>1000)])
histogram(bin_sizes)
f_list = nothing
using Dictionaries
using Combinatorics
using Random
using Arrow
using Tables

include("src/precursor.jl")
include("src/parseFASTA.jl")
include("src/PrecursorDatabase.jl")
include("src/applyMods.jl")


open("/Users/n.t.wamsley/Desktop/targets.csv", "w") do file
    write(file, "modified_sequence,collision_energy,precursor_charge\n")
    for (id, pep) in pairs(test_table.id_to_pep)
        sequence = replace(getSeq(pep), r"M\[MOx\]"=>"M(ox)")
        sequence = replace(sequence, r"C\[Carb\]"=>"C")

        if (occursin("[H", sequence)) | (occursin("U", sequence)) | (occursin("O", sequence)) |  (occursin("X", sequence)) | occursin("Z", getSeq(pep)) | occursin("B", getSeq(pep))
            continue
        end

        if (length(sequence) > 30) | (length(sequence) < 7)
            continue
        end 

        if (isDecoy(pep) == false)
            for charge in [2, 3, 4]
                write(file, "$sequence, 33, $charge \n")
            end
        end
    end
end

open("/Users/n.t.wamsley/Desktop/decoys.csv", "w") do file
    write(file, "modified_sequence,collision_energy,precursor_charge\n")
    for (id, pep) in pairs(test_table.id_to_pep)
        sequence = replace(getSeq(pep), r"M\[MOx\]"=>"M(ox)")
        sequence = replace(sequence, r"C\[Carb\]"=>"C")
        if (occursin("[H", sequence)) | (occursin("U", sequence)) | (occursin("O", sequence)) |  (occursin("X", sequence)) | occursin("Z", getSeq(pep)) | occursin("B", getSeq(pep))
            continue
        end
        if (length(sequence) > 30) | (length(sequence) < 7)
            continue
        end 
        if (isDecoy(pep) == true)
            for charge in [2, 3, 4]
                write(file, "$sequence, 33, $charge \n")
            end
        end
    end
end

for (id, pep) in pairs(test_table.id_to_pep)
    if occursin("Z", getSeq(pep))
        println(getSeq(pep))
    end
end

