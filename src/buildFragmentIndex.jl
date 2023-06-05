struct FragmentIon{T<:AbstractFloat}
    frag_mz::T
    pep_id::UInt32
    prec_mz::T
end

getFragMZ(f::FragmentIon) = f.frag_mz
getPepID(f::FragmentIon) = f.pep_id
getPrecMZ(f::FragmentIon) = f.prec_mz

import Base.<
import Base.>
<(y::FragmentIon, x::T) where {T<:Real} = getFragMZ(y) < x
>(y::FragmentIon, x::T) where {T<:Real} = getFragMZ(y) > x

function getSortedFragmentList(peptides::UnorderedDictionary{UInt32, Peptide}, mods_dict::Dict{String, T}; 
                                frag_charges::Vector{UInt8} = UInt8[1, 2], frag_isotopes::Vector{UInt8} = UInt8[0],
                                y_start::Int = 3, b_start::Int = 3, low_mz::Float64 = 200.0, high_mz::Float64 = 1700.0)::Vector{FragmentIon{T}} where {T<:AbstractFloat}
    fragment_list = Vector{FragmentIon{T}}()
    #fragment_list = Vector{Tuple{Float64, UInt32, Float64}}()
    for (id, peptide) in pairs(peptides)
        residues = getResidues(getSeq(peptide), mods_dict)
        prec_mz = getIonMZ(residues, UInt8(1))
        for (frag_charge, frag_isotope) in zip(frag_charges, frag_isotopes)
            for frag_mz in getFragIons(getResidues(getSeq(peptide), mods_dict), charge = frag_charge, isotope = frag_isotope, b_start = b_start, y_start = y_start)
                if (frag_mz > low_mz) & (frag_mz < high_mz)
                    push!(fragment_list, FragmentIon(frag_mz, id, prec_mz))
                end
            end
        end
    end
    sort!(fragment_list, by = x -> getFragMZ(x))
    return fragment_list
end

struct PrecursorBinItem{T<:AbstractFloat}
    prec_id::UInt32
    prec_mz::T
end

getPrecID(pbi::PrecursorBinItem) = pbi.prec_id
getPrecMZ(pbi::PrecursorBinItem) = pbi.prec_mz

struct PrecursorBin{T<:AbstractFloat}
    precs::Vector{PrecursorBinItem{T}}
end

getPrecursors(pb::PrecursorBin) = pb.precs

function setPrecursor!(pb::PrecursorBin, index::Int, pbi::PrecursorBinItem)
    pb.precs[index] = pbi
end

PrecursorBin(T::DataType, N::Int) = PrecursorBin(Vector{PrecursorBinItem{T}}(undef, N))

struct FragBin{T<:AbstractFloat}
    lb::T
    ub::T
    prec_bin::UInt32
end

getLB(fb::FragBin) = fb.lb
getUB(fb::FragBin) = fb.ub
FragBin() = FragBin(0.0, 0.0, UInt32(0))

struct FragmentIndex{T<:AbstractFloat}
    fragment_bins::Vector{FragBin{T}}
    precursor_bins::Vector{PrecursorBin{T}}
end

FragmentIndex(T::DataType, M::Int, N::Int) = FragmentIndex(fill(FragBin(), N), fill(PrecursorBin(T, M), N))

getFragmentBin(fi::FragmentIndex, bin::Int) = fi.fragment_bins[bin]

function setFragmentBin!(fi::FragmentIndex, bin::Int64, frag_bin::FragBin)
    fi.fragment_bins[bin] = frag_bin
end

getPrecursorBin(fi::FragmentIndex, bin::Int64) = fi.precursor_bins[bin]

function setPrecursorBinItem!(fi::FragmentIndex{T}, bin::Int64, index::Int64, prec_bin_item::PrecursorBinItem{T}) where {T<:AbstractFloat}
    setPrecursor!(fi.precursor_bins[bin], index, prec_bin_item)
end

function makeFragmentIndex!(ptable::PrecursorTable, frag_ions::Vector{FragmentIon{T}}, N::Int = 32, charges::Vector{UInt8} = UInt8[2, 3, 4]) where {T<:AbstractFloat}
   
    bin_count = length(frag_ions)÷N + 1
    bin_size = N*length(charges)
    frag_index = FragmentIndex(T, bin_size, bin_count) 

    function fillPrecursorBin!(frag_index::FragmentIndex, frag_ions::Vector{FragmentIon{T}}, bin::Int, test_start::Int, stop::Int)
        i = 1
        for ion_index in start:stop
            
            for charge in charges
                ion = frag_ions[ion_index]
                #Add precursor corresponding to the charge state
                setPrecursorBinItem!(frag_index, 
                                    bin, 
                                    i, 
                                    PrecursorBinItem(getPepID(ion), (getPrecMZ(ion)+ PROTON*charge)/charge)
                                    )
                i += 1
            end
        end
    end

    for bin in 1:(bin_count - 1)
        start = 1 + (bin - 1)*N
        stop = (bin)*N
        setFragmentBin!(frag_index, bin, FragBin(getFragMZ(frag_ions[start]),getFragMZ(frag_ions[stop]),UInt32(bin)));
        #prec_bin = getPrecursorBin(frag_index, bin)#frag_index.preursor_bins[bin]
        #Frag ions only includes precursor_mz for the +1 charge state.
        #We want the precursor bin to include the precursor_mzs for all charge states in `charges`
        fillPrecursorBin!(frag_index, frag_ions, bin, start, stop)
        #Sort fragment ions in the precursor bin by their precursor m/z
        sort!(getPrecursors(getPrecursorBin(frag_index, bin)), by = x->getPrecMZ(x));
    end

    #Last bin is special case
    bin = bin_count 
    start = 1 + (bin-1)*N
    stop = length(frag_ions)
    setFragmentBin!(frag_index, bin, FragBin(getFragMZ(frag_ions[start]),getFragMZ(frag_ions[stop]),UInt32(bin)));
    fillPrecursorBin!(frag_index, frag_ions, bin, start, stop)
    sort!(getPrecursors(getPrecursorBin(frag_index, bin)), by = x->getPrecMZ(x));

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

function query_frag_index(frag_index::Vector{FragBin{T}}, query::T) where {T<:AbstractFloat}
    lo, hi = 1, length(frag_index)
    while lo <= hi
        mid = (lo + hi) ÷ 2
        if getUB(frag_index[mid]) < (query)
             lo = mid + 1
        elseif getLB(frag_index[mid]) > (query)
            hi = mid - 1
        else
            return lo
        end
    end
end
BIN = f_index.fragment_bins
QUERY = 1500.0
FUNC = (t,x)->x.lb>t