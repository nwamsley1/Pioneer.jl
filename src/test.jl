
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
            #println("A")
            push!(frag_index.precursor_bins, PrecursorBin(Vector{PrecursorBinItem{T}}(undef, start - stop))) #PrecursorBin(T, start - stop))
            #println("B")
            fillPrecursorBin!(frag_index, frag_ions, bin, start, stop)

            sort!(getPrecursors(getPrecursorBin(frag_index, bin)), by = x->getPrecMZ(x));
            bin += 1
            start = i 
        end
    end

    #Last bin is special case
    #=bin = bin_count 
    start = 1 + (bin-1)*N
    stop = length(frag_ions)
    setFragmentBin!(frag_index, bin, FragBin(getFragMZ(frag_ions[start]),getFragMZ(frag_ions[stop]),UInt32(bin)));
    fillPrecursorBin!(frag_index, frag_ions, bin, start, stop)
    sort!(getPrecursors(getPrecursorBin(frag_index, bin)), by = x->getPrecMZ(x));=#

    return frag_index
end