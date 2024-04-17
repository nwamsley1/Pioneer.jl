function buildDesignMatrix!(H::SparseArray{UInt32,Float32}, matches::Vector{m},  misses::Vector{m}, nmatches::Int64, nmisses::Int64, precID_to_col::ArrayDict{UInt32, UInt16}; block_size = 10000) where {m<:MatchIon{Float32}}
    T = Float32
    #Number of rows equals the number of unique matched peaks
    #Remember "getPeakInd(x)" is hte index of the matched peak in the MS2 spectrum.
    M = 1
    for i in range(2, nmatches)
        if getPeakInd(matches[i])!= getPeakInd(matches[i - 1])
            M += 1
        end
    end
    M += nmisses

    #If M exceeds pre-allocated size
    if (nmatches + nmisses) >= length(H.colval) - 1
        block_size = max(block_size,nmatches + nmisses - length(H.colval))
        append!(H.colval, zeros(eltype(H.colval), block_size))
        append!(H.rowval, zeros(eltype(H.rowval), block_size))
        append!(H.nzval, zeros(eltype(H.nzval), block_size))
        append!(H.x, zeros(eltype(H.x), block_size))
        append!(H.matched, zeros(eltype(H.matched), block_size))
        #MUST BE ONES!!!
        append!(H.mask, ones(eltype(H.mask), block_size))
        append!(H.colptr, Vector{UInt32}(undef, block_size))
    end
    #Spectrum/empirical intensities for each peak. Zero by default (for unmatched/missed fragments)
    #X = zeros(T, M)
    #Maps a precursor id to a row of H. 

    #Current highest row encountered
    prec_col = zero(UInt32)
    row = 0
    #Number of unique peaks encountered. 
    last_peak_ind = 0
    for i in range(1, nmatches)#matches)
        match = matches[i]
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if iszero(precID_to_col[getPrecID(match)])
            prec_col += one(UInt32)
            update!(precID_to_col, getPrecID(match), UInt16(prec_col))
        end

        #If this peak has not been encountered yet, then start filling a new 
        if getPeakInd(match) != last_peak_ind
            row += 1
            last_peak_ind = getPeakInd(match)
            #X[row] = getIntensity(match)
        end

        H.colval[i] = precID_to_col[getPrecID(match)]
        H.rowval[i] = row
        H.nzval[i] = getPredictedIntenisty(match)
        H.x[i] = getIntensity(match)#X[row]
        H.matched[i] = true
        H.n_vals += 1

    end
    last_peak_ind = 0
    for i in range(nmatches + 1, nmatches + nmisses)
        miss = misses[i - nmatches]
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if iszero(precID_to_col[getPrecID(miss)])
            prec_col += one(UInt32)
            update!(precID_to_col, getPrecID(miss), UInt16(prec_col))
        end

        # if getPeakInd(miss) != last_peak_ind
        row += 1
        last_peak_ind = getPeakInd(miss)
        #end
        #H.row_col_nzval_x[i] = FRAG(row, Int64(first(precID_to_col[getPrecID(miss)])), getPredictedIntenisty(miss), zero(U))
        H.colval[i] = precID_to_col[getPrecID(miss)]
        #println("row $row")
        H.rowval[i] = row
        H.nzval[i] = getPredictedIntenisty(miss) #factor to reduce impact of unmatched ions 
        H.x[i] = zero(Float32)
        H.matched[i] = false

        H.n_vals += 1
    end
    #if i >= (nmatches - 1)
    #return X, sparse(vcat(H_COLS, U_COLS), vcat(H_ROWS, U_ROWS), vcat(H_VALS, U_VALS)), sparse(vcat(H_ROWS, U_ROWS), vcat(H_COLS, U_COLS), vcat(H_VALS, U_VALS)), precID_to_row, H_ncol
    #end
    #return X, sparse(H_COLS, H_ROWS, H_VALS), sparse(H_ROWS, H_COLS, H_VALS), precID_to_row, H_ncol
    #sortSparse!(H)
    #sort!(@view(matches[1:nmatches]), by = x->getPrecID(x), alg = PartialQuickSort(1:nmatches))
    sortSparse!(H)
end
