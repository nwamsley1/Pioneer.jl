function buildDesignMatrix(matches::Vector{FragmentMatch{Float32}},  misses::Vector{FragmentMatch{Float32}}, seed_size::Int) #where {T<:AbstractFloat}

    #Number of unique matched peaks.
    #Should be getPeakInd(matches[end])???
    #Number of rows equals the number of unique matched peaks
    M = length(unique([getPeakInd(x) for x in matches]))
    #Design matrix. One row for every precursor template. One column for every matched peak. 
    H = zeros(Float32, (seed_size, M))
    #Design matrix. One row for every precursor template. One column for every unmatched peak. 
    UNMATCHED = zeros(Float32, (seed_size, length(misses)))

    #Spectrum/empirical intensities for each peak. Zero by default (for unmatched/missed fragments)
    X = zeros(Float32, (1, M))

    #Maps a precursor id to a row of H. 
    precID_to_row = UnorderedDictionary{UInt32, UInt8}()

    #Current highest row encountered
    prec_row = UInt8(0)
    col = 0
    #Number of unique peaks encountered. 
    last_peak_ind = 0
    for match in matches
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(match))
            prec_row  += UInt8(1)
            insert!(precID_to_row, getPrecID(match), prec_row )
        end

        #If this peak has not been encountered yet, then start filling a new column
        if getPeakInd(match) != last_peak_ind
            col += 1
            last_peak_ind = getPeakInd(match)
        end

        row = precID_to_row[getPrecID(match)]
        H[row, col] = getPredictedIntenisty(match)
        X[1, col] = getIntensity(match)
    end

    col = 0
    last_peak_ind = 0
    for miss in misses
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(miss))
            prec_row  += UInt8(1)
            insert!(precID_to_row, getPrecID(miss), prec_row)
        end
        if getPeakInd(miss) != last_peak_ind
            col += 1
            last_peak_ind = getPeakInd(miss)
        end
        row = precID_to_row[getPrecID(miss)]
        UNMATCHED[row, col] = getPredictedIntenisty(miss)
    end

    return X, H, UNMATCHED, precID_to_row
end