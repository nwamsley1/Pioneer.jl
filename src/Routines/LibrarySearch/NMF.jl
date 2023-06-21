function buildDesignMatrix(matches::Vector{FragmentMatch{Float32}},  misses::Vector{FragmentMatch{Float32}}, seed_size::Int) #where {T<:AbstractFloat}

    #Number of unique matched peaks.
    matched_peaks = length(unique([getPeakInd(x) for x in matches]))
    #Number of rows equals the number of unique matched peaks + the number of expected fragments that 
    #failed to match a peak in the spectrm
    M = (matched_peaks + length(misses))
    #Design matrix. One row for every precursor template. One column for every matched + missed peak. 
    H = zeros(Float32, (seed_size, M))
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

    for miss in misses
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(miss))
            prec_row  += UInt8(1)
            insert!(precID_to_row, getPrecID(miss), prec_row)
        end
        col = getPeakInd(miss) + matched_peaks
        row = precID_to_row[getPrecID(miss)]
        H[row, col] = getPredictedIntenisty(miss)
    end

    return X, H, precID_to_row
end