function buildDesignMatrix(matches::Vector{FragmentMatch{Float32}},  misses::Vector{FragmentMatch{Float32}}) #where {T<:AbstractFloat}

    
    #Number of rows equals the number of unique matched peaks
    #Remember "getPeakInd(x)" is hte index of the matched peak in the MS2 spectrum.
    M = length(unique([getPeakInd(x) for x in matches]))

    #Design matrix. One row for every precursor template. One column for every matched peak. 
    H_COLS = zeros(Int64, length(matches))
    H_ROWS = zeros(Int64, length(matches))
    H_VALS = zeros(Float32, length(matches))

    #Design matrix. One row for every precursor template. One column for every unmatched peak. 
    U_COLS = zeros(Int64, length(misses))
    U_ROWS = zeros(Int64, length(misses))
    U_VALS = zeros(Float32, length(misses))

    #Spectrum/empirical intensities for each peak. Zero by default (for unmatched/missed fragments)
    X = zeros(Float32, M)

    #Maps a precursor id to a row of H. 
    precID_to_row = UnorderedDictionary{UInt32, UInt32}()

    #Current highest row encountered
    prec_row = zero(UInt32)
    col = 0
    #Number of unique peaks encountered. 
    last_peak_ind = 0
    for (i, match) in enumerate(matches)
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(match))
            prec_row += one(UInt32)
            insert!(precID_to_row, getPrecID(match), prec_row)
        end

        #If this peak has not been encountered yet, then start filling a new column
        if getPeakInd(match) != last_peak_ind
            col += 1
            last_peak_ind = getPeakInd(match)
            X[col] = getIntensity(match)
        end
        row = precID_to_row[getPrecID(match)]
        H_COLS[i] = col
        H_ROWS[i] = row
        H_VALS[i] = getPredictedIntenisty(match)
    end
    H_ncol = col
    #col = 0
    last_peak_ind = 0
    for (i, miss) in enumerate(misses)
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
        U_COLS[i] = col
        U_ROWS[i] = row
        U_VALS[i] = getPredictedIntenisty(miss)
    end

    return append!(X, zeros(Float32, length(U_VALS))), sparse(vcat(H_COLS, U_COLS), vcat(H_ROWS, U_ROWS), vcat(H_VALS, U_VALS)), sparse(vcat(H_ROWS, U_ROWS), vcat(H_COLS, U_COLS), vcat(H_VALS, U_VALS)), precID_to_row, H_ncol
end
#=

function buildDesignMatrix(matches::Vector{FragmentMatch{Float32}},  misses::Vector{FragmentMatch{Float32}}, seed_size::Int) #where {T<:AbstractFloat}

    #Number of unique matched peaks.
    #Should be getPeakInd(matches[end])???
    #Number of rows equals the number of unique matched peaks
    M = length(unique([getPeakInd(x) for x in matches]))
    #Design matrix. One row for every precursor template. One column for every matched peak. 
    H_COLS = zeros(Int64, length(matches))
    H_ROWS = zeros(Int64, length(matches))
    H_VALS = zeros(Float32, length(matches))
    #H = zeros(Float32, (seed_size, M))
    #Design matrix. One row for every precursor template. One column for every unmatched peak. 
    U_COLS = zeros(Int64, length(misses))
    U_ROWS = zeros(Int64, length(misses))
    U_VALS = zeros(Float32, length(misses))

    #Spectrum/empirical intensities for each peak. Zero by default (for unmatched/missed fragments)
    X = zeros(Float32, M)

    #Maps a precursor id to a row of H. 
    precID_to_row = UnorderedDictionary{UInt32, UInt32}()

    #Current highest row encountered
    prec_row = zero(UInt32)
    col = 0
    #Number of unique peaks encountered. 
    last_peak_ind = 0
    for (i, match) in enumerate(matches)
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(match))
            prec_row += one(UInt32)
            insert!(precID_to_row, getPrecID(match), prec_row)
        end

        #If this peak has not been encountered yet, then start filling a new column
        if getPeakInd(match) != last_peak_ind
            col += 1
            last_peak_ind = getPeakInd(match)
            X[col] = getIntensity(match)
        end
        row = precID_to_row[getPrecID(match)]
        H_COLS[i] = col
        H_ROWS[i] = row
        H_VALS[i] = getPredictedIntenisty(match)
    end
    H_ncol = col
    col = 0
    last_peak_ind = 0
    for (i, miss) in enumerate(misses)
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
        U_COLS[i] = col
        U_ROWS[i] = row
        U_VALS[i] = getPredictedIntenisty(miss)
    end
    U_ncol = col
    nrow = length(keys(precID_to_row))
    return X, SparseMatrix(nrow, H_ncol, H_COLS, H_ROWS, H_VALS), SparseMatrix(nrow, U_ncol, U_COLS, U_ROWS, U_VALS), precID_to_row
end


#=function buildDesignMatrix(matches::Vector{FragmentMatch{Float32}},  misses::Vector{FragmentMatch{Float32}}, seed_size::Int) #where {T<:AbstractFloat}

    #Number of unique matched peaks.
    #Should be getPeakInd(matches[end])???
    #Number of rows equals the number of unique matched peaks
    M = length(unique([getPeakInd(x) for x in matches]))
    #Design matrix. One row for every precursor template. One column for every matched peak. 
    H = zeros(Float32, (seed_size, M))
    #H = zeros(Float32, (seed_size, M))
    #Design matrix. One row for every precursor template. One column for every unmatched peak. 
    UNMATCHED = zeros(Float32, (seed_size, length(misses)))

    #Spectrum/empirical intensities for each peak. Zero by default (for unmatched/missed fragments)
    X = zeros(Float32, (1, M))

    #Maps a precursor id to a row of H. 
    precID_to_row = UnorderedDictionary{UInt32, UInt32}()

    #Current highest row encountered
    prec_row = zero(UInt32)
    col = 0
    #Number of unique peaks encountered. 
    last_peak_ind = 0
    for match in matches
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(match))
            prec_row += one(UInt32)
            insert!(precID_to_row, getPrecID(match), prec_row)
        end

        #If this peak has not been encountered yet, then start filling a new column
        if getPeakInd(match) != last_peak_ind
            col += 1
            last_peak_ind = getPeakInd(match)
            X[1, col] = getIntensity(match)
        end
        row = precID_to_row[getPrecID(match)]
        H[row, col] = getPredictedIntenisty(match)
        
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
struct SparseMatrix{Ti<:Integer, Tv<:AbstractFloat}
    m::Int                  # Number of rows
    n::Int                  # Number of columns
    colptr::Vector{Ti}      # Column indices of stored values
    rowval::Vector{Ti}      # Row indices of stored values
    nzval::Vector{Tv}       # Stored values, typically nonzeros
end
=#

function buildDesignMatrix(matches::Vector{FragmentMatch{Float32}}) #where {T<:AbstractFloat}

    #Number of unique matched peaks.
    #Should be getPeakInd(matches[end])???
    #Number of rows equals the number of unique matched peaks
    M = length(unique([getPeakInd(x) for x in matches]))
    #Design matrix. One row for every precursor template. One column for every matched peak. 
    H_COLS = zeros(Int64, length(matches))
    H_ROWS = zeros(Int64, length(matches))
    H_VALS = zeros(Float32, length(matches))


    #Spectrum/empirical intensities for each peak. Zero by default (for unmatched/missed fragments)
    X = zeros(Float32, M)

    #Maps a precursor id to a row of H. 
    precID_to_row = UnorderedDictionary{UInt32, UInt32}()

    #Current highest row encountered
    prec_row = zero(UInt32)
    col = 0
    #Number of unique peaks encountered. 
    last_peak_ind = 0
    for (i, match) in enumerate(matches)
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(match))
            prec_row += one(UInt32)
            insert!(precID_to_row, getPrecID(match), prec_row)
        end

        #If this peak has not been encountered yet, then start filling a new column
        if getPeakInd(match) != last_peak_ind
            col += 1
            last_peak_ind = getPeakInd(match)
            X[col] = getIntensity(match)
        end
        row = precID_to_row[getPrecID(match)]
        H_COLS[i] = col
        H_ROWS[i] = row
        H_VALS[i] = getPredictedIntenisty(match)
    end
    H_ncol = col
    nrow = length(keys(precID_to_row))
    return X, sparse(H_COLS, H_ROWS, H_VALS), sparse(H_ROWS, H_COLS, H_VALS), precID_to_row
end

function buildDesignMatrix(matches::Vector{FragmentMatch{Float32}}) #where {T<:AbstractFloat}

    #Number of unique matched peaks.
    #Should be getPeakInd(matches[end])???
    #Number of rows equals the number of unique matched peaks
    M = length(unique([getPeakInd(x) for x in matches]))
    N = length(unique([getPrecID(x) for x in matches]))
    #Design matrix. One row for every precursor template. One column for every matched peak. 
    H = zeros(Float32, (N, M))
    #Design matrix. One row for every precursor template. One column for every unmatched peak. 
    #UNMATCHED = zeros(Float32, (seed_size, length(misses)))

    #Spectrum/empirical intensities for each peak. Zero by default (for unmatched/missed fragments)
    X = zeros(Float32, (1, M))

    #Maps a precursor id to a row of H. 
    precID_to_row = UnorderedDictionary{UInt32, UInt32}()

    #Current highest row encountered
    prec_row = zero(UInt32)
    col = 0
    #Number of unique peaks encountered. 
    last_peak_ind = 0
    for match in matches
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(match))
            prec_row += one(UInt32)
            insert!(precID_to_row, getPrecID(match), prec_row)
        end

        #If this peak has not been encountered yet, then start filling a new column
        if getPeakInd(match) != last_peak_ind
            col += 1
            last_peak_ind = getPeakInd(match)
            X[1, col] = getIntensity(match)
        end
        row = precID_to_row[getPrecID(match)]
        H[row, col] = getPredictedIntenisty(match)
        
    end

    return X, H, precID_to_row
end=#

