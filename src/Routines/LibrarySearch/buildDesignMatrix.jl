function buildDesignMatrix!(H::SparseArray{Ti,U}, matches::Vector{m},  misses::Vector{m}, nmatches::Int64, nmisses::Int64, precID_to_col::ArrayDict{UInt32, UInt16}; block_size = 10000) where {m<:Match,Ti<:Integer,U<:AbstractFloat}
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
        append!(H.colptr, Vector{Ti}(undef, block_size))
    end
    #Spectrum/empirical intensities for each peak. Zero by default (for unmatched/missed fragments)
    X = zeros(T, M)
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
            X[row] = getIntensity(match)
        end

        #H.row_col_nzval_x[i] = FRAG(row, Int64(first(precID_to_col[getPrecID(match)])), getPredictedIntenisty(match), X[row])
        H.colval[i] = Int64(precID_to_col[getPrecID(match)])
        H.rowval[i] = row
        H.nzval[i] = getPredictedIntenisty(match)
        H.x[i] = X[row]
        H.matched[i] = true
        H.n_vals += 1
            #println("i ", i)
    end
    H_nrow = row
    #col = 0
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
        H.colval[i] = Int64(precID_to_col[getPrecID(miss)])
        H.rowval[i] = row
        H.nzval[i] = Float32(getPredictedIntenisty(miss)/2) #factor to reduce impact of unmatched ions 
        H.x[i] = zero(U)
        H.matched[i] = false

        H.n_vals += 1
    end
    #if i >= (nmatches - 1)
    #return X, sparse(vcat(H_COLS, U_COLS), vcat(H_ROWS, U_ROWS), vcat(H_VALS, U_VALS)), sparse(vcat(H_ROWS, U_ROWS), vcat(H_COLS, U_COLS), vcat(H_VALS, U_VALS)), precID_to_row, H_ncol
    #end
    #return X, sparse(H_COLS, H_ROWS, H_VALS), sparse(H_ROWS, H_COLS, H_VALS), precID_to_row, H_ncol
    #sortSparse!(H)
    sortSparse!(H)
end

function buildDesignMatrix(matches::Vector{m},  misses::Vector{m}, nmatches::Int64, nmisses::Int64, H_COLS::Vector{Int64}, H_ROWS::Vector{Int64}, H_VALS::Vector{Float32}; block_size = 10000) where {m<:Match}
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
    if (nmatches + nmisses) >= length(H_COLS) - 1
        H_COLS = zeros(Int64, length(H_COLS) + block_size)
        H_ROWS = zeros(Int64, length(H_ROWS) + block_size)
        H_VALS = zeros(Float32, length(H_VALS) + block_size)
    end
    #Spectrum/empirical intensities for each peak. Zero by default (for unmatched/missed fragments)
    X = zeros(T, M)
    #Maps a precursor id to a row of H. 
    precID_to_row = UnorderedDictionary{UInt32, Tuple{UInt32, UInt8}}()

    #Current highest row encountered
    prec_row = zero(UInt32)
    col = 0
    #Number of unique peaks encountered. 
    last_peak_ind = 0
    for i in range(1, nmatches)#matches)
        match = matches[i]
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(match))
            prec_row += one(UInt32)
            insert!(precID_to_row, getPrecID(match), (prec_row, zero(UInt8)))
        end
        if getRank(match) < 3
            precID_to_row[getPrecID(match)] = (precID_to_row[getPrecID(match)][1], precID_to_row[getPrecID(match)][2] + one(UInt8))
        end
        #If this peak has not been encountered yet, then start filling a new column
        if getPeakInd(match) != last_peak_ind
            col += 1
            last_peak_ind = getPeakInd(match)
            X[col] = getIntensity(match)
        end
        row = precID_to_row[getPrecID(match)][1]
        H_COLS[i] = col
        H_ROWS[i] = row
        H_VALS[i] = getPredictedIntenisty(match)
            #println("i ", i)
    end
    H_ncol = col
    #col = 0
    last_peak_ind = 0
    for i in range(nmatches + 1, nmatches + nmisses)
        miss = misses[i - nmatches]
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(miss))
            prec_row  += UInt8(1)
            insert!(precID_to_row, getPrecID(miss), (prec_row, false))
        end
       # if getPeakInd(miss) != last_peak_ind
            col += 1
            last_peak_ind = getPeakInd(miss)
        #end
        row = precID_to_row[getPrecID(miss)][1]
        H_COLS[i] = col
        H_ROWS[i] = row
        H_VALS[i] = getPredictedIntenisty(miss)
    end
    #if i >= (nmatches - 1)
    #return X, sparse(vcat(H_COLS, U_COLS), vcat(H_ROWS, U_ROWS), vcat(H_VALS, U_VALS)), sparse(vcat(H_ROWS, U_ROWS), vcat(H_COLS, U_COLS), vcat(H_VALS, U_VALS)), precID_to_row, H_ncol
    #end
    #return X, sparse(H_COLS, H_ROWS, H_VALS), sparse(H_ROWS, H_COLS, H_VALS), precID_to_row, H_ncol
    return X, sparse(@view(H_COLS[1:(nmatches + nmisses)]), @view(H_ROWS[1:(nmatches + nmisses)]), @view(H_VALS[1:(nmatches + nmisses)])), precID_to_row, H_ncol
    
end
#=
function buildDesignMatrix(matches::Vector{m},  misses::Vector{m}, nmatches::Int64, nmisses::Int64, H_COLS::Vector{Int64}, H_ROWS::Vector{Int64}, H_VALS::Vector{Float32}, H_MASK::Vector{Float32}; y_min_ind::Int64 = 4, b_min_ind::Int64 = 3, block_size = 10000) where {m<:Match}
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
    if (nmatches + nmisses) >= length(H_COLS) - 1
        H_COLS = zeros(Int64, length(H_COLS) + block_size)
        H_ROWS = zeros(Int64, length(H_ROWS) + block_size)
        H_VALS = zeros(Float32, length(H_VALS) + block_size)
        H_MASK = zeros(Bool, length(H_VALS) + block_size)
    end
    #Spectrum/empirical intensities for each peak. Zero by default (for unmatched/missed fragments)
    X = zeros(T, M)
    #Maps a precursor id to a row of H. 
    precID_to_row = UnorderedDictionary{UInt32, Tuple{UInt32, UInt8}}()

    #Current highest row encountered
    prec_row = zero(UInt32)
    col = 0
    #Number of unique peaks encountered. 
    last_peak_ind = 0
    for i in range(1, nmatches)#matches)
        match = matches[i]
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(match))
            prec_row += one(UInt32)
            insert!(precID_to_row, getPrecID(match), (prec_row, zero(UInt8)))
        end
        if getRank(match) < 3
            precID_to_row[getPrecID(match)] = (precID_to_row[getPrecID(match)][1], precID_to_row[getPrecID(match)][2] + one(UInt8))
        end
        #If this peak has not been encountered yet, then start filling a new column
        if getPeakInd(match) != last_peak_ind
            col += 1
            last_peak_ind = getPeakInd(match)
            X[col] = getIntensity(match)
        end
        row = precID_to_row[getPrecID(match)][1]
        H_COLS[i] = col
        H_ROWS[i] = row
        H_VALS[i] = getPredictedIntenisty(match)
        H_MASK[i] = getIonType(match) == 'y' ? Float32(getFragInd(match) >= y_min_ind) : Float32(getFragInd(match) >= b_min_ind)
        #println("i ", i)
    end
    H_ncol = col
    #col = 0
    last_peak_ind = 0
    for i in range(nmatches + 1, nmatches + nmisses)
        miss = misses[i - nmatches]
        #If a match for this precursor hasn't been encountered yet, then assign it an unused row of H
        if !haskey(precID_to_row,  getPrecID(miss))
            prec_row  += UInt8(1)
            insert!(precID_to_row, getPrecID(miss), (prec_row, false))
        end
       # if getPeakInd(miss) != last_peak_ind
            col += 1
            last_peak_ind = getPeakInd(miss)
        #end
        row = precID_to_row[getPrecID(miss)][1]
        H_COLS[i] = col
        H_ROWS[i] = row
        H_VALS[i] = getPredictedIntenisty(miss)
        H_MASK[i] = getIonType(miss) == 'y' ? Float32(getFragInd(miss) >= y_min_ind) : Float32(getFragInd(miss) >= b_min_ind)
    end
    #if i >= (nmatches - 1)
    #return X, sparse(vcat(H_COLS, U_COLS), vcat(H_ROWS, U_ROWS), vcat(H_VALS, U_VALS)), sparse(vcat(H_ROWS, U_ROWS), vcat(H_COLS, U_COLS), vcat(H_VALS, U_VALS)), precID_to_row, H_ncol
    #end
    #return X, sparse(H_COLS, H_ROWS, H_VALS), sparse(H_ROWS, H_COLS, H_VALS), precID_to_row, H_ncol
    return X, sparse(@view(H_COLS[1:(nmatches + nmisses)]), @view(H_ROWS[1:(nmatches + nmisses)]), @view(H_VALS[1:(nmatches + nmisses)])), sparse(@view(H_COLS[1:(nmatches + nmisses)]), @view(H_ROWS[1:(nmatches + nmisses)]), @view(H_MASK[1:(nmatches + nmisses)])) , precID_to_row, H_ncol
    
end
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

            #=
            last_open = 1
            filtered_nmisses = 0
            for i in range(1, nmatches)
                m = ionMatches[i]
                if haskey(IDtoROW_SUB, getPrecID(m))
                        filtered_nmatches += 1
                        ionMatches[last_open] = ionMatches[i]
                        last_open = filtered_nmatches + 1
                end
            end
            filtered_nmisses = 0
            for i in range(1, nmisses)
                m = ionMisses[i]
                if haskey(IDtoROW_SUB, getPrecID(m))
                        filtered_nmisses += 1
                        ionMisses[last_open] = ionMisses[i]
                        last_open = filtered_nmisses + 1
                end
            end

            prep_time += @elapsed X, Hs, IDtoROW, last_matched_col = buildDesignMatrix(ionMatches, ionMisses, filtered_nmatches, filtered_nmisses, H_COLS, H_ROWS, H_VALS)
            =#
            #=
            weights = zeros(eltype(Hs), Hs.n)#sparseNMF(Hs, X, λ, γ, regularize, max_iter=max_iter, tol=nmf_tol)[:]
            ids = zeros(UInt32, length(weights))
            for (id, row) in pairs(IDtoROW)
                ids[first(row)] = id
            end
            for (id, row) in pairs(IDtoROW)
                weights[first(row)] = precursor_weights[id]
            end
            #prep_time += @elapsed begin 
            col_mask = scores[:matched_ratio] .> 0.5
            Hs = Hs[:,col_mask]
            weights = weights[col_mask]
            ids = ids[col_mask]
            for (id, row) in pairs(IDtoROW)
                if !col_mask[first(row)]
                    delete!(IDtoROW, id)
                end
            end
            for (row, id) in enumerate(ids)
                IDtoROW[id] = (row, last(IDtoROW[id]))
            end
            =#