function getDistanceMetrics(H::SparseMatrixCSC{T, Int64}, Ht::SparseMatrixCSC{T, Int64}, X::Vector{T}, unmatched_col::Int) where {T<:AbstractFloat}
    function rowNormsAndSums(A::SparseMatrixCSC{T, Int64}, X::Vector{T}, unmatched_col::Int) where {T<:AbstractFloat}

        rownorms_A = zeros(T, (A.m,))
        rowsums_A = zeros(T, (A.m,))
        rowsums_C = zeros(T, (A.m,))
        rowsums_sqrt_A = zeros(T, (A.m,))
        rownorms_X = zeros(T, (A.m,))
        rowsums_sqrt_X = zeros(T, (A.m,))
        row_counts = zeros(Int64, (A.m,))
        row_dot = zeros(T, (A.m,))

        for col in 1:(unmatched_col)

            for i in range(A.colptr[col], (A.colptr[col+1]-1))
                rownorms_A[A.rowval[i]] += A.nzval[i]^2
                row_counts[A.rowval[i]] += 1
                rownorms_X[A.rowval[i]] += X[col]^2
                rowsums_sqrt_A[A.rowval[i]] += sqrt(A.nzval[i])
                rowsums_A[A.rowval[i]] += A.nzval[i]
                rowsums_sqrt_X[A.rowval[i]] += sqrt(X[col])
                row_dot[A.rowval[i]] +=  A.nzval[i]*X[col]
            end
        end

        rownorms_ALL = copy(rownorms_A)
        for col in unmatched_col:A.n 
            for i in A.colptr[col]:(A.colptr[col + 1] -1)
                rownorms_ALL[A.rowval[i]] +=  A.nzval[i]^2
                rowsums_C[A.rowval[i]] +=  A.nzval[i]
            end
        end

        return sqrt.(rownorms_A), rowsums_sqrt_A, sqrt.(rownorms_X), rowsums_sqrt_X, sqrt.(rownorms_ALL), row_dot, row_counts, rowsums_A, rowsums_C
    end

    rownorms_A, rowsums_sqrt_A, rownorms_X, rowsums_sqrt_X, rownorms_ALL, row_dot, row_counts, rowsums_MATCHED, rowsums_UNMATCHED = rowNormsAndSums(H, X, unmatched_col)

    function scribeScore(a::T, a_sum::T, b::T, b_sum::T) where {T<:AbstractFloat}

        return ((a/a_sum) - (b/b_sum))^2
     end

    function cityBlockDist(a::T, a_norm::T, b::T, b_norm::T) where {T<:AbstractFloat}
        abs(a/a_norm - b/b_norm)
    end

    N = H.m 
    scribe_squared_errors = zeros(T, (N,)) 
    city_block_dist = zeros(T,(N,))
    matched_ratio = zeros(T, (N,))
    spectral_contrast_all = zeros(T, (N,))
    for col in 1:(unmatched_col - 1)
        for i in range(H.colptr[col], H.colptr[col+1]-1)
                    scribe_squared_errors[H.rowval[i]] += scribeScore(sqrt(H.nzval[i]), 
                                                          rowsums_sqrt_A[H.rowval[i]], 
                                                          sqrt(X[col]),
                                                          rowsums_sqrt_X[H.rowval[i]]
                                                )
                    city_block_dist[H.rowval[i]] += cityBlockDist(H.nzval[i], 
                                            rownorms_A[H.rowval[i]], 
                                            X[col],
                                            rownorms_X[H.rowval[i]]
                                        )
        end
    end

    @turbo for (i, count) in enumerate(row_counts)
        scribe_squared_errors[i] = -1*log((scribe_squared_errors[i]^2)/count)
        city_block_dist[i] = log(city_block_dist[i]/count)
    end

    @turbo for (i, dot) in enumerate(row_dot)
        spectral_contrast_all[i] = dot/(rownorms_ALL[i]*rownorms_X[i])
    end

    entropy_sim = getEntropy(X, Ht)
    #=for i in range(1, H.m)
        A = collect(X[H[i,:].!=0.0])# A = X
        A, SA = weightedEntropy(allowmissing(A))
        B = collect(H[i,H[i,:].!=0.0])
        B, SB = weightedEntropy(allowmissing(B))
        A += B
        AB, SAB = weightedEntropy(allowmissing(A))
        entropy_sim[i] =  Float32(1 - (2*SAB - SA - SB)/(log(4)))
    end=#

    matched_ratio = rowsums_MATCHED

    return scribe_squared_errors, city_block_dist, matched_ratio, spectral_contrast_all, entropy_sim
end


#=function getDistanceMetrics(H::SparseMatrixCSC{T, Int64}, Ht::SparseMatrixCSC{T, Int64}, X::Vector{T}, unmatched_col::Int) where {T<:AbstractFloat}
    function rowNormsAndSums(A::SparseMatrixCSC{T, Int64}, X::Vector{T}, unmatched_col::Int) where {T<:AbstractFloat}

        rownorms_A = zeros(T, (A.m,))
        rowsums_A = zeros(T, (A.m,))
        rowsums_C = zeros(T, (A.m,))
        rowsums_sqrt_A = zeros(T, (A.m,))
        rownorms_X = zeros(T, (A.m,))
        rowsums_sqrt_X = zeros(T, (A.m,))
        row_counts = zeros(Int64, (A.m,))
        row_dot = zeros(T, (A.m,))

        for col in 1:(unmatched_col)

            for i in range(A.colptr[col], (A.colptr[col+1]-1))
                rownorms_A[A.rowval[i]] += A.nzval[i]^2
                row_counts[A.rowval[i]] += 1
                rownorms_X[A.rowval[i]] += X[col]^2
                rowsums_sqrt_A[A.rowval[i]] += sqrt(A.nzval[i])
                rowsums_A[A.rowval[i]] += A.nzval[i]
                rowsums_sqrt_X[A.rowval[i]] += sqrt(X[col])
                row_dot[A.rowval[i]] +=  A.nzval[i]*X[col]
            end
        end

        rownorms_ALL = copy(rownorms_A)
        for col in unmatched_col:A.n 
            for i in A.colptr[col]:(A.colptr[col + 1] -1)
                rownorms_ALL[A.rowval[i]] +=  A.nzval[i]^2
                rowsums_C[A.rowval[i]] +=  A.nzval[i]
            end
        end

        return sqrt.(rownorms_A), rowsums_sqrt_A, sqrt.(rownorms_X), rowsums_sqrt_X, sqrt.(rownorms_ALL), row_dot, row_counts, rowsums_A, rowsums_C
    end

    rownorms_A, rowsums_sqrt_A, rownorms_X, rowsums_sqrt_X, rownorms_ALL, row_dot, row_counts, rowsums_MATCHED, rowsums_UNMATCHED = rowNormsAndSums(H, X, unmatched_col)

    function scribeScore(a::T, a_sum::T, b::T, b_sum::T) where {T<:AbstractFloat}

        return ((a/a_sum) - (b/b_sum))^2
     end

    function cityBlockDist(a::T, a_norm::T, b::T, b_norm::T) where {T<:AbstractFloat}
        abs(a/a_norm - b/b_norm)
    end

    N = H.m 
    scribe_squared_errors = zeros(T, (N,)) 
    city_block_dist = zeros(T,(N,))
    matched_ratio = zeros(T, (N,))
    spectral_contrast_all = zeros(T, (N,))
    for col in 1:(unmatched_col - 1)
        for i in range(H.colptr[col], H.colptr[col+1]-1)
                    scribe_squared_errors[H.rowval[i]] += scribeScore(sqrt(H.nzval[i]), 
                                                          rowsums_sqrt_A[H.rowval[i]], 
                                                          sqrt(X[col]),
                                                          rowsums_sqrt_X[H.rowval[i]]
                                                )
                    city_block_dist[H.rowval[i]] += cityBlockDist(H.nzval[i], 
                                            rownorms_A[H.rowval[i]], 
                                            X[col],
                                            rownorms_X[H.rowval[i]]
                                        )
        end
    end

    @turbo for (i, count) in enumerate(row_counts)
        scribe_squared_errors[i] = -1*log((scribe_squared_errors[i]^2)/count)
        city_block_dist[i] = log(city_block_dist[i]/count)
    end

    @turbo for (i, dot) in enumerate(row_dot)
        spectral_contrast_all[i] = dot/(rownorms_ALL[i]*rownorms_X[i])
    end

    entropy_sim = getEntropy(X, Ht)
    #=for i in range(1, H.m)
        A = collect(X[H[i,:].!=0.0])# A = X
        A, SA = weightedEntropy(allowmissing(A))
        B = collect(H[i,H[i,:].!=0.0])
        B, SB = weightedEntropy(allowmissing(B))
        A += B
        AB, SAB = weightedEntropy(allowmissing(A))
        entropy_sim[i] =  Float32(1 - (2*SAB - SA - SB)/(log(4)))
    end=#

    matched_ratio = rowsums_MATCHED

    return scribe_squared_errors, city_block_dist, matched_ratio, spectral_contrast_all, entropy_sim
end=#


#=function test(H::SparseMatrix{Int64, T}, X::Vector{T}, UNMATCHED::SparseMatrix{Int64, T}) where {T<:AbstractFloat}
    for i in 1:1000
        getDistanceMetrics(H, X, UNMATCHED)
    end
end=#

#=
function getDistanceMetrics(H::SparseMatrixCSC{T, Int64}, X::Vector{T}, unmatched_col::Int) where {T<:AbstractFloat}
    #println(X)
    #println(UNMATCHED)
    function rowNormsAndSums(A:::SparseMatrixCSC{T, Int64}, X::Vector{T}, unmatched_col::Int)
        rownorms_A = zeros(T, A.m)
        rowsums_A = zeros(T, A.m)
        rowsums_C = zeros(T, A.m)
        rowsums_sqrt_A = zeros(T, A.m)
        rownorms_X = zeros(T, A.m)
        rowsums_sqrt_X = zeros(T, A.m)
        row_counts = zeros(Int64, A.m)
        row_dot = zeros(T, A.m)
        
        for (i, nzval) in enumerate(A.nzval)
            rownorms_A[A.rowval[i]] += nzval^2
            row_counts[A.rowval[i]] += 1
            rownorms_X[A.rowval[i]] += X[A.colptr[i]]^2
            rowsums_sqrt_A[A.rowval[i]] += sqrt(nzval)
            rowsums_A[A.rowval[i]] += nzval
            rowsums_sqrt_X[A.rowval[i]] += sqrt(X[A.colptr[i]])
            row_dot[A.rowval[i]] +=  nzval*X[A.colptr[i]]
        end

        rownorms_ALL = copy(rownorms_A)
        for (i, nzval) in enumerate(C.nzval)
            rownorms_ALL[C.rowval[i]] += nzval^2
            rowsums_C[C.rowval[i]] += nzval
        end

        return sqrt.(rownorms_A), rowsums_sqrt_A, sqrt.(rownorms_X), rowsums_sqrt_X, sqrt.(rownorms_ALL), row_dot, row_counts, rowsums_A, rowsums_C
    end

    rownorms_A, rowsums_sqrt_A, rownorms_X, rowsums_sqrt_X, rownorms_ALL, row_dot, row_counts, rowsums_MATCHED, rowsums_UNMATCHED = rowNormsAndSums(H, X, UNMATCHED)
    #println(rowsums_MATCHED[2])
    #println(rowsums_UNMATCHED[2])
    function scribeScore(a::T, a_sum::T, b::T, b_sum::T) where {T<:AbstractFloat}
        #-1*log(mean(((a/a_sum) .- (b/b_sum)).^2))
        ((a/a_sum) - (b/b_sum))^2
     end

    function cityBlockDist(a::T, a_norm::T, b::T, b_norm::T) where {T<:AbstractFloat}
        #-1*log(mean(((a/sum(a)) .- (b/sum(b))).^2))
        #log(sum(abs.(a./norm(a) .- b./norm(b))/length(a)))
        abs(a/a_norm - b/b_norm)
    end

    N = H.m 
    scribe_squared_errors = zeros(T, N) #Vector{T}(undef, N)
    city_block_dist = zeros(T, N)
    matched_ratio = zeros(T, N)
    spectral_contrast_matched = zeros(T, N)
    spectral_contrast_all = zeros(T, N)

    @turbo for (i, nzval) in enumerate(H.nzval)
        scribe_squared_errors[H.rowval[i]] += scribeScore(sqrt(nzval), 
                                                          rowsums_sqrt_A[H.rowval[i]], 
                                                          sqrt(X[H.colptr[i]]),
                                                          rowsums_sqrt_X[H.rowval[i]]
                                                )
        city_block_dist[H.rowval[i]] += cityBlockDist(nzval, 
                                            rownorms_A[H.rowval[i]], 
                                            X[H.colptr[i]],
                                            rownorms_X[H.rowval[i]]
                                        )
    end

    @turbo for (i, count) in enumerate(row_counts)
        scribe_squared_errors[i] = -1*log((scribe_squared_errors[i]^2)/count)
        city_block_dist[i] = log(city_block_dist[i]/count)
    end

    @turbo for (i, dot) in enumerate(row_dot)
        spectral_contrast_matched[i] = dot/(rownorms_A[i]*rownorms_X[i])
        spectral_contrast_all[i] = dot/(rownorms_ALL[i]*rownorms_X[i])
        rowsums_MATCHED[i] = rowsums_MATCHED[i]/rowsums_UNMATCHED[i]
    end

    matched_ratio = rowsums_MATCHED
    return scribe_squared_errors, city_block_dist, matched_ratio, spectral_contrast_matched, spectral_contrast_all
end
=#