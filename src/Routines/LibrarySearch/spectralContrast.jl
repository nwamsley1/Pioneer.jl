function getSpectralContrast(H::Matrix{T}, X::Matrix{T}) where {T<:AbstractFloat}
    function spectralContrast(a::Vector{T}, b::Vector{T}) where {T<:AbstractFloat}
        dot(a, b)/(norm(a)*norm(b))
    end

    N = size(H)[1]
    spectral_contrast = Vector{T}(undef, N)

    for row in range(1, N)
        non_zero = H[row,:].!=0
        spectral_contrast[row] = spectralContrast(H[row, non_zero], X[1, non_zero])
    end

    return spectral_contrast 
end