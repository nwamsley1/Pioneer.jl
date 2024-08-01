function getWittakerHendersonDesignMat(n::Int64, λ::AbstractFloat)
    n < 3 ? throw(ArgumentError("n must be greater than 2")) : nothing
    #λ <= 0 ? throw(ArgumentError("λ must be greater than zero ")) : nothing
    D = secondOrderFiniteDiffMatrix(n, dtype = eltype(λ))
    D2 = transpose(D)*D
    return λ.*D2 + I
end

function secondOrderFiniteDiffMatrix(n::Int64; dtype::DataType = Float64)
    n < 3 ? throw(ArgumentError("n must be greater than 2")) : nothing
    D = zeros(dtype, (n - 2, n))
    start = 1
    for i in range(1, n - 2)
        D[i, start] = 1
        D[i, start + 1] = -2
        D[i, start + 2] = 1
        start += 1
    end
    return D
end