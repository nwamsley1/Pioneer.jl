using Combinatorics, Polynomials, LinearAlgebra

function getLegendrePolynomials(d::Int64; dtype::DataType = Float64)
    l_coeff = Vector{Polynomial{dtype}}(undef, d+1)
    #recurrence start 
    l_coeff[1] = Polynomial([one(dtype)])
    l_coeff[2] = Polynomial([zero(dtype), one(dtype)])
    s = Polynomial([zero(dtype), one(dtype)])
    #Bonnet's recursion formula 
    for i in range(3, d+1)
        k = i - 1
        l_coeff[i] =((2*k - 1)/k)*s*l_coeff[i - 1] - ((k - 1)/k)*l_coeff[i - 2]
    end
    #Legendre polynomials from 0 up to degree d 
    return l_coeff
end

function getLegendreBasedSGKernel(N::Int64, d::Int64; dtype::DataType = Float64)
    P = getLegendrePolynomials(d, dtype = dtype)
    K = Vector{dtype}(undef, 2*N + 1)
    h = Polynomial([zero(dtype), one(dtype)])
    K = ((d)/2)*P[d - 1](zero(dtype))*(P[d]/h)
    #for (i, h) in enumerate(LinRange(-N, N, 2*N + 1))
    #    K[i] = (d/2)*P[d - 1](zero(dtype))*P[d](h/N)/h
    #end
    return K
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

function getWittakerHendersonDesignMat(n::Int64, λ::AbstractFloat)
    n < 3 ? throw(ArgumentError("n must be greater than 2")) : nothing
    λ <= 0 ? throw(ArgumentError("λ must be greater than zero ")) : nothing
    D = secondOrderFiniteDiffMatrix(n, dtype = eltype(λ))
    D2 = transpose(D)*D
    return λ.*D2 + I
end


subchrom = groupby(gchroms[(precursor_idx = best_precursors[N,1],)],:iso_rank)

y = subchrom[1][!,:intensity]
b = zeros(Float32, 200)
A = getWittakerHendersonDesignMat(length(b), 0.5f0)
prob = LinearProblem(A, b)
linsolve = init(prob)

fill!(linsolve.b, zero(eltype(linsolve.b)))
for i in range(1, length(y))
    linsolve.b[i] = y[i]
end
solve!(linsolve)

start, stop = argmax(y) - 20, argmax(y) + 20
plot(y, alpha = 0.5, xlim = (start, stop), seriestype=:scatter)
plot!(linsolve.u[1:length(y)], alpha = 0.5, seriestype=:scatter)
y2 = zeros(length(linsolve.u))
for i in range(2, length(y2) - 2)
    y2[i] = linsolve.u[i + 1] - 2*linsolve.u[i] + linsolve.u[i - 1]
end
plot!(y2, alpha = 0.5, seriestype=:scatter)

#plot!(savitzky_golay(y, 7, 3).y, alpha = 0.5, seriestype=:scatter)
N += 1