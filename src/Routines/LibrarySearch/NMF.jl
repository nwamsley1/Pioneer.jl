
function factorSpectrum(Wnew::Matrix{T}, Wold::Matrix{T}, HHt_diag::Vector{T}, WxHHt_VHt::Matrix{T}, HHt::SparseMatrixCSC{T, Int64}, λs::Vector{T}, max_iter::Int, tol::T) where {T<:AbstractFloat}
    #Fast Coordinate Descent Methods with Variable Selection for Non-negative Matrix Factorization
    #Cho-Jui Hsieh and Inderjit S. Dhillon
    a = Inf
    i = 1
    while (abs(a) > tol) & (i < max_iter)
        a = 0
        for r in 1:length(Wnew)
            
            #if iszero(Wnew[r]) #Is this legitimate? Several X speed boost from ignoring weights after they are set to zero. 
            #    continue
            #end
            Wnew[r] = max(zero(T), Wold[r] - (WxHHt_VHt[r] + (λs[r]))/HHt_diag[r])
            for i in HHt.colptr[r]:(HHt.colptr[r+1] -1)
                WxHHt_VHt[HHt.rowval[i]] += HHt.nzval[i]*(Wnew[r] - Wold[r])
            end
            a += abs(Wnew[r] - Wold[r])
            Wold[r] = Wnew[r]
        end
        i += 1
    end
end

function sparseNMF(H::SparseMatrixCSC{T, Int64}, Ht::SparseMatrixCSC{T, Int64}, X::Vector{T}; λ::T = zero(T), γ::T = zero(T)/2, max_iter::Int = 1000, tol::T = 100*one(T)) where {T<:AbstractFloat}

    Wnew = 100*ones(T, (1, H.m))
    Wold = copy(Wnew)
    λs = zeros(T, H.m)
    #initW!(Wnew, Ht, X)
    #Wnew, Wold = copy(W[:]), copy(W[:])

    HHt = H*Ht
    HHt_diag = collect(diag(HHt))
    VHt = X'*Ht
    WxHHt_VHt = collect(Wnew*HHt - VHt)

    ##OLS estimate with non-negative constraint since penalties are zero 
    factorSpectrum(Wnew, Wold, HHt_diag, WxHHt_VHt, HHt, λs, max_iter, tol);

    #Set adaptive weights 
    #setLambdas!(λs, Float32(λ*sqrt(H.m)), γ, Wnew)
    setLambdas!(λs, Float32(λ), γ, Wnew)

    #Addaptive LASSO estimation 
    factorSpectrum(Wnew, Wold, HHt_diag, WxHHt_VHt, HHt, λs, max_iter, tol);
    return Wnew

end

function setLambdas!(λs::Vector{T}, λ::T, γ::T, W::Matrix{T}, adaptive::Bool = true) where {T<:AbstractFloat}
    if adaptive
        for i in 1:length(λs)
            if !isinf(λ/(W[i]))
                λs[i] = λ/(W[i]^γ)
            end
        end
    else
        for i in 1:length(λs)
                λs[i] = λ
        end
    end
end


#mean(X.*sum(Hs_mat.>0, dims = 1)[:])*200*200
#=function factorSpectrum(Wnew::Vector{T}, Wold::Vector{T}, HHt_diag::Vector{T}, WxHHt_VHt::Matrix{T}, HHt::SparseMatrixCSC{T, Int64}, λ::T, max_iter::Int, tol::T) where {T<:AbstractFloat}
    a = Inf
    i = 1
    while (abs(a) > tol) & (i < max_iter)
        a = 0
        for r in 1:length(Wnew)
            Wnew[r] = max(zero(T), Wold[r] - (WxHHt_VHt[r] + λ)/HHt_diag[r])
            for i in HHt.colptr[r]:(HHt.colptr[r+1] -1)
                WxHHt_VHt[HHt.rowval[i]] += HHt.nzval[i]*(Wnew[r] - Wold[r])
            end
            a += abs(Wnew[r] - Wold[r])
            Wold[r] = Wnew[r]
        end
        i += 1
    end
end

non_zeros = Int64[]
test_lambdas = Float64[]
t_λ = 1e7
for i in 1:150 #λ in range(1e4, 3e6, step = 1.5e4)
    t_λ = t_λ*1.1
    testW =  sparseNMF(Hs, Hst, X, λ=Float32(t_λ), max_iter = 1000)
    push!(non_zeros, sum(testW.>100))
    push!(test_lambdas, t_λ)
end
plot(log.(test_lambdas), non_zeros, seriestype=:scatter)

function sparseNMF(H::SparseMatrixCSC{T, Int64}, Ht::SparseMatrixCSC{T, Int64}, X::Vector{T}, λ::T = zero(T), max_iter::Int = 1000, tol::T = 10*one(T)) where {T<:AbstractFloat}
    W = [abs(x) for x in randn(T, (1, H.m))] 
    Wnew, Wold = copy(W[:]), copy(W[:])
    #Initialize
    #HHt = H*H'
    HHt = H*Ht
    HHt_diag = collect(diag(HHt))
    VHt = X'*Ht
    WxHHt_VHt = collect(W*HHt - VHt)
    factorSpectrum(Wnew, Wold, HHt_diag, WxHHt_VHt, HHt, λ, max_iter, tol);
    return Wnew
end

non_zeros = Int64[]
test_lambdas = Float64[]
t_λ = 1e4
for i in 1:100 #λ in range(1e4, 3e6, step = 1.5e4)
    t_λ = t_λ*1.1
    testW = sparseNMF(Hs, X, λ=Float32(t_λ));
    push!(non_zeros, sum(testW.<100))
    push!(test_lambdas, t_λ)
end
plot(log2.(test_lambdas), non_zeros, seriestype=:scatter)
=#
######
#Get seperate regulrization coefficients for each template. Termed "Adaptive LASSO"
######

#=function initW!(W::Matrix{T}, Ht::SparseMatrixCSC{T, Int64}, X::Vector{T}) where {T<:AbstractFloat}
    for col in 1:(length(Ht.colptr)-1)
        n = 0
        for row in Ht.colptr[col]:(Ht.colptr[col+1] - 1)
            W[col] += X[Ht.rowval[row]]/Ht.nzval[row]
            n += 1
        end
        W[col] = W[col]/n
    end
end

sparseNMF(Hs, Hst, X, λ=one(Float32), max_iter = 1000)



norm_facs = zeros(Float32, 200)
for r in 1:(length(H_peak_pep.colptr)-1)
    n = 1
    for i in H_peak_pep.colptr[r]:(H_peak_pep.colptr[r+1] -1)
        norm_facs[r] += X[H_peak_pep.rowval[i]]
        n += 1
    end
    norm_facs[r]/n
end

Hs_mat_half = Hs_mat[1:100,:]
@benchmark (NMF.solve!(NMF.CoordinateDescent{Float32}(maxiter=1000, verbose = false, 
tol = 100, #Need a reasonable way to choose lambda?
update_H = false, #Important to keep H constant. 
regularization=:both,
l₁ratio = 1.0,
α = 1e5,
shuffle=true
), X_, W, Hs_mat).W)

W_half = W[:,1:100]
@benchmark (NMF.solve!(NMF.CoordinateDescent{Float32}(maxiter=1000, verbose = false, 
tol = 100, #Need a reasonable way to choose lambda?
update_H = false, #Important to keep H constant. 
regularization=:both,
l₁ratio = 1.0,
α = 1e5,
shuffle=true
), X_, W_half, Hs_mat_half).W[1,:])
#=@btime sparseNMF(Hs, X, λ=Float32(1e5))

@profview for i in 1:100
    sparseNMF(Hs, X, λ=Float32(1e5))
end
@save "/Users/n.t.wamsley/Projects/PROSIT/H.jld2"  H

Hs_mat = Matrix(HS)
X_ = Matrix(X')
weights = (NMF.solve!(NMF.CoordinateDescent{Float32}(maxiter=1000, verbose = false, 
    tol = 100, #Need a reasonable way to choose lambda?
    update_H = false, #Important to keep H constant. 
    regularization=:both,
    l₁ratio = 1.0,
    α = 1e5
    ), X_, W, Hs_mat).W[1,:])

    @benchmark (NMF.solve!(NMF.GreedyCD{Float32}(maxiter=1000, verbose = false, 
           tol = 100, #Need a reasonable way to choose lambda?
           update_H = false, #Important to keep H constant. 
           lambda_w = 1e5,
           lambda_h = 1e5
           ), X_, W, Hs_mat).W[1,:])

    @benchmark (NMF.solve!(NMF.CoordinateDescent{Float32}(maxiter=1000, verbose = false, 
           tol = 100, #Need a reasonable way to choose lambda?
           update_H = false, #Important to keep H constant. 
           regularization=:both,
           l₁ratio = 1.0,
           α = 1e5,
           shuffle=true
           ), X_, W, Hs_mat).W[1,:])
           =#


=#