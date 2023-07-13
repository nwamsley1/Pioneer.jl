

function factorSpectrum(Wnew::Vector{T}, Wold::Vector{T}, HHt_diag::Vector{T}, WxHHt_VHt::Matrix{T}, HHt::SparseMatrixCSC{T, Int64}, λ::T, max_iter::Int, tol::T) where {T<:AbstractFloat}
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

function sparseNMF(H::SparseMatrixCSC{T, Int64}, X::Vector{T}, λ::T = zero(T), max_iter::Int = 1000, tol::T = 100*one(T)) where {T<:AbstractFloat}
    W = [abs(x) for x in randn(T, (1, H.m))] 
    Wnew, Wold = copy(W[:]), copy(W[:])
    #Initialize
    HHt = H*H'
    HHt_diag = collect(diag(HHt))
    VHt = X'*H'
    WxHHt_VHt = collect(W*HHt - VHt)
    factorSpectrum(Wnew, Wold, HHt_diag, WxHHt_VHt, HHt, λ, max_iter, tol);
    return Wnew
end
Hs_mat_half = Hs_mat[1:100,:]
@benchmark (NMF.solve!(NMF.CoordinateDescent{Float32}(maxiter=1000, verbose = false, 
tol = 100, #Need a reasonable way to choose lambda?
update_H = false, #Important to keep H constant. 
regularization=:both,
l₁ratio = 1.0,
α = 1e5,
shuffle=true
), X_, W, Hs_mat).W[1,:])

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