function factorSpectrum(Wnew::Vector{Float32}, Wold::Vector{Float32}, HHt_diag::Vector{Float32}, WxHHt_VHt::Matrix{Float32}, HHt::SparseMatrixCSC{Float32, Int64}, λ::Float32 = zero(Float32))
    a = Inf
    i = 1
    #println("TEST")
    while (abs(a) > 100.0) & (i < 1000)
        a = 0
        for r in 1:length(Wnew)
            Wnew[r] = max(zero(Float32), Wold[r] - (WxHHt_VHt[r] + λ)/HHt_diag[r])
            #WxHHt_VHt = (Wnew*HHt - VHt)
           @turbo for i in HHt.colptr[r]:(HHt.colptr[r+1] -1)
                WxHHt_VHt[HHt.rowval[i]] += ((HHt.nzval[i]*Wnew[r]) - (HHt.nzval[i]*Wold[r]))
            end
            a += abs(Wnew[r] - Wold[r])
            Wold[r] = Wnew[r]
        end
        i += 1
    end
end

function sparseNMF(H::SparseMatrixCSC{Float32, Int64}, X::Vector{T}; λ::T = zero(T)) where {T<:AbstractFloat}
    W = [abs(x) for x in randn(Float32, (1, H.m))] 
    Wnew, Wold = copy(W[:]), copy(W[:])
    #Initialize
    HHt = H*H'
    HHt_diag = collect(diag(HHt))
    VHt = X'*Hs'
    WxHHt_VHt = collect(W*HHt - VHt)
    factorSpectrum(Wnew, Wold, HHt_diag, WxHHt_VHt, HHt, λ);
    return Wnew
end

function factorSpectrum(Wnew::Vector{Float32}, Wold::Vector{Float32}, HHt_diag::SparseVector{Float32, Int64}, WxHHt_VHt::Vector{Float32}, HHt::SparseMatrixCSC{Float32, Int64})
    a = Inf
    i = 1
    N = length(HHt.nzval)

    while (abs(a) > 100.0) & (i < 32)
        for r in 1:length(Wnew)
            Wnew[r] = max(zero(Float32), Wold[r] - (WxHHt_VHt[r] + λ)/HHt_diag[r])
            #WxHHt_VHt = (Wnew*HHt - VHt) #This is inefficient and is 99% of the computation. I don't need to rework the entire thing. 
            update!(WxHHt_VHt, HHt, Wold[r], Wnew[r], r)
            Wold[r] = Wnew[r]
        end
        #a = sum(Wnew .- W)
        #W = copy(Wnew)
        i += 1
    end
end

function update!(WxHHt::Vector{Float32}, HHt::SparseMatrixCSC{Float32, Int64}, Wold::Float32, Wnew::Float32, col::Int)
    for i in HHt.colptr[col]:min(HHt.colptr[col+1], length(HHt.nzval))
        WxHHt[HHt.rowval[i]] += ((HHt.nzval[i]*Wnew) - (HHt.nzval[i]*Wold))
        #Wnew[H.rowval[i]] = (HHt.nzval[i]*Wnew[col] - VHt[col]) - (HHt.nzval[i]*Wold[col] - VHt[col])
    end
end
Hs = sparse(H.rowval[sortperm(H.rowval)],H.colptr[sortperm(H.rowval)], H.nzval[sortperm(H.rowval)])
#=H_peak_pep = sparse(Hs')
norm_facs = zeros(Float32, 200)
for r in 1:(length(H_peak_pep.colptr)-1)
    n = 1
    for i in H_peak_pep.colptr[r]:(H_peak_pep.colptr[r+1] -1)
        norm_facs[r] += X[H_peak_pep.rowval[i]]
        n += 1
    end
    norm_facs[r]/n
end
for i in 1:length(Hs.rowval)
    Hs.nzval[i] = Hs.nzval[i]./norm_facs[Hs.rowval[i]]
end=#
test = [X, H, W, Hs, H]
@save "/Users/n.t.wamsley/Desktop/test.jld2" test
W = [abs(x) for x in 100*randn(Float32, (1, 200))] #1x200


Wnew = copy(W[:])
Wold = copy(W[:])
#Hs = Hs' #200x453
HHt = Hs*Hs'
HHt_diag = collect(diag(HHt))
VHt = X*Hs'
HHt = Hs*Hs'
WxHHt_VHt = collect(W*HHt - VHt)
@time factorSpectrum(Wnew, Wold, HHt_diag, WxHHt_VHt, HHt, Float32(1e-3));
#@profview factorSpectrum(Wnew, Wold, HHt_diag, WxHHt_VHt, HHt)
#sum(weights .- Wnew[:])

MHs =  Matrix(Hs)
weights = (NMF.solve!(NMF.CoordinateDescent{Float32}(maxiter=1000, verbose = false, 
    tol = 1, #Need a reasonable way to choose lambda?
    update_H = false, #Important to keep H constant. 
    regularization=:both,
    l₁ratio = 1.0,
    α = 1e5
    ), X, W, Matrix(Hs)).W[1,:])

weights[1:5]'
Wnew[1:5]'


W, H = NMF.randinit(X, 3)

X = Float32[1 2 3; 2 4 6; 4 8 12]
W, H = NMF.randinit(X, 3)

struct SparseMatrix{Ti<:Integer, Tv<:AbstractFloat}
    m::Int                  # Number of rows
    n::Int                  # Number of columns
    colptr::Vector{Ti}      # Column indices of stored values
    rowval::Vector{Ti}      # Row indices of stored values
    nzval::Vector{Tv}       # Stored values, typically nonzeros
end