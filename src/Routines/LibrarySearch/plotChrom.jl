include("src/Routines/LibrarySearch/searchRAW.jl")
test_fmatches = integrateMS2(MS_TABLE, 
    frag_list, 
    rt_index,
    UInt32(ms_file_idx), 
    frag_err_dist_dict[ms_file_idx],
    integrate_ms2_params, 
    scan_range = (49648, 52770)
#scan_range = (101357, 102357)
);

chromatogram = Dictionary{String, Vector{Tuple{Float32, Float32}}}()

for match in test_fmatches
    ion_type = match.ion_type#'y' ? 'y' : 'b'
    ion_name = ion_type*string(match.frag_index)*"+"*string(match.frag_charge)
    if !haskey(chromatogram, ion_name)
        println("ion_name ", ion_name, "rank ", match.predicted_rank)
        insert!(chromatogram, ion_name, Vector{Tuple{Float32, Float32}}())
    end
    push!(chromatogram[ion_name], (MS_TABLE[:retentionTime][match.scan_idx], match.intensity*2))
end

p = plot()
integratePrecursor(ms2_chroms, UInt32(best_psms_passing[N,:precursor_idx]),(best_psms_passing[N,:scan_idx]), isplot = true)
#p = plot()
for (color, key) in enumerate(sort(keys(chromatogram)))
    #if key == "b3+1"
    #    continue
    #end
    plot!(chromatogram[key], seriestype=:scatter, color = color, show = true, xlim = (34.0, 34.8), label = key)
    plot!(chromatogram[key], show = true, color = color, xlim = (34.0, 34.8), label = nothing)
end

PSMs[1][PSMs[1][:,:precursor_idx] .== best_psms_passing[N,:precursor_idx],[:sequence,:precursor_idx,:weight,:total_ions,:best_rank,:entropy_sim,:matched_ratio,:spectral_contrast,:scribe_score,:RT,:q_value,:prob]]

integratePrecursor(ms2_chroms, UInt32(best_psms_passing[N,:precursor_idx]),(best_psms_passing[N,:scan_idx]), isplot = true)

#integratePrecursor(ms2_chroms, UInt32(best_psms_passing[N,:precursor_idx]), 51807, isplot = true)

#plot!(ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:rt],
#ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:weight])
t = ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)]
t = t[t[:,:rank].>1,:]
plot(t[:,:rt],
    (t[:,:frag_count].^2).*t[:,:weight],
    seriestype=:scatter
)

plot(t[:,:rt],
    t[:,:weight],
    seriestype=:scatter
)


#=plot!(ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:rt],
    ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:weight]*36,
    seriestype=:scatter
)=#

N += 1


integratePrecursor(ms2_chroms, UInt32(best_psms_passing[N,:precursor_idx]),(best_psms_passing[N,:scan_idx]), isplot = true)
ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)]
N += 1

include("src/Routines/LibrarySearch/searchRAW.jl")
X, Hs, IDtoROW, last_matched_col = integrateMS2(MS_TABLE, 
    frag_list, 
    rt_index,
    UInt32(ms_file_idx), 
    frag_err_dist_dict[ms_file_idx],
    integrate_ms2_params, 
    scan_range = (49648, 52770)
#scan_range = (101357, 102357)
);

b = collect(X)
A = collect(Hs)

vX = 1000*ones(Hs.n)
vX_ = 1000*ones(Hs.n)
n = 1
for i in range(1, 1000)
    vX = vX .- (10000/sqrt(i)).*(A'*sign.(A*vX - X))
    vX = max.(vX, 0.0)
    n = sqrt(sum(abs.(vX .- vX_)))
    vX_ = vX
    #println("obj ", norm(A * vX - X, 1));
end


vX = 1000*ones(Hs.n)
vX_ = 1000*ones(Hs.n)
#vX = collect(w[:])
#vX = collect(w[:])
n = 1

δ = 1000000
Y = (A*vX .- X)/δ
function coorDescent!(A::Matrix{<:AbstractFloat}, vX::Vector{<:AbstractFloat}, X::Vector{<:AbstractFloat}; δ::Float64=10000.0, N::Int = 20)
    Y = (A*vX .- X)/δ
    for n in ProgressBar(range(1, N))
        for j in range(1, length(vX))
            #Y = (A*vX .- X)/δ
            L = 0.0
            for i in range(1, length(X))
            L += A[i, j]*Y[j]/((1 + Y[j]^2)^(1/2)) 
            end
            L = L*δ
            #println("L $L")
            #n = sum((δ*A')*(AXb./sqrt.(1 .+ (AXb./δ).^2)))
            vold = vX[j]
            vX[j] = max(vX[j] - L, 0.0)
            #vX[j] = vX[j] - 100L/sqrt(n)

            Y .+= (A[:,j]*vX[j] .- X)/δ - (A[:,j]*vold .- X)/δ 
        end
        #println(sum(abs.(vX_ .- vX)))
        #vX_ = vX[:]v
    end
    return vX
end
X[(A[:,5].!=0.0).&(X.!=0.0)]./A[(A[:,5].!=0.0).&(X.!=0.0),5]

N = 56
t = X[(A[:,N ].!=0.0).&(X.!=0.0)]./A[(A[:,N].!=0.0).&(X.!=0.0),N]
std(t)/mean(t)



δ = 1000000
x = LinRange(-1e6, 1e6, 10000)
y = (δ^2)*(sqrt.(1 .+ (x/δ).^2) .- 1)#(1/δ)*log.(exp.(δ*x) + exp.(-δ*x))
plot(x, y)


δ = 0.01
y = exp(δ*x)
y = exp(δ*x)


y = (1/δ)*log.(exp.(δ*x) + exp.(-δ*x))


#=function Lx′(AXb::Vector{<:AbstractFloat}, Aj::Vector{<:AbstractFloat}, δ::AbstractFloat)
    Aj'*((AXb)./(sqrt.(1 .+ (AXb./δ).^2)))
end 

function Lx′′(A::Matrix{<:AbstractFloat}, X::Vector{<:AbstractFloat}, b::Vector{<:AbstractFloat}, Aj::Vector{<:AbstractFloat}, δ::AbstractFloat)
    AXb = A*X .- b
    AXb2δ = (1 .+ (AXb./δ).^2)
    #println( (Aj'*b./δ)*(Aj'*((AXb./(AXb2δ.^(3/2))))))
    return Aj'*(Aj./sqrt.(AXb2δ)) .- ((Aj'*(A*X))./δ)*(Aj'*(AXb./(AXb2δ.^(3/2)))) .+ (
        (Aj'*b./δ)*(Aj'*((AXb./(AXb2δ.^(3/2)))))
    )
end=#

@load "/Users/n.t.wamsley/Desktop/X_test.jld2" X
@load "/Users/n.t.wamsley/Desktop/Hs_test.jld2" Hs
@load "/Users/n.t.wamsley/Desktop/IDtoROW.jld2" IDtoROW


δ = 100.0
weights = sparseNMF(Hs, X, zero(Float32), zero(Float32))
b = X[:]
X₁ = weight[:]#1000*ones(eltype(Hs), Hs.n)
X₀ = weights[:]#1000*ones(eltype(Hs), Hs.n)
r = Hs*X₁ .- b

Hinv = inv(Matrix(Hs'*Hs))
function solveHuber!(Hs::SparseMatrixCSC{Float32, Int64}, r::Vector{T}, X₀::Vector{T}, X₁::Vector{T}, δ::T; max_iter::Int = 1000, tol::AbstractFloat = 100.0) where {T<:AbstractFloat}
    L1 = zero(T)
    L2 = zero(T)
    #newton = zeros(T, Hs.n)
    ΔX = Inf
    i = 0
    while (ΔX > tol) & (i < max_iter)
        ΔX = 0.0
        for col in range(1, Hs.n)
            #println("L1 $L1")
            #println("L2 $L2")
            L1 = zero(T)
            L2 = zero(T)
            for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
                R = 1 + (r[Hs.rowval[i]]/δ)^2
                L1 += Hs.nzval[i]*r[Hs.rowval[i]]*(R)^(-1/2)
                L2 += (Hs.nzval[i]^2)*(R^(-1/2) - ((r[Hs.rowval[i]]^2)/δ)*(R^(-3/2)))
                #WOLFRAM
                #R = δ^2 + r[Hs.rowval[i]]^2
                #L1 += (Hs.nzval[i]*r[Hs.rowval[i]])/(δ*sqrt(R))
                #L2 += ((δ^2)*(Hs.nzval[i])^2)/(δ*(R)^(3/2))

                #SquaredError
                #R = δ^2 + r[Hs.rowval[i]]^2
                #L1 += 2*Hs.nzval[i]*r[Hs.rowval[i]]
                #L2 += 2*(Hs.nzval[i])^2
                #newtwon[col] +=( (Hs.nzval[i]*r[Hs.rowval[i]])/(δ*sqrt(R)) )/( 
                #    ((δ^2)*(Hs.nzval[i])^2)/(δ*(R)^(3/2)) )
            end
            #X₁[col] = max(X[col] - (L1/L2), 0.0)
            X₁[col] = max(X[col] - (L1/L2), 0.0)
            #X₁[col] = max(X[col] - (L1/L2))

            for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
                r[Hs.rowval[i]] += Hs.nzval[i]*(X₁[col] - X₀[col])
            end

            ΔX += abs(X₁[col] - X₀[col])
            X₀[col] = X₁[col]
        end
        i += 1
        println(log(ΔX))
    end
    println(i)
end
#solveHuber!(Hs, r, X₀, X₁, Float32(δ), max_iter = 10)
solveHuber!(Hs, r, X₀, X₁, Float32(δ))

dot(X₁, weights[:])/(norm(X₁)*norm(weights[:]))
plot(log10.(X₁), log10.(weights[:]), seriestype=:scatter)
sum(X₁'.==0.0)


plot!(log10.(X), log10.(weights[:]), seriestype=:scatter)

#=
@load "/Users/n.t.wamsley/Desktop/X_test.jld2" X
@load "/Users/n.t.wamsley/Desktop/Hs_test.jld2" Hs
@load "/Users/n.t.wamsley/Desktop/IDtoROW.jld2" IDtoROW
=#

function updateResiduals(r::Vector{T}, A::SparseMatrixCSC{T, Int64}, k::Int64, x₀::Vector{T}, x₁::Vector{T}) where {T<:AbstractFloat}
    for 

    end
end

δ = 100000.0
x = LinRange(-10000, 10000, 10000)
#δ = 100
y = (δ^2)*(sqrt.(1 .+ (x./δ).^2) .- 1)
plot(x, y)
plot!(x, x.^2)


y′ = (δ)*(x./(sqrt.(1 .+ (x./δ).^2)))
y[7500] .+ x.*(δ)*((500)./(sqrt.(1 .+ ((500)./δ).^2)));
plot!(x, y′)


#Slope at x = 500
(δ)*(500/(sqrt.(1 .+ (500/δ).^2)))
x2 = LinRange(0, 1000, 1000)
y2 = 0 .+ (δ)*(500/(sqrt.(1 .+ (500/δ).^2))).*x2
y2 = (-10000) .+ δ.*x2
plot!(x2, y2)
#plot!(x2, y[7500] .+ (x2.-500)*(δ)*((500)./(sqrt.(1 .+ ((500)./δ).^2))))#