abstract type Params end

struct HuberParams{T<:AbstractFloat} <: Params
    σ::T
    tᵣ::T
    τ::T
    H::T
end

import Base.-
function -(a::HuberParams{T}, b::HuberParams{U}) where {T,U<:AbstractFloat}
    return HuberParams(
                        a.σ - b.σ, 
                        a.tᵣ - b.tᵣ, 
                        a.τ - b.τ,
                        a.H - b.H)
end

import Base.+
function +(a::HuberParams{T}, b::HuberParams{U}) where {T,U<:AbstractFloat}
    return HuberParams(
                        a.σ + b.σ, 
                        a.tᵣ + b.tᵣ, 
                        a.τ + b.τ,
                        a.H + b.H)
end

function +(a::HuberParams{T}, d::R) where {T<:AbstractFloat,R<:Real}
    return HuberParams(
                        a.σ+d, 
                        a.tᵣ+d, 
                        a.τ+d,
                        a.H+d)
end

import Base.*
function *(a::HuberParams{T}, b::HuberParams{T}) where {T<:AbstractFloat}
    return HuberParams(
                        a.σ*b.σ, 
                        a.tᵣ*b.tᵣ, 
                        a.τ*b.τ,
                        a.H*b.H)
end

function *(a::HuberParams{T}, d::R) where {T<:AbstractFloat,R<:Real}
    return HuberParams(
                        a.σ*d, 
                        a.tᵣ*d, 
                        a.τ*d,
                        a.H*d)
end

*(d::R, a::HuberParams{T}) where {T<:AbstractFloat,R<:Real} = *(a, d)

import Base./
function /(a::HuberParams{T}, b::HuberParams{U}) where {T,U<:AbstractFloat}
    return HuberParams(
                        a.σ/b.σ, 
                        a.tᵣ/b.tᵣ, 
                        a.τ/b.τ,
                        a.H/b.H)
end

import Base./
function /(a::HuberParams{T}, d::R) where {T<:AbstractFloat,R<:Real}
    return HuberParams(
                        a.σ/d, 
                        a.tᵣ/d, 
                        a.τ/d,
                        a.H/d)
end

function norm(a::HuberParams{T}) where {T<:AbstractFloat}
    return sqrt(a.σ^2 + a.tᵣ^2 + a.τ^2 + a.H^2)
end

import Base.sqrt
function sqrt(a::HuberParams{T}) where {T<:AbstractFloat}
    return HuberParams(
        sqrt(a.σ),sqrt(a.tᵣ),sqrt(a.τ),sqrt(a.H))#sqrt(a.σ^2 + a.tᵣ^2 + a.τ^2 + a.H^2)
end

mutable struct GD_state{P<:Params, T<:Real, I,J<:Integer}
    params::P
    t::Vector{T}
    y::Vector{T}
    data::Vector{T}
    mask::BitVector
    n::I
    max_index::J
end

function F(state::GD_state{HuberParams{T}, U, I, J}, tᵢ::AbstractFloat) where {T,U<:AbstractFloat, I,J<:Integer}
    #Evaluate EGH function at eath time point tᵢ and store them in pre-allocated array 'f'. 
    #for (i, tᵢ) in enumerate(x)
    d = 2*state.params.σ + state.params.τ*(tᵢ - state.params.tᵣ)
    if real(d) > 0
        return state.params.H*exp((-(tᵢ - state.params.tᵣ)^2)/d)
    else
        return zero(T)
    end
end

function F!(state::GD_state{HuberParams{T}, U, I, J}) where {T,U<:AbstractFloat, I,J<:Integer}
    #Evaluate EGH function at eath time point tᵢ and store them in pre-allocated array 'f'. 
    #for (i, tᵢ) in enumerate(x)
    for i in range(1, state.max_index)
        state.mask[i] ? continue : nothing
        tᵢ = state.t[i]
        d = 2*state.params.σ + state.params.τ*(tᵢ - state.params.tᵣ)
        if real(d) > 0
            state.y[i] = state.params.H*exp((-(tᵢ - state.params.tᵣ)^2)/d)
        else
            state.y[i] = zero(T)
        end
    end
end

function Jacobian(state::GD_state{HuberParams{T}, U, I, J}, δ::V) where {T,U,V<:AbstractFloat, I,J<:Integer}

    #Initialize parameters
    H, τ, σ, tᵣ = state.params.H, state.params.τ, state.params.σ, state.params.tᵣ
    J_σ, J_tᵣ, J_τ, J_H = 0.0, 0.0, 0.0, 0.0

    for i in range(1, state.max_index)
        state.mask[i] ? continue : nothing #If data point is masked, skip

        #Compute repeating terms in the partial derivatives 
        tᵢ = state.t[i]
        DT = tᵢ - tᵣ
        T2 = DT^2
        D = 2*σ + DT*τ
        T2D2 = (T2/D^2)
        EXP = exp(-1.0*T2/D)
        N = H*EXP - state.data[i]
        Denom = sqrt(1 + (N/δ)^2)
        Common = N*EXP/Denom
        if 2*σ + τ*(DT) <= 0.0
            continue
        end

        J_σ += -2.0*H*Common*(-1.0*T2D2)
        J_tᵣ += H*(2*DT/D + τ*T2D2)*Common
        J_τ += -1.0*H*Common*DT*(-1.0*T2D2)
        J_H += Common

    end

    return HuberParams(T(J_σ), T(J_tᵣ), T(J_τ), T(J_H))
end

function reset!(state::GD_state{P,T,I,J}) where {P<:Params, T<:Real, I,J<:Integer} 
    for i in range(1, state.max_index)
        state.t[i], state.y[i], state.mask[i] = zero(T), zero(T), false
    end
    state.max_index = 0
    state.n = 0
    return 
end

function updateParams(params::HuberParams{T}, lower::HuberParams{T}, upper::HuberParams{T}) where {T,U<:AbstractFloat}
     return HuberParams(
         max(min(params.σ, upper.σ), lower.σ),
         max(min(params.tᵣ, upper.tᵣ), lower.tᵣ),
         max(min(params.τ, upper.τ), lower.τ),
         max(min(params.H, upper.H), lower.H)
     )
end

function GD(state::GD_state{P,T,I,J}, lower::P, upper::P; 
            tol::Float64 = 1e-3, 
            max_iter::Int64 = 1000, 
            δ::Float64 = 1000.0,
            α::Float64 = 0.1,
            β1 = 0.9,
            β2 = 0.999,
            ϵ = 1e-8) where {P<:Params, T<:Real, I,J<:Integer} 
    
    #Iteration counter
    state.n = 1

    #Initialize moment estimates 
    m1 = HuberParams(T(0.0), T(0.0),T(0.0),T(0.0))
    m2 = HuberParams(T(0.0), T(0.0),T(0.0),T(0.0))
    β1 = T(β1)
    β2 = T(β2)
    ϵ = T(ϵ)
    α = T(α)

    while state.n <= max_iter
        t = state.n #Time step
        F!(state) #Evaluate function with current parameters
        grad = Jacobian(state, δ) #Compute the gradient 
        m1 = β1*m1 + (1 - β1)*grad #First moment estimate
        m2 = β2*m2 + (1 - β2)*grad*grad #second moment estimate
        αₜ = α*sqrt(1 - β2^t)/(1 - β1^t) #Bias correction 
        diff = αₜ*m1/(sqrt(m2) + ϵ)
        #Get new parameters
        state.params = updateParams(state.params - diff, lower, upper)
        #If the euclidian norm of the update has fallen below a throeshold
        if norm(diff) <= tol
            return 
        end
        state.n += 1
    end    
end

function Integrate(state::GD_state{P,T,I,J}, x::Vector{U}, w::Vector{U}; α::AbstractFloat = 0.01) where {P<:Params, T,U<:AbstractFloat, I,J<:Integer} 

    #Use GuassLegendre Quadrature to integrate f on the integration bounds 
    #using FastGaussQuadrature
    #x, w = gausslegendre(n)

    function getBase(state::GD_state{HuberParams{T}, U, I, J}, α::AbstractFloat) where {T,U<:AbstractFloat, I,J<:Integer}
        τ = state.params.τ
        σ = state.params.σ
        B = (-1/2)*(sqrt(abs(log(α)*((τ^2)*log(α) - 8*σ) + τ*log(α))))
        A = (1/2)*( τ*log(α) - sqrt(abs(log(α)*((τ^2)*log(α) - 8*σ))))
        return state.params.tᵣ - abs(A), state.params.tᵣ + abs(B)
    end

    a, b = getBase(state, α)
    #a = state.params.tᵣ - 1.0
    #b = state.params.tᵣ + 1.0
    #println("a $a; b $b")
    #Quadrature rules are for integration bounds -1 to 1 but can shift
    #to arbitrary bounds a and b. 
    dotp = 0.0
    for i in eachindex(x)
        dotp += F(state, 
                    (x[i]*((b - a)/2) + (a + b)/2) #Evaluate at 
                )*w[i]
    end

    return ((b - a)/2)*dotp
end

function updateParams(params::HuberParams{T}, lower::HuberParams{T}, upper::HuberParams{T}) where {T<:AbstractFloat}
    return HuberParams(
        max(min(params.σ, upper.σ), lower.σ),
        max(min(params.tᵣ, upper.tᵣ), lower.tᵣ),
        max(min(params.τ, upper.τ), lower.τ),
        max(min(params.H, upper.H), lower.H)
    )
end


#σ=p[1], tᵣ=p[2], τ=p[3], H = p[4]
function getP0(α::T, B::T, A::T, tᵣ::T, H::T, lower::HuberParams{T}, upper::HuberParams{T}) where {T<:AbstractFloat}
    return HuberParams(
             max(min((-1/(2*log(α)))*(B*A), upper.σ), lower.σ),
             max(min(tᵣ,upper.tᵣ), lower.tᵣ),
             max(min((-1/log(α))*(B - A),upper.τ), lower.τ),
             max(min(H,upper.H), lower.H)
             )
end

function getP0(α::T, B::T, A::T, tᵣ::T, H::T,lower::Vector{T},upper::Vector{T}) where {T<:AbstractFloat}
    p0 = T[(-1/(2*log(α)))*(B*A),
             tᵣ,
             (-1/log(α))*(B - A),
             H]
    for (i, v) in enumerate(p0)
        if (v < lower[i]) | (v > upper[i]) 
            p0[i] = (lower[i] + upper[i])/2
        end
    end
    return p0
        
end


getP0(p0::NTuple{5, T}) where {T<:AbstractFloat} = getP0(p0[1], p0[2],p0[3],p0[4],p0[5])
getP0(p0::NTuple{5, T},lower::Vector{T},upper::Vector{T}) where {T<:AbstractFloat} = getP0(p0[1], p0[2],p0[3],p0[4],p0[5], lower, upper)

function getFWHM(α::AbstractFloat, τ::T, σ::T) where {T<:AbstractFloat}
    #FWHM given parameters for EGH function 

    #When is absolute value necessary. 
    B = (-1/2)*(sqrt(abs(log(α)*((τ^2)*log(α) - 8*σ) + τ*log(α))))
    A = (1/2)*( τ*log(α) - sqrt(abs(log(α)*((τ^2)*log(α) - 8*σ))))
    return abs(A) + abs(B)
end

function getFWHM(state::GD_state{HuberParams{T}, U, I, J}, α::AbstractFloat) where {T,U<:AbstractFloat, I,J<:Integer}
    #FWHM given parameters for EGH function 
    τ = state.params.τ
    σ = state.params.σ
    #When is absolute value necessary. 
    B = (-1/2)*(sqrt(abs(log(α)*((τ^2)*log(α) - 8*σ) + τ*log(α))))
    A = (1/2)*( τ*log(α) - sqrt(abs(log(α)*((τ^2)*log(α) - 8*σ))))
    return abs(A) + abs(B)
end


function getBestPSM(sdf::SubDataFrame)
    if any(sdf.q_value.<=0.01)
        return sdf[argmax((sdf.q_value.<=0.01).*(sdf.hyperscore)),:]
    else
        return sdf[argmax(sdf.prob),:]
    end
end

#=

function CrossCorrFFT(f::Matrix{Complex{T}}, g::Matrix{Complex{T}}, δt::T) where {T<:AbstractFloat}
    
    if length(f) != length(g)
        throw(ArgumentError("f and g must be of the same length"))
    end

    #Inplace fast fourier transform
    #Unfortunately rfft! is not implemented at this time so must use fft!. 
    fft!(@view(f[:,1]), 1)
    fft!(@view(g[:,1]), 1)

    #Elementwise multiplication in frequency domain. 
    #Overwrite output to f so we don't need to allocate
    #This is the cross correlation theorem. corr(f, g) = ifft(fft(f)*conj(fft(g)))
    #f⋆g = ℱ[conj(F(v))G(v)]
    for i in range(1, size(f)[1])
        f[i, 1] = conj((f[i, 1]))*(g[i, 1])
    end

    #Get Cross Correlation
    ifft!(f)
    cross_corr = f #Rename 

    #########
    #Get Offset with Maximum Cross Correlation
    max = 0.0
    offset = 0
    for i in eachindex(cross_corr)
        if real(cross_corr[i, 1]) > max
            max = real(cross_corr[i, 1])
            offset = i
        end
    end


    if offset > length(f)÷2
        return -abs(length(f) - offset)*δt
    else
        return (offset)*δt
    end
end

function pearson_corr(f::Matrix{Complex{T}}, g::Matrix{Complex{T}}) where {T<:AbstractFloat}
    mean1, mean2, = 0.0, 0.0 # = LinearAlgebra.norm(f)*LinearAlgebra.norm(g)#sqrt(sum((f .- mean(f)).^2))*sqrt(sum((g .- mean(g)).^2))
    for i in range(1, length(f))
        mean1 += real(f[i, 1])
        mean2 += real(g[i, 1])
    end

    mean1, mean2 = mean1/length(f), mean2/length(f)

    norm1, norm2, dot = 0.0, 0.0, 0.0

    for i in range(1, length(f))
        fᵢ, gᵢ = (real(f[i, 1]) - mean1), (real(g[i, 1]) - mean2)
        dot += fᵢ*gᵢ
        norm1 += (fᵢ)^2
        norm2 += (gᵢ)^2
    end

    return dot/sqrt(norm1*norm2)
end

function pearson_corr(ms1::Matrix{Complex{T}}, ms2::Matrix{Complex{T}}, ms1_params::NTuple{4, Union{T,Missing}}, ms2_params::NTuple{4, Union{T,Missing}}; N::Int64 = 500, width_t::Float64 = 2.0) where {T<:AbstractFloat}
   
    if any(ismissing.(ms1_params)) | any(ismissing.(ms2_params))
        return missing
    end
    #Points at which to evaluate the function. 
    #May need to re-evalute how we find the window center
    times = LinRange(T(ms2_params[2] - width_t), T(ms2_params[2] + width_t), N)

    #Fill ms1 and ms2 in place with values for teh exponential-gaussian hybrid function
    EGH!(ms1, times, ms1_params)
    EGH!(ms2, times, ms2_params)

    #Pearson Correlation
    return pearson_corr(ms1, ms2)
        
    #δt = CrossCorrFFT(ms1, ms2, T(width_t*2/N))

    #return δt, ρ

end

function pearson_corr!(psms::DataFrame; N::Int64 = 500, width_t::Float64 = 2.0) where {T<:AbstractFloat}
    ms1 = zeros(Complex{Float32}, N, 1)
    ms2 = zeros(Complex{Float32}, N, 1)
    psms[:,:ρ] = Vector{Union{Missing, Float32}}(undef, size(psms)[1])
    i = 0
    for psm in eachrow(psms)
        i += 1
        if ismissing(psm[:τ_ms1])
            psms[i,:ρ] = missing
            continue
        end
        ms1_params = (Float32(psm[:σ_ms1]),
        Float32(psm[:tᵣ_ms1]),
        Float32(psm[:τ_ms1]),
        Float32(psm[:H_ms1]))

        ms2_params = (Float32(psm[:σ]),
                        Float32(psm[:tᵣ]),
                        Float32(psm[:τ]),
                        Float32(psm[:H])
        )
        
        psms[i,:ρ] = pearson_corr(ms1, ms2, ms1_params, ms2_params, N = N, width_t = width_t)
    end
end

"""
Lan K, Jorgenson JW. A hybrid of exponential and gaussian functions as a simple model of asymmetric chromatographic peaks. J Chromatogr A. 2001 Apr 27;915(1-2):1-13. doi: 10.1016/s0021-9673(01)00594-5. PMID: 11358238.
"""
function EGH(t::Vector{T}, p::NTuple{4, T}) where {T<:AbstractFloat}
    #σ=p[1], tᵣ=p[2], τ=p[3], H = p[4]
     y = zeros(T, length(t))
     for (i, tᵢ) in enumerate(t)
        d = 2*p[1] + p[3]*(tᵢ - p[2])
        if d > 0
            y[i] = p[4]*exp((-(tᵢ - p[2])^2)/d)
        end
    end
    return y
end

function EGH(t::T, p::NTuple{4, T}) where {T<:AbstractFloat}
    #σ=p[1], tᵣ=p[2], τ=p[3], H = p[4]
    d = 2*p[1] + p[3]*(t - p[2])
    if d > 0
        return p[4]*exp((-(t - p[2])^2)/d)
    else
        return zero(T)
    end
end

"""
Lan K, Jorgenson JW. A hybrid of exponential and gaussian functions as a simple model of asymmetric chromatographic peaks. J Chromatogr A. 2001 Apr 27;915(1-2):1-13. doi: 10.1016/s0021-9673(01)00594-5. PMID: 11358238.
"""
function EGH!(f::Matrix{Complex{T}}, t::LinRange{T, Int64}, p::NTuple{4, T}) where {T<:AbstractFloat}
    #σ=p[1], tᵣ=p[2], τ=p[3], H = p[4]
    #Given parameters in 'p'
    #Evaluate EGH function at eath time point tᵢ and store them in pre-allocated array 'f'. 
     for (i, tᵢ) in enumerate(t)
        d = 2*p[1] + p[3]*(tᵢ - p[2])
        if real(d) > 0
            f[i] = p[4]*exp((-(tᵢ - p[2])^2)/d)
        else
            f[i] = zero(Complex{T})
        end
    end
    return nothing
end


function EGH_inplace(F::Vector{T}, x::Vector{T}, p::Vector{T}) where {T<:AbstractFloat}
    #σ=p[1], tᵣ=p[2], τ=p[3], H = p[4]
    #Given parameters in 'p'
    #Evaluate EGH function at eath time point tᵢ and store them in pre-allocated array 'f'. 
     for (i, tᵢ) in enumerate(x)
        d = 2*p[1] + p[3]*(tᵢ - p[2])
        if real(d) > 0
            F[i] = p[4]*exp((-(tᵢ - p[2])^2)/d)
        else
            F[i] = zero(T)
        end
    end
end

function JEGH_inplace(J::Matrix{T}, x::Vector{T}, p::Vector{T}) where {T<:AbstractFloat}
    for (i, tᵢ) in enumerate(x)
        δt = tᵢ - p[2]
        d = 2*p[1] + p[3]*δt
        q = δt/d
        f = exp((-δt^2)/(d))
        #f = exp((-δt^2)/(d))
        if d > 0
            J[i,1] = (2*(q)^2)*p[4]*f
            J[i,2] = p[4]*(2*q - (p[3]*(q^2)))*f
            J[i,3] = ((δt)*(q^2))*p[4]*f
            J[i,4] = f
        end
   end
end

"""
Lan K, Jorgenson JW. A hybrid of exponential and gaussian functions as a simple model of asymmetric chromatographic peaks. J Chromatogr A. 2001 Apr 27;915(1-2):1-13. doi: 10.1016/s0021-9673(01)00594-5. PMID: 11358238.
"""
function JEGH(t::Vector{T}, p::NTuple{4, T}) where {T<:AbstractFloat}
    #σ=p[1], tᵣ=p[2], τ=p[3], H = p[4]
    J = zeros(T, (length(t), length(p)))
    for (i, tᵢ) in enumerate(t)
        δt = tᵢ - p[2]
        d = 2*p[1] + p[3]*δt
        q = δt/d
        f = exp((-δt^2)/(d))
        #f = exp((-δt^2)/(d))
        if d > 0
            J[i,1] = (2*(q)^2)*p[4]*f
            J[i,2] = p[4]*(2*q - (p[3]*(q^2)))*f
            J[i,3] = ((δt)*(q^2))*p[4]*f
            J[i,4] = f
        end
   end
   return J
end

=#