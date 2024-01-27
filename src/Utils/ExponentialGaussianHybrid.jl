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

function updateParams(params::HuberParams{T}, lower::HuberParams{T}, upper::HuberParams{T}) where {T<:AbstractFloat}
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
