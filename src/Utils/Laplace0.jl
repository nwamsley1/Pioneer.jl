
struct Laplace0{T<:Real} <: ContinuousUnivariateDistributions
    θ::T
    Laplace0{T}(θ::T) where {T} = new{T}(θ)
end

function Laplace0(θ::T; check_args::Bool=true) where {T <: Real}
    @check_args Laplace0 (θ, θ > zero(θ))
    return Laplace0{T}(θ)
end

@distr_support Laplace0 -Inf Inf

function fit_mle(::Type{<:Laplace0}, x::AbstractArray{<:Real})
    xc = similar(x)
    copyto!(xc, x)
    m = median!(xc)
    xc .= abs.(x .- m)
    return Laplace0(mean(xc))
end

