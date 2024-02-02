
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

function getFragTolTest(mass::Float32, ppm_err::Float32, intensity::Float32, ppm_tol_param::Float32)
    mass -= Float32(ppm_err*mass/1e6)
    tol = (ppm_tol_param/sqrt(intensity))*mass/1e6
    println("tol $tol")
    println("ppm ", (ppm_tol_param/sqrt(intensity)))
    return Float32(mass - tol), Float32(mass + tol)
   #return Float32(mass*(1 - 16.1*mass/1e6)), Float32(mass*(1 + 16.1*mass/1e6))
end

function getPPMTest(mz::Float32, ppm_tol_param::Float32, intensity::Float32)
    tol = (ppm_tol_param/sqrt(intensity))*mz/1e6
    #tol = 40.0*mz/1e6
    #tol = (1/sqrt(intensity))*mz/1e6
    return mz - tol, mz + tol        
end

N = 100000
a = MS_TABLE[:masses][N]
sort(diff(a)./(a[2:end]./1e6))
N += 1