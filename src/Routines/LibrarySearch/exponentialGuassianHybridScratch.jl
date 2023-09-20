function EGH(t::Vector{T}, p::Vector{T}) where {T<:AbstractFloat}
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

function JEGH(t::Vector{T}, p::Vector{T}) where {T<:AbstractFloat}
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

function Integrate(f::Function, p::Vector{T}, bounds::Tuple{T, T}, n::Int64 = 1000) where {T<:AbstractFloat}
    x, w = gausslegendre(n)
    a, b = first(bounds), last(bounds)
    return ((b - a)/2)*dot(w, f(Float32.(x.*((b - a)/2) .+ (a + b)/2), p))
end

function getP0(α::T, B::T, A::T, tᵣ::T, H::T) where {T<:AbstractFloat}
    return T[(-1/(2*log(α)))*(B*A),
             tᵣ,
             (-1/log(α))*(B - A),
             H]
end

function getFWHM(α::AbstractFloat, τ::T, σ::T) where {T<:AbstractFloat}
    B = (-1/2)*(sqrt(log(α)*((τ^2)*log(α) - 8*σ)) + τ*log(α))
    A = (1/2)*( τ*log(α) - sqrt(log(α)*((τ^2)*log(α) - 8*σ)))
    return abs(A) + abs(B)
end

integratePrecursor(ms2_chroms, UInt32(best_psms_passing[N,:precursor_idx]),(best_psms_passing[N,:scan_idx]), isplot = true)
N  += 1

using FastGaussQuadrature
x, w = gausslegendre(100)
a = fit1.param[2] - 2.0
b = fit1.param[2] + 2.0
((b - a)/2)*dot(w, EGH(collect(x).*((b - a)/2) .+ (a + b)/2, fit1.param))

for prec in precursors_list
    if prec.sequence == "M(ox)YPIDFEK"
        println("Valine ", prec)
    end
    if prec.sequence == "M(ox)YPLDFEK"
        println("Isoleucine ", prec)
    end
end
VPAGLPDLK
NumericalIntegration.integrate(huber_loss[:,:rt], EGH(Float64.(collect(huber_loss[:,:rt])), fit1.param), TrapezoidalFast())
a = ts[1]
b = ts[end]
gausslegen
EGH(collect(ts), fit1.param)



pg = [0.1^2, 30.2, 1/16, 6e5]
@btime fit0 = curve_fit(EGH, Float64.(x[17:end]), Float64.(y[17:end]), pg)

pg = [0.1^2, 30.2, 1/16, 6e5]
@btime fit1 = curve_fit(EGH, JEGH, Float64.(x[17:end]), Float64.(y[17:end]), pg)

fit1 = curve_fit(EGH, JEGH, Float64.(x), Float64.(y), pg)

N = 9936

squared_error = ms2_chroms_square[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][1:35,:]
huber_loss = ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][1:35,:]
plot(squared_error[:,:rt],
squared_error[:,:weight], seriestype=:scatter,
alpha = 0.5)
plot!(huber_loss[:,:rt],
huber_loss[:,:weight], seriestype=:scatter,
alpha = 0.5)

x = huber_loss[1:30,:rt]
y = huber_loss[1:30,:weight]
fit1 = curve_fit(EGH, JEGH, Float64.(x), Float64.(y), pg)
Integrate(EGH, fit1.param, (fit1.param[2] - 1.0, fit1.param[2] + 1.0))
#plot!(collect(ts), EGH(collect(ts), fit0.param))
plot!(collect(ts), EGH(collect(ts), fit1.param))
N += 1



p = [6e5, 30.2, 0.15]
m(t, p) = p[1]*exp.((-1).*((t .- p[2]).^2)./(2*p[3]^2))
m(ts, p)
plot!(ts, m(ts, p))


using Splines2, GLM, Random
x = huber_loss[:,:rt]
y = huber_loss[:,:weight]
ns1 = Splines2.ns_(x, df = 15, intercept = true);
X = ns1(x);
fit1 = lm(X, y)

plot!(x,
GLM.predict(fit1, ns1(x)))