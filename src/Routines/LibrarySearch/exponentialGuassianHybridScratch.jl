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

function EGH!(f::Matrix{Complex{T}}, t::LinRange{T, Int64}, p::Vector{T}) where {T<:AbstractFloat}
    #σ=p[1], tᵣ=p[2], τ=p[3], H = p[4]
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

x = huber_loss[:,:rt]
y = huber_loss[:,:weight]
fit1 = curve_fit(EGH, JEGH, Float64.(x), Float64.(y), pg)
#Integrate(EGH, fit1.param, (fit1.param[2] - 1.0, fit1.param[2] + 1.0))
#plot!(collect(ts), EGH(collect(ts), fit0.param))
ts = LinRange(fit1.param[2] - 2.0, fit1.param[2] + 2.0, 2000)
plot(collect(ts), EGH(collect(ts), (fit1.param)))

f = EGH(collect(ts), (fit1.param))
fit1.param[2] = fit1.param[2] - 0.5
g = EGH(collect(ts), (fit1.param))
plot(collect(ts),f)
plot!(collect(ts),g)

@time EGH!(f, ts, fit1.param)

function getOffset(f::Matrix{Complex{T}}, g::Matrix{Complex{T}}, δt::T) where {T<:AbstractFloat}
    if length(f) != length(g)
        #Need to throw and error here
        return 
    end

    mean1, mean2, = 0.0, 0.0 # = LinearAlgebra.norm(f)*LinearAlgebra.norm(g)#sqrt(sum((f .- mean(f)).^2))*sqrt(sum((g .- mean(g)).^2))
    for i in range(1, length(f))
        mean1 += f[i, 1]
        mean2 += g[i, 1]
    end

    mean1, mean2 = mean1/length(f), mean2/length(f)

    norm1, norm2 = 0.0, 0.0
    for i in range(1, length(f))
        norm1 += (f[i, 1] - mean1)^2
        norm2 += (g[i, 1] - mean2)^2
    end

    norm1, norm2 = sqrt(norm1), sqrt(norm2)
    norm = norm1*norm2
    #norm = 1

    #Inplace fast fourier transform
    fft!(@view(f[:,1]), 1)
    fft!(@view(g[:,1]), 1)

    #Elementwise multiplication in frequency domain
    for i in range(1, size(f)[1])
        f[i, 1] = conj((f[i, 1]))*(g[i, 1])
    end

    #Get Cross Correlation
    ifft!(f)
    cross_corr = f
    #plot(real.(f)/norm, show = true)
    #########
    #Get Maximum Cross Correlation
    max = 0.0
    offset = 0
    for i in eachindex(cross_corr)
        if real(cross_corr[i, 1]) > max
            max = real(cross_corr[i, 1])
            offset = i
        end
    end

    #cross_corr = irfft(conj.(rfft(f)).*rfft(g), length(f))./norm

    if offset > length(f)÷2
        return -abs(length(f) - offset)*δt, max/norm
    else
        return offset*δt, max/norm
    end
end

x = huber_loss[:,:rt]
y = huber_loss[:,:weight]
fit1 = curve_fit(EGH, JEGH, Float64.(x), Float64.(y), pg)
N = 1000
ts = LinRange(fit1.param[2] - 2.0, fit1.param[2] + 2.0, N)
f = zeros(Complex{Float64}, N, 1)
g = zeros(Complex{Float64}, N, 1)
f[:,1] .= complex.(EGH(collect(ts), (fit1.param)))
fit1.param[2] = fit1.param[2] - 1.0
fit1.param[3] += 1.0
fit1.param[1] = 0.0
fit1.param[4] = fit1.param[4]/5000000
g[:,1] .= complex.(EGH(collect(ts), (fit1.param)))
#g = rand(Complex{Float64}, N, 1)
norm = LinearAlgebra.norm(real.(f))*LinearAlgebra.norm(real.(g))
f = f./norm
f = f .- mean(f)
g = g./norm
g = g .- mean(g)
plot(real.(g), seriestype=:scatter)
plot!(real.(f), seriestype=:scatter)
@time getOffset(f, g, 4.0/N)


#FFT in place  
fft!(@view(f[:,1]), 1)
fft!(@view(g[:,1]), 1)
for i in range(1, size(f)[1])
    f[i, 1] = conj(f[i, 1])*g[i, 1]
end
plot(real.(ifft!(f)))

p = plan_fft!(zeros(Float64, 2000); flags=FFTW.ESTIMATE, timelimit=Inf)


buf_FFT = rand(Complex{Float64}, 10,10)
#in-place FFT on the first column of the matrix
fft!( view(buf_FFT, 1:size(buf_FFT, 1),1), 1 ) 
#alternately, to save some typing:
ifft!( view(buf_FFT, Colon(), 1), 1 )
#Or, even simpler, use the @view macro
ifft!(@view(buf_FFT[:,1]), 1)

getOffset(f, g, 4.0/2000)
getOffset(g, f, 4.0/2000)
getOffset(f, f, 4.0/2000)


norm = sqrt(sum((f .- mean(f)).^2))*sqrt(sum((g .- mean(g)).^2))
cross_cor = irfft(conj.(rfft(f_pad)).*rfft(g_pad), 1999)./norm
(1999 - argmax(cross_cor))*((maximum(x) - minimum(x))/1000)
plot(cross_cor)

N = 2000
f = f[:,1]
g = g[:,1]
norm = sqrt(sum((f .- mean(f)).^2))*sqrt(sum((g .- mean(g)).^2))
cross_cor = irfft(conj.(rfft(f)).*rfft(g), N)./norm
(N - argmax(cross_cor))*((maximum(x) - minimum(x))/N)
plot(cross_cor)


norm = sqrt(sum((f .- mean(f)).^2))*sqrt(sum((f .- mean(f)).^2))
cross_cor = irfft(conj.(rfft(f)).*rfft(f), N)./norm
(N- argmax(cross_cor))*((maximum(x) - minimum(x))/N)
plot!(cross_cor)

norm = sqrt(sum((f .- mean(f)).^2))*sqrt(sum((g .- mean(g)).^2))
cross_cor = irfft(conj.(rfft(g)).*rfft(f), N)./norm
(N - argmax(cross_cor))*((maximum(x) - minimum(x))/N)
plot!(cross_cor)




norm = sqrt(sum((f_pad .- mean(f_pad)).^2))*sqrt(sum((f_pad .- mean(f_pad)).^2))
cross_cor = irfft(conj.(rfft(f_pad)).*rfft(f_pad), 1999)./norm
(1999 - argmax(cross_cor))*((maximum(x) - minimum(x))/1000)



norm = sqrt(sum((f .- mean(f)).^2))*sqrt(sum((g .- mean(g)).^2))
cross_cor = irfft(conj.(rfft(f)).*rfft(g), 1000)./norm
(1000 - argmax(cross_cor))*((maximum(x) - minimum(x))/1000)


norm = sqrt(sum((f_pad .- mean(f_pad)).^2))*sqrt(sum((f_pad .- mean(f_pad)).^2))
cross_cor = irfft(conj.(rfft(f_pad)).*rfft(f_pad), 1999)./norm
(1999 - argmax(cross_cor))*((maximum(x) - minimum(x))/1000)





cyclic_conv = ifft(fft(f).*fft(g))
plot(collect(ts), real.(cyclic_conv))

plot(collect(ts), imag.(cyclic_conv))
plot(collect(ts),f)
plot(collect(ts),g)
argmax(real.(cyclic_conv))


FastCyclicConv1D(f,g) = ifft(fft(f).*fft(g))
function FastLinearConvolution(f,g)
    N = length(f)
    M = length(g)

    f_pad = [ f; zeros(M-1) ]     
    g_pad = [ g; zeros(N-1) ]     

    return FastCyclicConv1D( f_pad, g_pad )
end

plot(collect(ts),f)
plot!(collect(ts),g)
(maximum(x) - minimum(x))/1000
(1999 - argmax(irfft(conj.(rfft(f_pad)).*rfft(g_pad), 1999)./norm))*((maximum(x) - minimum(x))/1000)

norm = sqrt(sum((f_pad .- mean(f_pad)).^2))*sqrt(sum((g_pad .- mean(g_pad)).^2))
cross_cor = irfft(conj.(rfft(f_pad)).*rfft(g_pad), 1999)./norm


(maximum(x) - minimum(x))/1000
(1999 - argmax(irfft(conj.(rfft(f_pad)).*rfft(g_pad), 1999)./norm))*((maximum(x) - minimum(x))/1000)

plot(LinRange(minimum(x), maximum(x), 1999),
cross_cor)

sqrt(sum((f_pad .- mean(f_pad)).^2))*sqrt(sum((g_pad .- mean(g_pad)).^2))

N = length(f)
M = length(g)

f_pad = [ f; zeros(M-1) ]     
g_pad = [ g; zeros(N-1) ]     

LinRange(minimum(x), maximum(x), 1000)[argmax(irfft(conj.(rfft(f_pad)).*rfft(g_pad), 1999))]
LinRange(minimum(x), maximum(x), 1000)[argmax(g)]


plot(1:1001, crosscor(f, g, -500:500))

plot!(1:1001, crosscor(f, f, -500:500))

plot(LinRange(minimum(x), maximum(x), 1999),
real.(FastLinearConvolution(f, f)))
plot!(LinRange(minimum(x), maximum(x), 1000), f*0.5e8)

plot!(collect(ts),f*0.5e8)
plot!(collect(ts),g*0.5e8)

plot(LinRange(minimum(x), maximum(x), 1000),
irfft(rfft(f).*rfft(g), 1000))
plot!(collect(ts),f*0.5e8)
plot!(collect(ts),g*0.5e8)

LinRange(minimum(x), maximum(x), 1999)[argmax(real.(FastLinearConvolution(f, g)))] - LinRange(minimum(x), maximum(x), 1999)[argmax([ f; zeros(1000-1) ] )]

LinRange(minimum(x), maximum(x), 1000)[argmax(irfft(rfft(f).*rfft(g), 1000))] - LinRange(minimum(x), maximum(x), 1000)[argmax(g)]
LinRange(minimum(x), maximum(x), 1000)[argmax(g)] - LinRange(minimum(x), maximum(x), 1000)[argmax(f)]


plot(LinRange(minimum(x), maximum(x), 1000),
irfft(rfft(f).*rfft(f), 1000))
plot!(collect(ts),f*0.5e8)
plot!(collect(ts),g*0.5e8)



N += 1
(1000 - 864)*((maximum(x) - minimum(x))/1000)


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