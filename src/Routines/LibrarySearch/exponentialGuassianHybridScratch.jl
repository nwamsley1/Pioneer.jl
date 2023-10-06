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

function Integrate(f::Function, p::NTuple{4, T}; α::AbstractFloat = 0.001, n::Int64 = 1000) where {T<:AbstractFloat}

    #Use GuassLegendre Quadrature to integrate f on the integration bounds 
    #using FastGaussQuadrature
    x, w = gausslegendre(n)

    function getBase(α::AbstractFloat, τ::T, σ::T) where {T<:AbstractFloat}
        B = (-1/2)*(sqrt(abs(log(α)*((τ^2)*log(α) - 8*σ) + τ*log(α))))
        A = (1/2)*( τ*log(α) - sqrt(abs(log(α)*((τ^2)*log(α) - 8*σ))))
        return abs(A), abs(B)
    end

    A, B = getBase(α, p[3], p[1])
    a = p[2] - A
    b = p[2] + B
    #Quadrature rules are for integration bounds -1 to 1 but can shift
    #to arbitrary bounds a and b. 
    return ((b - a)/2)*dot(w, f(Float32.(x.*((b - a)/2) .+ (a + b)/2), p))
end

function getP0(α::T, B::T, A::T, tᵣ::T, H::T) where {T<:AbstractFloat}
    return T[(-1/(2*log(α)))*(B*A),
             tᵣ,
             (-1/log(α))*(B - A),
             H]
end

getP0(p0::NTuple{5, T}) where {T<:AbstractFloat} = getP0(p0[1], p0[2],p0[3],p0[4],p0[5])

function getFWHM(α::AbstractFloat, τ::T, σ::T) where {T<:AbstractFloat}
    #FWHM given parameters for EGH function 

    #When is absolute value necessary. 
    B = (-1/2)*(sqrt(abs(log(α)*((τ^2)*log(α) - 8*σ) + τ*log(α))))
    A = (1/2)*( τ*log(α) - sqrt(abs(log(α)*((τ^2)*log(α) - 8*σ))))
    return abs(A) + abs(B)
end

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
        ms1_params = (psm[:σ_ms1],
                        psm[:tᵣ_ms1],
                        psm[:τ_ms1],
                        psm[:H_ms1])

        if any(ismissing.(ms1_params))
            psms[i,:ρ] = missing
            continue
        end

        ms2_params = (psm[:σ],
                        psm[:tᵣ],
                        psm[:τ],
                        psm[:H]
        )
        
        psms[i,:ρ] = pearson_corr(ms1, ms2, ms1_params, ms2_params, N = N, width_t = width_t)
    end
end

function getBestPSM(sdf::SubDataFrame)
    if any(sdf.q_value.<=0.01)
        return sdf[argmax((sdf.q_value.<=0.01).*(sdf.hyperscore)),:]
    else
        return sdf[argmax(sdf.prob),:]
    end
end

ms1 = zeros(Complex{Float32}, 500, 1)
ms2 = zeros(Complex{Float32}, 500, 1)
select!(best_psms, Not(:ρ))
select!(best_psms, Not(:δt))
best_psms[:,:δt] = Vector{Union{Missing, Float32}}(undef, size(best_psms)[1])#zeros(Float32, size(best_psms)[1])
best_psms[:,:ρ] = Vector{Union{Missing, Float32}}(undef, size(best_psms)[1])#zeros(Float32, size(best_psms)[1])





best_psms[best_psms[:,:GOF].>0.99,[:GOF,:points_above_FWHM,:points_above_FWHM_01,:sequence,:q_value,:decoy]]

value_counts(df, col) = combine(groupby(df, col), nrow)
value_counts(best_psms, :precursor_idx)

best_psms = combine(sdf -> sdf[argmax(sdf.prob),:], groupby(PSMs[ms_file_idx][PSMs[ms_file_idx][:,:q_value].<=1.0,:], [:sequence,:charge]));
  

grouped_df = groupby(best_psms, :sequence);
best_psms[:,:n_precs] = (combine(grouped_df) do sub_df
    #if size(sub_df)[1] == 28
    #    display(sub_df[:,[:precursor_idx,:sequence,:charge]])
    #end
    repeat([size(sub_df)[1]], size(sub_df)[1])
end)[:,:x1]

ms1 = zeros(Complex{Float32}, 500, 1)
ms2 = zeros(Complex{Float32}, 500, 1)
select!(best_psms, Not(:ρ))
select!(best_psms, Not(:δt))
best_psms[:,:δt] = Vector{Union{Missing, Float32}}(undef, size(best_psms)[1])#zeros(Float32, size(best_psms)[1])
best_psms[:,:ρ] = Vector{Union{Missing, Float32}}(undef, size(best_psms)[1])#zeros(Float32, size(best_psms)[1])
for i in range(1, size(best_psms)[1])


end

select!(best_psms, Not(:prec_ρ))
best_psms[:,:prec_ρ] = Vector{Union{Missing, Float32}}(undef, size(best_psms)[1])#zeros(Float32, size(best_psms)[1])
grouped_df = groupby(best_psms, :sequence);
best_psms[:,:prec_ρ] = (combine(grouped_df) do sub_df

    if size(sub_df)[1] == 1
        return repeat([missing], size(sub_df)[1])
    end

    precs = sortperm(sub_df[:,:peak_area], rev=true)[1:2]

    p1_params = [sub_df[precs[1],:σ],
                  sub_df[precs[1],:tᵣ],
                  sub_df[precs[1],:τ],
                  sub_df[precs[1],:H]
    ]

    p2_params = [sub_df[precs[2],:σ],
                  sub_df[precs[2],:tᵣ],
                  sub_df[precs[2],:τ],
                  sub_df[precs[2],:H]
    ]

    if (sum(ismissing.(p1_params))>0) |  (sum(ismissing.(p2_params))>0)
        return repeat([missing], size(sub_df)[1])
    end

    δt, ρ = CrossCorrMS1vMS2(ms1, ms2, p1_params, p2_params)

    return repeat([ρ], size(sub_df)[1])

end)[:,:x1]

sort(best_psms,:sequence)[:,[:precursor_idx,:sequence,:n_precs,:charge,:ρ]]
best_psms[(best_psms[:,:peak_area].>(10^7.9)),[:peak_area,:peak_area_ms1,:precursor_idx,:sequence,:ρ]] #.&(best_psms[:,:peak_area_ms1].<10^5),:]

sort(PSMs[1][PSMs[1][:,:precursor_idx].== 5455685,[:weight,:prob,:q_value,:RT]],:weight)
test = PSMs[(PSMs[:,:sequence].=="IWHHTFYNELR"),:] 
best_psms[best_psms[:,:GOF].>0.99,[:precursor_idx,:GOF,:points_above_FWHM,:points_above_FWHM_01,:sequence,:q_value,:decoy]]
#=
sort(PSMs[(PSMs[:,:precursor_idx].==UInt32(best_psms_passing[N,:precursor_idx])).&(PSMs[:,:q_value].<=0.01),[:sequence,:prob,:q_value,:hyperscore,:RT]],:hyperscore)

sort(PSMs[(PSMs[:,:sequence].=="IWHHTFYNELR").&(PSMs[:,:q_value].<=0.01),[:sequence,:prob,:q_value,:hyperscore,:RT]],:hyperscore)

best_psms_passing = best_psms[best_psms[:,:q_value].<=0.01,:]
integratePrecursor(ms1_chroms, UInt32(best_psms_passing[N,:precursor_idx]), (0.1f0, 0.15f0, 0.15f0, Float32(best_psms_passing[N,:RT]), best_psms_passing[N,:weight]), isplot = true)

integratePrecursor(ms2_chroms, UInt32(11328380, (0.1f0, 0.15f0, 0.15f0, Float32(best_psms_passing[N,:RT]), best_psms[N,:weight]), isplot = true)


integratePrecursor(ms1_chroms, UInt32(11328380), (0.1f0, 0.15f0, 0.15f0, Float32(66.2004), Float32(1e4)), isplot = true)
ms1_chroms[(precursor_idx=UInt32(11328380),)][:,:]


integratePrecursor(ms2_chroms, UInt32(11328380), (0.1f0, 0.15f0, 0.15f0, Float32(66.2004), Float32(1e4)), isplot = true)
ms2_chroms[(precursor_idx=UInt32(11328380),)][:,:]

N = 1000
integratePrecursor(ms2_chroms,UInt32(best_psms[N,:precursor_idx]), (0.1f0, 0.15f0, 0.15f0, Float32(66.2004), Float32(1e4)), isplot = true)
ms2_chroms[(precursor_idx=UInt32(best_psms[N,:precursor_idx]),)][:,:]
N += 1

integratePrecursor(ms1_chroms_huber, UInt32(11328380), (0.1f0, 0.15f0, 0.15f0, Float32(66.2004), Float32(1e4)), isplot = true)
ms1_chroms_huber[(precursor_idx=UInt32(11328380),)][:,:]

huber_loss = ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:]
N  += 1

PSMs[(PSMs[:,:precursor_idx].==UInt32(best_psms_passing[N,:precursor_idx])).&(PSMs[:,:q_value].<=0.01),[:sequence,:prob,:q_value,:hyperscore,:RT]]
best_psms[(best_psms[:,:decoy].==true).&(best_psms[:,:δt].==0.0),:]
N = findfirst(x->x== 5455685  , best_psms[:,:precursor_idx])
huber_loss = ms2_chroms[(precursor_idx=UInt32(best_psms_passing[N,:precursor_idx]),)][:,:]
ms1 = ms1_chroms[(precursor_idx=UInt32(best_psms[N,:precursor_idx]),)][:,:]
plot(huber_loss[:,:rt],
huber_loss[:,:weight], seriestype=:scatter,
alpha = 0.5)
plot!(ms1[:,:rt],
ms1[:,:weight], seriestype=:scatter,
alpha = 0.5)



integratePrecursor(ms2_chroms, UInt32(best_psms_passing[N,:precursor_idx]), (0.1f0, 0.15f0, 0.15f0, Float32(best_psms_passing[N,:RT]), best_psms[N,:weight]), isplot = true)
N  += 1

=#

bins = LinRange(0, 0.6, 100)
histogram(PSMs[(PSMs[:,:decoy]),:q_value], bins = bins, alpha = 0.5, normalize = :probability)
histogram!(PSMs[(PSMs[:,:decoy].==false).&(PSMs[:,:q_value].<=0.25),:q_value], bins = bins, alpha = 0.5, normalize = :probability)

frag_ppm_errs = [_getPPM(match.theoretical_mz, match.match_mz) for match in best_matches];
frag_intensities = [match.intensity for match in best_matches];
frag_mz =  [match.theoretical_mz for match in best_matches];

plot(sqrt.(frag_mz), frag_ppm_errs.*sqrt.(frag_intensities), seriestype=:scatter, alpha = 0.1)
plot((frag_mz), frag_ppm_errs.*sqrt.(frag_intensities), seriestype=:scatter, alpha = 0.1)

plot(sqrt.(frag_intensities), frag_ppm_errs, seriestype=:scatter, alpha = 0.1)

histogram2d(frag_ppm_errs, log2.(frag_intensities))

histogram(frag_ppm_errs.*sqrt.(frag_intensities).*sqrt.(frag_mz))

frag_err_dist = estimateErrorDistribution(frag_ppm_errs, Laplace{Float64}, 0.0, 3.0, 30.0,
f_out = "" );

frag_err_dist = estimateErrorDistribution((frag_ppm_errs .- 3.385).*sqrt.(frag_intensities), Laplace{Float64}, 0.0, 3.0, 30.0,
f_out = "" );

frag_err_dist = estimateErrorDistribution((frag_ppm_errs .- 3.385)./sqrt.(frag_intensities), Laplace{Float64}, 0.0, 3.0, 30.0,
f_out = "" );


frag_err_dist = estimateErrorDistribution((frag_ppm_errs .- 3.385).*sqrt.(1000*rand(length(frag_intensities))), Laplace{Float64}, 0.0, 3.0, 30.0,
f_out = "" );

plot(abs.((frag_ppm_errs .- 3.385)),
        sqrt.(frag_intensities), seriestype=:scatter, alpha = 0.1, ylim = (0, 1e3))

estimateErrorDistribution((frag_ppm_errs .- 3.385).*sqrt.(frag_intensities), Laplace{Float64}, 0.0, 3.0, 30.0,
f_out = "" );

bins = LinRange(-15, -5, 100)
histogram(best_psms[(best_psms[:,:decoy]),:err_likelihood_norm], bins = bins, alpha = 0.5, normalize = :probability)
histogram!(best_psms[(best_psms[:,:decoy].==false).&(best_psms[:,:q_value].<=0.01),:err_likelihood_norm], bins = bins, alpha = 0.5, normalize = :probability)


bins = LinRange(-30, 30, 50)
histogram(frag_ppm_errs[frag_intensities.<1e5], bins = bins, alpha = 0.5, normalize = :probability)
histogram!(frag_ppm_errs[frag_intensities.>1e6], bins = bins, alpha = 0.5, normalize = :probability)


bins = LinRange(-15, -5, 100)
histogram(PSMs[(PSMs[:,:decoy]),:err_likelihood_norm], bins = bins, alpha = 0.5, normalize = :probability)
histogram!(PSMs[(PSMs[:,:decoy].==false).&(PSMs[:,:q_value].<=0.01),:err_likelihood_norm], bins = bins, alpha = 0.5, normalize = :probability)

bins = LinRange(-5e8, 0, 100)
histogram(PSMs[(PSMs[:,:decoy]),:err_norm], bins = bins, alpha = 0.5, normalize = :probability)
histogram!(PSMs[(PSMs[:,:decoy].==false).&(PSMs[:,:q_value].<=0.01),:err_norm], bins = bins, alpha = 0.5, normalize = :probability)

bins = LinRange(-15e8, 0, 100)
histogram(PSMs[(PSMs[:,:decoy]),:error], bins = bins, alpha = 0.5, normalize = :probability)
histogram!(PSMs[(PSMs[:,:decoy].==false).&(PSMs[:,:q_value].<=0.01),:error], bins = bins, alpha = 0.5, normalize = :probability)
