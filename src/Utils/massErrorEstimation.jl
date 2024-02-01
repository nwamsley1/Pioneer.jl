struct MassErrorModel{D<:Distribution,T<:AbstractFloat}
    intensity_to_shape_param::AbstractExtrapolation
    distribution::Type{D}
    location::T
end

function (mem::MassErrorModel)(x::T)::T where {T<:AbstractFloat} 
    mem.intensity_to_shape_param(x) 
end

function getLocation(mem::MassErrorModel{D, T}) where {D<:Distribution,T<:AbstractFloat}
    return mem.location
end

function getMassErrorLogLik(model::MassErrorModel{Laplace{T}, T}, intensity::T, err::T) where {T<:AbstractFloat}
    b = model(intensity)
    return -log(b*2) - abs(err)/b
end

function EstimateMixtureWithUniformNoise(errs::AbstractVector{T}, #data
                                         err_model::Type{D}, #Distribution to model error
                                         μ::T, #Fixed/known location parameter
                                         w::T, #Known uniform distribution width
                                         b0::T, #Initial Scale estimate
                                         z::T, #Initial mixture estimate
                                         γ::AbstractVector{Bool}; #Latent state variable 
                                         max_iter::Int = 100,
                                         min_iter::Int = 10,
                                         atol::Float64 = 1e-6) where {T<:AbstractFloat,D<:Distribution}
    logpdf = zero(T)
    L = err_model(μ, b0)
    U = Distributions.Uniform(-w, w) #Fixed uniform distribution
    runif = Uniform(0, 1) 
    b = b0
    @inbounds @fastmath for err in errs #Calculate logpdf
        logpdf += log(pdf(L, err)*z + (1 - z)*pdf(U, err))
    end
    n = length(errs)
    for iter in range(1, max_iter)
        L = err_model(μ, b) #Current err distibution
        Threads.@threads for i in range(1, n)
            err = errs[i]
            p0, p1 = z*pdf(L, err),  (1 - z)*pdf(U, err)
            if rand(runif) < p0/(p0 + p1)
                γ[i] = true
            else
                γ[i] = false
            end
        end

        #Update params
        SAD = 0.0 #Sum of absolute deviations
        N = 0
        @inbounds @fastmath for (i, err) in enumerate(errs)
            SAD += γ[i]*abs(err - μ)
            N += γ[i]
        end
        b = SAD/N
        z = sum(γ)/length(γ)

        oldlogpdf = logpdf
        logpdf = zero(T)
        @inbounds @fastmath for err in errs #Calculate logpdf
            logpdf += log(pdf(L, err)*z + (1 - z)*pdf(U, err))
        end

        if (abs((oldlogpdf - logpdf)/oldlogpdf) < atol) & (iter > min_iter)
            
            break
        end
    end
    
    return L, z

end

function ModelMassErrs(intensities::Vector{T},
                       ppm_errs::Vector{U},
                       frag_tol::U;
                       n_intensity_bins::Int = 10) where {T,U<:AbstractFloat}
    log2_intensities = log2.(intensities)

    bins = cut(log2_intensities, n_intensity_bins)
    err_df = DataFrame(Dict(
            :ppm_errs => ppm_errs, 
            :log2_intensities => log2_intensities, 
            :bins => bins)
            )
    err_df[!,:γ] .= false
    median_intensities = zeros(T, n_intensity_bins)
    shape_estimates = zeros(T, n_intensity_bins)
    μ = median(err_df[!,:ppm_errs]) #global location estimate
    bin_idx = 0
    for (int_bin, subdf) in pairs(groupby(err_df, :bins))
        bin_idx += 1 #Intensity bin counter
        median_intensities[bin_idx] = median(subdf[!,:log2_intensities])
        b = mean(abs.(subdf[!,:ppm_errs] .- μ)) #Mean absolute deviation estimate
        L, z = EstimateMixtureWithUniformNoise(
            subdf[!,:ppm_errs],
            Laplace{Float64},
            μ,
            frag_tol,
            b,
            0.5, #mixture estimate
            subdf[!,:γ] 
        )
        #println("L.θ ", L.θ, " z ", z)
        #println("median_intensities[bin_idx] ", 2^median_intensities[bin_idx] )
        #println("cutof, ", L.θ*log(1*z*2*L.θ*80/(1-z)))
        println("quantile ", quantile(Laplace(0.0, L.θ), 0.975))
        shape_estimates[bin_idx] =  quantile(Laplace(0.0, L.θ), 0.975)# L.θ
    end


    intensities = 2 .^ median_intensities;
    println(intensities)
    #=
    return extrapolate(
            Interpolations.scale(
                interpolate(hcat(intensities, shape_estimates), 
                            #(BSpline(Quadratic(Natural(OnGrid()))), NoInterp())
                            BSpline(Quadratic(InPlace(OnGrid())))
                            ),
                LinRange(0, intensities, 100), 1:2), Flat()
            )
            intensities = Float32[11034.614, 15646.866, 20811.754, 27601.855, 36833.156, 49842.383, 70764.81, 106766.65, 192568.69, 661740.94]
                plot([x for x in intensities],
        [mass_err_model(x) for x in intensities])
                        plot!([x for x in intensities],
        [2931/sqrt(x) for x in intensities]
    )
(1 ./sqrt.(intensities))\errs
    errs = [mass_err_model(x) for x in intensities]
    sqr
    vline!([intensities[5]])
    vline!([intensities[6]])
    =##


    return MassErrorModel(
        LinearInterpolation(intensities, shape_estimates, extrapolation_bc = Flat()),
        Laplace{T},
        Float32(μ)
    )
    #return LinearInterpolation(intensities, shape_estimates, extrapolation_bc = Flat())
    #return Interpolations.linear_interpolation(median_intensities, shape_estimates, extrapolation_bc = Flat())#median_intensities, shape_estimates
end


#=
@time intensity_to_scale = ModelMassErrs(frag_ppm_intensities,
       frag_ppm_errs,
       40.0)
getMassErrorLogLik(testmodel, Float32(1e6), 1.0f0)



=#

