struct MassErrorModel{T<:AbstractFloat}
    power::T
    factor::T
    location::T
end

function getMassCorrection(mem::MassErrorModel{T}) where {T<:AbstractFloat}
    return mem.location
end

function getLocation(mem::MassErrorModel{T}) where {T<:AbstractFloat}
    return mem.location
end

function (mem::MassErrorModel{T})(intensity::T) where {T<:AbstractFloat}
    return mem.factor*(intensity^(mem.power))
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
                       n_intensity_bins::Int = 10,
                       frag_err_quantile::Float64 = 0.999,
                       out_fdir::String = "./",
                       out_fname = "mass_err_estimate") where {T,U<:AbstractFloat}
    log2_intensities = log2.(intensities)

    bins = cut(log2_intensities, n_intensity_bins)
    err_df = DataFrame(Dict(
            :ppm_errs => ppm_errs, 
            :log2_intensities => log2_intensities, 
            :bins => bins)
            )
    err_df[!,:γ] .= false
    median_intensities = zeros(T, n_intensity_bins)
    shape_estimates = Vector{Float32}(undef, n_intensity_bins)#zeros(T, n_intensity_bins)
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
        shape_estimates[bin_idx] = quantile(Laplace(0.0, L.θ), frag_err_quantile)# L.θ
    end


    intensities = 2 .^ median_intensities;

    intensities = hcat(log.(intensities), ones(length(intensities)))
    m, b = intensities\log.(shape_estimates)

    p = Plots.plot(exp.(intensities[:,1]), 
                    shape_estimates,
                    seriestype=:scatter,
                    title = out_fname,
                    xlabel = "Median intensity in Bin",
                    ylabel = "$frag_err_quantile quantile of laplace \n distributed mass errors",
                    label = nothing)
    bins = LinRange(0, 
                        maximum(exp.(intensities[:,1])), 
                        1000
                    )
    Plots.plot!(p, 
    bins,
    [exp(b)*(x^m) for x in bins]
    )

    savefig(p, joinpath(out_fdir, out_fname)*".pdf")

    return MassErrorModel(Float32(m), Float32(exp(b)), Float32(μ))
end