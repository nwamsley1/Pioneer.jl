#Estimate Mass Errors
function estimateErrorDistribution(errs::Vector{T}, err_model::Type{D}, frag_err::T, frag_tol::T, window::T; f_out::String = "./") where {T<:AbstractFloat,D<:Distribution}
    #Mixture model with initial parameter guesses 
    #mix_guess = MixtureModel[err_model(frag_err, frag_tol), Distributions.Uniform(-window, window), [0.9, 1-0.9]]
    mix_guess = ExpectationMaximization.MixtureModel(
                        [
                        err_model(frag_err, frag_tol), 
                        #TDist(1),
                        Distributions.Uniform(-window, window)], 
                        [0.9, 1-0.9]
                        )

    #Estimte mixture model parameters by EM algorithm
    mix_mle = fit_mle(mix_guess, errs; atol = 1e-5, robust = true, infos = false)

    println(Distributions.params(mix_mle))

    err_dist = err_model(Distributions.params(mix_mle)[1][1][1], Distributions.params(mix_mle)[1][1][2])
    background_dist = Distributions.Uniform(Distributions.params(mix_mle)[1][2][1], Distributions.params(mix_mle)[1][2][2])
    bernoulli = Distributions.Bernoulli(Distributions.params(mix_mle)[2][1])

    #Sample from fitted model 
    dist = zeros(Float64, 200000)#)length(errs))
    for i in eachindex(dist)
        w = rand(bernoulli)*1.0
        dist[i] = w*rand(err_dist) + (1 - w)*rand(background_dist)
    end

    bins = LinRange(quantile(errs, 0.01), quantile(errs, 0.99), 100)
    p = Plots.histogram(dist, normalize = :probability, bins = bins, alpha = 0.5, label = "Samples from Fitted Dist.",
                        xlabel = "Mass Error (ppm)", ylabel = "Probability",
                        size = 100*[13.3, 7.5],
                        fontsize = 24,
                        titlefontsize = 24,
                        legendfontsize = 18,
                        tickfontsize = 24,
                        guidefontsize = 24,
                        margin = 10Plots.mm, legend = :topleft,
                        dpi = 300)

    Plots.histogram!(p, (errs), normalize = :probability, bins = bins, alpha = 0.5, label = "Observed Mass Errors (ppm)")
    #Plots.vline!([median(frag_ppm_errs)])
    Plots.vline!(p, [Distributions.params(mix_mle)[1][1][1]], lw = 6.0, color = :black, label = "Estimated Mass Error (ppm)")
    Plots.savefig(p, f_out*"mass_error.pdf")
    
    return err_dist
end

#Estimate Mass Errors
μ = median(frag_ppm_errs)
errs = gdf[10][!,:ppm_errs]
n = length(errs)
b = mean(abs.(errs.- μ))
U = Distributions.Uniform(-40.0, 40.0)
L = Distributions.Laplace(μ, b)
γ = zeros(Bool, length(errs))
z = 0.5
for iter in range(1, 50)
    #Estimate latent variable
    logpdf = 0.0
    for (i, err) in enumerate(errs)
        logpdf += log(pdf(L, err)*z + (1 - z)*pdf(U, err))
    end
    println("iter $iter, z $z, b $b, logpdf $logpdf")

    L = Distributions.Laplace(μ, b)
    for (i, err) in enumerate(errs)
        if rand(Uniform(0, 1)) < z*pdf(L, err)/(z*pdf(L, err) + (1 - z)*pdf(U, err))
            γ[i] = true
        else
            γ[i] = false
        end
    end

    #Update params
    SAD = 0.0 #Sum of absolute deviations
    N = 0
    for (i, err) in enumerate(errs)
        SAD += γ[i]*abs(err - μ)
        N += γ[i]
    end
    b = SAD/N
    z = sum(γ)/length(γ)
end


function EstimateMixtureWithUniformNoise(errs::AbstractVector{T}, #data
                                         err_model::Type{D}, #Distribution to model error
                                         μ::T, #Fixed/known location parameter
                                         w::T, #Known uniform distribution width
                                         b0::T, #Initial Scale estimate
                                         z::T, #Initial mixture estimate
                                         γ::Vector{Bool}; #Latent state variable 
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

@time EstimateMixtureWithUniformNoise(
    errs,
    Laplace{Float64},
    μ,
    40.0,
    b,
    0.5,
    γ
)

samples = zeros(Float64, 100000)
for i in range(1, length(samples))
    if rand(Uniform(0, 1)) < z
        samples[i] = rand(L)
    else
        samples[i] = rand(U)
    end
end

histogram(samples, alpha = 0.5, bins = LinRange(-40, 40, 100), normalize=:probability)
histogram!(errs, alpha = 0.5, bins = LinRange(-40, 40, 100), normalize=:probability)

#Estimate Mass Errors
#=
μ = median(frag_ppm_errs)
errs = gdf[1][!,:ppm_errs]
n = length(errs)
b = mean(abs.(errs.- μ))
U = Distributions.Uniform(-40.0, 40.0)
L = Distributions.Laplace(μ, b)
γ = zeros(Bool, length(errs))
z = 0.5
for iter in range(1, 50)
    #Estimate latent variable
    logpdf = 0.0
    for (i, err) in enumerate(errs)
        logpdf += log(pdf(L, err)*z + (1 - z)*pdf(U, err))
    end
    println("iter $iter, z $z, b $b, logpdf $logpdf")

    L = Distributions.Laplace(μ, b)
    for (i, err) in enumerate(errs)
        if rand(Uniform(0, 1)) < z*pdf(L, err)/(z*pdf(L, err) + (1 - z)*pdf(U, err))
            γ[i] = true
        else
            γ[i] = false
        end
    end

    #Update params
    SAD = 0.0 #Sum of absolute deviations
    N = 0

    pairs = 0
    while pairs < n*(n + 1)/2

    end
    #for (i, err) in enumerate(errs)
    #    SAD += γ[i]*abs(err - μ)
    #    N += γ[i]
    #end
    b = SAD/N
    z = sum(γ)/length(γ)
end


samples = zeros(Float64, 100000)
for i in range(1, length(samples))
    if rand(Uniform(0, 1)) < z
        samples[i] = rand(L)
    else
        samples[i] = rand(U)
    end
end
=#
histogram(samples, alpha = 0.5, bins = LinRange(-40, 40, 100), normalize=:probability)
histogram!(errs, alpha = 0.5, bins = LinRange(-40, 40, 100), normalize=:probability)
