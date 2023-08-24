#Estimate Mass Errors
function estimateErrorDistribution(errs::Vector{T}, err_model::Type{D}, frag_err::T, frag_tol::T, window::T; f_out::String = "./") where {T<:AbstractFloat,D<:Distribution}
    #Mixture model with initial parameter guesses 
    #mix_guess = MixtureModel[err_model(frag_err, frag_tol), Distributions.Uniform(-window, window), [0.9, 1-0.9]]
    mix_guess = ExpectationMaximization.MixtureModel(
                        [err_model(frag_err, frag_tol), 
                        Distributions.Uniform(-window, window)], 
                        [0.9, 1-0.9]
                        )

    #Estimte mixture model parameters by EM algorithm
    mix_mle = fit_mle(mix_guess, errs; atol = 1e-5, robust = true, infos = false)

    err_dist = err_model(Distributions.params(mix_mle)[1][1][1], Distributions.params(mix_mle)[1][1][2])
    background_dist = Distributions.Uniform(Distributions.params(mix_mle)[1][2][1], Distributions.params(mix_mle)[1][2][2])
    bernoulli = Distributions.Bernoulli(Distributions.params(mix_mle)[2][1])

    #Sample from fitted model 
    dist = zeros(Float64, length(errs))
    for i in eachindex(dist)
        w = rand(bernoulli)*1.0
        dist[i] = w*rand(err_dist) + (1 - w)*rand(background_dist)
    end

    
    p = Plots.histogram(dist, normalize = :probability, bins = 100, alpha = 0.5, label = "Samples from Fitted Dist.",
                        xlabel = "Mass Error (ppm)", ylabel = "Probability",
                        size = 100*[13.3, 7.5],
                        fontsize = 24,
                        titlefontsize = 24,
                        legendfontsize = 18,
                        tickfontsize = 24,
                        guidefontsize = 24,
                        margin = 10Plots.mm, legend = :topleft,
                        dpi = 300)

    Plots.histogram!(p, (errs), normalize = :probability, bins = 100, alpha = 0.5, label = "Observed Mass Errors (ppm)")
    #Plots.vline!([median(frag_ppm_errs)])
    Plots.vline!(p, [Distributions.params(mix_mle)[1][1][1]], lw = 6.0, color = :black, label = "Estimated Mass Error (ppm)")
    Plots.savefig(p, f_out*"mass_error.pdf")
    
    return err_dist
end