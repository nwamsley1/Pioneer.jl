
struct MassErrorModel{T<:AbstractFloat}
    mass_offset::T#UniformSpline{N, T}
    mass_tolerance::Tuple{T, T}#UniformSpline{N, T}
end

getRightTol(mem::MassErrorModel) = last(mem.mass_tolerance)
getLeftTol(mem::MassErrorModel) = first(mem.mass_tolerance)
getMassOffset(mem::MassErrorModel) = last(mem.mass_offset)

function getMassCorrection(mem::MassErrorModel)
    return mem.mass_offset
end

function getLocation(mem::MassErrorModel)
    return mem.location
end

function (mem::MassErrorModel)(mass::Float32, log_intensity::Float32, quantile::Float32)
    ppm_norm = Float32(1e6)
    ppm = mass/ppm_norm
    mass -= getMassOffset(mem)*ppm
    ppm = mass/ppm_norm
    r_tol = getRightTol(mem)*ppm
    l_tol = getLeftTol(mem)*ppm
    return Float32(mass - r_tol), Float32(mass - r_tol), Float32(mass + l_tol)
end

function ModelMassErrs(ppm_errs::Vector{Float32};
                       frag_err_quantile::Float32 = 0.01f0,
                       out_fdir::String = "./",
                       out_fname::String = "mass_err_estimate")

    bins = LinRange(minimum(ppm_errs), maximum(ppm_errs), 100)
    mass_err = median(ppm_errs)
    ppm_errs = ppm_errs .- mass_err
    r_bound = quantile(ppm_errs, 1 - frag_err_quantile) 
    l_bound = quantile(ppm_errs, frag_err_quantile)


    errs = ppm_errs .+ mass_err
    plot_title = ""
    n = 0
    for i in range(1, length(out_fname))
        n += 1
        if n > 24
            n = 1
            plot_title *= "\n"
        end
        plot_title *= out_fname[i]
    end
    n = length(errs)
    p = Plots.histogram(errs,
                    orientation = :h, 
                    yflip = true,
                    #seriestype=:scatter,
                    title = plot_title*"\n n = $n",
                    xlabel = "Count",
                    ylabel = "Mass Error (ppm)",
                    label = nothing,
                    bins = bins,
                    ylim = (minimum(ppm_errs)-2, maximum(ppm_errs)+2),
                    topmargin =15mm,
                    #bottommargin = 10mm,
                    )

    Plots.hline!([l_bound + mass_err, mass_err, r_bound + mass_err], label = nothing, color = :black, lw = 2)
    l_err = l_bound + mass_err
    Plots.annotate!(last(xlims(p)), l_err, text("$l_err", :black, :right, :bottom, 12))
    r_err = r_bound + mass_err
    Plots.annotate!(last(xlims(p)), r_err, text("$r_err", :black, :right, :bottom, 12))
    Plots.annotate!(last(xlims(p)), mass_err, text("$mass_err", :black, :right, :bottom, 12))
    savefig(p, joinpath(out_fdir, out_fname)*".pdf")

    MassErrorModel(
                   Float32(mass_err),
                (Float32(abs(l_bound)), Float32(abs(r_bound)))
                    )
end

#=
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

function getIntensityBinIds(log2_intensities::Vector{T},
                            max_n_bins::Int64,
                            min_bin_size::Int64)::Vector{UInt32} where {T<:AbstractFloat}

    #Width of intensity bins 
    bin_width = (maximum(log2_intensities) - minimum(log2_intensities))/max_n_bins
    
    #Assumes log2_intensities is sorted in descenging order 
    #Allocate bin ids and staring indices 
    N = length(log2_intensities)
    bins = zeros(UInt32, N)
    bin_idx = 1
    #Intensity of stoping and starting bins 
    start_intensity, stop_intensity = first(log2_intensities), first(log2_intensities)
    start_idx, stop_idx = 1, 1
    #Construct intensity bins 
    for i in range(1, N)
        stop_intensity = log2_intensities[i]
        stop_idx = i
        #If the bin exceeds the minimum width and number of data points,
        #then start a new bin 
        if (abs(stop_intensity - start_intensity) > bin_width) & ((stop_idx - start_idx) > min_bin_size)
            bin_idx += 1
            start_idx, stop_idx = i, i
            start_intensity = stop_intensity
        end
        bins[i] = bin_idx
    end

    #If last bin has fwewer than the minimum number of fragments, merge
    #it with the second to last bin 
    if abs(start_idx - stop_idx) < min_bin_size
        bins[start_idx:stop_idx] .= max(bin_idx - 1, 1)
    end
    return bins      
end

function ModelMassErrs(intensities::Vector{T},
                       ppm_errs::Vector{U},
                       frag_tol::U;
                       min_ppm::Float32 = zero(Float32),
                       max_ppm::Float32 = typemax(Float32),
                       max_n_bins::Int = 30,
                       min_bin_size::Int = 300,
                       frag_err_quantile::Float64 = 0.999,
                       out_fdir::String = "./",
                       out_fname = "mass_err_estimate") where {T,U<:AbstractFloat}

    #sort intensities/ppm_errs in increasing order of intensity 
    log2_intensities = log2.(intensities)
    new_indices = sortperm(log2_intensities, rev = true)
    log2_intensities = log2_intensities[new_indices]
    ppm_errs = ppm_errs[new_indices]

    #bins = getIntensityBinIds(log2_intensities, max_n_bins, min_bin_size)
    bins = getIntensityBinIds(log2_intensities, max_n_bins, min_bin_size)

    err_df = DataFrame(Dict(
            :ppm_errs => ppm_errs, 
            :log2_intensities => log2_intensities, 
            :bins => bins,
            :γ => zeros(Bool, length(ppm_errs)))
            )
    #return err_df
    #Pre-allocate outpues 
    median_intensities = zeros(T, bins[end])
    shape_estimates = Vector{Float32}(undef, bins[end])#zeros(T, n_intensity_bins)
    μ_estimates = Vector{Float32}(undef, bins[end])
    #Estimate Laplace-uniform mixture distribution for each intensity bin 
    bin_idx = 0
    for (int_bin, subdf) in pairs(groupby(err_df, :bins))
        bin_idx += 1 #Intensity bin counter
        median_intensities[bin_idx] = median(subdf[!,:log2_intensities])
        println("median intensity ", median_intensities[bin_idx])
        bin_μ = median(subdf[!,:ppm_errs])
        b = mean(abs.(subdf[!,:ppm_errs] .- bin_μ)) #Mean absolute deviation estimate
        L, z = EstimateMixtureWithUniformNoise(
            subdf[!,:ppm_errs],
            Laplace{Float64},
            bin_μ,
            frag_tol,
            b,
            0.5, #mixture estimate
            subdf[!,:γ] 
        )
        shape_estimates[bin_idx] = quantile(abs.(subdf[!,:ppm_errs].-median(subdf[!,:ppm_errs])), 0.98)#L.θ#3*mad(subdf[!,:ppm_errs], normalize = true)#quantile(abs.(subdf[!,:ppm_errs] .- bin_μ), 0.99)#quantile(Normal(0.0, L.σ), 0.99)#L.θ
        μ_estimates[bin_idx] = bin_μ
    end


    intensities = median_intensities;
    new_order = sortperm(intensities)
    intensities = intensities[new_order]
    shape_estimates = shape_estimates[new_order]
    μ_estimates = μ_estimates[new_order]
    #return intensities, shape_estimates
    #Build Splines 
    bins = LinRange(minimum(intensities), maximum(intensities), 1000)
    #=
    log_intensity_to_shape = BSplineApprox(shape_estimates,
                                            intensities, 3, 4, :Uniform, :Uniform, extrapolate = true)
    log_intensity_to_μ = BSplineApprox(μ_estimates,
                                            intensities, 3, 4, :Uniform, :Uniform, extrapolate = true)                                        
    =#
    log_intensity_to_shape = UniformSpline(shape_estimates, intensities, 3, 4)
    log_intensity_to_μ = UniformSpline(μ_estimates, intensities, 3, 4)

    p = Plots.plot((intensities), 
                    shape_estimates,#.*(-1.0*log(2.0f0*(1.0f0 - frag_err_quantile))),
                    #shape_estimates,
                    seriestype=:scatter,
                    title = out_fname,
                    xlabel = "Median intensity in Bin",
                    ylabel = "$frag_err_quantile quantile of laplace \n distributed mass errors",
                    label = nothing)


    Plots.plot!(p, 
   # bins, [log_intensity_to_shape(x)*.*(-1.0*log(2.0f0*(1.0f0 - frag_err_quantile))) for x in bins]
   bins, [log_intensity_to_shape(x) for x in bins]
    )

    savefig(p, joinpath(out_fdir, out_fname)*".pdf")

    p = Plots.plot((intensities), 
    μ_estimates,
    #shape_estimates,
    seriestype=:scatter,
    title = out_fname,
    xlabel = "Median intensity in Bin",
    ylabel = "mu estimate",
    label = nothing)

    Plots.plot!(p, 
    bins, [ log_intensity_to_μ(x) for x in bins]
    )

    savefig(p, joinpath(out_fdir, out_fname*"_mu")*".pdf")

    MassErrorModel(
                    log_intensity_to_μ,
                    log_intensity_to_shape,
                    min_ppm,
                    max_ppm
                    )
end
=#