function integratePrecursorMS2(chroms::GroupedDataFrame{DataFrame}, gauss_weights::Vector{Float64}, gauss_x::Vector{Float64}, precursor_idx::UInt32, p0::NTuple{5, T}; max_smoothing_window::Int = 15, min_smoothing_order::Int = 3, min_scans::Int = 5, min_width::AbstractFloat = 1.0/6.0, integration_width::AbstractFloat = 4.0, integration_points::Int = 1000, isplot::Bool = false) where {T<:AbstractFloat}

    if !((precursor_idx=precursor_idx,) in keys(chroms)) #If the precursor is not found
        return (missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing)
    end

    #Chromatogram for the precursor. 
    #Has columns "weight" and "rt". 
    chrom = chroms[(precursor_idx=precursor_idx,)]
    return integratePrecursorMS2(chrom, gauss_weights, gauss_x, p0, eltype(chrom[:,:weight]), max_smoothing_window = max_smoothing_window, min_smoothing_order = min_smoothing_order, min_scans = min_scans, min_width = min_width, integration_width = integration_width, integration_points = integration_points, isplot = isplot)
end

function integratePrecursorMS2(chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, x::Vector{Float64}, w::Vector{Float64}, p0::NTuple{5, U}, T::Type; max_smoothing_window::Int = 15, min_smoothing_order::Int = 3, min_scans::Int = 5, min_width::AbstractFloat = 1.0/6.0, integration_width::AbstractFloat = 4.0, integration_points::Int = 1000, isplot::Bool = false) where {U<:AbstractFloat}
    
    function getBestPSM(filter::BitVector, weights::AbstractVector{<:AbstractFloat}, total_ions::AbstractVector{UInt32}, q_values::AbstractVector{<:AbstractFloat})
        best_scan_idx = 1
        best_scan_score = zero(Float32)
        has_reached_fdr = false
        for i in range(1, length(weights))

            #Could be a better hueristic?
            score = weights[i]*total_ions[i]
            #Don't consider because deconvolusion set weights to zero
            #Or because it is filtered out 
            if iszero(score) | filter[i]
                continue
            end
            #If passing a threshold based on the logistic regression model
            #Automatically give priority. Could choose a better heuristic.
            #Possibly take into account retention time accuracy
            if q_values[i] <= 0.01
                #println("score $score, i $i, best_scan_score $best_scan_score ")
                #If this is the first scan reached to pass the 
                #logistic regression model threshold
                #then reset the best score to zero
                if has_reached_fdr == false
                    best_scan_score = zero(Float32)
                    has_reached_fdr = true
                end

                #Is best score? If yes then set new best score and scan index
                if score > best_scan_score
                    best_scan_score = score
                    best_scan_idx = i
                end

            elseif has_reached_fdr == false

                #Is best score? If yes then set new best score and scan index
                if score > best_scan_score
                    best_scan_score = score
                    best_scan_idx = i
                end

            end
        end
        #println("best_scan_idx $best_scan_idx")
        return best_scan_idx
    end

    #Same precursor may be isolated multiple times within a single cycle_idx
    #Retain only most abundant within each cycle .
    function setFilter!(filter::BitVector, weights::AbstractVector{T}, scan_idxs::AbstractVector{I}) where {T<:AbstractFloat, I<:Integer}
        for i in range(1, length(weights)-1)
            if abs(scan_idxs[i]-scan_idxs[i + 1]) == 1
                #May not be appropriate criterion to choose between competing scans
                if weights[i] > weights[i + 1]
                    filter[i + 1] = true
                else
                    filter[i] = true
                end
            end
        end
    end

    function filterLowIntensity!(filter::BitVector, min_intensity::T, weights::AbstractVector{T}) where {T<:AbstractFloat}
        for i in range(1, length(weights))
            if weights[i] <= min_intensity
                filter[i] = true
            end
        end
    end

    function fillIntensityandRT!(intensity::Vector{<:AbstractFloat}, rt::Vector{<:AbstractFloat},filter::BitVector, chrom_weight::AbstractVector{<:AbstractFloat}, chrom_rt::AbstractVector{<:AbstractFloat})
        n = 1
        for i in range(1, length(chrom_rt))
            if !filter[i]
                rt[n + 1] = chrom_rt[i]
                intensity[n + 1] = chrom_weight[i]
                n += 1
            end
        end 
    end

    function fillLsqFitWeights!(lsq_fit_weight::Vector{T}, intensity::Vector{T}) where {T<:AbstractFloat}
        for i in range(1, length(intensity) - 2)
            if (intensity[i + 1] < intensity[i]) & (intensity[i + 1] < intensity[i + 2])
                lsq_fit_weight[i + 1] = zero(Float32)
            end
            if (intensity[i + 1] < intensity[i]) & (intensity[i + 1] < intensity[min(i + 3, length(intensity))])
                lsq_fit_weight[i + 1] = zero(Float32)
                if (intensity[i + 2] < intensity[i]) & (intensity[i + 2] < intensity[min(i + 3, length(intensity))])#intensity[min(i + 1, length(intensity))])
                    lsq_fit_weight[i + 2] = zero(Float32)
                end
            end
        end
    end

    function filterOnRT!(filter::BitVector, best_rt::T, rts::AbstractVector{T}) where {T<:AbstractFloat}
        for i in eachindex(rts)
            if (rts[i] > (best_rt - 1.0)) & (rts[i] < (best_rt + 1.0))
                continue
            else
                filter[i] = true
            end
        end
    end
    #display(chrom)
    #@with chrom begin 
        filter = falses(size(chrom)[1])
        setFilter!(filter, chrom.weight, chrom.scan_idx)
        best_scan = getBestPSM(filter, chrom.weight, chrom.total_ions, chrom.q_value)
        best_rt, height = chrom.RT[best_scan], chrom.weight[best_scan]
        filterOnRT!(filter, best_rt, chrom.RT)
        #Needs to be setable parametfer. 1% is currently hard-coded
        filterLowIntensity!(filter, height/Float32(100), chrom.weight)
    #end
    #Filter out samples below intensity threshold. 
    intensity, rt, lsq_fit_weight = zeros(Float32, (size(chrom)[1] - sum(filter)) + 2),  zeros(Float32, (size(chrom)[1] - sum(filter)) + 2),  ones(Float32, (size(chrom)[1] - sum(filter)) + 2)
    #Offset needs to be a parameter/optional argument

    fillIntensityandRT!(intensity, rt, filter, chrom.weight, chrom.RT)
    rt[1], rt[end] = rt[2] - 0.25f0, rt[end-1] + 0.25f0
    intensity[1], intensity[end] = zero(Float32), zero(Float32)
    fillLsqFitWeights!(lsq_fit_weight, intensity)

    #=
    right_w = rt[end - 1] - rt[argmax(intensity)]
    left_w = rt[argmax(intensity)] - rt[2]
    w = 0.2
    if left_w < w
        if right_w > w
            rt[1] = rt[argmax(intensity)] - right_w
            intensity[1] = intensity[end - 1]
        else
            rt[1] = rt[2] - 0.25f0
        end
    else
        lsq_fit_weight[1] = 0.0f0
    end

    if right_w < w
        if left_w > w
            rt[end] = rt[argmax(intensity)] + left_w
            intensity[end - 1] = intensity[2]
        else
            rt[end] = rt[end - 1] + 0.25f0
        end
    else
        lsq_fit_weight[end] = 0.0f0
    end
    =#
    #=
    rt = rt[2:end - 1]
    intensity = intensity[2:end - 1]
    lsq_fit_weight = lsq_fit_weight[2:end - 1]
    
    println("intensity $intensity")
    println("rt $rt")
    println("lsq_fit_weight $lsq_fit_weight")
    =#
    #Scan with the highest score 
    #Too few points to attempt integration
    if (size(chrom)[1] - sum(filter)) < 2
        return (missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing)
    end
    #########
    #Probably don't need to allocate all of this memory. 
    #=
    min_intensity = maximum(chrom[:,:weight])/100.0
    rt = Float32.(chrom[chrom[:,:weight].>=min_intensity,:RT])
    pushfirst!(rt, Float32(chrom[1,:RT] - 0.25))
    push!(rt, Float32(chrom[end,:RT] + 0.25))
    intensity = collect(chrom[chrom[:,:weight].>=min_intensity,:weight])
    pushfirst!(intensity, Float32(0.0))
    push!(intensity, Float32(0.0))
    #Use smoothed first derivative crossings to identify peak apex and left/right boundaries
    start, stop = 1, length(intensity)# getIntegrationBounds(Int64.(frag_counts), rt, intensity, window_size, order, min_width, min_scans)
    =#
    ########
    #Fit EGH Curve 
    EGH_FIT = nothing
    try
        EGH_FIT = LsqFit.curve_fit(EGH_inplace, 
                                    JEGH_inplace, 
                                    rt,
                                    intensity,
                                    lsq_fit_weight,
                                    getP0((0.1f0, 
                                            0.15f0, 
                                            0.15f0, 
                                            Float32(best_rt), 
                                            Float32(height)
                                        ));
                                    x_tol = 1e-3,
                                    maxIter = 100,
                                    inplace = true)
        #=
        LORENZ_FIT = LsqFit.curve_fit(LORENZ_inplace, 
                                    #JGAUSS_inplace, 
                                    rt[mask],
                                    #intensity[start:stop],
                                    Float32.(smooth[mask]),
                                    #Float32.(smooth),
                                    #w,
                                    #(chrom[start:stop, :rank] .+ 1)./4,
                                    #[Float32(height)/Float32(2*π), Float32(best_rt), 0.2f0];
                                    [Float32(height)*Float32(sqrt(2*π)), Float32(best_rt), 1.0f0];
                                    inplace = true,
                                    #show_trace = true,
                                    #lower = Float32[1e3, 0, 0.01],
                                    #upper = Float32[1e10, 1000, 1],
                                    autodiff=:finiteforward)
        println([Float32(height), Float32(best_rt), 0.2f0])
        println("GAUSS_FIT ", LORENZ_FIT.param)
        =#
    catch
        return (missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing)
    end

    ##########
    #Plots                                           
    if isplot
        #Plot Data
        Plots.plot(rt, intensity, show = true, seriestype=:scatter)
        Plots.plot!(rt[lsq_fit_weight.!=0.0], intensity[lsq_fit_weight.!=0.0], show = true, alpha = 0.5, color = :green, seriestype=:scatter)
        #Plots.plot!(chrom[:,:RT], chrom[:,:weight], show = true, alpha = 0.5, seriestype=:scatter)

        X = LinRange(T(best_rt - 1.0), T(best_rt + 1.0), integration_points)
        Plots.plot!(X,  
                    EGH(T.(collect(X)), Tuple(EGH_FIT.param)), 
                    fillrange = [0.0 for x in 1:integration_points], 
                    alpha = 0.25, color = :grey, show = true
                    ); 
    
        Plots.vline!([best_rt], color = :blue);
        Plots.vline!([rt[1]], color = :red);
        Plots.vline!([rt[end]], color = :red);
    end

    ############
    #Calculate Features
    MODEL_INTENSITY = EGH(rt, Tuple(EGH_FIT.param))

    peak_area = Integrate(EGH, x, w, Tuple(EGH_FIT.param), n = integration_points)
    #return (missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing)

    RESID = abs.(intensity .- MODEL_INTENSITY)
    GOF_IDX = (MODEL_INTENSITY.>(EGH_FIT.param[end]*0.1)) .& (intensity.>(EGH_FIT.param[end]*0.1))
    GOF_IDX = (MODEL_INTENSITY.>(EGH_FIT.param[end]*0.0)) .& (intensity.>(EGH_FIT.param[end]*0.0))

    GOF = 1 - sum(RESID[GOF_IDX])/sum(intensity[GOF_IDX])
    FWHM = getFWHM(0.5, EGH_FIT.param[3], EGH_FIT.param[1])
    FWHM_01 = getFWHM(0.01, EGH_FIT.param[3], EGH_FIT.param[1])
    asymmetry = atan(abs(EGH_FIT.param[3])/sqrt(abs(EGH_FIT.param[1]/2)))
    points_above_FWHM = sum(intensity.>(EGH_FIT.param[end]*0.5))
    points_above_FWHM_01 = sum(intensity.>(EGH_FIT.param[end]*0.01))

    ############
    #Model Parameters
    σ = EGH_FIT.param[1]
    tᵣ = EGH_FIT.param[2]
    τ = EGH_FIT.param[3]
    H = EGH_FIT.param[4]

    ############
    #Summary Statistics 
    sum_of_weights = sum(chrom[:,:weight])
    mean_spectral_contrast = mean(chrom[:,:spectral_contrast])
    entropy_sum = sum(chrom[:,:entropy_sim])
    mean_log_probability = mean(log2.(chrom[:,:prob]))
    ion_sum = sum(chrom[:,:total_ions])
    data_points = Int64(sum(lsq_fit_weight) - 2) #May need to change
    mean_ratio = mean(chrom[:,:matched_ratio])
    base_width = rt[end] - rt[1]

    return peak_area, GOF, FWHM, FWHM_01, asymmetry, points_above_FWHM, points_above_FWHM_01, σ, tᵣ, τ, H, sum_of_weights, mean_spectral_contrast, entropy_sum, mean_log_probability, ion_sum, data_points, mean_ratio, base_width
end