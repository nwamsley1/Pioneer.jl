function integratePrecursorMS2(chroms::GroupedDataFrame{DataFrame}, precursor_idx::UInt32, p0::NTuple{5, T}; max_smoothing_window::Int = 15, min_smoothing_order::Int = 3, min_scans::Int = 5, min_width::AbstractFloat = 1.0/6.0, integration_width::AbstractFloat = 4.0, integration_points::Int = 1000, isplot::Bool = false) where {T<:AbstractFloat}

    if !((precursor_idx=precursor_idx,) in keys(chroms)) #If the precursor is not found
        return (missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing)
    end

    #Chromatogram for the precursor. 
    #Has columns "weight" and "rt". 
    chrom = chroms[(precursor_idx=precursor_idx,)]
    return integratePrecursorMS2(chrom, p0, eltype(chrom[:,:weight]), max_smoothing_window = max_smoothing_window, min_smoothing_order = min_smoothing_order, min_scans = min_scans, min_width = min_width, integration_width = integration_width, integration_points = integration_points, isplot = isplot)
end

function integratePrecursorMS2(chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, p0::NTuple{5, U}, T::Type; max_smoothing_window::Int = 15, min_smoothing_order::Int = 3, min_scans::Int = 5, min_width::AbstractFloat = 1.0/6.0, integration_width::AbstractFloat = 4.0, integration_points::Int = 1000, isplot::Bool = false) where {U<:AbstractFloat}
    
    function getBestPSM(sdf::DataFrame)
        if any(sdf.q_value.<=0.01)
            #return argmax((sdf.q_value.<=0.01).*(sdf.hyperscore))
            return argmax((sdf.q_value.<=0.01).*(sdf.weight).*(sdf.total_ions))
        else
            return argmax(sdf.weight.*(sdf.total_ions))
        end
    end
    #display(chrom)
    non_zero = BitVector(undef, size(chrom)[1])
    fillNonZero!(non_zero, chrom[:,:weight])
    #Same precursor may be isolated multiple times within a single cycle_idx
    #Retain only most abundant within each cycle .
    for i in range(1, size(chrom)[1]-1)
        if abs(chrom[i,:scan_idx]-chrom[i+1,:scan_idx]) == 1
            #May not be appropriate criterion to choose between competing scans
            if chrom[i,:weight] > chrom[i + 1,:weight]
                non_zero[i + 1] = false
            else
                non_zero[i] = false
            end
        end
    end
    display(chrom)
    chrom = chrom[non_zero,:]
    
    best_scan = getBestPSM(chrom)
    best_rt = chrom[best_scan,:RT]
    height = chrom[best_scan,:weight]
    #println(best_rt)
    chrom = chrom[(chrom[:,:RT].>(best_rt - 0.5)).&((chrom[:,:RT]).<(best_rt + 0.5)),:]
    display(chrom)
    #return chrom
    #display(chrom)
    #Scan with the highest score 
    #Too few points to attempt integration
    if size(chrom)[1] < 2
        return (missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing)
    end
    #########
    #Probably don't need to allocate all of this memory. 
    min_intensity = maximum(chrom[:,:weight])/100.0
    rt = Float32.(chrom[chrom[:,:weight].>=min_intensity,:RT])
    pushfirst!(rt, Float32(chrom[1,:RT] - 0.25))
    push!(rt, Float32(chrom[end,:RT] + 0.25))
    intensity = collect(chrom[chrom[:,:weight].>=min_intensity,:weight])
    pushfirst!(intensity, Float32(0.0))
    push!(intensity, Float32(0.0))
    #Use smoothed first derivative crossings to identify peak apex and left/right boundaries
    start, stop = 1, length(intensity)# getIntegrationBounds(Int64.(frag_counts), rt, intensity, window_size, order, min_width, min_scans)

    ########
    #Fit EGH Curve 
    EGH_FIT = nothing
    mask = ones(Bool, length(intensity))
    #try
       
        for i in range(1, length(intensity) - 2)
            if (intensity[i + 1] < intensity[i]) & (intensity[i + 1] < intensity[i + 2])
                mask[i + 1] = 0
                #intensity[i + 1] = (intensity[i] + intensity[i + 2])/2
            end
            if (intensity[i + 1] < intensity[i]) & (intensity[i + 1] < intensity[min(i + 3, length(intensity))])
                #a = (intensity[min(i + 3, length(intensity))])/2
                #intensity[i + 1] = a 
                #intensity[ i + 2] = a
                mask[i + 1] = 0
                if (intensity[i + 2] < intensity[min(i + 3, length(intensity))])
                    mask[i + 2] = 0
                end
            end
        end
        #println(mask)
        #mask = [true, true, false, true, false, false, true]
        p0 = (0.1f0, 0.15f0, 0.15f0, Float32(best_rt), Float32(height))
        smooth = intensity
        EGH_FIT = LsqFit.curve_fit(EGH_inplace, 
                                    JEGH_inplace, 
                                    rt[mask],
                                    #intensity[start:stop],
                                    Float32.(smooth[mask]),
                                    #Float32.(smooth),
                                    #w,
                                    #(chrom[start:stop, :rank] .+ 1)./4,
                                    getP0(p0);
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
    #catch
    #    return (missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing)
    #end
    peak_area = Integrate(EGH, Tuple(EGH_FIT.param), n = integration_points)
    ##########
    #Plots                                           
    if isplot
        #Plot Data
        Plots.plot(rt, intensity, show = true, seriestype=:scatter)
        Plots.plot!(rt[mask], intensity[mask], show = true, seriestype=:scatter)
        Plots.plot!(chrom[:,:RT], chrom[:,:weight], show = true, alpha = 0.5, seriestype=:scatter)
        #Plots.plot!(rt, savitzky_golay(intensity, 7, 3).y, show = true, seriestype=:scatter)
        #Plots.plot!(rt, filter, show = true, seriestype=:scatter)
        #Plots Fitted EGH Curve
        #X = LinRange(T(EGH_FIT.param[2] - 1.0), T(EGH_FIT.param[2] + 1.0), integration_points)
        X = LinRange(T(best_rt - 1.0), T(best_rt + 1.0), integration_points)
        Plots.plot!(X,  
                    EGH(T.(collect(X)), Tuple(EGH_FIT.param)), 
                    fillrange = [0.0 for x in 1:integration_points], 
                    alpha = 0.25, color = :grey, show = true
                    ); 
        
        #=
                    Plots.plot!(X,  
                    LORENZ(T.(collect(X)), Tuple(LORENZ_FIT.param)), 
                    fillrange = [0.0 for x in 1:integration_points], 
                    alpha = 0.25, color = :blue, show = true
                    ); 
        =#
        #println(Tuple(GAUSS_FIT.param))
        #Plots.plot!(X,  
        #            GAUSS(T.(collect(X)), Tuple(GAUSS_FIT.param)), 
        #            fillrange = [0.0 for x in 1:integration_points], 
        #            alpha = 0.25, color = :red, show = true
        #            ); 
        Plots.vline!([best_rt], color = :blue);
        Plots.vline!([rt[start]], color = :red);
        Plots.vline!([rt[stop]], color = :red);
    end

    ############
    #Calculate Features
    MODEL_INTENSITY = EGH(rt, Tuple(EGH_FIT.param))
    RESID = abs.(intensity .- MODEL_INTENSITY)
    GOF_IDX = (MODEL_INTENSITY.>(EGH_FIT.param[end]*0.1)) .& (intensity.>(EGH_FIT.param[end]*0.1))
    GOF_IDX = (MODEL_INTENSITY.>(EGH_FIT.param[end]*0.0)) .& (intensity.>(EGH_FIT.param[end]*0.0))

    var = std(RESID)/EGH_FIT.param[4]
    GOF = 1 - sum(RESID[GOF_IDX])/sum(intensity[GOF_IDX])
    FWHM = getFWHM(0.5, EGH_FIT.param[3], EGH_FIT.param[1])
    FWHM_01 = getFWHM(0.01, EGH_FIT.param[3], EGH_FIT.param[1])
    asymmetry = atan(abs(EGH_FIT.param[3])/sqrt(abs(EGH_FIT.param[1]/2)))
    points_above_FWHM = sum(intensity[start:stop].>(EGH_FIT.param[end]*0.5))
    points_above_FWHM_01 = sum(intensity[start:stop].>(EGH_FIT.param[end]*0.01))

    σ = EGH_FIT.param[1]
    tᵣ = EGH_FIT.param[2]
    τ = EGH_FIT.param[3]
    H = EGH_FIT.param[4]
    return peak_area, GOF, FWHM, FWHM_01, asymmetry, points_above_FWHM, points_above_FWHM_01, σ, tᵣ, τ, H, sum(chrom[:,:weight]), sum(log2.(chrom[:,:spectral_contrast])), sum(chrom[:,:entropy_sim]), sum(log2.(chrom[:,:prob])), sum(chrom[:,:total_ions]), size(chrom)[1], mean(chrom[:,:matched_ratio]), rt[end] - rt[1]
end