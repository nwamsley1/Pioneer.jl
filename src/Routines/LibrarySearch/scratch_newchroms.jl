function integratePrecursorMS2(chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, gauss_quad_x::Vector{Float64}, gauss_quad_w::Vector{Float64}; intensity_filter_fraction::Float32 = 0.01f0, α::Float32 = 0.01f0, half_width_at_α::Float32 = 0.15f0, LsqFit_tol::Float64 = 1e-3, Lsq_max_iter::Int = 100, tail_distance::Float32 = 0.25f0, isplot::Bool = false)
    
    function getBestPSM(filter::BitVector, weights::AbstractVector{<:AbstractFloat}, total_ions::AbstractVector{<:Integer}, q_values::AbstractVector{<:AbstractFloat})
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
   
    T = eltype(chrom.weight)

    filter = falses(size(chrom)[1])
    setFilter!(filter, chrom.weight, chrom.scan_idx)
    best_scan = getBestPSM(filter, chrom.weight, chrom.total_ions, chrom.q_value)
    best_rt, height = chrom.RT[best_scan], chrom.weight[best_scan]
    filterOnRT!(filter, best_rt, chrom.RT)
    #Needs to be setable parametfer. 1% is currently hard-coded
    filterLowIntensity!(filter, height*intensity_filter_fraction, chrom.weight)
    #Filter out samples below intensity threshold. 
    intensity, rt, lsq_fit_weight = zeros(Float32, (size(chrom)[1] - sum(filter)) + 2),  zeros(Float32, (size(chrom)[1] - sum(filter)) + 2),  ones(Float32, (size(chrom)[1] - sum(filter)) + 2)
    #Offset needs to be a parameter/optional argument

    fillIntensityandRT!(intensity, rt, filter, chrom.weight, chrom.RT)
    rt[1], rt[end] = rt[2] - tail_distance, rt[end-1] + tail_distance
    intensity[1], intensity[end] = zero(Float32), zero(Float32)
    fillLsqFitWeights!(lsq_fit_weight, intensity)

    #Scan with the highest score 
    #Too few points to attempt integration
    if (size(chrom)[1] - sum(filter)) < 2
        return nothing
    end

    ########
    #Fit EGH Curve 
    EGH_FIT = nothing
    try
        EGH_FIT = LsqFit.curve_fit(EGH_inplace, 
                                    JEGH_inplace, 
                                    rt,
                                    intensity,
                                    lsq_fit_weight,
                                    getP0((α, 
                                            half_width_at_α, 
                                            half_width_at_α, 
                                            Float32(best_rt), 
                                            Float32(height)
                                        ));
                                    x_tol = LsqFit_tol,
                                    maxIter = Lsq_max_iter,
                                    inplace = true)
    catch
        return nothing
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
                    fillrange = [0.0 for x in 1:length(x)], 
                    alpha = 0.25, color = :grey, show = true
                    ); 
    
        Plots.vline!([best_rt], color = :blue);
        Plots.vline!([rt[1]], color = :red);
        Plots.vline!([rt[end]], color = :red);
    end

    ############
    #Calculate Features
    
    MODEL_INTENSITY = EGH(rt, Tuple(EGH_FIT.param))
    
    peak_area = Integrate(EGH, gauss_quad_x, gauss_quad_w, Tuple(EGH_FIT.param))
   
    sum_of_residuals = zero(Float32)
    points_above_FWHM = zero(Int32)
    points_above_FWHM_01 = zero(Int32)
    for (i, I) in enumerate(intensity)
        if I  > (EGH_FIT.param[end]*0.5)
            points_above_FWHM += one(Int32)
            points_above_FWHM_01 += one(Int32)
        elseif I > (EGH_FIT.param[end]*0.01)
            points_above_FWHM_01 += one(Int32)
        end 
        sum_of_residuals += abs.(I - MODEL_INTENSITY[i])
    end

    GOF = 1 - sum_of_residuals/sum(intensity)
    FWHM = getFWHM(0.5, EGH_FIT.param[3], EGH_FIT.param[1])
    FWHM_01 = getFWHM(0.01, EGH_FIT.param[3], EGH_FIT.param[1])
    #asymmetry = atan(abs(EGH_FIT.param[3])/sqrt(abs(EGH_FIT.param[1]/2)))



    ############
    #Model Parameters
    σ = EGH_FIT.param[1]
    tᵣ = EGH_FIT.param[2]
    τ = EGH_FIT.param[3]
    H = EGH_FIT.param[4]
    
    ############
    #Summary Statistics 
    log_sum_of_weights = log2(sum(chrom.weight))
    mean_log_spectral_contrast = mean(log2.(chrom.spectral_contrast))
    mean_log_entropy = mean((max.(log2.(chrom.entropy_sim), Float16(1e-3))))
    mean_scribe_score = mean(chrom.scribe_score)
    mean_log_probability = mean(log2.(chrom.prob))
    ions_sum = sum(chrom.total_ions)
    data_points = Int64(sum(lsq_fit_weight) - 2) #May need to change
    mean_ratio = mean(chrom.matched_ratio)
    base_width = rt[end] - rt[1]
    
    chrom.peak_area[best_scan] = peak_area
    chrom.GOF[best_scan] = GOF
    chrom.FWHM[best_scan] = FWHM
    chrom.FWHM_01[best_scan] = FWHM_01
    #chrom.asymmetry[best_scan] = asymmetry
    chrom.points_above_FWHM[best_scan] = points_above_FWHM
    chrom.points_above_FWHM_01[best_scan] = points_above_FWHM_01
    chrom.σ[best_scan] = σ
    chrom.tᵣ[best_scan] = tᵣ
    chrom.τ[best_scan] = τ
    chrom.H[best_scan] = H
    chrom.log_sum_of_weights[best_scan] = log_sum_of_weights
    chrom.mean_log_spectral_contrast[best_scan] = mean_log_spectral_contrast
    chrom.mean_log_entropy[best_scan] = mean_log_entropy
    chrom.mean_scribe_score[best_scan] = mean_scribe_score
    chrom.mean_log_probability[best_scan] = mean_log_probability
    chrom.ions_sum[best_scan] = ions_sum
    chrom.data_points[best_scan] = data_points
    chrom.mean_matched_ratio[best_scan] = mean_ratio
    chrom.base_width_min[best_scan] = base_width
    chrom.best_scan[best_scan] = true
    

    return nothing
end

function integratePrecursors(grouped_precursor_df::GroupedDataFrame{DataFrame}; n_quadrature_nodes::Int64 = 100, intensity_filter_fraction::Float32 = 0.01f0, α::Float32 = 0.01f0, half_width_at_α::Float32 = 0.15f0, LsqFit_tol::Float64 = 1e-3, Lsq_max_iter::Int = 100, tail_distance::Float32 = 0.25f0, isplot::Bool = false)

    gx, gw = gausslegendre(n_quadrature_nodes)

    for i in ProgressBar(range(1, length(grouped_precursor_df)))
        integratePrecursorMS2(grouped_precursor_df[i]::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}},
                                gx::Vector{Float64},
                                gw::Vector{Float64},
                                intensity_filter_fraction = intensity_filter_fraction,
                                α = α,
                                half_width_at_α = half_width_at_α,
                                LsqFit_tol = LsqFit_tol,
                                Lsq_max_iter = Lsq_max_iter,
                                tail_distance = tail_distance,
                                isplot = isplot
                                )
    end

    ############
    #Clean

end
#=
function integratePrecursorMS2(chroms::GroupedDataFrame{DataFrame}, chroms_keys::Set{UInt32}, gauss_weights::Vector{Float64}, gauss_x::Vector{Float64}, precursor_idx::UInt32; max_smoothing_window::Int = 15, min_smoothing_order::Int = 3, min_scans::Int = 5, min_width::AbstractFloat = 1.0/6.0, integration_width::AbstractFloat = 4.0, integration_points::Int = 1000, isplot::Bool = false) where {T<:AbstractFloat}
   
    if (precursor_idx ∉ chroms_keys) #If the precursor is not found
        return (missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing)
    end

    #Chromatogram for the precursor. 
    #Has columns "weight" and "rt". 
    #chrom = chroms[(precursor_idx=precursor_idx,)]
    return integratePrecursorMS2(chroms[(precursor_idx=precursor_idx,)], gauss_weights, gauss_x, max_smoothing_window = max_smoothing_window, min_smoothing_order = min_smoothing_order, min_scans = min_scans, min_width = min_width, integration_width = integration_width, integration_points = integration_points, isplot = isplot)
end
=#