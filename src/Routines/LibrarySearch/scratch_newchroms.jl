function integratePrecursorMS2(chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, state::GD_state{HuberParams{U}, V, I, J}, gauss_quad_x::Vector{Float64}, gauss_quad_w::Vector{Float64}; intensity_filter_fraction::Float32 = 0.1f0, α::Float32 = 0.01f0, half_width_at_α::Float32 = 0.15f0, LsqFit_tol::Float64 = 1e-3, Lsq_max_iter::Int = 100, tail_distance::Float32 = 0.05f0, isplot::Bool = false) where {U,V<:AbstractFloat, I,J<:Integer}
    
    function getBestPSM(filter::BitVector,  hyperscore::AbstractVector{<:AbstractFloat}, weights::AbstractVector{<:AbstractFloat}, total_ions::AbstractVector{<:Integer}, q_values::AbstractVector{<:AbstractFloat}, RT_error::AbstractVector{<:AbstractFloat})
        best_scan_idx = 1
        best_scan_score = zero(Float32)
        has_reached_fdr = false
        for i in range(1, length(weights))

            #Could be a better hueristic?
            score = #hyperscore[i]*sqrt(weights[i])#weights[i]*total_ions[i]/RT_error[i]
            #Don't consider because deconvolusion set weights to zero
            #Or because it is filtered out 
            if iszero(weights[i]) #| filter[i]
                continue
            end
            #If passing a threshold based on the logistic regression model
            #Automatically give priority. Could choose a better heuristic.
            #Possibly take into account retention time accuracy
            if q_values[i] <= 0.1
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
    function setFilter!(state::GD_state{HuberParams{T}, U, I, J}, weights::AbstractVector{T}, scan_idxs::AbstractVector{I}) where {T,U<:AbstractFloat, I,J<:Integer}
        for i in range(1, length(weights)-1)
            if abs(scan_idxs[i]-scan_idxs[i + 1]) == 1
                #May not be appropriate criterion to choose between competing scans
                if weights[i] > weights[i + 1]
                    if state.mask[i] == false
                        state.mask[i + 1] = true
                    end
                else
                    state.mask[i] = true
                end
            end
        end
    end

    function filterLowIntensity!(state::GD_state{HuberParams{T}, U, I, J}, min_intensity::T, weights::AbstractVector{T}) where {T,U<:AbstractFloat, I,J<:Integer}
        for i in range(1, length(weights))
            if weights[i] <= min_intensity
                state.mask[i] = true
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

    function fillLsqFitWeights!(state::GD_state{HuberParams{T}, U, I, J}) where {T,U<:AbstractFloat, I,J<:Integer}

        intensity = state.data
        for i in range(1, length(intensity) - 2)
            if (intensity[i + 1] < intensity[i]) & (intensity[i + 1] < intensity[i + 2])
                state.mask[i + 1] = true
            end
            if (intensity[i + 1] < intensity[i]) & (intensity[i + 1] < intensity[min(i + 3, length(intensity))])
                state.mask[i + 1] = true
                if (intensity[i + 2] < intensity[i]) & (intensity[i + 2] < intensity[min(i + 3, length(intensity))])#intensity[min(i + 1, length(intensity))])
                    state.mask[i + 2] = true
                end
            end
        end
    end

    function filterOnRT!(state::GD_state{HuberParams{T}, U, I, J}, best_rt::T, rts::AbstractVector{T}) where {T,U<:AbstractFloat, I,J<:Integer}
        for i in eachindex(rts)
            if (rts[i] > (best_rt - 1.0)) & (rts[i] < (best_rt + 1.0))
                continue
            else
                state.mask[i] = true
            end
        end
    end
   
    function truncateAfterSkip!(state::GD_state{HuberParams{T}, U, I, J}, best_scan::Int64, rts::AbstractVector{<:AbstractFloat}) where {T,U<:AbstractFloat, I,J<:Integer}
        
        
        for i in range(best_scan, length(rts)-1)
            if (rts[i+ 1] - rts[i]) > 0.3
                for n in range(i + 1, length(rts))
                    state.mask[n] = true
                end
            end
        end
        
        for i in range(1, best_scan - 1)
            if (rts[best_scan - i + 1] - rts[best_scan - i]) > 0.3
                for n in range(1, best_scan-i)
                    state.mask[n] = true
                end
            end
        end

        return 

    end

    function filterOnMatchedRatio!(state::GD_state{HuberParams{T}, U, I, J}, best_scan_idx::Int64, matched_ratios::AbstractVector{U}) where {T,U<:AbstractFloat, I,J<:Integer}
        for i in range(1, state.max_index)
            if matched_Ratios[i] < (matched_ratios[best_scan_idx] - 1)
                state.mask[i] = true
            end
        end
    end

    function fillState!(state::GD_state{HuberParams{T}, U, I, J}, chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}) where {T,U<:AbstractFloat, I,J<:Integer}
        for i in range(1, length(chrom.weight))
            state.t[i] = chrom.RT[i]
            state.data[i] = chrom.weight[i]
        end
        state.max_index = length(chrom.weight)
    end

    T = eltype(chrom.weight)

    #Initialize state with observed data 
    fillState!(state, chrom)
    #best_scan = getBestPSM(filter, chrom.matched_ratio, chrom.weight, chrom.total_ions, chrom.q_value, chrom.RT_error)
    best_scan = argmax(chrom.weight)
    truncateAfterSkip!(state, best_scan, chrom.RT) #Only allow a fixed ammount of time without a scan
    setFilter!(state, chrom.weight, chrom.scan_idx) #For adjacent scans, choose only the most abundant
    best_rt, best_height = chrom.RT[best_scan], chrom.weight[best_scan] #Height and rt of scan with largest coefficient
    filterOnRT!(state, best_rt, chrom.RT) #Remove scans greater than a certain distance from the best scan 
    #Maybe don't need this?
    #filterLowIntensity!(state, best_height*intensity_filter_fraction, chrom.weight) #Remove scans with intensity less than a fixed fraction of the most intense scan 
    filterOnMatchedRatio!(state, best_scan, chrom.matched_ratio)

    fillLsqFitWeights!(state)

    best_rt, best_height = 0.0, 0.0
    for i in range(1, state.max_index)
        if state.data[i]*state.mask[i] > best_height
            best_height = state.data[i]
            best_rt = state.t[i]
        end
    end

    data_points = Int64(state.max_index) #May need to change
    ########
    #Fit EGH Curve to data 

    #Parameter constraints
    lower, upper = HuberParams(0.001, 0, -1, 0.05), HuberParams(1, Inf, 1, 2.0);
    #Initial parameters
    state.params = getP0(α, 
                            half_width_at_α, 
                            half_width_at_α,
                            best_rt,
                            best_height,
                            lower, upper)
    #end

    #Estimate model parameters in `state` by gradint descent using Adam algorithm
    GD(state,
                lower,
                upper,
                tol = 5e-3, 
                max_iter = 100, 
                δ = 1e-2,
                α=0.1,
                β1 = 0.9,
                β2 = 0.999,
                ϵ = 1e-8)

    ##########
    #Plots                                           
    if isplot
        display(chrom)
        p = Plots.plot(state.t, #Plot observed data points 
                        state.data, 
                        show = true, seriestype=:scatter, reuse = false)

        Plots.plot!(p, #Plot fitted model curve
                    state.t[state.mask.==false], 
                    state.y[state.mask.==false], 
                    show = true, alpha = 0.5, color = :green, seriestype=:scatter)

        X = LinRange(T(best_rt - 1.0), T(best_rt + 1.0), length(gauss_quad_w))
        Plots.plot!(p, X,  
                    [F(state, x) for x in X], #Evaluate fitted model at each point. 
                    fillrange = [0.0 for x in 1:length(gauss_quad_w)], 
                    alpha = 0.25, color = :grey, show = true
                    ); 

        Plots.vline!(p, [best_rt], color = :blue);
        Plots.vline!(p,[state.t[1]], color = :red);
        Plots.vline!(p, [state.t[end][end]], color = :red);
        display(p)
        peak_area = Integrate(state, gauss_quad_x, gauss_quad_w)
        println("peak_area $peak_area")
    end

    ############
    #Calculate Features
    
    #MODEL_INTENSITY = EGH(rt, Tuple(EGH_FIT.param))
    peak_area = Integrate(state, gauss_quad_x, gauss_quad_w)
   
    sum_of_residuals = zero(Float32)
    points_above_FWHM = zero(Int32)
    points_above_FWHM_01 = zero(Int32)
    
    for (i, I) in enumerate(state.data)
        if I > (state.params.H*0.5)
            points_above_FWHM += one(Int32)
        elseif I > (state.params.H*0.01)
            points_above_FWHM_01 += one(Int32)
        end 
        sum_of_residuals += abs.(I - state.y[i])
    end
    
    GOF = 1 - sum_of_residuals/sum(intensity)
    FWHM = getFWHM(state, 0.5)
    FWHM_01 = getFWHM(state, 0.01)
    #asymmetry = atan(abs(EGH_FIT.param[3])/sqrt(abs(EGH_FIT.param[1]/2)))

    ############
    #Model Parameters
    σ = state.params.σ
    tᵣ = state.params.tᵣ
    τ = state.params.τ
    H = state.params.H
    
    ############
    #Summary Statistics 
    log_sum_of_weights = 0.0
    mean_scribe_score = 0.0
    max_weight = 0.0
    mean_log_entropy = -10.0
    mean_log_probability = 0.0
    mean_log_spectral_contrast = -10.0
    count = 0
    ions_sum = 0
    for i in range(1, length(state.mask))
        if !state.mask[i]
            if chrom.scribe[i]>mean_scribe_score
                mean_scribe_score = chrom.scribe[i]
            end
            if log2(max(chrom.entropy_score[i], 0.001))>mean_log_entropy
                mean_log_entropy=log2(max(chrom.entropy_score[i], 0.001))
            end
            if chrom.scribe[i]>mean_scribe_score
                mean_scribe_score = chrom.scribe[i]
            end
            if log2(chrom.spectral_contrast[i])>mean_log_spectral_contrast
                mean_log_spectral_contrast= log2(chrom.spectral_contrast[i])
            end
            if chrom.weight[i]>log_sum_of_weights 
                log_sum_of_weights  = chrom.weight[i]
                max_weight = chrom.weight[i]
            end
            if chrom.total_ions[i] > ions_sum
                ions_sum = chrom.total_ions[i]
            end

            mean_log_probability += log2(chrom.prob[i])

            count += 1
        end
    end    
    log_sum_of_weights = log2(log_sum_of_weights)
    mean_log_probability = mean_log_probability/count
   # mean_log_entropy = mean_log_entropy/count
   # mean_scribe_score = mean_scribe_score/count
   # mean_log_spectral_contrast = mean_log_spectral_contrast/count
    #mean_scribe_score = mean_scribe_score/count   

    #mean_scribe_score = sum(chrom.scribe_score)
    #mean_log_probability = mean(log2.(chrom.prob))
    #ions_sum = sum(chrom.total_ions)
    #data_points = Int64(sum(lsq_fit_weight) - 2) #May need to change
    mean_ratio = mean(chrom.matched_ratio)
    base_width = rt[end] - rt[1]

    best_scan = argmax(chrom.weight)#argmax(chrom.prob.*(chrom.weight.!=0.0))
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
    chrom.max_weight[best_scan] = max_weight
    

    return nothing
end

function integratePrecursors(grouped_precursor_df::GroupedDataFrame{DataFrame}; n_quadrature_nodes::Int64 = 100, intensity_filter_fraction::Float32 = 0.01f0, α::Float32 = 0.01f0, half_width_at_α::Float32 = 0.15f0, LsqFit_tol::Float64 = 1e-3, Lsq_max_iter::Int = 100, tail_distance::Float32 = 0.25f0, isplot::Bool = false)

    gx, gw = gausslegendre(n_quadrature_nodes)
    N = 500
    state = GD_state(
                    HuberParams(0.0, 0.0, 0.0), #Initial params
                    zeros(Float64, N), #t
                    zeros(Float64, N), #y
                    zeros(Float64, N), #data
                    falses(N), #mask
                    0, #number of iterations
                    N #max index
                    )
    #Threads.@threads for i in ProgressBar(range(1, length(grouped_precursor_df)))
    for i in range(1, length(grouped_precursor_df))
        integratePrecursorMS2(grouped_precursor_df[i]::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}},
                                state,
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
        reset!(state)
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