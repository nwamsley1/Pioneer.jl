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
    function setFilter!(state::GD_state{HuberParams{T}, U, I, J}, weights::SubArray{Float32, 1, Vector{Float32}, Tuple{Vector{Int64}}, false}, scan_idxs::SubArray{UInt32, 1, Vector{UInt32}, Tuple{Vector{Int64}}, false}) where {T,U<:AbstractFloat, I,J<:Integer}
        for i in range(1, length(weights)-1)
            #println("scan_idxs[i] ", scan_idxs[i])
            #println("scan_idxs[i + 1] ", scan_idxs[i + 1])
            #println(" abs(scan_idxs[i]-scan_idxs[i + 1]) ", abs(scan_idxs[i]-scan_idxs[i + 1]))
            if abs(scan_idxs[i+1]-scan_idxs[i]) == 1
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

    function getDataPoints(state::GD_state{HuberParams{T}, U, I, J}) where {T,U<:AbstractFloat, I,J<:Integer}
        data_points = 0
        for i in range(1, state.max_index)
            if state.mask[i]==false
                data_points += 1
            end
        end
        return data_points
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

        #println("state.max_index ", state.max_index)            
        #println("state.mask[1:state.max_index] ", state.mask[1:state.max_index])
        for i in range(1, state.max_index - 2)

            if (intensity[i + 1] < intensity[i]) & (intensity[i + 1] < intensity[i + 2])
                state.mask[i + 1] = true
            end
            if (intensity[i + 1] < intensity[i]) & (intensity[i + 1] < intensity[min(i + 3, state.max_index)])
                #println("i $i")
                state.mask[i + 1] = true
                if (intensity[i + 2] < intensity[i]) & (intensity[i + 2] < intensity[min(i + 3, state.max_index)])#intensity[min(i + 1, length(intensity))])
                    state.mask[i + 2] = true
                end
            end
            #println("state.mask[1:state.max_index] ", state.mask[1:state.max_index])
        end
    end

    function filterOnRT!(state::GD_state{HuberParams{T}, U, I, J}, best_rt::T, rts::SubArray{Float32, 1, Vector{Float32}, Tuple{Vector{Int64}}, false}) where {T,U<:AbstractFloat, I,J<:Integer}
        for i in eachindex(rts)
            if (rts[i] > (best_rt - 1.0)) & (rts[i] < (best_rt + 1.0))
                continue
            else
                state.mask[i] = true
            end
        end
    end
   
    function truncateAfterSkip!(state::GD_state{HuberParams{T}, U, I, J}, best_scan::Int64, rts::SubArray{Float32, 1, Vector{Float32}, Tuple{Vector{Int64}}, false}) where {T,U<:AbstractFloat, I,J<:Integer}
        
        
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

    function filterOnMatchedRatio!(state::GD_state{HuberParams{T}, U, I, J}, best_scan_idx::Int64, matched_ratios::SubArray{V, 1, Vector{V}, Tuple{Vector{Int64}}, false}) where {T,U,V<:AbstractFloat, I,J<:Integer}
        for i in range(1, state.max_index)
            if (matched_ratios[i] < (matched_ratios[best_scan_idx] - 1)) & (matched_ratios[i] < 0)
                state.mask[i] = true
            end
            if matched_ratios[i] < -1
                state.mask[i] = true
            end
        end
    end

    function fillState!(state::GD_state{HuberParams{T}, U, I, J}, chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}) where {T,U<:AbstractFloat, I,J<:Integer}
        max_weight = maximum(chrom.weight)
        for i in range(1, length(chrom.weight))
            state.t[i] = chrom.RT[i]
            state.data[i] = chrom.weight[i]/max_weight
        end
        state.max_index = length(chrom.weight)
        return max_weight
    end

    T = eltype(chrom.weight)

    #Initialize state with observed data 
    norm_factor = fillState!(state, chrom)
    #best_scan = getBestPSM(filter, chrom.matched_ratio, chrom.weight, chrom.total_ions, chrom.q_value, chrom.RT_error)
    best_scan = argmax(chrom.weight)
    truncateAfterSkip!(state, best_scan, chrom.RT) #Only allow a fixed ammount of time without a scan
    setFilter!(state, chrom.weight, chrom.scan_idx) #For adjacent scans, choose only the most abundant
    best_rt, best_height = chrom.RT[best_scan], chrom.weight[best_scan] #Height and rt of scan with largest coefficient
    filterOnRT!(state, best_rt, chrom.RT) #Remove scans greater than a certain distance from the best scan 
    #Maybe don't need this?
    filterLowIntensity!(state, best_height*intensity_filter_fraction, chrom.weight) #Remove scans with intensity less than a fixed fraction of the most intense scan 
    filterOnMatchedRatio!(state, best_scan, chrom.matched_ratio) #Remove scans with matched ratios significantly lower than that of the best scan
    fillLsqFitWeights!(state) #Enforce smoothness 
    best_rt, best_height = 0.0, 0.0
    for i in range(1, state.max_index)
        if state.data[i]*(1 - state.mask[i]) >= best_height
            best_height = state.data[i]
            best_rt = state.t[i]
        end
    end

    data_points = zero(UInt32)#Int64(state.max_index) #May need to change
    for i in range(1, state.max_index)
        if state.mask[i]==false
            data_points += one(UInt32)
        end
    end

    #println("state.mask ", state.mask[1:state.max_index])
    ########
    #Fit EGH Curve to data 
    #println("data_points $data_points")
    #Parameter constraints
    #lower, upper = HuberParams(T(0.001), T(0), T(-1), T(0.05)), HuberParams(T(1), T(Inf), T(1), T(2.0));
    if data_points > 6
        lower, upper = HuberParams(T(0.001), T(0), T(-1), T(1.0)), HuberParams(T(1), T(Inf), T(1), T(1.0));
    else
        lower, upper = HuberParams(T(0.01), T(0), T(-0.01), T(1.0)), HuberParams(T(0.025), T(Inf), T(0.05), T(1.0));
    end
    #Initial parameters
    state.params = getP0(T(α), 
                            T(half_width_at_α), 
                            T(half_width_at_α),
                            T(best_rt),
                            T(best_height),
                            lower, upper)

    state.t[state.max_index + 1], state.t[state.max_index + 2] = T(best_rt - 0.5), T(best_rt + 0.5)
    state.data[state.max_index + 1], state.data[state.max_index + 2] = zero(T), zero(T)
    state.max_index += 2
    GD(state,
                lower,
                upper,
                tol = 1e-3, 
                max_iter = 200, 
                δ = 1e-3,
                #δ = 100.0,
                α=0.01,
                β1 = 0.9,
                β2 = 0.999,
                ϵ = 1e-8)
    #println("area first ", Integrate(state, gauss_quad_x, gauss_quad_w, α = α))
    ##########
    #Plots                                           
    if isplot
        #display(chrom)
        println("best_height $best_height")
        p = Plots.plot(state.t[1:state.max_index], #Plot observed data points 
                        state.data[1:state.max_index].*norm_factor, 
                        show = true, seriestype=:scatter, reuse = false)
        t = [state.t[i] for i in range(1, state.max_index) if state.mask[i] == false]
        d = [state.data[i] for i in range(1, state.max_index) if state.mask[i] == false]
        Plots.plot!(p, #Plot fitted model curve
                    t,d.*norm_factor,
                    show = true, alpha = 0.5, color = :green, seriestype=:scatter)

        
        X = LinRange(T(best_rt - 1.0), T(best_rt + 1.0), length(gauss_quad_w))
        Plots.plot!(p, X,  
                    [F(state, x) for x in X].*norm_factor, #Evaluate fitted model at each point. 
                    fillrange = [0.0 for x in 1:length(gauss_quad_w)], 
                    alpha = 0.25, color = :grey, show = true
                    ); 

        #lower, upper = HuberParams(T(0.01), T(0), T(-0.01), T(0.05)), HuberParams(T(0.025), T(Inf), T(0.05), T(2.0));
        lower, upper = HuberParams(T(0.01), T(0), T(-0.01), T(1.0)), HuberParams(T(0.025), T(Inf), T(0.05), T(1.0));
        #Initial parameters
        state.params = getP0(T(α), 
                                T(half_width_at_α), 
                                T(half_width_at_α),
                                T(best_rt),
                                T(best_height),
                                lower, upper)
       
        state.t[state.max_index + 1], state.t[state.max_index + 2] = T(best_rt - 0.5), T(best_rt + 0.5)
        state.data[state.max_index + 1], state.data[state.max_index + 2] = zero(T), zero(T)
        state.max_index += 2
        println("state.params.H ", state.params.H)
        GD(state,
                    lower,
                    upper,
                    tol = 1e-3, 
                    max_iter = 200, 
                    δ = 1e-3,
                    #δ = 100.0,
                    α=0.01,
                    β1 = 0.9,
                    β2 = 0.999,
                    ϵ = 1e-8)
        println("state.params.H ", state.params.H)
        println("area second ", Integrate(state, gauss_quad_x, gauss_quad_w, α = α))
        println("state.n ", state.n)
        X = LinRange(T(best_rt - 1.0), T(best_rt + 1.0), length(gauss_quad_w))
        Plots.plot!(p, X,  
                    [F(state, x) for x in X].*norm_factor, #Evaluate fitted model at each point. 
                    fillrange = [0.0 for x in 1:length(gauss_quad_w)], 
                    alpha = 0.25, color = :red, show = true
                    ); 


        Plots.vline!(p, [best_rt], color = :blue);
        Plots.vline!(p,[state.t[1]], color = :red);
        Plots.vline!(p, [state.t[state.max_index]], color = :red);
        
        display(p)
        peak_area = Integrate(state, gauss_quad_x, gauss_quad_w, α = α)
        println("peak_area $peak_area")
    end

    ############
    #Calculate Features
    
    #MODEL_INTENSITY = EGH(rt, Tuple(EGH_FIT.param))
    peak_area = Integrate(state, gauss_quad_x, gauss_quad_w, α = α)
   
    sum_of_residuals = zero(Float32)
    points_above_FWHM = zero(Int32)
    points_above_FWHM_01 = zero(Int32)
    
    total_intensity = zero(Float32)
    for i in range(1, state.max_index)
        intensity = state.data[i]
        if intensity > (state.params.H*0.5)
            points_above_FWHM += one(Int32)
        elseif intensity > (state.params.H*0.01)
            points_above_FWHM_01 += one(Int32)
        end 
        if state.mask[i] == false
            sum_of_residuals += abs.(intensity - state.y[i])
            total_intensity += intensity
        end
    end
    GOF = 1 - sum_of_residuals/total_intensity
    if isinf(GOF)
        GOF = missing
    end
    FWHM = getFWHM(state, 0.5)
    FWHM_01 = getFWHM(state, 0.01)
   

    ############
    #Model Parameters
    σ = state.params.σ
    tᵣ = state.params.tᵣ
    τ = state.params.τ
    H = state.params.H
    
    ############
    #Summary Statistics 
    function getSummaryScores!(state::GD_state{HuberParams{T}, U, I, J}, chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}) where {T,U<:AbstractFloat, I,J<:Integer}
        max_scribe_score = -100.0
        max_matched_ratio = -100.0
        max_entropy = -100.0
        max_score = -100.0
        mean_score = 0.0
        max_city_fitted = -100.0
        mean_city_fitted = 0.0
        count = 0
        y_ions_sum = 0
        max_y_ions = 0

        start = max(1, best_scan - 2)
        stop = min(length(chrom.weight), best_scan + 2)

        for i in range(1, length(chrom.weight))#range(start, stop)
            if !state.mask[i]
                if chrom.scribe[i]>max_scribe_score
                    max_scribe_score = chrom.scribe[i]
                end

                if chrom.matched_ratio[i]>max_matched_ratio
                    max_matched_ratio = chrom.matched_ratio[i]
                end

                if chrom.entropy_score[i]>max_entropy
                    max_entropy=chrom.entropy_score[i]
                end

                if chrom.city_block_fitted[i]>max_city_fitted
                    max_city_fitted=chrom.max_city_fitted[i]
                end

                if chrom.score[i]>max_score
                    max_score = chrom.score[i]
                end

                mean_score += chrom.score[i]

                y_ion_count = (chrom.b_count[i] + chrom.y_count[i])
                y_ions_sum += y_ion_count

                if y_ion_count > max_y_ions
                    max_y_ions = y_ion_count
                end
                #sum_log_probability += log2(chrom.prob[i])
                mean_city_fitted += chrom.city_block_fitted[i]
                count += 1
            end
        end    

        chrom.max_scribe_score[best_scan] = max_scribe_score
        chrom.max_matched_ratio[best_scan] = max_matched_ratio
        chrom.max_entropy[best_scan] = max_entropy
        chrom.max_city_fitted[best_scan] = max_city_fitted
        chrom.mean_city_fitted[best_scan] = mean_city_fitted/count
        chrom.max_score[best_scan] = max_score
        chrom.mean_score[best_scan] = mean_score/count
        chrom.y_ions_sum[best_scan] = y_ions_sum
        chrom.max_y_ions[best_scan] = max_y_ions

    end

    getSummaryScores!(state, chrom);

    best_scan = argmax(chrom.weight)#argmax(chrom.prob.*(chrom.weight.!=0.0))
    chrom.peak_area[best_scan] = peak_area*norm_factor
    chrom.GOF[best_scan] = GOF
    chrom.FWHM[best_scan] = FWHM
    chrom.FWHM_01[best_scan] = FWHM_01
    chrom.assymetry[best_scan] = atan(abs(state.params.τ)/sqrt(abs(state.params.σ/2)))
    chrom.points_above_FWHM[best_scan] = points_above_FWHM
    chrom.points_above_FWHM_01[best_scan] = points_above_FWHM_01
    chrom.σ[best_scan] = σ
    chrom.tᵣ[best_scan] = tᵣ
    chrom.τ[best_scan] = τ
    chrom.H[best_scan] = H*norm_factor
    chrom.data_points[best_scan] = data_points#getDataPoints(state) - 4
    chrom.fraction_censored[best_scan] =  Float16(((state.max_index - 2) - chrom.data_points[best_scan])/(state.max_index - 2))

    chrom.base_width_min[best_scan] = state.t[state.max_index-2] - state.t[1]
    chrom.best_scan[best_scan] = true

    

    return nothing
end

function integratePrecursors(grouped_precursor_df::GroupedDataFrame{DataFrame}; n_quadrature_nodes::Int64 = 100, intensity_filter_fraction::Float32 = 0.01f0, α::Float32 = 0.01f0, half_width_at_α::Float32 = 0.15f0, LsqFit_tol::Float64 = 1e-3, Lsq_max_iter::Int = 100, tail_distance::Float32 = 0.25f0, isplot::Bool = false)

    gx, gw = gausslegendre(n_quadrature_nodes)
    N = 500
    dtype = eltype(grouped_precursor_df[1].weight)
    tasks_per_thread = 10
    n_chroms = length(grouped_precursor_df)
    chunk_size = max(1, n_chroms ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:n_chroms, chunk_size) # partition your data into chunks that
    #for i in ProgressBar(range(1, length(grouped_precursor_df)))
    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            state = GD_state(
                HuberParams(zero(dtype), zero(dtype),zero(dtype),zero(dtype)), #Initial params
                zeros(dtype, N), #t
                zeros(dtype, N), #y
                zeros(dtype, N), #data
                falses(N), #mask
                0, #number of iterations
                N #max index
                )
            for i in chunk
        #for i in range(1, length(grouped_precursor_df))
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
        end
    end
    fetch.(tasks)
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