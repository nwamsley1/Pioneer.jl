function integratePrecursorMS2(chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, 
                                state::GD_state{HuberParams{U}, V, I, J}, 
                                gauss_quad_x::Vector{Float64}, 
                                gauss_quad_w::Vector{Float64}; 
                                intensity_filter_fraction::Float32 = 0.1f0, 
                                α::Float32 = 0.01f0, 
                                half_width_at_α::Float32 = 0.15f0,
                                min_prob = 0f0,
                                max_scan_gap = 0.2f0,
                                max_peak_width = 1.0f0) where {U,V<:AbstractFloat, I,J<:Integer}
    
    #########
    #Helper Functions  
    #########
    function getBestPSM(weights::AbstractVector{<:AbstractFloat}, probs::AbstractVector{<:AbstractFloat}, min_prob::Float32)
        best_scan_idx = 1
        max_weight = zero(Float32)
        best_scan_under_prob_thresh_idx = 1 
        max_weight_over_prob_thresh = zero(Float32)
        for i in range(1, length(weights))
            #If passing a threshold based on the logistic regression model
            #Automatically give priority. Could choose a better heuristic.
            if probs[i] <= min_prob
                if weights[i] > max_weight_over_prob_thresh
                    max_weight_over_prob_thresh = weights[i]
                    best_scan_under_prob_thresh_idx = i
                end
            else
                if weights[i] > max_weight
                    max_weight = weights[i]
                    best_scan_idx = i
                end
            end
        end

        if iszero(max_weight)#max_weight > max_weight_over_prob_thresh
            return max_weight_over_prob_thresh, best_scan_under_prob_thresh_idx
        else
            return max_weight, best_scan_idx
        end
    end

    function fillState!(state::GD_state{HuberParams{T}, U, I, J}, chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, max_weight::Float32) where {T,U<:AbstractFloat, I,J<:Integer}
        for i in range(1, length(chrom.weight))
            state.t[i] = chrom.RT[i]
            state.data[i] = chrom.weight[i]/max_weight
        end
        state.max_index = length(chrom.weight)
        return nothing
    end

    function truncateAfterSkip!(state::GD_state{HuberParams{T}, U, I, J}, 
                                best_scan::Int64, 
                                rts::SubArray{Float32, 1, Vector{Float32}, Tuple{Vector{Int64}}, false},
                                max_scan_gap::Float32) where {T,U<:AbstractFloat, I,J<:Integer}
        
        
        for i in range(best_scan, length(rts)-1)
            if (rts[i+ 1] - rts[i]) > max_scan_gap
                for n in range(i + 1, length(rts))
                    state.mask[n] = true
                end
            end
        end
        
        for i in range(1, best_scan - 1)
            if (rts[best_scan - i + 1] - rts[best_scan - i]) > max_scan_gap
                for n in range(1, best_scan-i)
                    state.mask[n] = true
                end
            end
        end

        return 

    end

    function filterOnRT!(state::GD_state{HuberParams{T}, U, I, J}, 
                            best_rt::T, 
                            max_peak_width::T,
                            rts::SubArray{Float32, 1, Vector{Float32}, Tuple{Vector{Int64}}, false}) where {T,U<:AbstractFloat, I,J<:Integer}
        for i in eachindex(rts)
            if (rts[i] > (best_rt - max_peak_width)) & (rts[i] < (best_rt + max_peak_width))
                continue
            else
                state.mask[i] = true
            end
        end
    end
   
    function getDataPoints(state::GD_state{HuberParams{T}, U, I, J}) where {T,U<:AbstractFloat, I,J<:Integer}
        data_points = zero(UInt32)
        for i in range(1, state.max_index)
            if state.mask[i]==false
                data_points += one(UInt32)
            end
        end
        return data_points
    end

    function filterLowIntensity!(state::GD_state{HuberParams{T}, U, I, J}, 
                                min_intensity::T, 
                                weights::AbstractVector{T}) where {T,U<:AbstractFloat, I,J<:Integer}
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

    function filterOnMatchedRatio!(state::GD_state{HuberParams{T}, U, I, J}, best_scan_idx::Int64, matched_ratios::SubArray{V, 1, Vector{V}, Tuple{Vector{Int64}}, false}) where {T,U,V<:AbstractFloat, I,J<:Integer}
        for i in range(1, state.max_index)
            if (matched_ratios[i] < (matched_ratios[best_scan_idx] - 1)) & (matched_ratios[i] < 0)
                state.mask[i] = true
            end
        end
    end

    function fitEGH(state::GD_state{HuberParams{T}, U, I, J}, 
                    lower_bounds::HuberParams{T}, 
                    upper_bounds::HuberParams{T},
                    max_peak_width::T,
                    α::T,
                    half_width_at_α::T,
                    best_rt::T) where {T,U<:AbstractFloat, I,J<:Integer}

        #Initial Parameter Guesses
        state.params = getP0(T(α), 
                            T(half_width_at_α), 
                            T(half_width_at_α),
                            T(best_rt),
                            T(best_height),
                            lower_bounds, upper_bounds)

        #Zero at boundaries 
        state.t[state.max_index + 1], state.t[state.max_index + 2] = T(best_rt - max_peak_width), T(best_rt + max_peak_width)
        state.data[state.max_index + 1], state.data[state.max_index + 2] = zero(T), zero(T)
        state.max_index += 2
        GD(state,
                lower_bounds,
                upper_bounds,
                tol = 1e-3, 
                max_iter = 200, 
                δ = 1e-3, #Huber loss parameter. 
                α=0.01,
                β1 = 0.9,
                β2 = 0.999,
                ϵ = 1e-8)
        
    end

    function sumResiduals(state::GD_state{HuberParams{T}, U, I, J}) where {T,U<:AbstractFloat, I,J<:Integer}
        sum_of_residuals = zero(Float32) #sum of residuals 
        points_above_FWHM = zero(Int32) #full width above 50% of height
        points_above_FWHM_01 = zero(Int32) #full width above 1% of height
        
        total_intensity = zero(Float32)
        for i in range(1, state.max_index)
            state.mask[i] ? continue : nothing

            intensity = state.data[i]
            if intensity > (state.params.H*0.5)
                points_above_FWHM += one(Int32)
            end
            if intensity > (state.params.H*0.01)
                points_above_FWHM_01 += one(Int32)
            end 

            sum_of_residuals += abs.(intensity - state.y[i])
            total_intensity += intensity

        end
        GOF = 1 - sum_of_residuals/total_intensity

        if isinf(GOF)
            GOF = missing
        end

        FWHM = getFWHM(state, 0.5)
        FWHM_01 = getFWHM(state, 0.01)
        return GOF, FWHM, FWHM_01, points_above_FWHM, points_above_FWHM_01
    end

    function getSummaryScores!(chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}},
                                best_scan::Int64,
                                state::GD_state{HuberParams{T}, U, I, J}, 
                                scribe::AbstractVector{Float16},
                                matched_ratio::AbstractVector{Float16},
                                entropy::AbstractVector{Float16},
                                city_block_fitted::AbstractVector{Float16},
                                combined_score::AbstractVector{Float32},
                                y_count::AbstractVector{UInt8}) where {T,U<:AbstractFloat, I,J<:Integer}
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
        stop = min(length(scribe), best_scan + 2)

        for i in range(start, stop)
            if !state.mask[i]
                if scribe[i]>max_scribe_score
                    max_scribe_score =scribe[i]
                end

                if matched_ratio[i]>max_matched_ratio
                    max_matched_ratio = matched_ratio[i]
                end

                if entropy[i]>max_entropy
                    max_entropy=entropy[i]
                end

                if city_block_fitted[i]>max_city_fitted
                    max_city_fitted = city_block_fitted[i]
                end

                if combined_score[i]>max_score
                    max_score = combined_score[i]
                end

                
                y_ions_sum += y_count[i]
                if y_count[i] > max_y_ions
                    max_y_ions = y_count[i]
                end

                mean_score += combined_score[i]
                mean_city_fitted += city_block_fitted[i]
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

    
    T = eltype(chrom.weight)

    ##########
    #Initialize state with observed data 
    max_peak_width = 1.0f0
    best_height, best_scan = getBestPSM(chrom.weight, chrom.prob, min_prob)
    norm_factor = best_height 
    fillState!(state, chrom, best_height)
    truncateAfterSkip!(state, best_scan, chrom.RT, max_scan_gap) #Only allow a fixed ammount of time without a scan
    filterOnRT!(state, chrom.RT[best_scan], 
        max_peak_width, 
        chrom.RT) #Remove scans greater than a certain distance from the best scan 
    filterLowIntensity!(state, best_height*intensity_filter_fraction, chrom.weight) #Remove scans with intensity less than a fixed fraction of the most intense scan 
    filterOnMatchedRatio!(state, best_scan, chrom.matched_ratio) #Remove scans with matched ratios significantly lower than that of the best scan
    data_points = getDataPoints(state)
    ##########
    #Fit EGH Curve to data 
    fitEGH(state, 
            HuberParams(T(0.001), T(0), T(-1), T(1.0)),
            HuberParams(T(1), T(Inf), T(1), T(1.0)),
            max_peak_width,
            α,
            half_width_at_α,
            chrom.RT[best_scan]
            )

    ##########
    #Fit EGH Curve to data
    peak_area = Integrate(state, gauss_quad_x, gauss_quad_w, α = α)
    GOF, FWHM, FWHM_01, points_above_FWHM, points_above_FWHM_01 = sumResiduals(state)

   
    ############
    #Summary Statistics 
    getSummaryScores!(chrom,
                        best_scan,
                        state,
                        chrom.scribe,
                        chrom.matched_ratio,
                        chrom.entropy_score,
                        chrom.city_block_fitted,
                        chrom.score,
                        chrom.y_count)

    ##########
    #Fill Features 
    chrom.peak_area[best_scan] = peak_area*norm_factor
    chrom.GOF[best_scan] = GOF
    chrom.FWHM[best_scan] = FWHM
    chrom.FWHM_01[best_scan] = FWHM_01
    chrom.assymetry[best_scan] = atan(abs(state.params.τ)/sqrt(abs(state.params.σ/2)))
    chrom.points_above_FWHM[best_scan] = points_above_FWHM
    chrom.points_above_FWHM_01[best_scan] = points_above_FWHM_01
    chrom.σ[best_scan] = state.params.σ
    chrom.tᵣ[best_scan] = state.params.tᵣ
    chrom.τ[best_scan] = state.params.τ
    chrom.H[best_scan] = state.params.H*norm_factor
    chrom.data_points[best_scan] = data_points#getDataPoints(state) - 4
    chrom.fraction_censored[best_scan] =  Float16(((state.max_index - 2) - data_points)/(state.max_index - 2))
    chrom.base_width_min[best_scan] = state.t[state.max_index-2] - state.t[1]
    chrom.best_scan[best_scan] = true

    

    return nothing
end

function integratePrecursors(grouped_precursor_df::GroupedDataFrame{DataFrame}; 
                                n_quadrature_nodes::Int64 = 100, 
                                intensity_filter_fraction::Float32 = 0.01f0, 
                                α::Float32 = 0.01f0, 
                                half_width_at_α::Float32 = 0.15f0                                                                    )

    gx, gw = gausslegendre(n_quadrature_nodes)
    N = 500
    dtype = eltype(grouped_precursor_df[1].weight)
    thread_tasks = partitionThreadTasks(length(grouped_precursor_df), 10, Threads.nthreads())
    #for i in ProgressBar(range(1, length(grouped_precursor_df)))
    println("thread_tasks $thread_tasks")
    tasks = map(thread_tasks) do chunk
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
                integratePrecursorMS2(grouped_precursor_df[i],
                                        state,
                                        gx::Vector{Float64},
                                        gw::Vector{Float64},
                                        intensity_filter_fraction = intensity_filter_fraction,
                                        α = α,
                                        half_width_at_α = half_width_at_α
                                        )
                reset!(state)
            end
        end
    end
    fetch.(tasks)
end
