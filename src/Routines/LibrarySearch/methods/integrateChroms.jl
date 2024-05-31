function integratePrecursorMS2(chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, 
                                state::GD_state{HuberParams{U}, V, I, J}, 
                                gauss_quad_x::Vector{Float64}, 
                                gauss_quad_w::Vector{Float64}; 
                                intensity_filter_fraction::Float32 = 0.1f0, 
                                α::Float32 = 0.01f0, 
                                half_width_at_α::Float32 = 0.15f0,
                                min_prob = 0f0,
                                max_scan_gap = 0.2f0, #Needs to be between 3 and4 cycles. 
                                max_peak_width = 1.0f0,
                                isplot::Bool = false) where {U,V<:AbstractFloat, I,J<:Integer}
    
    #########
    #Helper Functions  
    #########
    function getBestPSM(weights::AbstractVector{<:AbstractFloat}, topn::AbstractVector{UInt8}, y_count::AbstractVector{UInt8}, min_prob::Float32)
        best_scan_idx = 1
        max_weight = zero(Float32)
        best_scan_under_prob_thresh_idx = 1 
        max_weight_over_prob_thresh = zero(Float32)
        for i in range(1, length(weights))
            #If passing a threshold based on the logistic regression model
            #Automatically give priority. Could choose a better heuristic.
            if true==true #(y_count[i] < 3) | (topn[i] < 2)# (probs[i] <= min_prob) | 
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
        #start_rt, stop_rt = chrom.rt[1], chrom.rt[end]

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

        #half_width_at_α = 0.15
        #Initial Parameter Guesses
        state.params = getP0(T(α), 
                            T(half_width_at_α), 
                            T(half_width_at_α),
                            T(best_rt),
                            T(1.0),
                            #T(best_height),
                            lower_bounds, upper_bounds)

        #Zero at boundaries 
        #start = state.t[1] - (state.t[2] - state.t[1])
        #start = state.t[max_index] - (state.t[state.max_index - 1] - state.t[state.max_index - 1])
        #max_peak_width = 0.062*3 #cal curve
        max_peak_width = 0.0482*3
        start, stop = 0, state.max_index
        found_start = false
        for i in range(1, state.max_index)
            if (!state.mask[i])
                if !found_start
                    start = state.t[i]
                    found_start = true
                end
                stop = state.t[i]
            end
        end
        #state.t[state.max_index + 1], state.t[state.max_index + 2] = T(start - max_peak_width), T(stop + max_peak_width)
        #state.data[state.max_index + 1], state.data[state.max_index + 2] = zero(T), zero(T)
        #state.mask[state.max_index + 1], state.mask[state.max_index + 2] = true, true
        #state.max_index += 2
        GD(state,
                lower_bounds,
                upper_bounds,
                tol = 1e-3, 
                max_iter = 200, 
                δ = 1e-1,#1e-3, #Huber loss parameter. 
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
    best_height, best_scan = getBestPSM(chrom.weight, chrom.topn, chrom.y_count, min_prob)
    norm_factor = best_height 
    fillState!(state, chrom, best_height)
    truncateAfterSkip!(state, best_scan, chrom.RT, max_scan_gap) #Only allow a fixed ammount of time without a scan
    filterOnRT!(state, chrom.RT[best_scan], 
        max_peak_width, 
        chrom.RT) #Remove scans greater than a certain distance from the best scan 
    filterLowIntensity!(state, best_height*intensity_filter_fraction, chrom.weight) #Remove scans with intensity less than a fixed fraction of the most intense scan 
    #filterOnMatchedRatio!(state, best_scan, chrom.matched_ratio) #Remove scans with matched ratios significantly lower than that of the best scan
    data_points = getDataPoints(state)
    ##########
    #Fit EGH Curve to data 
    points_above_half_max = 0
    for i in range(1, state.max_index)
        if !state.mask[i]
            if state.data[i] >= 0.1
                points_above_half_max += 1
            end
        end
    end
    if points_above_half_max < 0
    fitEGH(state, 
            HuberParams(T(0.001), T(0), T(-1), T(1.0)),
            HuberParams(T(1), T(Inf), T(1), T(1.0)),
            max_peak_width,
            α,
            half_width_at_α,
            chrom.RT[best_scan]
            )
    else
        fitEGH(state, 
        #HuberParams(T(0.001), T(0), T(-1), T(0.75)),
        #HuberParams(T(1), T(Inf), T(1), T(1.25)),
        HuberParams(T(-Inf), T(0), T(-Inf), T(0.75)),
        HuberParams(T(-Inf), T(Inf), T(Inf), T(1.25)),
        max_peak_width,
        α,
        half_width_at_α,
        chrom.RT[best_scan]
        )
    end
    if isplot
        println("TEST2")
        #p = plot()
        mi = state.max_index
        #plot(state.t[1:mi], state.data[1:mi], seriestype=:scatter, show = true)

        plot(state.t[1:mi][state.mask[1:mi].==false], norm_factor .* state.data[1:mi][state.mask[1:mi].==false], seriestype=:scatter, alpha = 0.5, show = true)
        #plot!(state.t[1:mi][state.mask[1:mi]], state.data[1:mi][state.mask[1:mi]], seriestype=:scatter, alpha = 0.5)
        xbins = LinRange(state.t[1] - 0.5, state.t[state.max_index] + .5, 100)
        plot!(xbins, [norm_factor*F(state, x) for x in xbins])
        #plot!(chrom.RT, chrom.weight, seriestype=:scatter, alpha = 0.5, show = true)
        println("area = ", norm_factor*Integrate(state, gauss_quad_x, gauss_quad_w, α = α))
        println("norm_factor ", norm_factor)
        println("H ", norm_factor*state.params.H)
    end
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

    start = max(1, best_scan - 1)
    stop = min(length(chrom.scribe), best_scan + 1)
    total = 0.0
    n = 0
    for i in start:stop
        if !state.mask[i]
            total += chrom.weight[i]
            n += 1
        end
    end
    chrom.H[best_scan] = state.params.H#total/n#state.params.H*norm_factor

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
    #chrom.H[best_scan] = state.params.H*norm_factor
    chrom.data_points[best_scan] = data_points#getDataPoints(state) - 4
    chrom.fraction_censored[best_scan] =  Float16(((state.max_index - 2) - data_points)/(state.max_index - 2))
    chrom.base_width_min[best_scan] = state.t[state.max_index-2] - state.t[1]
    chrom.best_scan[best_scan] = true

    

    return nothing
end

function integrateChrom(chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, 
                                linsolve::LinearSolve.LinearCache,
                                u2::Vector{Float32},
                                state::GD_state{HuberParams{U}, V, I, J},
                                gw::Vector{Float64},
                                gx::Vector{Float64}; 
                                α::Float32 = 0.01f0, 
                                height_at_integration_width::Float32 = 0.001f0,
                                n_pad::Int64 = 0,
                                isplot::Bool = false) where {U,V<:AbstractFloat, I,J<:Integer}
    
    #########
    #Helper Functions  
    #########
    function WHSmooth!( linsolve::LinearSolve.LinearCache, 
                        intensities::AbstractVector{Float32},
                        n_pad::Int64)
        #Reset linsolve and second derivative 
        @inbounds for i in range(1, length(linsolve.b))
            linsolve.b[i] = zero(Float32)
            linsolve.u[i] = zero(Float32)
            u2[i] = zero(Float32)
        end
        #Copy data to linsolve
        best_scan, max_intensity = 0, typemin(Float32)
        @inbounds for i in range(1, size(chrom, 1))
            if intensities[i]>max_intensity
                max_intensity = intensities[i]
                best_scan = i
            end
            linsolve.b[i+n_pad] = intensities[i]
        end

        #WH smoothing 

        solve!(linsolve)


        #Best scan is the most intense 
        #best_scan = argmax(linsolve.b)
        
        return best_scan, linsolve.b[best_scan+n_pad]
    end

    function fillU2!(
        u2::Vector{Float32},
        u::Vector{Float32})
        #Get second-order descrete derivative 
        u2[1], u2[end] = zero(Float32), zero(Float32)
        @inbounds @fastmath for i in range(2, length(linsolve.b) - 1)
            u2[i] = u[i + 1] - 2*u[i] + u[i - 1]
        end
    end

    function getIntegrationBounds!(u2::Vector{Float32},
                                   u::Vector{Float32},
                                   N::Int64,
                                   best_scan::Int64,
                                   n_pad::Int64)
        max_stop = N
        start_search, stop_search = best_scan + n_pad - 1, best_scan + n_pad + 1
        start, stop = start_search, stop_search
        N += n_pad

        #get RH boundary
        @inbounds @fastmath begin 

            for i in range(stop_search, N-1)
                if (u2[i-1] < u2[i]) & (u2[i+1]<u2[i])
                    stop = min(i, N)
                    break
                end
            end

            for i in range(stop, N-1)
                if u[i + 1] > u[i]
                    break
                else
                    stop = i
                end
            end

            #get LH boundary 
            for i in reverse(range(2, start_search))
                if (u2[i] > u2[i - 1]) & (u2[i+1]<u2[i])
                    start = max(i, 1)
                    break
                end
            end

            for i in reverse(range(2, start))
                if u[i - 1] > u[i]
                    break
                else
                    start = i
                end
            end
        end
        return range(max(start-n_pad, 1), min(stop-n_pad,  max_stop))#range(max(start-n_pad, 1), max(stop-n_pad, 1))#range( min(best_scan-3, start), max(stop,best_scan+3))
    end

    function fillState!(state::GD_state{HuberParams{Float32}, Float32, Int64, Int64},
                        u::Vector{Float32},
                        rt::AbstractVector{Float16},
                        start::Int64, 
                        stop::Int64,
                        best_scan::Int64,
                        n_pad::Int64
                        )

        start_rt = rt[start]
        best_rt = rt[best_scan]
        #start_rt, best_rt = rt[start], rt[best_scan]
        rt_width = rt[stop] - start_rt

        norm_factor = u[best_scan+n_pad]

        #Write data to state
        #Normalize so that maximum intensity is 1 
        #And time difference from start to finish is 1. 
        @inbounds @fastmath for i in range(1, stop - start + 1)
            n = start + i - 1
            state.t[i] = (chrom[n,:rt] - start_rt)/rt_width
            state.data[i] = u[n+n_pad]/norm_factor
        end

        state.max_index = stop - start + 1
        best_rt = Float32((best_rt - start_rt)/rt_width)
        return norm_factor, start_rt, rt_width, best_rt
    end

    function fitEGH(state::GD_state{HuberParams{T}, U, I, J}, 
                    lower_bounds::HuberParams{T}, 
                    upper_bounds::HuberParams{T},
                    α::T,
                    half_width_at_α::T,
                    best_rt::T) where {T,U<:AbstractFloat, I,J<:Integer}

        #half_width_at_α = 0.15
        #Initial Parameter Guesses
        state.params = getP0(T(α), 
                            T(half_width_at_α), 
                            T(half_width_at_α),
                            T(best_rt),
                            T(0.6),
                            lower_bounds, upper_bounds)

        GD(state,
                lower_bounds,
                upper_bounds,
                tol = 1e-4, 
                max_iter = 300, 
                δ = 1e-5,#1e-3, #Huber loss parameter. 
                α=Float64(α),
                β1 = 0.9,
                β2 = 0.999,
                ϵ = 1e-8)
        
    end

    function getPeakProperties(state::GD_state{HuberParams{T}, U, I, J}) where {T,U<:AbstractFloat, I,J<:Integer}
        points_above_FWHM = zero(Int32)
        points_above_FWHM_01 = zero(Int32)
        for i in range(1, state.max_index)
            intensity = state.data[i]
            if intensity > (state.params.H*0.5)
                points_above_FWHM += one(Int32)
            end
            if intensity > (state.params.H*0.01)
                points_above_FWHM_01 += one(Int32)
            end 
        end

        FWHM = getFWHM(state, 0.5)
        FWHM_01 = getFWHM(state, 0.01)

        return FWHM, FWHM_01, points_above_FWHM, points_above_FWHM_01
    end
    #Baseline subtraction?
    function subtractBaseline!(
        u::Vector{Float32}, #smoothed data
        best_scan::Int64, #peak apex
        scan_range::UnitRange{Int64},
        n_pad::Int64) #start and stop of integration bounds 
        
        best_scan  = best_scan + n_pad
        scan_range = (first(scan_range) + n_pad, last(scan_range) + n_pad)
        #Fine LH baseline 
        lmin,li = typemax(Float32),first(scan_range)
        @inbounds @fastmath for i in range(first(scan_range), best_scan)
            if u[i] < lmin
                lmin = u[i]
                li = i
            end
        end

        #Find RH baseline 
        rmin,ri = typemax(Float32),last(scan_range)
        @inbounds @fastmath for i in range(best_scan, last(scan_range))
            if u[i] < rmin
                rmin = u[i]
                ri = i
            end
        end


        #= Another option is to just use the extrema of the integration boundary 
        lmin = linsolve.u[first(scan_range)]
        li = first(scan_range)
        ri = last(scan_range)
        rmin = linsolve.u[last(scan_range)]
        =#
        #Subtract the baseline 
        
        h = (rmin - lmin)/(ri - li)
        @inbounds @fastmath for i in scan_range
            u[i] = u[i]-(lmin + (i - li)*h)
        end

    end

    function integrateTrapezoidal(state::GD_state{HuberParams{T}, U, I, J}) where {T,U<:AbstractFloat, I,J<:Integer}
        retval = state.data[2]
        #Assumption that state.max_index is at least 4. 
        for i in range(3, state.max_index - 1)
            @inbounds retval = (state.t[2] - state.t[1])*(state.data[1] + state.data[2])
            @inbounds @fastmath for i in 2:(state.max_index - 1)
                retval += (state.t[i + 1] - state.t[i])*(state.data[i] + state.data[i + 1])
            end
        end
        return (1//2)*retval
    end
    #Whittaker Henderson Smoothing
    best_scan, max_intensity = WHSmooth!(
        linsolve,
        chrom[!,:intensity],
        n_pad
    )
    #Second discrete derivative of smoothed data
    fillU2!(
        u2,
        linsolve.u,
    )
    
    #Integration boundaries based on smoothed second derivative 
    scan_range = getIntegrationBounds!(
        u2,
        linsolve.u,
        size(chrom, 1),
        best_scan,
        n_pad
    )

    subtractBaseline!(
        linsolve.u,
        best_scan,
        scan_range,
        n_pad
    )

    #File `state` to fit EGH function. Get the inensity, and rt normalization factors 
    norm_factor, start_rt, rt_norm, best_rt = fillState!(
        state,
        linsolve.u,
        chrom[!,:rt],
        first(scan_range),
        last(scan_range),
        best_scan,
        n_pad
    )
    
    #Initial estimate for FWHM
    a = max(best_scan - 1, 1)
    b = min(best_scan + 1, size(chrom, 1))

    half_width_at_α = Float32(((chrom[b,:rt]-start_rt)/rt_norm - (chrom[a,:rt] - start_rt)/rt_norm)/2)
    T = eltype(chrom.intensity)
    ##########
    #Fit EGH to data. 
    fitEGH(state, 
            HuberParams(T(0.001), T(0), T(-1), T(0.95)),
            HuberParams(T(1),  T(Inf), T(1), T(1.05)),
            α,
            half_width_at_α,
            best_rt
            )

    if isplot
        mi = state.max_index
        start = max(best_scan - 18, 1)
        stop = min(best_scan + 18, length(chrom.rt))
        plot(chrom.rt[start:stop], chrom.intensity[start:stop], seriestype=:scatter, alpha = 0.5, show = true)
        vline!([chrom.rt[first(scan_range)], chrom.rt[last(scan_range)]])
        plot!(state.t[1:mi].*rt_norm .+ start_rt, norm_factor.*state.data[1:mi], seriestype=:scatter, alpha = 0.5, show = true)
        xbins = LinRange(state.t[1]-0.5, state.t[state.max_index]+0.5, 100)
        plot!(xbins.*rt_norm .+ start_rt, [norm_factor*F(state, x) for x in xbins])
        plot!(chrom.rt[start:stop], u2[start+n_pad:stop+n_pad])
        hline!([norm_factor*0.95])
    end

    peak_area = rt_norm*norm_factor*Integrate(state, 
                                        gx,
                                        gw, 
                                        α = height_at_integration_width)


    trapezoid_area = rt_norm*norm_factor*integrateTrapezoidal(state)

    #trapezoid_area = 0.0f0
    FWHM, FWHM_01, points_above_FWHM, points_above_FWHM_01 = getPeakProperties(state)
    best_scan_idx = chrom[best_scan,:scan_idx]
    return best_scan_idx, peak_area, trapezoid_area, max_intensity,  FWHM, FWHM_01, points_above_FWHM, points_above_FWHM_01 
end

    function getSummaryScores!(
                                psms_integrated::DataFrame,
                                best_scan::Int64,
                                psms_idx::Int64,
                                scribe::AbstractVector{Float16},
                                matched_ratio::AbstractVector{Float16},
                                entropy::AbstractVector{Float16},
                                city_block_fitted::AbstractVector{Float16},
                                combined_score::AbstractVector{Float32},
                                y_count::AbstractVector{UInt8})

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

        psms_integrated.max_scribe_score[psms_idx] = max_scribe_score
        psms_integrated.max_matched_ratio[psms_idx] = max_matched_ratio
        psms_integrated.max_entropy[psms_idx] = max_entropy
        psms_integrated.max_city_fitted[psms_idx] = max_city_fitted
        psms_integrated.mean_city_fitted[psms_idx] = mean_city_fitted/count
        psms_integrated.max_score[psms_idx] = max_score
        psms_integrated.mean_score[psms_idx] = mean_score/count
        psms_integrated.y_ions_sum[psms_idx] = y_ions_sum
        psms_integrated.max_y_ions[psms_idx] = max_y_ions

    end
function integratePrecursors(spectral_scores::GroupedDataFrame{DataFrame},
                             chromatograms::GroupedDataFrame{DataFrame},
                             precs_to_integrate::Vector{DataFrames.GroupKey{GroupedDataFrame{DataFrame}}},
                             psms_integrated::DataFrame; 
                             λ::Float32 = 1.0f0,
                             α::Float32 = 0.01f0,
                             height_at_integration_width::Float32 = 0.001f0,
                             n_quadrature_nodes::Int64 = 100)

    gx, gw = gausslegendre(n_quadrature_nodes)
    dtype = eltype(spectral_scores[1].weight)
    thread_tasks = partitionThreadTasks(length(precs_to_integrate), 10, Threads.nthreads())

    N = 0
    for i in range(1, length(precs_to_integrate))
        if size(chromatograms[i], 1) > N
            N =  size(chromatograms[i], 1)
        end
    end
    N += 20

    tasks = map(thread_tasks) do chunk
        Threads.@spawn begin

            b = zeros(Float32, N);
            A = getWittakerHendersonDesignMat(length(b), λ);
            prob = LinearProblem(A, b);
            linsolve = init(prob);
            u2 = zeros(Float32, length(linsolve.b));

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
                prec_key = precs_to_integrate[i]

                #Chromatograms must be sorted by retention time 
                chroms = chromatograms[(precursor_idx = prec_key[:precursor_idx], isotopes_captured = prec_key[:isotopes_captured])]
                sort!(chroms,:rt, alg = QuickSort)
                best_scan_idx, peak_area, trapezoid_area, max_intensity, FWHM, FWHM_01, points_above_FWHM, points_above_FWHM_01 = integrateChrom(
                                chroms,
                                linsolve,
                                u2,
                                state,
                                gw,gx,
                                α = α,
                                height_at_integration_width = height_at_integration_width,
                                n_pad = 10,
                                isplot = false
                                );
                #println("peak_area $peak_area i $i")
                psms_integrated[i,:scan_idx] = best_scan_idx
                psms_integrated[i,:peak_area] = peak_area
                psms_integrated[i,:trapezoid_area] = trapezoid_area
                psms_integrated[i,:weight] = max_intensity
                psms_integrated[i,:FWHM] = FWHM
                psms_integrated[i,:FWHM_01] = FWHM_01
                psms_integrated[i,:points_above_FWHM] = points_above_FWHM
                psms_integrated[i,:points_above_FWHM_01] = points_above_FWHM_01
                psms_integrated[i,:precursor_idx] = prec_key[:precursor_idx]
                psms_integrated[i,:isotopes_captured] = prec_key[:isotopes_captured]
                reset!(state)


                psms = spectral_scores[(precursor_idx = prec_key[:precursor_idx], isotopes_captured = prec_key[:isotopes_captured])]
                sort!(psms, :RT, alg = QuickSort)
                
                best_scan = missing 
                for i in range(1, size(psms, 1))
                    if psms[i,:scan_idx] == best_scan_idx
                        best_scan = i
                        break
                    end
                end
                if ismissing(best_scan)
                    psms_integrated[i,:scan_idx] = missing
                    continue
                end

                psms_integrated[i,:best_rank] = psms[best_scan,:best_rank]
                psms_integrated[i,:topn] = psms[best_scan,:topn]
                psms_integrated[i,:longest_y] = psms[best_scan,:longest_y]
                psms_integrated[i,:b_count] = psms[best_scan,:b_count]
                psms_integrated[i,:y_count] = psms[best_scan,:y_count]
                psms_integrated[i,:p_count] = psms[best_scan,:p_count]
                psms_integrated[i,:isotope_count] = psms[best_scan,:isotope_count]
                psms_integrated[i,:total_ions] = psms[best_scan,:total_ions]
                psms_integrated[i,:poisson] = psms[best_scan,:poisson]
                psms_integrated[i,:hyperscore] = psms[best_scan,:hyperscore]
                psms_integrated[i,:log2_intensity_explained] = psms[best_scan,:log2_intensity_explained]
                psms_integrated[i,:error] = psms[best_scan,:error]
                psms_integrated[i,:error_norm] = psms[best_scan,:err_norm]

                psms_integrated[i,:scribe] = psms[best_scan,:scribe]
                psms_integrated[i,:scribe_fitted] = psms[best_scan,:scribe_fitted]
                psms_integrated[i,:city_block] = psms[best_scan,:city_block]
                psms_integrated[i,:city_block_fitted] = psms[best_scan,:city_block_fitted]
                psms_integrated[i,:spectral_contrast] = psms[best_scan,:spectral_contrast]
                psms_integrated[i,:matched_ratio] = psms[best_scan,:matched_ratio]
                psms_integrated[i,:entropy_score] = psms[best_scan,:entropy_score]
                psms_integrated[i,:RT] = psms[best_scan,:RT]
                psms_integrated[i,:target] = psms[best_scan,:target]
                psms_integrated[i,:cv_fold] = psms[best_scan,:cv_fold]
                getSummaryScores!(
                    psms_integrated,
                    best_scan,
                    i,
                    psms[!,:scribe],
                    psms[!,:matched_ratio],
                    psms[!,:entropy_score],
                    psms[!,:city_block_fitted],
                    psms[!,:score],
                    psms[!,:y_count]
                )
            end
        end
    end
    fetch.(tasks)
end
