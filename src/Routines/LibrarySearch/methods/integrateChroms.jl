mutable struct Chromatogram{T<:Real, J<:Integer}
    t::Vector{T}
    data::Vector{T}
    max_index::J
end

function reset!(state::Chromatogram)
    for i in range(1, state.max_index)
        state.t[i], state.data[i] = zero(eltype(state.t)), zero(eltype(state.data))
    end
    state.max_index = 0
    return 
end



function integrateChrom(chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, 
                                linsolve::LinearSolve.LinearCache,
                                u2::Vector{Float32},
                                state::Chromatogram,
                                gw::Vector{Float64},
                                gx::Vector{Float64}; 
                                height_at_integration_width::Float32 = 0.001f0,
                                n_pad::Int64 = 0,
                                isplot::Bool = false)
    
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

    function fillState!(state::Chromatogram,
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

    function getPeakProperties(state::Chromatogram)
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

    function integrateTrapezoidal(state::Chromatogram)
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

    
    peak_area = 0.0f0


    trapezoid_area = rt_norm*norm_factor*integrateTrapezoidal(state)

    #trapezoid_area = 0.0f0
    FWHM, FWHM_01, points_above_FWHM, points_above_FWHM_01 = 0.0f0, 0.0f0, 0.0f0, 0.0f0#getPeakProperties(state)
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

            state = Chromatogram(
                zeros(dtype, N), #t
                zeros(dtype, N), #data
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
