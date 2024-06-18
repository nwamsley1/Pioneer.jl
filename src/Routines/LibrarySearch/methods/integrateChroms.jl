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
                                best_scan::UInt32,
                                linsolve::LinearSolve.LinearCache,
                                u2::Vector{Float32},
                                state::Chromatogram; 
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
        @inbounds for i in range(1, size(chrom, 1))
            #if 
            #    max_intensity = intensities[i]
            #    best_scan = i
            #end
            linsolve.b[i+n_pad] = intensities[i]
        end
        #WH smoothing 
        solve!(linsolve)
        return
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
    WHSmooth!(
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
    trapezoid_area = rt_norm*norm_factor*integrateTrapezoidal(state)

    #trapezoid_area = 0.0f0
    return trapezoid_area, chrom[!,:scan_idx][best_scan]
end

function getSummaryScores!(
                            psms::SubDataFrame,
                            weight::AbstractVector{Float32},
                            scribe::AbstractVector{Float16},
                            matched_ratio::AbstractVector{Float16},
                            entropy::AbstractVector{Float16},
                            city_block_fitted::AbstractVector{Float16},
                            y_count::AbstractVector{UInt8},

                        )

    max_scribe_score = -100.0
    max_matched_ratio = -100.0
    max_entropy = -100.0
    max_city_fitted = -100.0
    mean_city_fitted = 0.0
    count = 0
    y_ions_sum = 0
    max_y_ions = 0

    best_scan = argmax(psms[!,:weight])
    #Need to make sure there is not a big gap. 
    start = max(1, best_scan - 2)
    stop = min(length(weight), best_scan + 2)

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
    
        y_ions_sum += y_count[i]
        if y_count[i] > max_y_ions
            max_y_ions = y_count[i]
        end

        mean_city_fitted += city_block_fitted[i]
        count += 1
    end    

    psms.max_scribe_score[best_scan] = max_scribe_score
    psms.max_matched_ratio[best_scan] = max_matched_ratio
    psms.max_entropy[best_scan] = max_entropy
    psms.max_city_fitted[best_scan] = max_city_fitted
    psms.mean_city_fitted[best_scan] = mean_city_fitted/count
    psms.y_ions_sum[best_scan] = y_ions_sum
    psms.max_y_ions[best_scan] = max_y_ions
    psms.best_scan[best_scan] = true

end

function integratePrecursors(chromatograms::GroupedDataFrame{DataFrame},
                             precursor_idx::AbstractVector{UInt32},
                             isotopes_captured::AbstractVector{Tuple{Int8, Int8}},
                             apex_scan_idx::AbstractVector{UInt32},
                             peak_area::AbstractVector{Float32},
                             new_best_scan::AbstractVector{UInt32}; 
                             λ::Float32 = 1.0f0)

    dtype = Float32
    thread_tasks = partitionThreadTasks(length(precursor_idx), 10, Threads.nthreads())

    #Maximal size of a chromatogram
    N = 0
    for (chrom_id, chrom) in pairs(chromatograms)
        if size(chrom, 1) > N
            N = size(chrom, 1)
        end
    end
    N += 20
    group_keys = keys(chromatograms)
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
                prec_id = precursor_idx[i]
                iso_set = isotopes_captured[i]
                apex_scan = apex_scan_idx[i]
                #Chromatograms must be sorted by retention time 
                (precursor_idx = prec_id, isotopes_captured = iso_set) ∉ group_keys ? continue : nothing
                chrom = chromatograms[(precursor_idx = prec_id, isotopes_captured = iso_set)]
                apex_scan = findfirst(x->x==apex_scan,chrom[!,:scan_idx]::AbstractVector{UInt32}) 
                if isnothing(apex_scan) ? continue : nothing
                sort!(chrom,:rt, alg = QuickSort)
                peak_area[i], new_best_scan[i] = integrateChrom(
                                chrom,
                                apex_scan,
                                linsolve,
                                u2,
                                state,
                                n_pad = 10,
                                isplot = false
                                );
                                
                #println("peak_area $peak_area i $i")
                reset!(state)
            end
        end
    end
    fetch.(tasks)
end
