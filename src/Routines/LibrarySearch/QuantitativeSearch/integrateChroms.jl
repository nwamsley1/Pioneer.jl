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
                                apex_scan::Int64,
                                linsolve::LinearSolve.LinearCache,
                                u2::Vector{Float32},
                                state::Chromatogram; 
                                max_apex_offset = 2,
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

    function getApexScan(
        apex_scan::Int64,
        max_offset::Int64,
        intensities::AbstractVector{<:AbstractFloat})
        N = length(intensities)
        max_intensity = zero(Float32)
        for i in range(max(1, apex_scan - max_offset), min(N, apex_scan+max_offset))
            if intensities[i] > max_intensity#intensities[apex_scan]
                apex_scan = i
                max_intensity = intensities[i]
            end
        end
        return apex_scan 
    end

    function getIntegrationBounds!(u2::Vector{Float32},
                                   u::Vector{Float32},
                                   N::Int64,
                                   apex_scan::Int64,
                                   n_pad::Int64)
        max_stop = N
        start_search, stop_search = apex_scan + n_pad - 1, apex_scan + n_pad + 1
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
                        apex_scan::Int64,
                        n_pad::Int64
                        )

        start_rt = rt[start]
        best_rt = rt[apex_scan]
        #start_rt, best_rt = rt[start], rt[best_scan]
        rt_width = rt[stop] - start_rt

        norm_factor = u[apex_scan+n_pad]

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

    #Baseline subtraction?
    function subtractBaseline!(
        u::Vector{Float32}, #smoothed data
        apex_scan::Int64, #peak apex
        scan_range::UnitRange{Int64},
        n_pad::Int64) #start and stop of integration bounds 
        
        apex_scan  = apex_scan + n_pad
        scan_range = (first(scan_range) + n_pad, last(scan_range) + n_pad)
        #Fine LH baseline 
        lmin,li = typemax(Float32),first(scan_range)
        @inbounds @fastmath for i in range(first(scan_range), apex_scan)
            if u[i] < lmin
                lmin = u[i]
                li = i
            end
        end

        #Find RH baseline 
        rmin,ri = typemax(Float32),last(scan_range)
        @inbounds @fastmath for i in range(apex_scan, last(scan_range))
            if u[i] < rmin
                rmin = u[i]
                ri = i
            end
        end

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
    
    apex_scan = getApexScan(
        apex_scan,
        max_apex_offset,
        chrom[!,:intensity]
    )
    #Integration boundaries based on smoothed second derivative 
    scan_range = getIntegrationBounds!(
        u2,
        linsolve.u,
        size(chrom, 1),
        apex_scan,
        n_pad
    )

    subtractBaseline!(
        linsolve.u,
        apex_scan,
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
        apex_scan,
        n_pad
    )
    
    if isplot
        mi = state.max_index
        start = max(apex_scan - 18, 1)
        stop = min(apex_scan + 18, length(chrom.rt))
        plot(chrom.rt[start:stop], chrom.intensity[start:stop], seriestype=:scatter, alpha = 0.5, show = true, label = "raw")
        vline!([chrom.rt[first(scan_range)], chrom.rt[last(scan_range)]], label = nothing)
        plot!(state.t[1:mi].*rt_norm .+ start_rt, norm_factor.*state.data[1:mi], seriestype=:scatter, alpha = 0.5, show = true, label = "smooth")
        #xbins = LinRange(state.t[1]-0.5, state.t[state.max_index]+0.5, 100)
        #plot!(xbins.*rt_norm .+ start_rt, [norm_factor*F(state, x) for x in xbins])
        #plot!(chrom.rt[start:stop], u2[start+n_pad:stop+n_pad])
        #hline!([norm_factor*0.95])
    end

    trapezoid_area = rt_norm*norm_factor*integrateTrapezoidal(state)

    #trapezoid_area = 0.0f0
    return trapezoid_area, chrom[!,:scan_idx][apex_scan]
end

function getSummaryScores!(
                            psms::SubDataFrame,
                            weight::AbstractVector{Float32},
                            gof::AbstractVector{Float16},
                            matched_ratio::AbstractVector{Float16},
                            #entropy::AbstractVector{Float16},
                            fitted_manhattan_distance::AbstractVector{Float16},
                            fitted_spectral_contrast::AbstractVector{Float16},
                            y_count::AbstractVector{UInt8},

                        )

    max_gof = -100.0
    max_matched_ratio = -100.0
   # max_entropy = -100.0
    max_fitted_manhattan_distance = -100.0
    max_fitted_spectral_contrast= -100
    count = 0
    y_ions_sum = 0
    max_y_ions = 0

    apex_scan = argmax(psms[!,:weight])
    #Need to make sure there is not a big gap. 
    start = max(1, apex_scan - 2)
    stop = min(length(weight), apex_scan + 2)

    for i in range(start, stop)
        if gof[i]>max_gof
            max_gof =gof[i]
        end

        if matched_ratio[i]>max_matched_ratio
            max_matched_ratio = matched_ratio[i]
        end

        #if entropy[i]>max_entropy
        #    max_entropy=entropy[i]
        #end

        if fitted_manhattan_distance[i]>max_fitted_manhattan_distance
            max_fitted_manhattan_distance = fitted_manhattan_distance[i]
        end

        if fitted_spectral_contrast[i]>max_fitted_spectral_contrast
            max_fitted_spectral_contrast = fitted_spectral_contrast[i]
        end
    
        y_ions_sum += y_count[i]
        if y_count[i] > max_y_ions
            max_y_ions = y_count[i]
        end

        count += 1
    end    

    psms.max_gof[apex_scan] = max_gof
    psms.max_matched_ratio[apex_scan] = max_matched_ratio
   # psms.max_entropy[apex_scan] = max_entropy
    psms.max_fitted_manhattan_distance[apex_scan] = max_fitted_manhattan_distance
    psms.max_fitted_spectral_contrast[apex_scan] = max_fitted_spectral_contrast
    psms.y_ions_sum[apex_scan] = y_ions_sum
    psms.max_y_ions[apex_scan] = max_y_ions
    psms.best_scan[apex_scan] = true

end


function integratePrecursors(chromatograms::DataFrame,
                             isotope_trace_type::IsotopeTraceType,
                             precursor_idx::AbstractVector{UInt32},
                             isotopes_captured::AbstractVector{Tuple{Int8, Int8}},
                             apex_scan_idx::AbstractVector{UInt32},
                             peak_area::AbstractVector{Float32},
                             new_best_scan::AbstractVector{UInt32}; 
                             λ::Float32 = 1.0f0,
                             n_pad::Int64 = 20,
                             max_apex_offset::Int64 = 2,
                             )
    chromatogram_keys = [:precursor_idx]
    if seperateTraces(isotope_trace_type)
        chromatogram_keys = [:precursor_idx,:isotopes_captured]
    end
    grouped_chroms = groupby(chromatograms, chromatogram_keys)
    dtype = Float32
    thread_tasks = partitionThreadTasks(length(precursor_idx), 10, Threads.nthreads())

    #Maximal size of a chromatogram
    N = 0
    for (chrom_id, chrom) in pairs(grouped_chroms)
        if size(chrom, 1) > N
            N = size(chrom, 1)
        end
    end
    N += n_pad*2
    group_keys = keys(grouped_chroms)
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
                #grouped_chroms must be sorted by retention time 
                if seperateTraces(isotope_trace_type)
                    (precursor_idx = prec_id, isotopes_captured = iso_set) ∉ group_keys ? continue : nothing
                    chrom = grouped_chroms[(precursor_idx = prec_id, isotopes_captured = iso_set)]
                else
                    (precursor_idx = prec_id,) ∉ group_keys ? continue : nothing
                    chrom = grouped_chroms[(precursor_idx = prec_id,)]
                end
                sort!(chrom,
                        :scan_idx, 
                        alg = QuickSort) #Could alternatively sort by :rt
                apex_scan = findfirst(x->x==apex_scan,chrom[!,:scan_idx]::AbstractVector{UInt32}) 
                isnothing(apex_scan) ? continue : nothing
                
                peak_area[i], new_best_scan[i] = integrateChrom(
                                chrom,
                                apex_scan,
                                linsolve,
                                u2,
                                state,
                                n_pad = n_pad,
                                max_apex_offset = max_apex_offset,
                                isplot = false
                                );
                #println("peak_area $peak_area i $i")
                reset!(state)
            end
        end
    end
    fetch.(tasks)
end


