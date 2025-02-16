"""
    mutable struct Chromatogram{T<:Real, J<:Integer}
        t::Vector{T}
        data::Vector{T}
        max_index::J
    end

Container for chromatogram data during integration.

Holds time points, intensities, and current array size.
"""
mutable struct Chromatogram{T<:Real, J<:Integer}
    t::Vector{T}
    data::Vector{T}
    max_index::J
end

"""
    reset!(state::Chromatogram)

Reset chromatogram state by zeroing arrays and index.
"""
function reset!(state::Chromatogram)
    for i in range(1, state.max_index)
        state.t[i], state.data[i] = zero(eltype(state.t)), zero(eltype(state.data))
    end
    state.max_index = 0
    return 
end


"""
    integrate_chrom(chrom::SubDataFrame, apex_scan::Int64,
                   linsolve::LinearSolve.LinearCache, u2::Vector{Float32},
                   state::Chromatogram; max_apex_offset=2, n_pad::Int64=0,
                   isplot::Bool=false) -> Tuple{Float32, UInt32}

Integrate a single chromatographic peak.

# Process
1. Applies Whittaker-Henderson smoothing
2. Calculates second derivative
3. Finds true apex within allowed offset
4. Determines integration bounds
5. Subtracts baseline
6. Normalizes and fills state
7. Performs trapezoidal integration

# Returns
- Peak area
- Updated apex scan index

#Internal Chromatogram Processing Functions:

- `WHSmooth!`: Apply Whittaker-Henderson smoothing
- `fillU2!`: Calculate second derivatives
- `getApexScan`: Find true apex within allowed offset
- `getIntegrationBounds!`: Determine integration boundaries
- `subtractBaseline!`: Remove linear baseline
- `fillState!`: Normalize and fill chromatogram state
- `integrateTrapezoidal`: Perform trapezoidal integration
"""
function integrate_chrom(chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, 
                                apex_scan::Int64,
                                b::Vector{Float32},
                                u2::Vector{Float32},
                                state::Chromatogram,
                                λ::Float32;
                                max_apex_offset = 2,
                                n_pad::Int64 = 0,
                                isplot::Bool = false)
    
    #########
    #Helper Functions  
    #########
    function WHSmooth!( b::Vector{Float32}, 
                        intensities::AbstractVector{Float32},
                        n_pad::Int64)
        
        n = length(b)

        # Reset b and second derivative 
        @inbounds for i in range(1, n)
            b[i] = zero(Float32)
            u2[i] = zero(Float32)
        end
       
        # Copy data to b
        @inbounds for i in range(1, size(chrom, 1))
            b[i+n_pad] = intensities[i]
        end

        # Create weight matrix of precursor isolated percentage
        max_isolation = maximum(chrom[!, :precursor_fraction_transmitted])
        
        w = max_isolation*ones(Float32, n)
        @inbounds for i in range(1, size(chrom, 1))
            w[i+n_pad] = chrom[!, :precursor_fraction_transmitted][i]
        end
        
        # Get RT spacing 
        rts = chrom[!, :rt]
        start_rt = rts[1]
        last_rt = rts[end]
        default_spacing = size(chrom, 1) < 2 ? 1.0f0 : rts[2] - start_rt

        # Fill RT differences
        x = zeros(Float32, n)

        # left padding
        for i in range(1, n_pad)
            x[i] = start_rt - (default_spacing * ((n_pad - i) + 1))
        end
        # real values
        for i in range(1, size(chrom, 1)-1 )
            x[i + n_pad] = rts[i]
        end
        # right padding
        for i in range(n_pad + size(chrom, 1), n-1)
            x[i] = last_rt + (default_spacing * (i - (n_pad + size(chrom, 1))))
        end

        # normalize by RT width so the same smoothing parameter is about optimal for all cases
        rt_width = last_rt-start_rt
        x = x / rt_width

        # solve
        z = whitsmddw(x, b, w, λ)

        return z
    end

    function fillU2!(
        u2::Vector{Float32},
        u::Vector{Float32})
        #Get second-order descrete derivative 
        u2[1], u2[end] = zero(Float32), zero(Float32)
        @inbounds @fastmath for i in range(2, length(b) - 1)
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
                    stop = i
                    break
                end
            end

            for i in range(stop, N-1)
                if u[i + 1] >= u[i]
                    break
                else
                    stop = i
                end
            end

            #get LH boundary 
            for i in reverse(range(2, start_search))
                if (u2[i] > u2[i - 1]) & (u2[i+1]<u2[i])
                    start = i
                    break
                end
            end

            for i in reverse(range(2, start))
                if u[i - 1] >= u[i]
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
    z = WHSmooth!(
        b,
        chrom[!,:intensity],
        n_pad
    )
    #Second discrete derivative of smoothed data
    fillU2!(
        u2,
        z,
    )
    
    apex_scan = getApexScan(
        apex_scan,
        max_apex_offset,
        chrom[!,:intensity]
    )
    #Integration boundaries based on smoothed second derivative 
    scan_range = getIntegrationBounds!(
        u2,
        z,
        size(chrom, 1),
        apex_scan,
        n_pad
    )

    subtractBaseline!(
        z,
        apex_scan,
        scan_range,
        n_pad
    )

    #File `state` to fit EGH function. Get the inensity, and rt normalization factors 
    norm_factor, start_rt, rt_norm, best_rt = fillState!(
        state,
        z,
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
