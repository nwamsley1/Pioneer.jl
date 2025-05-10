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
                                avg_cycle_time::Float32,
                                λ::Float32;
                                max_apex_offset = 2,
                                n_pad::Int64 = 0,
                                isplot::Bool = false)
    
    #########
    #Helper Functions  
    #########
    function WHSmooth!( b::Vector{Float32}, 
                        chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}},
                        n_pad::Int64,
                        λ::Float32)
        
        n = length(b)
        m = size(chrom, 1)

        # Reset b and second derivative 
        @inbounds for i in range(1, n)
            b[i] = zero(Float32)
            u2[i] = zero(Float32)
        end
       
        # Copy data to b
        @inbounds for i in range(1, m)
            b[i+n_pad] = chrom[!, :intensity][i]
        end

        # Create weight matrix of precursor isolated percentage
        max_isolation = maximum(chrom[!, :precursor_fraction_transmitted])
        
        w = max_isolation*ones(Float32, n)
        @inbounds for i in range(1, m)
            w[i+n_pad] = chrom[!, :precursor_fraction_transmitted][i]
        end
        
        # Get RT spacing 
        rts = chrom[!, :rt]
        start_rt = rts[1]
        last_rt = rts[end]
        default_spacing = m < 2 ? 1.0f0 : rts[2] - start_rt

        # Fill RT differences
        x = zeros(Float32, n)

        # left padding
        for i in range(1, n_pad)
            x[i] = start_rt - (default_spacing * ((n_pad - i) + 1))
        end
        # real values
        for i in range(1, m)
            x[i + n_pad] = rts[i]
        end
        # right padding
        for i in range(n_pad+m+1, n)
            x[i] = last_rt + (default_spacing * (i - (n_pad + m)))
        end

        # normalize by RT width so the same smoothing parameter is about optimal for all cases
        rt_width = last_rt-start_rt
        x = x / rt_width

        if m <= 1
            return b, x
        end

        active_length = m + (2*n_pad)
        z = whitsmddw(x, b, w, active_length, λ)

        return z, x
    end

    function fillU2!(
        u2::Vector{Float32},
        u::Vector{Float32},
        t::Vector{Float32}  # time points corresponding to u
    )
        u2[1], u2[end] = 0.0f0, 0.0f0
        @inbounds @fastmath for i in 2:(length(u) - 1)
            dt1 = t[i] - t[i - 1]
            dt2 = t[i + 1] - t[i]
            dt_total = t[i + 1] - t[i - 1]
            a = (u[i + 1] - u[i]) / dt2
            b = (u[i] - u[i - 1]) / dt1
            u2[i] = 2f0 / dt_total * (a - b)
        end
    end

    function getApexScan(
        apex_scan::Int64,
        n_pad::Int64,
        intensities::AbstractVector{<:AbstractFloat})

        N = length(intensities)
        max_intensity = zero(Float32)
        apex_scan += n_pad
        #for i in range(n_pad+1, length(intensities))
        for i in range(max(1, apex_scan - 2), min(N, apex_scan+2))
            if intensities[i] > max_intensity#intensities[apex_scan]
                apex_scan = i
                max_intensity = intensities[i]
            end
        end
        return apex_scan - n_pad
    end

    """
    getIntegrationBounds!(u2, u, N, apex_scan, n_pad) -> UnitRange

    Find the start/stop scan indices of a chromatographic peak whose apex
    (in the *unpadded* region) is at `apex_scan`.

    * `u2` - second-derivative-like array (length = n_pad + N + n_pad)
    * `u`  - intensity array (same length as `u2`)
    * `N`  - length of the **central** window (without padding)
    * `apex_scan` - 1-based apex position inside the central window
    * `n_pad` - number of padded samples on each side

    Returns a `UnitRange{Int}` with indices **in the un-padded domain** (`1:N`).
    """
    function getIntegrationBounds!(u2::Vector{Float32},
                                u::Vector{Float32},
                                N::Int,
                                apex_scan::Int,
                                n_pad::Int)::UnitRange{Int}

        # indices in the *padded* coordinate system
        pad_start   = n_pad + 1                   # first index of the real window
        pad_end     = n_pad + N
        apex_padded = n_pad + apex_scan           # apex index in padded coords
        
        #return pad_start:pad_end

        # initialise search bounds
        start = apex_padded - 1
        stop  = apex_padded + 1

        # ──────────────── search to the right (RH boundary) ────────────────
        # 1. advance to first local maximum of u2  (peak of d²/dt² < 0)
        @inbounds @fastmath for i in stop:pad_end-1
            if (u2[i-1] < u2[i]) && (u2[i+1] < u2[i])
                stop = i
                break
            end
        end
        # 2. keep going while the intensity still decreases
        @inbounds @fastmath for i in stop:pad_end-1
            if u[i+1] >= u[i]
                break
            else
                stop = i+1
            end
        end

        # ──────────────── search to the left  (LH boundary) ────────────────
        # 1. first local maximum of u2 when scanning left
        @inbounds @fastmath for i in reverse(pad_start+1:start)
            if (u2[i] > u2[i-1]) && (u2[i+1] < u2[i])
                start = i
                break
            end
        end
        # 2. keep going left while intensity decreases
        @inbounds @fastmath for i in reverse(pad_start+1:start)
            if u[i-1] >= u[i]
                break
            else
                start = i-1
            end
        end

        # convert from padded indices back to 1…N
        start_mid = max(start - n_pad, 1)
        stop_mid  = min(stop  - n_pad, N)

        return start_mid:stop_mid
    end


    function fillState!(state::Chromatogram,
                        u::AbstractVector{Float32},
                        rt::AbstractVector{Float32},
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

    function subtractBaseline!(
        x::Vector{Float32},  # time (or x-axis values)
        u::Vector{Float32},  # smoothed signal values
        apex_scan::Int,      # peak apex index
        scan_range::UnitRange{Int},
        n_pad::Int
    )
        # Apply pad offset
        apex_idx  = apex_scan + n_pad
        scan_start = first(scan_range) + n_pad
        scan_stop  = last(scan_range) + n_pad
    
        # Find left baseline: minimum value between scan_start and apex_idx
        lmin, li = typemax(Float32), scan_start
        @inbounds @fastmath for i in scan_start:apex_idx
            if u[i] < lmin
                lmin = u[i]
                li = i
            end
        end
    
        # Find right baseline: minimum value between apex_idx and scan_stop
        rmin, ri = typemax(Float32), scan_stop
        @inbounds @fastmath for i in apex_idx:scan_stop
            if u[i] < rmin
                rmin = u[i]
                ri = i
            end
        end
    
        # Handle special case where li == apex_idx or ri == apex_idx
        if li == apex_idx
            lmin = rmin
        end
        if ri == apex_idx
            rmin = lmin
        end

        # Calculate slope based on actual x values, not index distance
        x_left = x[li]
        x_right = x[ri]
        dx = x_right - x_left
        if dx == 0
            return nothing
        end

        slope = (rmin - lmin) / dx
    
        # Subtract interpolated baseline
        @inbounds @fastmath for i in scan_start:scan_stop
            xi = x[i]
            baseline = lmin + (xi - x_left) * slope
            u[i] -= baseline
        end
    
        return nothing
    end


    function integrateTrapezoidal(state::Chromatogram, avg_cycle_time::Float32)
        if state.max_index == 1
            # Special case: 1 point only, treat like a triangle on each side
            height = state.data[1]
            return avg_cycle_time * height
        elseif state.max_index >= 2
            retval = 0.0f0
            retval += avg_cycle_time * (state.data[1] + state.data[state.max_index]) # triangle area on each side
            @inbounds @fastmath for i in 1:(state.max_index - 1)
                dt = state.t[i + 1] - state.t[i]
                retval += dt * (state.data[i] + state.data[i + 1])
            end
            return 0.5f0 * retval
        else
            # No points, no area
            return 0.0f0
        end
    end

    #Whittaker Henderson Smoothing
    z, x = WHSmooth!(
        b,
        chrom,
        n_pad,
        λ
    )

    #Second discrete derivative of smoothed data
    fillU2!(
        u2,
        z,
        x
    )
    
    apex_scan = getApexScan(
        apex_scan,
        n_pad,
        z
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
        x,
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
    #if chrom.precursor_idx[1] == 4910216
        mi = state.max_index
        start = max(apex_scan - 18, 1)
        stop = min(apex_scan + 18, length(chrom.rt))
        plot(chrom.rt[start:stop], chrom.intensity[start:stop], seriestype=:scatter, alpha = 0.5, show = true, label = "raw")
        vline!([chrom.rt[first(scan_range)], chrom.rt[last(scan_range)]], label = nothing)
        plot!(state.t[1:mi].*rt_norm .+ start_rt, norm_factor.*state.data[1:mi], seriestype=:scatter, alpha = 1.0, show = true, label = "smooth", color = "purple")
        #savefig("/Users/dennisgoldfarb/Downloads/" * string(chrom.precursor_idx[1]) * " " * string(chrom.intensity[1]) * ".png")
        #xbins = LinRange(state.t[1]-0.5, state.t[state.max_index]+0.5, 100)
        #plot!(xbins.*rt_norm .+ start_rt, [norm_factor*F(state, x) for x in xbins])
        #plot!(chrom.rt[start:stop], u2[start+n_pad:stop+n_pad])
        #hline!([norm_factor*0.95])
    end

    trapezoid_area = rt_norm * norm_factor * integrateTrapezoidal(state, avg_cycle_time)

    #trapezoid_area = 0.0f0
    return trapezoid_area, chrom[!,:scan_idx][apex_scan]
end
