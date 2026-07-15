# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
    mutable struct Chromatogram{T<:Real, J<:Integer}
        t::Vector{T}
        data::Vector{T}
        max_index::J
    end

Container for chromatogram data during integration.

Holds the normalized, baseline-subtracted peak trace used by trapezoidal
integration.
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
    debug_save_chromatogram_integration_plot(...)

Save the final debug view for one integrated chromatogram. The figure shows raw
points, WH-smoothed signal, the baseline-subtracted segment that was integrated,
the expected PSM apex, the selected bounds, and precursor transmission.
"""
function debug_save_chromatogram_integration_plot(
    plot_path::AbstractString,
    rt_col::AbstractVector{<:AbstractFloat},
    scan_idx_col::AbstractVector{<:Integer},
    intensity_col::AbstractVector{<:AbstractFloat},
    fraction_col::AbstractVector{<:AbstractFloat},
    wh_smoothed_col::AbstractVector{<:AbstractFloat},
    baseline_subtracted_col::AbstractVector{<:AbstractFloat},
    scan_range::UnitRange{Int},
    apex_scan::Int,
    peak_area::Float32,
    points_integrated::Integer,
    min_fraction_transmitted::Float32,
    status::AbstractString,
    title::AbstractString,
)
    isempty(plot_path) && return nothing
    isempty(rt_col) && return nothing
    ensure_directory_exists(String(plot_path))

    withenv("GKSwstype" => "100", "GKS_WSTYPE" => "100") do
        Plots.gr()

        rt_vals = Float64.(rt_col)
        raw_vals = Float64.(intensity_col)
        wh_vals = Float64.(wh_smoothed_col)
        baseline_vals = Float64.(baseline_subtracted_col)
        fraction_vals = Float64.(fraction_col)
        safe_scan_range = max(first(scan_range), 1):min(last(scan_range), length(rt_vals))

        plot_title = isempty(title) ? "Chromatogram integration debug" : String(title)
        area_text = @sprintf("%.6g", Float64(peak_area))
        p_signal = Plots.plot(
            rt_vals,
            raw_vals,
            seriestype = :scatter,
            alpha = 0.55,
            ms = 4,
            label = "raw",
            xlabel = "retention time",
            ylabel = "intensity",
            title = plot_title * " status=$(status) area=$(area_text) points=$(points_integrated)",
            size = (1100, 800),
            dpi = 150,
        )
        Plots.plot!(p_signal, rt_vals, wh_vals, lw = 2.0, color = :purple, label = "WH smoothed")
        if !isempty(safe_scan_range)
            Plots.plot!(
                p_signal,
                rt_vals[safe_scan_range],
                baseline_vals[safe_scan_range],
                lw = 2.0,
                color = :orange,
                ls = :dash,
                label = "baseline-subtracted",
            )
        end
        Plots.vline!(p_signal, [rt_vals[apex_scan]], color = :red, lw = 1.5, ls = :dash, label = "apex")
        Plots.vline!(
            p_signal,
            [rt_vals[first(safe_scan_range)], rt_vals[last(safe_scan_range)]],
            color = :black,
            lw = 1.0,
            ls = :dot,
            label = "bounds",
        )
        p_fraction = Plots.plot(
            rt_vals,
            fraction_vals,
            seriestype = :scatter,
            alpha = 0.7,
            ms = 4,
            label = "fraction transmitted",
            xlabel = "retention time",
            ylabel = "precursor fraction transmitted",
        )
        Plots.hline!(
            p_fraction,
            [Float64(min_fraction_transmitted)],
            color = :red,
            lw = 1.0,
            ls = :dash,
            label = "minimum",
        )

        plot_obj = Plots.plot(p_signal, p_fraction, layout = (2, 1), size = (1100, 800), dpi = 150)
        Plots.savefig(plot_obj, String(plot_path))
    end

    return nothing
end



"""
WHSmooth! — set up and solve the Whittaker-Henderson smoothing problem.
All I/O goes through ws fields: ws.b (raw signal), ws.z (smoothed output),
ws.u2 (zeroed), ws.w_tmp (weights), ws.x_tmp (RT positions).
"""
function WHSmooth!(ws::WHWorkspace,
                    intensity_col::AbstractVector,
                    fraction_col::AbstractVector,
                    rt_col::AbstractVector,
                    m::Int,
                    n_pad::Int64,
                    min_fraction_transmitted::Float32,
                    λ::Float32)

    b  = ws.b
    u2 = ws.u2
    n  = ws.n_max

    # Reset b and second derivative
    @inbounds for i in range(1, n)
        b[i] = zero(Float32)
        u2[i] = zero(Float32)
    end

    # Copy transmission-corrected data to b
    max_isolation = 0.0f0
    @inbounds for i in range(1, m)
        frac = fraction_col[i]
        if frac >= min_fraction_transmitted
            b[i + n_pad] = intensity_col[i] / frac
            max_isolation = max(max_isolation, frac)
        else
            b[i + n_pad] = 0.0f0
        end
    end

    # Create weight vector using workspace (no allocation)
    w = ws.w_tmp
    @inbounds for i in range(1, n)
        w[i] = max(max_isolation, 1.0f0)
    end
    @inbounds for i in range(1, m)
        frac = fraction_col[i]
        w[i+n_pad] = frac >= min_fraction_transmitted ? frac : 0.0f0
    end

    # Get RT spacing using workspace (no allocation)
    x = ws.x_tmp
    start_rt = rt_col[1]
    last_rt = rt_col[end]
    default_spacing = m < 2 ? 1.0f0 : rt_col[2] - start_rt

    # left padding
    for i in range(1, n_pad)
        x[i] = start_rt - (default_spacing * ((n_pad - i) + 1))
    end
    # real values
    for i in range(1, m)
        x[i + n_pad] = rt_col[i]
    end
    # right padding
    for i in range(n_pad+m+1, n)
        x[i] = last_rt + (default_spacing * (i - (n_pad + m)))
    end

    # normalize by RT width in-place
    rt_width = last_rt-start_rt
    @inbounds for i in range(1, n)
        x[i] /= rt_width
    end

    active_length = m + (2*n_pad)

    if m <= 1 || active_length < 3
        # No smoothing possible — copy raw data into ws.z
        @inbounds for i in 1:active_length
            ws.z[i] = b[i]
        end
        return nothing
    end

    # Solves in-place; result lives in ws.z[1:active_length]
    whitsmddw!(ws, x, b, w, active_length, λ)
    return nothing
end

function fillU2!(
    u2::Vector{Float32},
    u::Vector{Float32},
    t::Vector{Float32},
    n::Int
)
    u2[1] = 0.0f0
    u2[n] = 0.0f0
    @inbounds @fastmath for i in 2:(n - 1)
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
    intensities::Vector{Float32},
    N::Int)

    apex_idx = apex_scan + n_pad

    # Walk right: climb while next scan is strictly higher
    right_apex = apex_idx
    while right_apex < N && intensities[right_apex + 1] > intensities[right_apex]
        right_apex += 1
    end

    # Walk left: climb while next scan is strictly higher
    left_apex = apex_idx
    while left_apex > 1 && intensities[left_apex - 1] > intensities[left_apex]
        left_apex -= 1
    end

    # Pick whichever local maximum is greater
    best = intensities[right_apex] >= intensities[left_apex] ? right_apex : left_apex
    return best - n_pad
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

    # initialise search bounds (clamp to valid padded range)
    start = max(apex_padded - 1, pad_start)
    stop  = min(apex_padded + 1, pad_end)

    # ──────────────── search to the right (RH boundary) ────────────────
    # 1. advance to first local maximum of u2  (peak of d²/dt² < 0)
    @inbounds @fastmath for i in stop:pad_end-1
        if (u2[i-1] < u2[i]) && (u2[i+1] < u2[i])
            stop = i
            break
        end
    end
    # 2. keep going while the intensity is still near its running minimum
    #    (tolerate small noise bumps up to 15% above the min seen so far)
    #    Always snap the boundary to the running min position.
    running_min = u[stop]
    running_min_idx = stop
    @inbounds @fastmath for i in stop:pad_end-1
        if u[i] < running_min
            running_min = u[i]
            running_min_idx = i
        end
        if u[i+1] > running_min * 1.15f0
            break
        end
    end
    stop = running_min_idx

    # ──────────────── search to the left  (LH boundary) ────────────────
    # 1. first local maximum of u2 when scanning left
    @inbounds @fastmath for i in reverse(pad_start+1:start)
        if (u2[i] > u2[i-1]) && (u2[i+1] < u2[i])
            start = i
            break
        end
    end
    # 2. keep going left while intensity is near its running minimum
    running_min = u[start]
    running_min_idx = start
    @inbounds @fastmath for i in reverse(pad_start+1:start)
        if u[i] < running_min
            running_min = u[i]
            running_min_idx = i
        end
        if u[i-1] > running_min * 1.15f0
            break
        end
    end
    start = running_min_idx

    # convert from padded indices back to 1…N
    start_mid = max(start - n_pad, 1)
    stop_mid  = min(stop  - n_pad, N)

    return start_mid:stop_mid
end

function ensureMinimumScanRange(scan_range::UnitRange{Int}, apex_scan::Int, N::Int; min_points::Int = 3)
    N <= 0 && return scan_range
    target_points = min(min_points, N)
    start = clamp(first(scan_range), 1, N)
    stop = clamp(last(scan_range), 1, N)
    start, stop = minmax(start, stop)

    if apex_scan < start
        start = apex_scan
    elseif apex_scan > stop
        stop = apex_scan
    end
    start = clamp(start, 1, N)
    stop = clamp(stop, 1, N)

    while stop - start + 1 < target_points
        left_span = apex_scan - start
        right_span = stop - apex_scan
        if start > 1 && (left_span <= right_span || stop == N)
            start -= 1
        elseif stop < N
            stop += 1
        elseif start > 1
            start -= 1
        else
            break
        end
    end

    return start:stop
end

function getForcedBoundaryRange(
    scan_idx_col::AbstractVector{<:Integer},
    start_scan::UInt32,
    stop_scan::UInt32,
)
    (start_scan == 0 || stop_scan == 0) && return nothing
    start_idx = find_nearest_scan(scan_idx_col, start_scan)
    stop_idx = find_nearest_scan(scan_idx_col, stop_scan)
    start_idx, stop_idx = minmax(start_idx, stop_idx)
    return start_idx:stop_idx
end


function fillState!(state::Chromatogram,
                    u::Vector{Float32},
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
        state.t[i] = (rt[n] - start_rt)/rt_width
        state.data[i] = u[n+n_pad]/norm_factor
    end

    state.max_index = stop - start + 1
    best_rt = Float32((best_rt - start_rt)/rt_width)
    return norm_factor, start_rt, rt_width, best_rt
end

function subtractBaseline!(
    x::Vector{Float32},  # time (or x-axis values)
    u::Vector{Float32},  # smoothed signal values
    _apex_scan::Int,     # retained for call-site compatibility
    scan_range::UnitRange{Int},
    n_pad::Int
)
    # Apply pad offset
    scan_start = first(scan_range) + n_pad
    scan_stop  = last(scan_range) + n_pad

    lmin = u[scan_start]
    rmin = u[scan_stop]
    x_left = x[scan_start]
    x_right = x[scan_stop]
    dx = x_right - x_left
    if dx == 0
        u[scan_start] = 0.0f0
        return nothing
    end

    slope = (rmin - lmin) / dx

    @inbounds @fastmath for i in scan_start:scan_stop
        xi = x[i]
        baseline = lmin + (xi - x_left) * slope
        u[i] = max(0, u[i] - baseline)
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


"""
    integrate_chrom(chrom::SubDataFrame, apex_scan::Int64,
                   ws::WHWorkspace, state::Chromatogram,
                   avg_cycle_time::Float32, λ::Float32;
                   n_pad::Int64=0,
                   isplot::Bool=false) -> Tuple{Float32, UInt32, Int, UInt32, UInt32}

Integrate a single chromatographic peak.

All scratch arrays (b, u2, z, x) live inside `ws` (one WHWorkspace per thread).
No SubArray views escape into the caller — every downstream function operates
on concrete `Vector{Float32}` fields of `ws`, eliminating GC-root / view-lifetime issues.

# Process
1. Applies transmission-corrected Whittaker-Henderson smoothing
2. Calculates the second derivative of the smoothed trace
3. Moves from the expected PSM apex to the adjacent local maximum when needed
4. Determines bounds from second-derivative extrema and nearby valleys
5. Ensures at least three points are covered when possible
6. Subtracts a linear baseline through the selected endpoint values
7. Normalizes the selected segment and performs trapezoidal integration

# Returns
- Peak area
- Updated apex scan index
- Number of points integrated
- Integration start scan index
- Integration stop scan index

# Internal Chromatogram Processing Functions

- `WHSmooth!`: Apply Whittaker-Henderson smoothing (result in ws.z)
- `fillU2!`: Calculate second derivatives (result in ws.u2)
- `getApexScan`: Move to the local maximum reachable from the expected apex
- `getIntegrationBounds!`: Determine second-derivative/valley boundaries
- `subtractBaseline!`: Remove a linear baseline anchored at selected endpoints
- `fillState!`: Normalize and fill chromatogram state
- `integrateTrapezoidal`: Perform trapezoidal integration
"""
function integrate_chrom(rt_col::AbstractVector{<:AbstractFloat},
                                scan_idx_col::AbstractVector{<:Integer},
                                intensity_col::AbstractVector{<:AbstractFloat},
                                fraction_col::AbstractVector{<:AbstractFloat},
                                apex_scan::Int64,
                                ws::WHWorkspace,
                                state::Chromatogram,
                                avg_cycle_time::Float32,
                                λ::Float32;
                                min_fraction_transmitted::Float32 = 0.0f0,
                                n_pad::Int64 = 0,
                                isplot::Bool = false,
                                debug_plot_path::Union{Nothing, AbstractString} = nothing,
                                debug_plot_title::AbstractString = "",
                                debug_plot_data::Union{Nothing, Base.RefValue} = nothing,
                                debug_apex_scan::Union{Nothing, Int64} = nothing,
                                forced_boundary_start_scan::UInt32 = UInt32(0),
                                forced_boundary_stop_scan::UInt32 = UInt32(0))

    m = length(rt_col)
    # Helper functions (WHSmooth!, fillU2!, getApexScan, getIntegrationBounds!,
    # ensureMinimumScanRange, getForcedBoundaryRange, fillState!, subtractBaseline!,
    # integrateTrapezoidal) are defined at module scope above to avoid closure boxing.

    # Whittaker-Henderson smoothing. `fraction_col` corrects transmitted
    # precursor signal and controls observation weights.
    WHSmooth!(ws, intensity_col, fraction_col, rt_col, m, n_pad, min_fraction_transmitted, λ)

    # Local aliases to workspace vectors (concrete Vector{Float32}, no SubArray views)
    n_active = m + 2*n_pad
    z  = ws.z
    x  = ws.x_tmp
    u2 = ws.u2

    # Second discrete derivative of smoothed data.
    fillU2!(u2, z, x, n_active)

    apex_scan = clamp(getApexScan(
        apex_scan,
        n_pad,
        z,
        n_active
    ), 1, m)

    # Integration boundaries based on the smoothed second derivative.
    scan_range = getIntegrationBounds!(
        u2,
        z,
        m,
        apex_scan,
        n_pad
    )
    fallback_scan_range = ensureMinimumScanRange(scan_range, apex_scan, m)
    forced_boundary_range = getForcedBoundaryRange(
        scan_idx_col,
        forced_boundary_start_scan,
        forced_boundary_stop_scan,
    )
    debug_display_apex_scan = debug_apex_scan === nothing ?
        apex_scan :
        clamp(Int(debug_apex_scan), 1, m)
    scan_range = forced_boundary_range === nothing ? fallback_scan_range : forced_boundary_range
    scan_range = ensureMinimumScanRange(scan_range, apex_scan, m)

    debug_enabled = debug_plot_path !== nothing || debug_plot_data !== nothing
    wh_smoothed_debug = debug_enabled ? copy(@view z[(n_pad + 1):(n_pad + m)]) : Float32[]

    subtractBaseline!(
        x,
        z,
        apex_scan,
        scan_range,
        n_pad
    )

    baseline_subtracted_debug = debug_enabled ? copy(@view z[(n_pad + 1):(n_pad + m)]) : Float32[]

    # Fill integration state and get the intensity/RT normalization factors.
    # Guard: if smoothed apex is zero or NaN after baseline subtraction, skip integration
    _apex_val = z[apex_scan + n_pad]
    if _apex_val <= 0f0 || isnan(_apex_val)
        if debug_enabled
            status = "skipped_zero_apex"
            debug_plot_data !== nothing && (debug_plot_data[] = (
                scan_range = scan_range,
                apex_scan_idx = UInt32(scan_idx_col[debug_display_apex_scan]),
                wh_smoothed = wh_smoothed_debug,
                baseline_subtracted = baseline_subtracted_debug,
                peak_area = 0.0f0,
                points_integrated = UInt32(0),
                status = status,
            ))
            if debug_plot_path !== nothing
                debug_save_chromatogram_integration_plot(
                    debug_plot_path,
                    rt_col,
                    scan_idx_col,
                    intensity_col,
                    fraction_col,
                    wh_smoothed_debug,
                    baseline_subtracted_debug,
                    scan_range,
                    debug_display_apex_scan,
                    0.0f0,
                    0,
                    min_fraction_transmitted,
                    status,
                    debug_plot_title,
                )
            end
        end
        return 0f0, scan_idx_col[apex_scan], Int(0), UInt32(0), UInt32(0)
    end

    norm_factor, start_rt, rt_norm, best_rt = fillState!(
        state,
        z,
        rt_col,
        first(scan_range),
        last(scan_range),
        apex_scan,
        n_pad
    )

    if isplot
    #if chrom.precursor_idx[1] == 2098
        mi = state.max_index
        start = max(apex_scan - 18, 1)
        stop = min(apex_scan + 18, length(rt_col))
        plot(rt_col[start:stop], intensity_col[start:stop], seriestype=:scatter, alpha = 0.5, show = true, label = "raw")
        vline!([rt_col[first(scan_range)], rt_col[last(scan_range)]], label = nothing)
        plot!(state.t[1:mi].*rt_norm .+ start_rt, norm_factor.*state.data[1:mi], seriestype=:scatter, alpha = 1.0, show = true, label = "smooth", color = "purple")
        #savefig("/Users/dennisgoldfarb/Downloads/" * string(chrom.precursor_idx[1]) * " " * string(chrom.intensity[1]) * ".png")
        #xbins = LinRange(state.t[1]-0.5, state.t[state.max_index]+0.5, 100)
        #plot!(xbins.*rt_norm .+ start_rt, [norm_factor*F(state, x) for x in xbins])
        #plot!(chrom.rt[start:stop], u2[start+n_pad:stop+n_pad])
        #hline!([norm_factor*0.95])
    end

    raw_trap = integrateTrapezoidal(state, avg_cycle_time)
    trapezoid_area = rt_norm * norm_factor * raw_trap

    # Log NaN diagnosis (limited to avoid spam)
    if isnan(trapezoid_area) || isnan(raw_trap)
        z_has_nan = any(isnan, @view z[1:n_active])
        z_all_zero = all(==(0f0), @view z[1:n_active])
        @warn "[NaN area] norm_factor=$norm_factor rt_norm=$rt_norm raw_trap=$raw_trap max_idx=$(state.max_index) m=$m n_pad=$n_pad z_all_zero=$z_all_zero z_has_nan=$z_has_nan scan_range=$scan_range" maxlog=10
    end

    # Count points within the full width at 20% maximum of the smoothed signal
    # Restrict evaluation to the scan range actually used for integration
    min_abundance = 0.1f0 * norm_factor
    scan_start = first(scan_range) + n_pad
    scan_stop  = last(scan_range) + n_pad
    num_points_integrated = count(i -> z[i] >= min_abundance, scan_start:scan_stop)

    if debug_enabled
        status = "integrated"
        debug_plot_data !== nothing && (debug_plot_data[] = (
            scan_range = scan_range,
            apex_scan_idx = UInt32(scan_idx_col[debug_display_apex_scan]),
            wh_smoothed = wh_smoothed_debug,
            baseline_subtracted = baseline_subtracted_debug,
            peak_area = trapezoid_area,
            points_integrated = UInt32(num_points_integrated),
            status = status,
        ))
        if debug_plot_path !== nothing
            debug_save_chromatogram_integration_plot(
                debug_plot_path,
                rt_col,
                scan_idx_col,
                intensity_col,
                fraction_col,
                wh_smoothed_debug,
                baseline_subtracted_debug,
                scan_range,
                debug_display_apex_scan,
                trapezoid_area,
                num_points_integrated,
                min_fraction_transmitted,
                status,
                debug_plot_title,
            )
        end
    end

    #trapezoid_area = 0.0f0
    return (
        trapezoid_area,
        scan_idx_col[apex_scan],
        num_points_integrated,
        UInt32(scan_idx_col[first(scan_range)]),
        UInt32(scan_idx_col[last(scan_range)]),
    )
end

function integrate_chrom(chrom::SubDataFrame,
                                apex_scan::Int64,
                                ws::WHWorkspace,
                                state::Chromatogram,
                                avg_cycle_time::Float32,
                                λ::Float32;
                                min_fraction_transmitted::Float32 = 0.0f0,
                                n_pad::Int64 = 0,
                                isplot::Bool = false,
                                debug_plot_path::Union{Nothing, AbstractString} = nothing,
                                debug_plot_title::AbstractString = "",
                                debug_plot_data::Union{Nothing, Base.RefValue} = nothing,
                                debug_apex_scan::Union{Nothing, Int64} = nothing,
                                forced_boundary_start_scan::UInt32 = UInt32(0),
                                forced_boundary_stop_scan::UInt32 = UInt32(0))
    return integrate_chrom(
        chrom[!, :rt],
        chrom[!, :scan_idx],
        chrom[!, :intensity],
        chrom[!, :precursor_fraction_transmitted],
        apex_scan,
        ws,
        state,
        avg_cycle_time,
        λ;
        min_fraction_transmitted = min_fraction_transmitted,
        n_pad = n_pad,
        isplot = isplot,
        debug_plot_path = debug_plot_path,
        debug_plot_title = debug_plot_title,
        debug_plot_data = debug_plot_data,
        debug_apex_scan = debug_apex_scan,
        forced_boundary_start_scan = forced_boundary_start_scan,
        forced_boundary_stop_scan = forced_boundary_stop_scan,
    )
end
