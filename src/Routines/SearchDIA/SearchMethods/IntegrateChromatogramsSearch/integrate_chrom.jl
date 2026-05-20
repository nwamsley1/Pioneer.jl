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

const BOUNDARY_CANDIDATE_CATEGORY_LABELS = (
    "fallback",
    "valley_combination",
)
const BOUNDARY_OUTSIDE_FEATURE_SCAN_WINDOW = 4
const MAX_BOUNDARY_VALLEY_ENDPOINTS_PER_SIDE = 3

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
    mainsearch_evidence_range::Union{Nothing,UnitRange{Int}} = nothing,
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
        if mainsearch_evidence_range !== nothing
            evidence_range = max(first(mainsearch_evidence_range), 1):min(last(mainsearch_evidence_range), length(rt_vals))
            if !isempty(evidence_range)
                Plots.vline!(
                    p_signal,
                    unique([rt_vals[first(evidence_range)], rt_vals[last(evidence_range)]]),
                    color = :green,
                    lw = 1.2,
                    ls = :dashdot,
                    label = "main search q interval",
                )
            end
        end

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

    debug_l1("chromatogram_debug_plot path=$(plot_path)")
    return nothing
end


"""
    integrate_chrom(chrom::SubDataFrame, apex_scan::Int64,
                   ws::WHWorkspace, state::Chromatogram,
                   avg_cycle_time::Float32, λ::Float32;
                   n_pad::Int64=0,
                   isplot::Bool=false) -> Tuple{Float32, UInt32, Int}

Integrate a single chromatographic peak.

All scratch arrays (b, u2, z, x) live inside `ws` (one WHWorkspace per thread).
No SubArray views escape into the caller — every downstream function operates
on concrete `Vector{Float32}` fields of `ws`, eliminating GC-root / view-lifetime issues.

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
- Number of points integrated

#Internal Chromatogram Processing Functions:

- `WHSmooth!`: Apply Whittaker-Henderson smoothing (result in ws.z)
- `fillU2!`: Calculate second derivatives (result in ws.u2)
- `getApexScan`: Find true apex within allowed offset
- `getIntegrationBounds!`: Determine integration boundaries
- `subtractBaseline!`: Remove linear baseline
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
                                mainsearch_1pct_start_scan::UInt32 = UInt32(0),
                                mainsearch_1pct_stop_scan::UInt32 = UInt32(0),
                                rt_to_irt_model::RtConversionModel = IdentityModel(),
                                forced_boundary_start_scan::UInt32 = UInt32(0),
                                forced_boundary_stop_scan::UInt32 = UInt32(0),
                                boundary_candidate_data::Union{Nothing, Base.RefValue} = nothing)

    m = length(rt_col)

    #########
    #Helper Functions
    #########

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

    function pushBoundaryCandidate!(
        candidates::Vector,
        seen::Set{Tuple{Int,Int}},
        start::Int,
        stop::Int,
        apex_scan::Int,
        N::Int,
        category::AbstractString,
    )
        start = clamp(start, 1, N)
        stop = clamp(stop, 1, N)
        start <= stop || return nothing
        start <= apex_scan <= stop || return nothing
        scan_range = ensureMinimumScanRange(start:stop, apex_scan, N)
        key = (first(scan_range), last(scan_range))
        if key ∉ seen
            push!(seen, key)
            push!(candidates, (scan_range = scan_range, category = String(category)))
        end
        return nothing
    end

    function isLocalMinimumEndpoint(
        u::Vector{Float32},
        idx::Int,
        N::Int,
        n_pad::Int,
    )
        padded_idx = idx + n_pad
        val = u[padded_idx]
        left_ok = idx == 1 || val <= u[padded_idx - 1]
        right_ok = idx == N || val <= u[padded_idx + 1]
        return left_ok && right_ok
    end

    function collectValleyEndpoints(
        u::Vector{Float32},
        scan_range::UnitRange{Int},
        apex_scan::Int,
        N::Int,
        n_pad::Int;
        left_side::Bool,
    )
        valleys = Int[]
        start = left_side ? first(scan_range) : max(apex_scan + 1, first(scan_range))
        stop = left_side ? min(apex_scan - 1, last(scan_range)) : last(scan_range)
        start <= stop || return valleys

        if left_side
            for i in stop:-1:start
                if isLocalMinimumEndpoint(u, i, N, n_pad)
                    push!(valleys, i)
                    length(valleys) >= MAX_BOUNDARY_VALLEY_ENDPOINTS_PER_SIDE && break
                end
            end
        else
            for i in start:stop
                if isLocalMinimumEndpoint(u, i, N, n_pad)
                    push!(valleys, i)
                    length(valleys) >= MAX_BOUNDARY_VALLEY_ENDPOINTS_PER_SIDE && break
                end
            end
        end
        return valleys
    end

    function addValleyCombinationCandidates!(
        candidates::Vector,
        seen::Set{Tuple{Int,Int}},
        u::Vector{Float32},
        current_range::UnitRange{Int},
        N::Int,
        apex_scan::Int,
        n_pad::Int,
    )
        left_valleys = collectValleyEndpoints(
            u,
            current_range,
            apex_scan,
            N,
            n_pad;
            left_side = true,
        )
        right_valleys = collectValleyEndpoints(
            u,
            current_range,
            apex_scan,
            N,
            n_pad;
            left_side = false,
        )

        for left in left_valleys, right in right_valleys
            pushBoundaryCandidate!(
                candidates,
                seen,
                left,
                right,
                apex_scan,
                N,
                "valley_combination",
            )
        end
        return nothing
    end

    function generateBoundaryCandidates(
        current_range::UnitRange{Int},
        u::Vector{Float32},
        N::Int,
        apex_scan::Int,
        n_pad::Int,
    )
        candidates = NamedTuple{(:scan_range, :category), Tuple{UnitRange{Int}, String}}[]
        seen = Set{Tuple{Int,Int}}()

        pushBoundaryCandidate!(
            candidates,
            seen,
            first(current_range),
            last(current_range),
            apex_scan,
            N,
            "fallback",
        )
        addValleyCombinationCandidates!(candidates, seen, u, 1:N, N, apex_scan, n_pad)

        return candidates
    end

    function localMinimumNear(u::Vector{Float32}, idx::Int, n_pad::Int, N::Int; radius::Int = 1)
        start = max(idx - radius, 1) + n_pad
        stop = min(idx + radius, N) + n_pad
        local_min = typemax(Float32)
        @inbounds for i in start:stop
            local_min = min(local_min, u[i])
        end
        return local_min
    end

    function secondaryPeakPenalty(u::Vector{Float32}, scan_range::UnitRange{Int}, apex_scan::Int, n_pad::Int, apex_height::Float32)
        penalty = 0.0f0

        min_seen = apex_height
        @inbounds for i in (apex_scan - 1):-1:first(scan_range)
            val = u[i + n_pad]
            if val > min_seen
                penalty = max(penalty, (val - min_seen) / apex_height)
            else
                min_seen = val
            end
        end

        min_seen = apex_height
        @inbounds for i in (apex_scan + 1):last(scan_range)
            val = u[i + n_pad]
            if val > min_seen
                penalty = max(penalty, (val - min_seen) / apex_height)
            else
                min_seen = val
            end
        end

        return penalty
    end

    function outsideBoundarySignalFeatures(
        u::Vector{Float32},
        scan_range::UnitRange{Int},
        N::Int,
        n_pad::Int,
        apex_height::Float32;
        left_side::Bool,
        window::Int = BOUNDARY_OUTSIDE_FEATURE_SCAN_WINDOW,
    )
        if left_side
            stop = first(scan_range) - 1
            start = max(1, stop - window + 1)
            boundary_idx = first(scan_range)
        else
            start = last(scan_range) + 1
            stop = min(N, start + window - 1)
            boundary_idx = last(scan_range)
        end
        start <= stop || return (0.0f0, 0.0f0, 0.0f0)

        signal_sum = 0.0f0
        outside_peak = 0.0f0
        n = 0
        @inbounds for i in start:stop
            val = max(0.0f0, u[i + n_pad])
            signal_sum += val
            outside_peak = max(outside_peak, val)
            n += 1
        end

        boundary_val = max(0.0f0, u[boundary_idx + n_pad])
        excluded_signal_fraction = signal_sum / (apex_height * Float32(n))
        boundary_recovery_fraction = max(0.0f0, outside_peak - boundary_val) / apex_height
        outside_peak_fraction = outside_peak / apex_height
        return (
            excluded_signal_fraction,
            boundary_recovery_fraction,
            outside_peak_fraction,
        )
    end

    function internalDipRecoveryScore(
        u::Vector{Float32},
        scan_range::UnitRange{Int},
        apex_scan::Int,
        n_pad::Int,
        apex_height::Float32,
    )
        score = 0.0f0

        min_seen = apex_height
        min_idx = apex_scan
        @inbounds for i in (apex_scan - 1):-1:first(scan_range)
            val = max(0.0f0, u[i + n_pad])
            if val < min_seen
                min_seen = val
                min_idx = i
            else
                distance = max(indexIrtDistance(i, min_idx), eps(Float32))
                score = max(score, (val - min_seen) / (apex_height * sqrt(distance)))
            end
        end

        min_seen = apex_height
        min_idx = apex_scan
        @inbounds for i in (apex_scan + 1):last(scan_range)
            val = max(0.0f0, u[i + n_pad])
            if val < min_seen
                min_seen = val
                min_idx = i
            else
                distance = max(indexIrtDistance(i, min_idx), eps(Float32))
                score = max(score, (val - min_seen) / (apex_height * sqrt(distance)))
            end
        end

        return score
    end

    function edgeValleyLog2Ratio(
        u::Vector{Float32},
        edge_idx::Int,
        apex_scan::Int,
        N::Int,
        n_pad::Int,
        apex_height::Float32;
        left_side::Bool,
    )
        apex_height <= 0.0f0 && return 0.0f0

        start = left_side ? edge_idx + 1 : apex_scan + 1
        stop = left_side ? apex_scan - 1 : edge_idx - 1
        start <= stop || return 0.0f0

        internal_valley = typemax(Float32)
        found_internal_valley = false
        @inbounds for i in start:stop
            if isLocalMinimumEndpoint(u, i, N, n_pad)
                internal_valley = min(internal_valley, max(0.0f0, u[i + n_pad]))
                found_internal_valley = true
            end
        end
        found_internal_valley || return 0.0f0

        edge_fraction = max(0.0f0, u[edge_idx + n_pad]) / apex_height
        internal_fraction = internal_valley / apex_height
        floor_fraction = 1.0f-6
        ratio = log2((edge_fraction + floor_fraction) / (internal_fraction + floor_fraction))
        return isfinite(ratio) ? Float32(ratio) : 0.0f0
    end

    function includedNonApexMaxFeatures(
        u::Vector{Float32},
        scan_range::UnitRange{Int},
        apex_scan::Int,
        n_pad::Int,
        apex_height::Float32,
    )
        apex_height <= 0.0f0 && return (0.0f0, 0.0f0)

        max_nonapex = 0.0f0
        max_idx = apex_scan
        found_nonapex = false
        @inbounds for i in scan_range
            i == apex_scan && continue
            val = max(0.0f0, u[i + n_pad])
            val <= 0.0f0 && continue
            distance = indexIrtDistance(apex_scan, i)
            if !found_nonapex ||
                    val > max_nonapex * (1.0f0 + 1.0f-6) ||
                    (abs(val - max_nonapex) <= max(max_nonapex, eps(Float32)) * 1.0f-6 &&
                        distance < indexIrtDistance(apex_scan, max_idx))
                max_nonapex = val
                max_idx = i
                found_nonapex = true
            end
        end
        found_nonapex || return (0.0f0, 0.0f0)

        floor_fraction = 1.0f-6
        ratio = log2(max(max_nonapex / apex_height, floor_fraction))
        log2_ratio = isfinite(ratio) ? Float32(ratio) : 0.0f0
        return (log2_ratio, indexIrtDistance(apex_scan, max_idx))
    end

    function mainsearchEvidenceExclusionPenalty(
        scan_range::UnitRange{Int},
        mainsearch_evidence_range::Union{Nothing,UnitRange{Int}},
    )
        mainsearch_evidence_range === nothing && return 0.0f0
        evidence_width = candidateWidthIrt(mainsearch_evidence_range)
        if evidence_width <= eps(Float32)
            evidence_idx = first(mainsearch_evidence_range)
            return first(scan_range) <= evidence_idx <= last(scan_range) ? 0.0f0 : 1.0f0
        end
        overlap_start = max(first(scan_range), first(mainsearch_evidence_range))
        overlap_stop = min(last(scan_range), last(mainsearch_evidence_range))
        overlap_width = overlap_start <= overlap_stop ?
            candidateWidthIrt(overlap_start:overlap_stop) :
            0.0f0
        return clamp(1.0f0 - overlap_width / evidence_width, 0.0f0, 1.0f0)
    end

    function mainsearchEvidenceBoundDeltas(
        scan_range::UnitRange{Int},
        mainsearch_evidence_range::Union{Nothing,UnitRange{Int}},
    )
        mainsearch_evidence_range === nothing && return (0.0f0, 0.0f0)
        return (
            indexIrt(first(scan_range)) - indexIrt(first(mainsearch_evidence_range)),
            indexIrt(last(scan_range)) - indexIrt(last(mainsearch_evidence_range)),
        )
    end

    function indexIrt(idx::Int)
        irt = Float32(rt_to_irt_model(rt_col[idx]))
        return isfinite(irt) ? irt : 0.0f0
    end

    function indexIrtDistance(left_idx::Int, right_idx::Int)
        return abs(indexIrt(right_idx) - indexIrt(left_idx))
    end

    function candidateWidthIrt(scan_range::UnitRange{Int})
        width = indexIrtDistance(first(scan_range), last(scan_range))
        return isfinite(width) ? width : 0.0f0
    end

    function irtAsymmetryDelta(scan_range::UnitRange{Int}, apex_idx::Int)
        left_span = indexIrt(apex_idx) - indexIrt(first(scan_range))
        right_span = indexIrt(last(scan_range)) - indexIrt(apex_idx)
        delta = right_span - left_span
        return isfinite(delta) ? Float32(delta) : 0.0f0
    end

    function positiveSmoothedAreaIrt(
        u::Vector{Float32},
        start_idx::Int,
        stop_idx::Int,
        n_pad::Int,
    )
        start_idx >= stop_idx && return 0.0f0

        area = 0.0f0
        @inbounds for i in start_idx:(stop_idx - 1)
            left_val = max(0.0f0, u[i + n_pad])
            right_val = max(0.0f0, u[i + 1 + n_pad])
            area += 0.5f0 * (left_val + right_val) * indexIrtDistance(i, i + 1)
        end
        return area
    end

    function areaAsymmetryPenalty(
        u::Vector{Float32},
        scan_range::UnitRange{Int},
        apex_scan::Int,
        n_pad::Int,
    )
        left_area = positiveSmoothedAreaIrt(u, first(scan_range), apex_scan, n_pad)
        right_area = positiveSmoothedAreaIrt(u, apex_scan, last(scan_range), n_pad)
        ratio = log2(max(right_area, eps(Float32)) / max(left_area, eps(Float32)))
        return isfinite(ratio) ? Float32(ratio) : 0.0f0
    end

    function residualLobeAreaIrt(
        residual_u::Vector{Float32},
        lobe_start::Int,
        lobe_stop::Int,
        scan_range::UnitRange{Int},
        n_pad::Int,
    )
        lobe_start > lobe_stop && return 0.0f0

        area = 0.0f0
        if lobe_start == lobe_stop
            left_width = lobe_start > first(scan_range) ?
                indexIrtDistance(lobe_start - 1, lobe_start) :
                0.0f0
            right_width = lobe_stop < last(scan_range) ?
                indexIrtDistance(lobe_stop, lobe_stop + 1) :
                0.0f0
            width = max(0.5f0 * (left_width + right_width), eps(Float32))
            return max(0.0f0, residual_u[lobe_start + n_pad]) * width
        end

        @inbounds for i in lobe_start:(lobe_stop - 1)
            left_val = max(0.0f0, residual_u[i + n_pad])
            right_val = max(0.0f0, residual_u[i + 1 + n_pad])
            area += 0.5f0 * (left_val + right_val) * indexIrtDistance(i, i + 1)
        end
        if lobe_start > first(scan_range)
            area += 0.5f0 *
                max(0.0f0, residual_u[lobe_start + n_pad]) *
                indexIrtDistance(lobe_start - 1, lobe_start)
        end
        if lobe_stop < last(scan_range)
            area += 0.5f0 *
                max(0.0f0, residual_u[lobe_stop + n_pad]) *
                indexIrtDistance(lobe_stop, lobe_stop + 1)
        end
        return area
    end

    function baselineDisconnectedLobeFeatures(
        residual_u::Vector{Float32},
        scan_range::UnitRange{Int},
        apex_scan::Int,
        n_pad::Int,
        apex_height::Float32,
    )
        threshold = max(eps(Float32), 1.0f-6 * max(apex_height, 0.0f0))
        total_positive_area = 0.0f0
        nonapex_area = 0.0f0
        apex_lobe_area = 0.0f0
        largest_nonapex_area = 0.0f0

        i = first(scan_range)
        while i <= last(scan_range)
            if residual_u[i + n_pad] <= threshold
                i += 1
                continue
            end

            lobe_start = i
            while i <= last(scan_range) && residual_u[i + n_pad] > threshold
                i += 1
            end
            lobe_stop = i - 1
            lobe_area = residualLobeAreaIrt(
                residual_u,
                lobe_start,
                lobe_stop,
                scan_range,
                n_pad,
            )
            total_positive_area += lobe_area
            if lobe_start <= apex_scan <= lobe_stop
                apex_lobe_area += lobe_area
            else
                nonapex_area += lobe_area
                largest_nonapex_area = max(largest_nonapex_area, lobe_area)
            end
        end

        if total_positive_area <= eps(Float32)
            return (0.0f0, -10.0f0)
        end

        disconnected_fraction = Float32(nonapex_area / total_positive_area)
        if largest_nonapex_area <= eps(Float32) || apex_lobe_area <= eps(Float32)
            return (disconnected_fraction, -10.0f0)
        end
        ratio = log2(largest_nonapex_area / apex_lobe_area)
        return (
            disconnected_fraction,
            isfinite(ratio) ? Float32(clamp(ratio, -10.0f0, 10.0f0)) : -10.0f0,
        )
    end

    function baselineInternalTroughFeatures(
        residual_u::Vector{Float32},
        scan_range::UnitRange{Int},
        n_pad::Int,
        apex_height::Float32,
    )
        first(scan_range) + 1 <= last(scan_range) - 1 || return (0.0f0, 0.0f0)

        candidate_peak = 0.0f0
        @inbounds for i in scan_range
            candidate_peak = max(candidate_peak, max(0.0f0, residual_u[i + n_pad]))
        end
        floor_value = max(eps(Float32), 1.0f-6 * max(apex_height, candidate_peak, 0.0f0))
        candidate_peak > floor_value || return (0.0f0, 0.0f0)

        best_score = 0.0f0
        best_log2_ratio = 0.0f0
        @inbounds for i in (first(scan_range) + 1):(last(scan_range) - 1)
            left_peak = 0.0f0
            for j in first(scan_range):(i - 1)
                left_peak = max(left_peak, max(0.0f0, residual_u[j + n_pad]))
            end

            right_peak = 0.0f0
            for j in (i + 1):last(scan_range)
                right_peak = max(right_peak, max(0.0f0, residual_u[j + n_pad]))
            end

            support = min(left_peak, right_peak)
            support > floor_value || continue
            trough = max(0.0f0, residual_u[i + n_pad])
            depth_fraction = max(0.0f0, 1.0f0 - trough / support)
            support_fraction = min(1.0f0, support / candidate_peak)
            score = depth_fraction * support_fraction
            if score > best_score
                best_score = score
                ratio = log2((trough + floor_value) / (support + floor_value))
                best_log2_ratio = isfinite(ratio) ? Float32(clamp(ratio, -20.0f0, 0.0f0)) : 0.0f0
            end
        end

        return (best_score, best_log2_ratio)
    end

    function boundaryCandidateFeatures(
        scan_range::UnitRange{Int},
        u::Vector{Float32},
        baseline_subtracted::Vector{Float32},
        N::Int,
        apex_scan::Int,
        n_pad::Int,
        mainsearch_evidence_range::Union{Nothing,UnitRange{Int}},
        fallback_range::UnitRange{Int},
    )
        apex_height = u[apex_scan + n_pad]
        log2_smoothed_apex_weight = isfinite(apex_height) && apex_height > 0.0f0 ?
            log2(apex_height) : 0.0f0
        if apex_height <= 0.0f0 || isnan(apex_height)
            mainsearch_left_bound_delta, mainsearch_right_bound_delta = mainsearchEvidenceBoundDeltas(
                scan_range,
                mainsearch_evidence_range,
            )
            return (
                candidate_width = candidateWidthIrt(scan_range),
                endpoint_height_fraction = 1.0f0,
                peak_prominence_score = 0.0f0,
                endpoint_valley_score = 0.0f0,
                log2_smoothed_apex_weight = log2_smoothed_apex_weight,
                secondary_peak_penalty = 1.0f0,
                asymmetry_penalty = 1.0f0,
                irt_asymmetry_delta = irtAsymmetryDelta(scan_range, apex_scan),
                baseline_disconnected_signal_fraction = 0.0f0,
                baseline_largest_nonapex_lobe_log2_ratio = -10.0f0,
                baseline_internal_trough_score = 0.0f0,
                baseline_internal_trough_log2_ratio = 0.0f0,
                left_excluded_signal_fraction = 0.0f0,
                right_excluded_signal_fraction = 0.0f0,
                left_boundary_recovery_fraction = 0.0f0,
                right_boundary_recovery_fraction = 0.0f0,
                left_outside_peak_fraction = 0.0f0,
                right_outside_peak_fraction = 0.0f0,
                internal_dip_recovery_score = 0.0f0,
                left_edge_valley_log2_ratio = 0.0f0,
                right_edge_valley_log2_ratio = 0.0f0,
                included_nonapex_max_log2_ratio = 0.0f0,
                included_nonapex_max_irt_distance = 0.0f0,
                mainsearch_exclusion_penalty = mainsearchEvidenceExclusionPenalty(scan_range, mainsearch_evidence_range),
                mainsearch_left_bound_delta = mainsearch_left_bound_delta,
                mainsearch_right_bound_delta = mainsearch_right_bound_delta,
            )
        end

        left = first(scan_range)
        right = last(scan_range)
        width = candidateWidthIrt(scan_range)
        left_val = max(0.0f0, u[left + n_pad])
        right_val = max(0.0f0, u[right + n_pad])

        endpoint_height_fraction = (left_val + right_val) / (2.0f0 * apex_height)
        peak_prominence_fraction = (apex_height - max(left_val, right_val)) / apex_height
        peak_prominence_score = peak_prominence_fraction / sqrt(max(width, eps(Float32)))

        left_local_min = localMinimumNear(u, left, n_pad, N)
        right_local_min = localMinimumNear(u, right, n_pad, N)
        left_valley = 1.0f0 - min(1.0f0, max(0.0f0, left_val - left_local_min) / apex_height)
        right_valley = 1.0f0 - min(1.0f0, max(0.0f0, right_val - right_local_min) / apex_height)
        endpoint_valley_score = (left_valley + right_valley) / 2.0f0

        secondary_peak_penalty = secondaryPeakPenalty(u, scan_range, apex_scan, n_pad, apex_height)
        left_excluded_signal_fraction, left_boundary_recovery_fraction, left_outside_peak_fraction =
            outsideBoundarySignalFeatures(
                u,
                scan_range,
                N,
                n_pad,
                apex_height;
                left_side = true,
            )
        right_excluded_signal_fraction, right_boundary_recovery_fraction, right_outside_peak_fraction =
            outsideBoundarySignalFeatures(
                u,
                scan_range,
                N,
                n_pad,
                apex_height;
                left_side = false,
            )
        internal_dip_recovery_score = internalDipRecoveryScore(
            u,
            scan_range,
            apex_scan,
            n_pad,
            apex_height,
        )
        left_edge_valley_log2_ratio = edgeValleyLog2Ratio(
            u,
            left,
            apex_scan,
            N,
            n_pad,
            apex_height;
            left_side = true,
        )
        right_edge_valley_log2_ratio = edgeValleyLog2Ratio(
            u,
            right,
            apex_scan,
            N,
            n_pad,
            apex_height;
            left_side = false,
        )
        included_nonapex_max_log2_ratio, included_nonapex_max_irt_distance =
            includedNonApexMaxFeatures(
                u,
                scan_range,
                apex_scan,
                n_pad,
                apex_height,
        )

        asymmetry_penalty = areaAsymmetryPenalty(u, scan_range, apex_scan, n_pad)
        irt_asymmetry_delta = irtAsymmetryDelta(scan_range, apex_scan)
        baseline_disconnected_signal_fraction, baseline_largest_nonapex_lobe_log2_ratio =
            baselineDisconnectedLobeFeatures(
                baseline_subtracted,
                scan_range,
                apex_scan,
                n_pad,
                apex_height,
            )
        baseline_internal_trough_score, baseline_internal_trough_log2_ratio =
            baselineInternalTroughFeatures(
                baseline_subtracted,
                scan_range,
                n_pad,
                apex_height,
            )
        mainsearch_exclusion_penalty = mainsearchEvidenceExclusionPenalty(
            scan_range,
            mainsearch_evidence_range,
        )
        mainsearch_left_bound_delta, mainsearch_right_bound_delta = mainsearchEvidenceBoundDeltas(
            scan_range,
            mainsearch_evidence_range,
        )

        return (
            candidate_width = width,
            endpoint_height_fraction = endpoint_height_fraction,
            peak_prominence_score = peak_prominence_score,
            endpoint_valley_score = endpoint_valley_score,
            log2_smoothed_apex_weight = log2_smoothed_apex_weight,
            secondary_peak_penalty = secondary_peak_penalty,
            asymmetry_penalty = asymmetry_penalty,
            irt_asymmetry_delta = irt_asymmetry_delta,
            baseline_disconnected_signal_fraction = baseline_disconnected_signal_fraction,
            baseline_largest_nonapex_lobe_log2_ratio = baseline_largest_nonapex_lobe_log2_ratio,
            baseline_internal_trough_score = baseline_internal_trough_score,
            baseline_internal_trough_log2_ratio = baseline_internal_trough_log2_ratio,
            left_excluded_signal_fraction = left_excluded_signal_fraction,
            right_excluded_signal_fraction = right_excluded_signal_fraction,
            left_boundary_recovery_fraction = left_boundary_recovery_fraction,
            right_boundary_recovery_fraction = right_boundary_recovery_fraction,
            left_outside_peak_fraction = left_outside_peak_fraction,
            right_outside_peak_fraction = right_outside_peak_fraction,
            internal_dip_recovery_score = internal_dip_recovery_score,
            left_edge_valley_log2_ratio = left_edge_valley_log2_ratio,
            right_edge_valley_log2_ratio = right_edge_valley_log2_ratio,
            included_nonapex_max_log2_ratio = included_nonapex_max_log2_ratio,
            included_nonapex_max_irt_distance = included_nonapex_max_irt_distance,
            mainsearch_exclusion_penalty = mainsearch_exclusion_penalty,
            mainsearch_left_bound_delta = mainsearch_left_bound_delta,
            mainsearch_right_bound_delta = mainsearch_right_bound_delta,
        )
    end

    function getMainsearchEvidenceRange(
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

        baseline_start = min(li, ri)
        baseline_stop = max(li, ri)

        # Subtract interpolated baseline only between the baseline anchors.
        # Points outside those anchors are beyond the local valleys and should
        # not be amplified by line extrapolation.
        @inbounds @fastmath for i in scan_start:scan_stop
            if i < baseline_start || i > baseline_stop
                u[i] = 0.0f0
                continue
            end
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

    function candidateBaselineSubtractedTrace(
        scan_range::UnitRange{Int},
        smoothed::Vector{Float32},
        n_active::Int,
    )
        candidate_u = copy(@view smoothed[1:n_active])
        subtractBaseline!(
            x,
            candidate_u,
            apex_scan,
            scan_range,
            n_pad,
        )
        return candidate_u
    end

    function integrateCandidateRange(
        scan_range::UnitRange{Int},
        candidate_u::Vector{Float32},
    )
        apex_val = candidate_u[apex_scan + n_pad]
        if apex_val <= 0.0f0 || isnan(apex_val)
            return 0.0f0, UInt32(scan_idx_col[apex_scan]), UInt32(0)
        end

        candidate_state = Chromatogram(zeros(Float32, m), zeros(Float32, m), 0)
        norm_factor, _, rt_norm, _ = fillState!(
            candidate_state,
            candidate_u,
            rt_col,
            first(scan_range),
            last(scan_range),
            apex_scan,
            n_pad,
        )
        raw_trap = integrateTrapezoidal(candidate_state, avg_cycle_time)
        area = Float32(rt_norm * norm_factor * raw_trap)
        if isnan(area) || area < 0.0f0
            area = 0.0f0
        end

        min_abundance = 0.1f0 * norm_factor
        scan_start = first(scan_range) + n_pad
        scan_stop = last(scan_range) + n_pad
        points = count(i -> candidate_u[i] >= min_abundance, scan_start:scan_stop)
        return area, UInt32(scan_idx_col[apex_scan]), UInt32(points)
    end

    function collectBoundaryCandidateData(
        candidates::Vector,
        fallback_range::UnitRange{Int},
        smoothed::Vector{Float32},
        n_active::Int,
        mainsearch_evidence_range::Union{Nothing,UnitRange{Int}},
    )
        rows = NamedTuple[]
        sizehint!(rows, length(candidates))

        for (candidate_index, candidate) in enumerate(candidates)
            scan_range = candidate.scan_range
            baseline_subtracted = candidateBaselineSubtractedTrace(scan_range, smoothed, n_active)
            features = boundaryCandidateFeatures(
                scan_range,
                smoothed,
                baseline_subtracted,
                m,
                apex_scan,
                n_pad,
                mainsearch_evidence_range,
                fallback_range,
            )
            area, best_scan, points = integrateCandidateRange(scan_range, baseline_subtracted)
            push!(rows, (;
                candidate_index = UInt16(candidate_index),
                candidate_start_idx = UInt16(first(scan_range)),
                candidate_stop_idx = UInt16(last(scan_range)),
                candidate_start_scan = UInt32(scan_idx_col[first(scan_range)]),
                candidate_stop_scan = UInt32(scan_idx_col[last(scan_range)]),
                candidate_category = candidate.category,
                peak_area = Float32(area),
                new_best_scan = UInt32(best_scan),
                points_integrated = UInt32(points),
                is_fallback = first(scan_range) == first(fallback_range) &&
                    last(scan_range) == last(fallback_range),
                features...
            ))
        end

        return rows
    end

    #Whittaker Henderson Smoothing — result written into ws.z
    WHSmooth!(ws, intensity_col, fraction_col, rt_col, m, n_pad, min_fraction_transmitted, λ)

    # Local aliases to workspace vectors (concrete Vector{Float32}, no SubArray views)
    n_active = m + 2*n_pad
    z  = ws.z
    x  = ws.x_tmp
    u2 = ws.u2

    #Second discrete derivative of smoothed data
    fillU2!(u2, z, x, n_active)

    apex_scan = clamp(getApexScan(
        apex_scan,
        n_pad,
        z,
        n_active
    ), 1, m)

    mainsearch_evidence_range = getMainsearchEvidenceRange(
        scan_idx_col,
        mainsearch_1pct_start_scan,
        mainsearch_1pct_stop_scan,
    )

    #Integration boundaries based on smoothed second derivative
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
    if boundary_candidate_data !== nothing
        candidates = generateBoundaryCandidates(
            fallback_scan_range,
            z,
            m,
            apex_scan,
            n_pad,
        )
        boundary_candidate_data[] = collectBoundaryCandidateData(
            candidates,
            fallback_scan_range,
            z,
            n_active,
            mainsearch_evidence_range,
        )
    end
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

    #File `state` to fit EGH function. Get the inensity, and rt normalization factors
    # Guard: if smoothed apex is zero or NaN after baseline subtraction, skip integration
    _apex_val = z[apex_scan + n_pad]
    if _apex_val <= 0f0 || isnan(_apex_val)
        if debug_enabled
            status = "skipped_zero_apex"
            debug_plot_data !== nothing && (debug_plot_data[] = (
                scan_range = scan_range,
                wh_smoothed = wh_smoothed_debug,
                baseline_subtracted = baseline_subtracted_debug,
                peak_area = 0.0f0,
                points_integrated = UInt32(0),
                status = status,
                mainsearch_evidence_range = mainsearch_evidence_range,
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
                    apex_scan,
                    0.0f0,
                    0,
                    min_fraction_transmitted,
                    status,
                    debug_plot_title,
                    mainsearch_evidence_range,
                )
            end
        end
        return 0f0, scan_idx_col[apex_scan], Int(0)
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
            wh_smoothed = wh_smoothed_debug,
            baseline_subtracted = baseline_subtracted_debug,
            peak_area = trapezoid_area,
            points_integrated = UInt32(num_points_integrated),
            status = status,
            mainsearch_evidence_range = mainsearch_evidence_range,
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
                apex_scan,
                trapezoid_area,
                num_points_integrated,
                min_fraction_transmitted,
                status,
                debug_plot_title,
                mainsearch_evidence_range,
            )
        end
    end

    #trapezoid_area = 0.0f0
    return trapezoid_area, scan_idx_col[apex_scan], num_points_integrated
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
                                mainsearch_1pct_start_scan::UInt32 = UInt32(0),
                                mainsearch_1pct_stop_scan::UInt32 = UInt32(0),
                                rt_to_irt_model::RtConversionModel = IdentityModel(),
                                forced_boundary_start_scan::UInt32 = UInt32(0),
                                forced_boundary_stop_scan::UInt32 = UInt32(0),
                                boundary_candidate_data::Union{Nothing, Base.RefValue} = nothing)
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
        mainsearch_1pct_start_scan = mainsearch_1pct_start_scan,
        mainsearch_1pct_stop_scan = mainsearch_1pct_stop_scan,
        rt_to_irt_model = rt_to_irt_model,
        forced_boundary_start_scan = forced_boundary_start_scan,
        forced_boundary_stop_scan = forced_boundary_stop_scan,
        boundary_candidate_data = boundary_candidate_data,
    )
end
