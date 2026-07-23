# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# Sciex ZT scanning-DIA: estimate the effective quadrupole transmission from the
# overlapping-scan fragment envelope, as an EmpiricalQuadModel (uniform LUT knots).
#
# The physical ZT quad is a ~5 Da window swept in ~1 Da steps, so a precursor's ions
# appear across many neighboring Q1 bins within ONE cycle. For each confident precursor
# whose m/z sits within `inbin` Da of a bin center (so its center bin ≈ peak transmission),
# we take each of its top-3 LIBRARY fragments' monoisotopic intensity in each neighbor bin
# and normalize to that SAME fragment's intensity in the SAME cycle's center bin. Pooling
# these ratios over precursors gives the transmission vs offset; the per-1-Da-knot median
# (symmetrized, peak-normalized) is the LUT. Normalization is strictly within-cycle.

@inline function _zt_frag_peak(mz::AbstractVector, it::AbstractVector, target::Float32, tol_ppm::Float32)
    tol = target * tol_ppm * 1f-6
    lo = target - tol; hi = target + tol
    best = 0.0f0
    @inbounds for i in eachindex(mz)
        m = mz[i]
        ismissing(m) && continue
        if (m >= lo) & (m <= hi)
            v = it[i]
            (!ismissing(v) && v > best) && (best = v)
        end
    end
    return best
end

"""
    estimate_zt_quad_transmission(psms, spectra, search_context, ms_file_idx; kwargs...)
      -> (EmpiricalQuadModel | nothing, knot_samples::Vector{Vector{Float32}})

Estimate the ZT effective quad transmission. `psms` must carry :precursor_idx and
:scan_idx (center-bin candidates from QuadTuning). Returns the fitted LUT model (or
`nothing` if too little signal) plus the per-knot ratio samples for the QC plot.
"""
function estimate_zt_quad_transmission(
        psms::DataFrame,
        spectra::MassSpecData,
        search_context::SearchContext,
        ms_file_idx::Int64;
        knot_step::Float32 = 1.0f0,
        knot_half::Float32 = 7.0f0,
        inbin::Float32 = 0.10f0,
        floor_int::Float32 = 200.0f0,
        tol_ppm::Float32 = 15.0f0,
        min_per_knot::Int = 40,
        neighbor_span::Int = 9,
        symmetric::Bool = true)

    lib = getSpecLib(search_context)
    lookup = getFragmentLookupTable(lib)
    prec_mz = getMz(getPrecursors(lib))
    nprec = length(prec_mz)

    nknot = round(Int, 2 * knot_half / knot_step) + 1
    mid = (nknot + 1) ÷ 2
    samples = [Float32[] for _ in 1:nknot]

    top3 = Dict{UInt32, NTuple{3, Float32}}()
    function frags(pid::UInt32)
        get!(top3, pid) do
            rr = getPrecFragRange(lookup, pid)
            n = length(rr)
            n == 0 && return (0.0f0, 0.0f0, 0.0f0)
            m1 = Float32(getMz(lookup.frags[rr[1]]))
            m2 = n >= 2 ? Float32(getMz(lookup.frags[rr[2]])) : 0.0f0
            m3 = n >= 3 ? Float32(getMz(lookup.frags[rr[3]])) : 0.0f0
            (m1, m2, m3)
        end
    end

    pids = psms[!, :precursor_idx]::AbstractVector
    scns = psms[!, :scan_idx]::AbstractVector
    seen = Set{Tuple{UInt32, UInt32}}()
    nscan = length(spectra)

    for r in eachindex(pids)
        pid = UInt32(pids[r]); c = Int(scns[r])
        (pid < 1 || pid > nprec) && continue
        (c < 1 || c > nscan) && continue
        (pid, UInt32(c)) in seen && continue
        push!(seen, (pid, UInt32(c)))

        P = prec_mz[pid]
        cmz = getCenterMz(spectra, c)
        (ismissing(cmz) || abs(Float32(cmz) - P) > inbin) && continue   # near bin center only
        fm = frags(pid)
        fm[1] == 0.0f0 && continue
        cycle = getCycleIdx(spectra, c)

        cmzc = getMzArray(spectra, c); citc = getIntensityArray(spectra, c)
        cent = (_zt_frag_peak(cmzc, citc, fm[1], tol_ppm),
                _zt_frag_peak(cmzc, citc, fm[2], tol_ppm),
                _zt_frag_peak(cmzc, citc, fm[3], tol_ppm))
        (cent[1] <= floor_int && cent[2] <= floor_int && cent[3] <= floor_int) && continue

        for d in -neighbor_span:neighbor_span
            s = c + d
            (s < 1 || s > nscan) && continue
            getMsOrder(spectra, s) == 2 || continue
            getCycleIdx(spectra, s) == cycle || continue      # STRICTLY within-cycle
            scm = getCenterMz(spectra, s); ismissing(scm) && continue
            off = Float32(scm) - P
            abs(off) >= knot_half && continue
            ki = round(Int, off / knot_step) + mid
            (ki < 1 || ki > nknot) && continue
            smz = getMzArray(spectra, s); sit = getIntensityArray(spectra, s)
            @inbounds for j in 1:3
                cent[j] <= floor_int && continue
                v = _zt_frag_peak(smz, sit, fm[j], tol_ppm)
                push!(samples[ki], v / cent[j])
            end
        end
    end

    vals = zeros(Float32, nknot)
    for i in 1:nknot
        vals[i] = length(samples[i]) >= min_per_knot ? Float32(median(samples[i])) : 0.0f0
    end
    maximum(vals) <= 0.0f0 && return (nothing, samples)

    if symmetric
        @inbounds for i in 1:(mid - 1)
            a = (vals[i] + vals[nknot + 1 - i]) / 2
            vals[i] = a; vals[nknot + 1 - i] = a
        end
    end
    vals ./= maximum(vals)
    return (EmpiricalQuadModel(vals, knot_step, knot_half), samples)
end

"""
    plot_zt_quad_transmission(model, samples, fname; knot_step, knot_half) -> Plots.Plot

QC plot for the ZT empirical transmission: per-knot median markers with IQR ribbon and
the interpolated model curve.
"""
function plot_zt_quad_transmission(model::EmpiricalQuadModel, samples::Vector{Vector{Float32}},
                                   fname::String; knot_step::Float32 = 1.0f0, knot_half::Float32 = 7.0f0)
    nknot = length(model.values)
    mid = (nknot + 1) ÷ 2
    kx = [(i - mid) * knot_step for i in 1:nknot]
    kmed = Float32[]; kq1 = Float32[]; kq3 = Float32[]; kxx = Float32[]
    for i in 1:nknot
        length(samples[i]) < 5 && continue
        v = samples[i]
        push!(kxx, kx[i]); push!(kmed, median(v))
        push!(kq1, quantile(v, 0.25)); push!(kq3, quantile(v, 0.75))
    end
    # peak-normalize the displayed medians to match the model
    pk = isempty(kmed) ? 1.0f0 : maximum(kmed)
    pk <= 0 && (pk = 1.0f0)
    p = Plots.plot(title = "ZT quad transmission — $fname", xlabel = "Q1 offset from precursor m/z (Da)",
                   ylabel = "transmission (÷ center bin)", legend = :topright,
                   size = (700, 450), dpi = 100)
    if !isempty(kxx)
        Plots.plot!(p, kxx, kmed ./ pk; ribbon = (max.(kmed .- kq1, 0f0) ./ pk, max.(kq3 .- kmed, 0f0) ./ pk),
                    seriestype = :scatter, ms = 4, alpha = 0.8, label = "median ± IQR", color = :black)
    end
    dx = collect(-knot_half:0.1f0:knot_half)
    f = getQuadTransmissionFunction(model, 0.0f0, 1.0f0)
    Plots.plot!(p, dx, [Float64(f(Float32(x))) for x in dx]; lw = 2.5, color = :dodgerblue,
                label = "EmpiricalQuadModel (LUT)")
    return p
end
