# MS1 mass-error diagnostic.
#
# For each PSM accepted by the MS2 mass-error fit, look up the nearest MS1 scan
# and find the closest peak to the precursor's predicted M+0/M+1/M+2 m/z within
# ±MS1_DIAG_WINDOW_PPM. Record observed-vs-theoretical ppm residuals and write a
# per-file histogram. Diagnostic only — does not fit a model or affect downstream
# search behavior.

const MS1_DIAG_WINDOW_PPM = 100.0f0
const MS1_DIAG_TOL_K_MAD  = 5.0f0

# Inlined minimal nearest-peak finder. Uses Julia's searchsortedfirst rather
# than the perf-tuned bsearch_hybrid since this only runs once per file on a
# bounded PSM set.
function _ms1_diag_find_peak(mz::Vector{Float32}, mz_target::Float32, tol_ppm::Float32)
    half_tol = mz_target * tol_ppm * 1f-6
    lo = mz_target - half_tol
    hi = mz_target + half_tol
    n = length(mz)
    n == 0 && return (false, 0f0)
    i0 = searchsortedfirst(mz, lo)
    best_diff = Inf32
    best_mz   = 0f0
    found = false
    @inbounds for j in i0:n
        m = mz[j]
        m > hi && break
        d = abs(m - mz_target)
        if d < best_diff
            best_diff = d
            best_mz   = m
            found = true
        end
    end
    return found, best_mz
end

# Strip Union{Missing, Float32} into a dense sorted-by-m/z Vector{Float32}.
function _ms1_diag_clean_mz(raw_mz)
    n_raw = length(raw_mz)
    out = Vector{Float32}()
    sizehint!(out, n_raw)
    @inbounds for j in 1:n_raw
        m = raw_mz[j]
        ismissing(m) || push!(out, Float32(m))
    end
    return out
end

"""
    collect_ms1_residuals(spectra, psms, search_context, ms_file_idx) -> Vector{Float32}

For each (scan_idx, precursor_idx) row in `psms` (typically the MS2-accepted
PSMs from ParameterTuningSearch), find the nearest MS1 scan by RT and record
the ppm residual of the closest peak to the precursor's M+0, M+1, M+2 isotopes
within ±MS1_DIAG_WINDOW_PPM. Returns a flat vector of signed ppm residuals
(positive = observed peak above theoretical m/z).
"""
function collect_ms1_residuals(spectra, psms::DataFrame, search_context, ms_file_idx::Integer)
    n = nrow(psms)
    n == 0 && return Float32[]

    n_scans = length(spectra)
    ms1_scan_idxs = Int[]
    ms1_scan_rts  = Float32[]
    for s in 1:n_scans
        if getMsOrder(spectra, s) == 1
            push!(ms1_scan_idxs, s)
            push!(ms1_scan_rts,  Float32(getRetentionTime(spectra, s)))
        end
    end
    isempty(ms1_scan_idxs) && return Float32[]

    n_ms1 = length(ms1_scan_rts)
    scan_to_ms1 = Vector{Int32}(undef, n_scans)
    @inbounds for s in 1:n_scans
        scan_rt = Float32(getRetentionTime(spectra, s))
        pos = searchsortedfirst(ms1_scan_rts, scan_rt)
        scan_to_ms1[s] = if pos == 1
            Int32(ms1_scan_idxs[1])
        elseif pos > n_ms1
            Int32(ms1_scan_idxs[end])
        else
            d_after  = abs(ms1_scan_rts[pos]   - scan_rt)
            d_before = abs(ms1_scan_rts[pos-1] - scan_rt)
            d_before <= d_after ? Int32(ms1_scan_idxs[pos-1]) : Int32(ms1_scan_idxs[pos])
        end
    end

    precursors  = getPrecursors(getSpecLib(search_context))
    prec_mzs     = getMz(precursors)
    prec_charges = getCharge(precursors)

    residuals = Float32[]
    sizehint!(residuals, n * 3)

    cached_mz::Vector{Float32} = Float32[]
    cached_ms1_idx::Int = -1

    @inbounds for i in 1:n
        ms1_idx = Int(scan_to_ms1[Int(psms.scan_idx[i])])
        if ms1_idx != cached_ms1_idx
            cached_ms1_idx = ms1_idx
            cached_mz = _ms1_diag_clean_mz(getMzArray(spectra, ms1_idx))
        end
        pid = UInt32(psms.precursor_idx[i])
        prec_mz  = Float32(prec_mzs[pid])
        prec_chg = Int(prec_charges[pid]); prec_chg == 0 && (prec_chg = 1)

        for iso in 0:2
            target = prec_mz + Float32(iso) * C13_C12_MASS_DIFF_F32 / Float32(prec_chg)
            hit, obs_mz = _ms1_diag_find_peak(cached_mz, target, MS1_DIAG_WINDOW_PPM)
            if hit
                push!(residuals, (obs_mz - target) / target * 1f6)
            end
        end
    end
    return residuals
end

"""
    fit_ms1_model_from_residuals(residuals) -> (model::MassErrorModel, median_ppm, mad_ppm) or nothing

Fit `MassErrorModel(median, (k·MAD, k·MAD))` with k = MS1_DIAG_TOL_K_MAD.
Returns `nothing` if fewer than 10 residuals are provided.
"""
function fit_ms1_model_from_residuals(residuals::AbstractVector{<:Real})
    n = length(residuals)
    n < 10 && return nothing
    med  = Float32(median(residuals))
    mad_ = Float32(mad(residuals; normalize=true))
    tol  = MS1_DIAG_TOL_K_MAD * mad_
    return MassErrorModel(med, (tol, tol)), med, mad_
end

"""
    generate_ms1_residual_histogram(residuals, fname, out_path)

Save a histogram of MS1 ppm residuals. Black dashed lines mark the
median ± k·MAD bounds installed as the MS1 mass error tolerance. No-op
if `length(residuals) < 10`.
"""
function generate_ms1_residual_histogram(residuals::AbstractVector{<:Real},
                                         fname::AbstractString,
                                         out_path::AbstractString)
    n = length(residuals)
    n < 10 && return nothing
    med  = median(residuals)
    mad_ = mad(residuals; normalize=true)
    tol  = MS1_DIAG_TOL_K_MAD * mad_
    rng  = max(quantile(abs.(residuals), 0.999), abs(med) + (MS1_DIAG_TOL_K_MAD + 1) * mad_)
    bins = collect(-rng:max(rng/50, 0.1):rng)
    p = Plots.histogram(residuals;
        bins = bins, label = nothing,
        title = _split_title(fname, "MS1 ppm residuals") *
                "\nn=$n, median=$(round(med, digits=2)) ppm, " *
                "MAD=$(round(mad_, digits=2)) ppm, tol=±$(round(tol, digits=2)) ppm",
        xlabel = "Observed − theoretical (ppm)",
        ylabel = "Count",
        size = (900, 600))
    Plots.vline!(p, [med]; label = nothing, color = :red, lw = 2)
    Plots.vline!(p, [med - tol, med + tol]; label = nothing, color = :black, lw = 1, ls = :dash)
    Plots.savefig(p, out_path)
    return p
end
