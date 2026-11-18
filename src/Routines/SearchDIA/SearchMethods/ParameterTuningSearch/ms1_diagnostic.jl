# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# MS1 mass-error diagnostic.
#
# For each PSM accepted by the MS2 mass-error fit, look up the nearest MS1 scan
# and find the most intense peak within ±MS1_DIAG_WINDOW_PPM of the precursor's
# predicted M+0/M+1/M+2 m/z. (Most intense, not nearest: the window is centred on
# zero offset, and on timsTOF data strong peaks carry weaker satellite centroids
# ~40 ppm away, so once the true offset passes ~20 ppm the nearest peak is a
# satellite and the median lands between the two clusters.) Record observed-vs-theoretical ppm residuals and write a
# per-file histogram. These observations also fit the production MS1 model.

const MS1_DIAG_WINDOW_PPM = 100.0f0
const MS1_DIAG_TOL_K_MAD  = 5.0f0

# Inlined minimal most-intense-peak finder. Uses Julia's searchsortedfirst rather
# than the perf-tuned bsearch_hybrid since this only runs once per file on a
# bounded PSM set. Ties keep the first (lowest m/z) peak.
# Index of the most intense peak within ±tol_ppm of `mz_target`, or 0.
function _ms1_diag_find_peak_idx(mz::Vector{Float32}, intensity::Vector{Float32}, mz_target::Float32, tol_ppm::Float32)
    half_tol = mz_target * tol_ppm * 1f-6
    lo = mz_target - half_tol
    hi = mz_target + half_tol
    n = length(mz)
    n == 0 && return 0
    i0 = searchsortedfirst(mz, lo)
    best_int = -Inf32
    best = 0
    @inbounds for j in i0:n
        m = mz[j]
        m > hi && break
        if intensity[j] > best_int
            best_int = intensity[j]
            best = j
        end
    end
    return best
end

# Strip Union{Missing, Float32} into dense sorted-by-m/z vectors (m/z, intensity; missing intensity -> 0).
function _ms1_diag_clean_peaks!(out_mz::Vector{Float32}, out_int::Vector{Float32}, raw_mz, raw_int)
    empty!(out_mz); empty!(out_int)
    n_raw = length(raw_mz)
    sizehint!(out_mz, n_raw); sizehint!(out_int, n_raw)
    @inbounds for j in 1:n_raw
        m = raw_mz[j]
        ismissing(m) && continue
        push!(out_mz, Float32(m))
        it = raw_int[j]
        push!(out_int, ismissing(it) ? 0f0 : Float32(it))
    end
    return out_mz, out_int
end

"""
    collect_ms1_samples(spectra, psms, search_context, ms_file_idx) -> Vector{MassErrSample}

For each (scan_idx, precursor_idx) row in `psms`, find the nearest MS1 scan (nearest MS1 frame at the row's own
mobility on ion-mobility data) and, for the precursor's M+0, M+1, M+2 isotopes, the most intense peak within
±MS1_DIAG_WINDOW_PPM. One sample per matched isotope: theoretical m/z, observed m/z, peak intensity (0 when the
spectrum has none), MS1 scan RT and precursor. The input of `collect_ms1_residuals` and of the intensity-aware
MS1 model fitted for DDA MS1 quantification.
"""
function collect_ms1_samples(spectra, psms::DataFrame, search_context, ms_file_idx::Integer)
    n = nrow(psms)
    samples = MassErrSample[]
    n == 0 && return samples

    # Nearest MS1 scan per scan (nearest MS1 frame at the row's own mobility on ion-mobility data);
    # shared with MainSearch's add_ms1_lookup_features!.
    scan_to_ms1 = build_scan_to_ms1(spectra)
    all(iszero, scan_to_ms1) && return samples

    precursors  = getPrecursors(getSpecLib(search_context))
    prec_mzs     = getMz(precursors)
    prec_charges = getCharge(precursors)
    sizehint!(samples, n * 3)

    cached_mz::Vector{Float32} = Float32[]
    cached_int::Vector{Float32} = Float32[]
    cached_ms1_idx::Int = -1
    cached_rt::Float32 = 0f0
    decode_buf = PeakDecodeBuffer()

    @inbounds for i in 1:n
        ms1_idx = Int(scan_to_ms1[Int(psms.scan_idx[i])])
        if ms1_idx != cached_ms1_idx
            cached_ms1_idx = ms1_idx
            _ms1_diag_clean_peaks!(cached_mz, cached_int, getPeaks!(decode_buf, spectra, ms1_idx)...)
            cached_rt = Float32(getRetentionTime(spectra, ms1_idx))
        end
        pid = UInt32(psms.precursor_idx[i])
        prec_mz  = Float32(prec_mzs[pid])
        prec_chg = Int(prec_charges[pid]); prec_chg == 0 && (prec_chg = 1)

        for iso in 0:2
            target = prec_mz + Float32(iso) * C13_C12_MASS_DIFF_F32 / Float32(prec_chg)
            j = _ms1_diag_find_peak_idx(cached_mz, cached_int, target, MS1_DIAG_WINDOW_PPM)
            j == 0 && continue
            push!(samples, MassErrSample(target, cached_mz[j], cached_int[j], cached_rt, pid))
        end
    end
    return samples
end

"""
    collect_ms1_residuals(spectra, psms, search_context, ms_file_idx) -> Vector{Float32}

For each (scan_idx, precursor_idx) row in `psms` (typically the MS2-accepted
PSMs from ParameterTuningSearch), find the nearest MS1 scan by RT and record
the ppm residual of the most intense peak near the precursor's M+0, M+1, M+2 isotopes
within ±MS1_DIAG_WINDOW_PPM (`collect_ms1_samples`). Returns a flat vector of signed ppm residuals
(positive = observed peak above theoretical m/z). Optional `qc_coordinates` receives
matching theoretical m/z and scan RT values for within-file residual checks;
`qc_precursor_ids` receives precursor IDs for grouping whole peptide sequences.
"""
function collect_ms1_residuals(spectra, psms::DataFrame, search_context, ms_file_idx::Integer;
    qc_coordinates::Union{Nothing,Tuple{Vector{Float32},Vector{Float32}}}=nothing,
    qc_precursor_ids::Union{Nothing,Vector{UInt32}}=nothing)
    samples = collect_ms1_samples(spectra, psms, search_context, ms_file_idx)
    residuals = Float32[(s.observed_mz - s.theoretical_mz) / s.theoretical_mz * 1f6 for s in samples]
    if qc_coordinates !== nothing
        append!(qc_coordinates[1], (s.theoretical_mz for s in samples))
        append!(qc_coordinates[2], (s.rt for s in samples))
    end
    qc_precursor_ids !== nothing && append!(qc_precursor_ids, (s.precursor_idx for s in samples))
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
