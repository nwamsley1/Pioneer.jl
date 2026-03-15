"""
    estimate_precursor_fwhm(irt_obs, weights)

Estimate FWHM for a single precursor from its multi-scan weight profile.
Walks outward from the max-weight scan to find half-maximum boundaries.
Returns `missing` if either boundary cannot be determined.
"""
function estimate_precursor_fwhm(irt_obs::AbstractVector{Float32},
                                  weights::AbstractVector{Float32})
    max_idx = argmax(weights)
    half_max = weights[max_idx] / 2.0f0

    # Walk left from max
    left_irt = missing
    for i in (max_idx - 1):-1:1
        if weights[i] > half_max
            left_irt = irt_obs[i]
        else
            break
        end
    end

    # Walk right from max
    right_irt = missing
    for i in (max_idx + 1):length(weights)
        if weights[i] > half_max
            right_irt = irt_obs[i]
        else
            break
        end
    end

    (ismissing(left_irt) || ismissing(right_irt)) && return missing
    return right_irt - left_irt
end

"""
    estimate_file_fwhm(psms; max_pep=0.1f0)

Estimate per-file median and MAD of precursor FWHMs from multi-scan PSMs.
Only includes precursors with ≥3 scans. Returns named tuple (median_fwhm, mad_fwhm).
Falls back to (0.2, 0.2) if no valid FWHMs found.
"""
function estimate_file_fwhm(psms::DataFrame; max_pep::Float32 = 0.1f0)
    t_start = time()
    sorted = sort(psms, [:precursor_idx, :irt_obs])
    gdf = groupby(sorted, :precursor_idx)

    fwhms = Float32[]
    n_precursors_evaluated = 0
    for g in gdf
        nrow(g) < 3 && continue
        n_precursors_evaluated += 1
        fwhm = estimate_precursor_fwhm(g[!, :irt_obs], g[!, :weight])
        ismissing(fwhm) && continue
        fwhm > 0 && push!(fwhms, fwhm)
    end

    elapsed = round(time() - t_start, digits=3)
    if isempty(fwhms)
        @info "    FWHM estimation: 0/$(n_precursors_evaluated) precursors had valid FWHM ($(elapsed)s), using defaults"
        return (median_fwhm=0.2f0, mad_fwhm=0.2f0)
    end

    med = Float32(median(fwhms))
    md = Float32(mad(fwhms, normalize=true))
    @info "    FWHM estimation: $(length(fwhms))/$(n_precursors_evaluated) precursors with valid FWHM, " *
          "median=$(round(med, digits=4)), MAD=$(round(md, digits=4)) ($(elapsed)s)"
    return (median_fwhm=med, mad_fwhm=md)
end

"""
    compute_chromatographic_tolerance!(search_context, median_irt_range, ms_data, n_files)

Set chromatographic iRT tolerance for all files using the median iRT range
(max - min observed iRT across scans) from high-confidence precursors (q ≤ 0.1%).

The tolerance is applied as ±irt_tol around each scan's iRT, so we use half
the full range to get the one-sided tolerance.

Sets the result into `getIrtErrors(search_context)[ms_file_idx]` for each file.
"""
function compute_chromatographic_tolerance!(
    search_context::SearchContext,
    median_irt_range::Float32,
    ms_data,
    n_files::Int
)
    irt_tol = median_irt_range / 2.0f0
    for ms_file_idx in 1:n_files
        getFailedIndicator(ms_data, ms_file_idx) && continue
        getIrtErrors(search_context)[ms_file_idx] = irt_tol
    end
    @info "Chromatographic tolerance: median_irt_range=$(round(median_irt_range, digits=4)), irt_tol=±$(round(irt_tol, digits=4))"
end
