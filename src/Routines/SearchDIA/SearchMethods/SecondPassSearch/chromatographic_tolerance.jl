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
    compute_chromatographic_tolerance!(search_context, file_fwhms, ms_data, n_files;
                                       fwhm_nstd=2.0f0, irt_nstd=3.0f0)

Compute chromatographic iRT tolerance per file, combining:
- Peak width: median_FWHM + fwhm_nstd * MAD_FWHM
- Cross-run iRT variance: irt_nstd * median(per-precursor iRT std across runs)

Sets the result into `getIrtErrors(search_context)[ms_file_idx]` for each file.
"""
function compute_chromatographic_tolerance!(
    search_context::SearchContext,
    file_fwhms::Dict{Int, @NamedTuple{median_fwhm::Float32, mad_fwhm::Float32}},
    ms_data,
    n_files::Int;
    fwhm_nstd::Float32 = 2.0f0,
    irt_nstd::Float32 = 3.0f0
)
    t_start = time()

    # Collect per-precursor iRT across files from fold Arrow files
    prec_irts = Dictionary{UInt32, Vector{Float32}}()
    n_psms_read = 0
    for ms_file_idx in 1:n_files
        getFailedIndicator(ms_data, ms_file_idx) && continue
        base_path = getSecondPassPsms(ms_data, ms_file_idx)
        isempty(base_path) && continue
        for fold in UInt8[0, 1]
            fold_path = "$(base_path)_fold$(fold).arrow"
            isfile(fold_path) || continue
            tbl = Arrow.Table(fold_path)
            hasproperty(tbl, :irt_obs) || continue
            for i in eachindex(tbl[:precursor_idx])
                pid = tbl[:precursor_idx][i]
                irt = Float32(tbl[:irt_obs][i])
                n_psms_read += 1
                if haskey(prec_irts, pid)
                    push!(prec_irts[pid], irt)
                else
                    insert!(prec_irts, pid, Float32[irt])
                end
            end
        end
    end
    t_read = time()

    # Cross-run iRT std (prefer >2 runs, fallback to 2, then 0)
    stds_gt2 = Float32[Float32(std(irts)) for (_, irts) in pairs(prec_irts) if length(irts) > 2]
    if !isempty(stds_gt2)
        irt_std = Float32(median(stds_gt2))
    else
        stds_eq2 = Float32[Float32(abs(irts[1] - irts[2]) / sqrt(2.0f0))
                           for (_, irts) in pairs(prec_irts) if length(irts) == 2]
        irt_std = isempty(stds_eq2) ? 0.0f0 : Float32(median(stds_eq2))
    end

    # Set per-file tolerance
    for ms_file_idx in 1:n_files
        getFailedIndicator(ms_data, ms_file_idx) && continue
        fwhm_stats = get(file_fwhms, ms_file_idx, (median_fwhm=0.2f0, mad_fwhm=0.2f0))
        peak_width = fwhm_stats.median_fwhm + fwhm_nstd * fwhm_stats.mad_fwhm
        chrom_tol = peak_width + irt_nstd * irt_std
        getIrtErrors(search_context)[ms_file_idx] = chrom_tol
    end

    elapsed = round(time() - t_start, digits=2)
    t_read_elapsed = round(t_read - t_start, digits=2)
    n_precs_multi = count(p -> length(p) > 1 for (_, p) in pairs(prec_irts))

    @info "Chromatographic iRT tolerance computed in $(elapsed)s (read=$(t_read_elapsed)s)"
    @info "  $(n_psms_read) PSMs read, $(length(prec_irts)) unique precursors, $(n_precs_multi) seen in >1 file"
    @info "  cross_run_irt_std=$(round(irt_std, digits=4)) ($(length(stds_gt2)) precs with >2 files)"
    for (fidx, stats) in sort(collect(file_fwhms), by=first)
        fname = getFileIdToName(ms_data, fidx)
        tol = getIrtErrors(search_context)[fidx]
        @info "  $(fname): FWHM=$(round(stats.median_fwhm, digits=4))±$(round(stats.mad_fwhm, digits=4)), chrom_tol=$(round(tol, digits=4))"
    end
end
