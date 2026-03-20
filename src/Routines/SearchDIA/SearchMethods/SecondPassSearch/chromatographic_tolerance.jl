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
