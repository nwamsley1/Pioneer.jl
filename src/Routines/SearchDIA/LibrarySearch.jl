"""
    LibrarySearch(spectra::Arrow.Table, params::NamedTuple; kwargs...)

Perform library search with required parameters specified in kwargs.

Required kwargs:
- frag_index: Fragment index for search
- precursors: Precursor information 
- fragment_lookup_table: Fragment lookup table
- rt_to_irt_spline: RT to iRT conversion function
- ms_file_idx: MS file index
- irt_tol: iRT tolerance
- ion_matches: Vector of fragment matches
- ion_misses: Vector of fragment misses
- fmatches: Fragment matches
- id_to_col: ID to column mapping
- ion_templates: Ion templates
- iso_splines: Isotope splines
- scored_psms: Scored PSMs
- unscored_psms: Unscored PSMs
- spectral_scores: Spectral scores
- prec_to_score: Precursor to score mapping
- mass_err_model: Mass error model
- sample_rate: Sample rate
- params: Search parameters
- isotope_err_bounds: Isotope error bounds
- quad_transmission_model: Quadrupole transmission model
"""
function LibrarySearch(
    spectra::Arrow.Table,
    params::NamedTuple;
    frag_index,
    precursors,
    fragment_lookup_table,
    rt_to_irt_spline,
    ms_file_idx::UInt32,
    irt_tol::Float64,
    ion_matches::Vector,
    ion_misses::Vector,
    fmatches,
    id_to_col::Vector,
    ion_templates::Vector,
    iso_splines,
    scored_psms::Vector,
    unscored_psms::Vector,
    spectral_scores::Vector,
    prec_to_score::Vector,
    mass_err_model,
    sample_rate::Float64,
    params::Dict,
    isotope_err_bounds::Tuple,
    quad_transmission_model
)
    # Partition scans across threads
    thread_tasks, _ = partitionScansToThreads(
        spectra[:mz_array],
        spectra[:retentionTime],
        spectra[:centerMz],
        spectra[:msOrder],
        Threads.nthreads(),
        1
    )

    # Initialize scan to precursor mapping
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra[:msOrder]))

    # Fragment index search
    precursors_passed_scoring = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            searchFragmentIndex(
                spectra,
                last(thread_task),
                frag_index,
                scan_to_prec_idx,
                rt_to_irt_spline,
                mass_err_model,
                searchScan!,
                prec_to_score[thread_id],
                Tuple(Int64.(isotope_err_bounds)),
                quad_transmission_model,
                UInt8(params["min_index_search_score"]),
                Float64(irt_tol),
                sample_rate,
                Set(2)
            )
        end |> fetch
    end

    # Determine precursor estimation strategy
    prec_estimation = params["abreviate_precursor_calc"] ? PartialPrecCapture() : FullPrecCapture()

    # PSM generation
    map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            getPSMS(
                spectra,
                last(thread_task),
                precursors,
                scan_to_prec_idx,
                precursors_passed_scoring[thread_id],
                fragment_lookup_table,
                rt_to_irt_spline,
                ms_file_idx,
                mass_err_model,
                quad_transmission_model,
                ion_matches[thread_id],
                ion_misses[thread_id],
                id_to_col[thread_id],
                ion_templates[thread_id],
                iso_splines,
                scored_psms[thread_id],
                unscored_psms[thread_id],
                spectral_scores[thread_id],
                Tuple(Int64.(isotope_err_bounds)),
                params["min_frag_count"],
                Float32(params["min_spectral_contrast"]),
                Float32(params["min_log2_matched_ratio"]),
                Tuple(Int64.(params["min_topn_of_m"])),
                params["max_best_rank"],
                Int64(params["n_frag_isotopes"]),
                UInt8(params["max_frag_rank"]),
                prec_estimation,
                Float32(irt_tol),
                Set(2)
            )
        end |> fetch
    end
end