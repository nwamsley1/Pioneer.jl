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
    FirstPassSearch

First pass deconvolution and LightGBM prescore search.

Pipeline:
1. process_file!: Load fragment_index_matches, deconvolve with prescore fragment settings
2. process_search_results!: Compute prescore features, train LightGBM, select best scan per precursor
3. summarize_results!: Global prescore aggregation, filter fragment_index_matches to passing precursors,
   write filtered_fragment_matches for SecondPassSearch
"""
struct FirstPassSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for first pass prescore search.
"""
struct FirstPassSearchResults <: SearchResults
    psms::Base.Ref{DataFrame}
end

#==========================================================
Interface Implementation
==========================================================#

get_parameters(::FirstPassSearch, params::Any) = SecondPassSearchParameters(params)

function init_search_results(::FirstPassSearch, ::P, search_context::SearchContext) where {P<:SecondPassSearchParameters}
    prescore_dir = joinpath(getDataOutDir(search_context), "temp_data", "prescore_scores")
    mkpath(prescore_dir)
    filtered_dir = joinpath(getDataOutDir(search_context), "temp_data", "filtered_fragment_matches")
    mkpath(filtered_dir)
    return FirstPassSearchResults(
        DataFrame()
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single file: load fragment index matches and deconvolve with prescore settings.
"""
function process_file!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:SecondPassSearchParameters}

    if check_and_skip_failed_file(search_context, ms_file_idx, "FirstPassSearch")
        return results
    end

    try
        t_start = time()

        # Load fragment index mapping from FragmentIndexSearch
        frag_match_path = getFragmentIndexMatches(getMSData(search_context), ms_file_idx)
        scan_to_prec_idx, precursors_passed = load_fragment_index_matches(
            frag_match_path, length(spectra)
        )
        t_load = time()

        # Deconvolve with Phase 1 prescore fragment settings
        # Profile the deconvolution step (first file only)
        profile_dir = "/Users/nathanwamsley/Desktop"
        do_profile = (ms_file_idx == 1)
        if do_profile
            Profile.clear()
            Profile.init(n = 10_000_000, delay = 0.0001)  # 100μs sampling
        end

        if do_profile
            Profile.@profile search_result = perform_second_pass_search(
                spectra,
                scan_to_prec_idx,
                precursors_passed,
                search_context,
                params,
                ms_file_idx,
                MS2CHROM();
                n_frag_isotopes = params.prescore_n_frag_isotopes,
                max_frag_rank = params.prescore_max_frag_rank,
                min_frag_count = params.prescore_min_frag_count,
                min_spectral_contrast = params.prescore_min_spectral_contrast,
                min_log2_matched_ratio = params.prescore_min_log2_matched_ratio,
                min_topn_of_m = params.prescore_min_topn_of_m
            )
        else
            search_result = perform_second_pass_search(
                spectra,
                scan_to_prec_idx,
                precursors_passed,
                search_context,
                params,
                ms_file_idx,
                MS2CHROM();
                n_frag_isotopes = params.prescore_n_frag_isotopes,
                max_frag_rank = params.prescore_max_frag_rank,
                min_frag_count = params.prescore_min_frag_count,
                min_spectral_contrast = params.prescore_min_spectral_contrast,
                min_log2_matched_ratio = params.prescore_min_log2_matched_ratio,
                min_topn_of_m = params.prescore_min_topn_of_m
            )
        end

        if do_profile
            # Save profile to Desktop
            pprof_path = joinpath(profile_dir, "pioneer_firstpass_profile.pb.gz")
            PProf.pprof(out = pprof_path, web = false)
            @info "PProf profile saved to $pprof_path"
            @info "View with: julia -e 'using PProf; PProf.pprof(\"$pprof_path\"; web=true)'"

            # Also save a text summary for analysis
            txt_path = joinpath(profile_dir, "pioneer_firstpass_profile.txt")
            open(txt_path, "w") do io
                Profile.print(io; mincount=10, noisefloor=2)
            end
            @info "Profile text summary saved to $txt_path"

            # Save flat profile (sorted by self time)
            flat_path = joinpath(profile_dir, "pioneer_firstpass_profile_flat.txt")
            open(flat_path, "w") do io
                Profile.print(io; format=:flat, sortedby=:count, mincount=10)
            end
            @info "Flat profile saved to $flat_path"
        end

        psms = search_result.psms
        t_deconv = time()

        results.psms[] = psms

        # Save precursors-per-scan histogram to QC folder
        try
            if nrow(psms) > 0
                counts_per_scan = combine(groupby(psms, :scan_idx), nrow => :n_precursors)
                file_name = getParsedFileName(search_context, ms_file_idx)
                qc_dir = joinpath(getDataOutDir(search_context), "qc_plots")
                mkpath(qc_dir)
                p = Plots.histogram(counts_per_scan.n_precursors,
                    xlabel = "Precursors per scan",
                    ylabel = "Number of scans",
                    title = "FirstPass: $file_name",
                    legend = false,
                    bins = range(0, maximum(counts_per_scan.n_precursors) + 1, step=1),
                    color = :steelblue
                )
                Plots.savefig(p, joinpath(qc_dir, "firstpass_precs_per_scan_$(file_name).pdf"))
            end
        catch e
            @warn "Failed to save precursors-per-scan histogram" exception=(e, catch_backtrace())
        end

        # Save solveOLS iteration QC plots
        try
            iter_counts = search_result.iter_counts
            col_counts = search_result.col_counts
            if !isempty(iter_counts)
                file_name = getParsedFileName(search_context, ms_file_idx)
                qc_dir = joinpath(getDataOutDir(search_context), "qc_plots")
                mkpath(qc_dir)

                p_hist = Plots.histogram(iter_counts;
                    xlabel = "Outer iterations",
                    ylabel = "Count",
                    title = "solveOLS iterations ($(file_name), n=$(length(iter_counts)))",
                    legend = false,
                    bins = min(100, max(20, div(length(iter_counts), 50)))
                )
                Plots.savefig(p_hist, joinpath(qc_dir, "solveOLS_iter_hist_$(file_name).pdf"))

                p_scatter = Plots.scatter(col_counts, iter_counts;
                    xlabel = "Number of columns (precursors)",
                    ylabel = "Outer iterations",
                    title = "solveOLS cols vs iters ($(file_name))",
                    legend = false,
                    markersize = 2,
                    markeralpha = 0.3
                )
                Plots.savefig(p_scatter, joinpath(qc_dir, "solveOLS_cols_vs_iters_$(file_name).pdf"))
            end
        catch e
            @warn "Failed to save solveOLS QC plots" exception=(e, catch_backtrace())
        end

        println()
        @info "FirstPassSearch deconvolution: $(nrow(psms)) PSMs, load=$(round(t_load - t_start, digits=2))s, deconv=$(round(t_deconv - t_load, digits=2))s"
        println()

    catch e
        handle_search_error!(search_context, ms_file_idx, "FirstPassSearch", e, createFirstPassFallbackResults!, results)
    end

    return results
end

"""
Create empty PSM results for a failed file.
"""
function createFirstPassFallbackResults!(results::FirstPassSearchResults, ms_file_idx::Int64)
    results.psms[] = DataFrame(
        ms_file_idx = UInt32[],
        scan_idx = UInt32[],
        precursor_idx = UInt32[],
        rt = Float32[],
        q_value = Float32[],
        score = Float32[],
        prob = Float32[]
    )
end

"""
Per-file scoring: compute prescore features, train LightGBM, select best scan per precursor.
"""
function process_search_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:SecondPassSearchParameters}

    if check_and_skip_failed_file(search_context, ms_file_idx, "FirstPassSearch results processing")
        return nothing
    end

    try
        t_start = time()
        psms = results.psms[]
        file_name = getParsedFileName(search_context, ms_file_idx)

        # Compute only prescore features (skip columns recomputed in Phase 2)
        prepare_psm_features!(psms, params, search_context, ms_file_idx, spectra, prescore_only=true)
        t_features = time()

        if nrow(psms) == 0
            @debug_l2 "No PSMs for file $ms_file_idx after feature computation"
            return nothing
        end

        # Train LightGBM on ALL PSMs, select best scan per precursor
        best_psms, scores, q_values, lgbm_timings = train_lgbm_and_select_best(psms; min_frag_count = params.prescore_min_frag_count)
        t_lgbm = time()

        # RT recalibration: refit iRT spline from high-confidence PSMs
        recalibrate_rt!(search_context, ms_file_idx, best_psms, scores)

        # Write prescore table (scores only, NO fold Arrow files yet)
        prescore_dir = joinpath(getDataOutDir(search_context), "temp_data", "prescore_scores")
        mkpath(prescore_dir)
        score_df = DataFrame(
            precursor_idx = best_psms[!, :precursor_idx],
            lgbm_prob = scores,
            target = best_psms[!, :target],
            irt_obs = best_psms[!, :irt_obs],
            irt_pred = best_psms[!, :irt_pred],
            rt = best_psms[!, :rt],
            scan_idx = best_psms[!, :scan_idx]
        )
        writeArrow(joinpath(prescore_dir, "$(file_name).arrow"), score_df)
        t_write = time()

        # Timing summary
        n_pass_1 = count((q_values .<= 0.01) .& best_psms[!, :target])
        n_pass_5 = count((q_values .<= 0.05) .& best_psms[!, :target])
        t_total = t_write - t_start
        r = s -> round(s, digits=2)
        println()
        @info "FirstPassSearch scoring: $(nrow(psms)) PSMs → $(nrow(best_psms)) precursors ($n_pass_1 @ 1%, $n_pass_5 @ 5% FDR)\n" *
              "  features=$(r(t_features - t_start))s, write=$(r(t_write - t_lgbm))s\n" *
              "  lgbm: matrix=$(r(lgbm_timings.matrix))s, train=$(r(lgbm_timings.train))s, predict=$(r(lgbm_timings.predict))s, select=$(r(lgbm_timings.select))s\n" *
              "  total=$(r(t_total))s"
        println()

    catch e
        file_name = try
            getMassSpecData(search_context).file_id_to_name[ms_file_idx]
        catch
            "file_$ms_file_idx"
        end
        reason = "FirstPassSearch scoring failed: $(typeof(e))"
        markFileFailed!(search_context, ms_file_idx, reason)
        setFailedIndicator!(getMSData(search_context), ms_file_idx, true)
        @user_warn "First pass scoring failed for MS data file: $file_name. Error: $(sprint(showerror, e))\n$(sprint(Base.show_backtrace, catch_backtrace()))"
    end

    return nothing
end

"""
Reset results containers.
"""
function reset_results!(results::FirstPassSearchResults)
    results.psms[] = DataFrame()
end

"""
Global aggregation: aggregate prescore across files, filter fragment_index_matches
to passing precursors, write filtered_fragment_matches for SecondPassSearch.
"""
function summarize_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:SecondPassSearchParameters}

    @info "=== FirstPassSearch: Global prescore aggregation + fragment index filtering ==="

    # Step 1: Global prescore aggregation → passing precursor set + Phase 1 scan lookup
    passing_precs, prec_best_scan = aggregate_prescore_globally!(search_context, params.global_prescore_qvalue_threshold, params.prescore_aggregation)

    # Store results for SecondPassSearch to retrieve Phase 1 best scans
    store_results!(search_context, FirstPassSearch, (passing_precs=passing_precs, prec_best_scan=prec_best_scan))

    # Step 2: Filter fragment_index_matches to only passing precursors, write filtered versions
    ms_data = getMSData(search_context)
    n_files = length(ms_data)
    filtered_dir = joinpath(getDataOutDir(search_context), "temp_data", "filtered_fragment_matches")
    mkpath(filtered_dir)

    n_filtered_files = 0
    n_total_entries = 0
    n_kept_entries = 0

    for ms_file_idx in 1:n_files
        getFailedIndicator(ms_data, ms_file_idx) && continue

        frag_match_path = getFragmentIndexMatches(ms_data, ms_file_idx)
        isempty(frag_match_path) && continue
        !isfile(frag_match_path) && continue

        file_name = getParsedFileName(ms_data, ms_file_idx)

        # Load original fragment index matches
        tbl = DataFrame(Tables.columntable(Arrow.Table(frag_match_path)))
        n_before = nrow(tbl)
        n_total_entries += n_before

        # Filter to only globally-passing precursors
        filter!(row -> row.precursor_idx in passing_precs, tbl)
        n_after = nrow(tbl)
        n_kept_entries += n_after

        # Write filtered version
        filtered_path = joinpath(filtered_dir, "$(file_name).arrow")
        writeArrow(filtered_path, tbl)
        setFilteredFragmentMatches!(ms_data, ms_file_idx, filtered_path)

        n_filtered_files += 1
        pct = round(100.0 * n_after / max(1, n_before), digits=1)
        @debug "  $file_name: $n_after / $n_before entries kept ($pct%)"
    end

    overall_pct = round(100.0 * n_kept_entries / max(1, n_total_entries), digits=1)
    @info "Filtered fragment matches: $n_kept_entries / $n_total_entries entries ($overall_pct%) across $n_filtered_files files"

    return nothing
end
