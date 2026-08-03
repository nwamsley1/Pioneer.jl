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

#==========================================================
PSM scoring (experiment-wide LightGBM, identical to MainSearch's per-file
classifier). Both stages call `train_psm_classifier_with_fallback`
(MainSearch/scoring.jl), so they share hyperparameters, the 2-fold CV split,
and the low-data probit fallback.
==========================================================#

"""
    score_precursor_isotope_traces(second_pass_folder, file_paths, precursors,
                                    fragment_lookup, max_psms_in_memory,
                                    q_value_threshold, force_oom;
                                    match_between_runs=true)

Score PSMs experiment-wide with the same LightGBM classifier the per-file
MainSearch uses (`SHARED_LGBM_HP` in `MainSearch/scoring.jl`). Loads PSMs
into memory, calls the shared training helper, and writes `:trace_prob` back
to the per-file Arrow files.

# Arguments
- `second_pass_folder`: Folder containing the per-file second-pass PSM
  Arrow files (output of MainSearch's prescore filter step).
- `file_paths`: Vector of per-file (or fold-split per-file) Arrow paths to
  write the scored output back to.
- `precursors`: Library precursor metadata. Retained in the calling convention; the
  `:accession_numbers` column is added later by `process_final_psms!`, not here.
- `fragment_lookup`: Retained in the internal calling convention; integrated
  MBR evidence is built later during chromatogram integration.
- `max_psms_in_memory`: Memory budget; surfaced for backward compatibility,
  currently unused (in-memory always; per-fold sub-sampling inside the
  shared helper handles large datasets).
- `q_value_threshold`: Surfaced for backward compatibility; currently unused.
- `force_oom`: Surfaced for backward compatibility; ignored.

# Returns
`nothing` (scores are written to file).
"""
function score_precursor_isotope_traces(
    second_pass_folder::String,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    ::LibraryFragmentLookup,
    ::Int64,                       # max_psms_in_memory (unused)
    ::Float32 = 0.01f0,            # q_value_threshold (unused)
    ::Bool = false;                # force_oom (unused)
    match_between_runs::Bool = true,
)
    # MBR-on freezes OOF and in-fold Pass-1 scores for later integrated
    # transfer rescoring. MBR-off keeps the legacy in-memory path.
    if match_between_runs
        return _score_precursor_isotope_traces_mbr(file_paths, precursors)
    else
        return _score_precursor_isotope_traces_no_mbr(
            second_pass_folder, file_paths, precursors,
        )
    end
end

# MBR-on Pass-1 path. Candidate selection, donor comparisons, and transfer
# rescoring now happen after the initial q-value filters and chromatogram
# integration. This stage only produces the frozen pre-MBR OOF and in-fold
# scores needed by those later steps.
function _score_precursor_isotope_traces_mbr(
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
)
    features = copy(ADVANCED_FEATURE_SET)
    pass1 = train_and_predict_pass1_oom!(
        file_paths;
        features        = features,
        compute_infold  = true,
        lgbm_hp         = SCORING_LGBM_HP,
        semisupervised  = true,
    )

    if pass1.last_classifier !== nothing
        lgbm_model = LightGBMModel(pass1.last_classifier, pass1.available_features, nothing)
        imp = importance(lgbm_model)
        if imp !== nothing
            sorted_imp = sort(imp, by = x -> -x[2])
            lines = ["ScoringSearch Pass-1 LGBM feature gains (all $(length(sorted_imp))):"]
            for (fname, gain) in sorted_imp
                push!(lines, "    $(rpad(string(fname), 40)) $(round(Int, gain))")
            end
            @debug_l1 join(lines, "\n")
        end
    end

    _merge_pass1_into_main!(file_paths, precursors)
    return nothing
end

# Legacy in-memory path used only when match_between_runs = false. Kept for
# small / non-MBR runs where the full best_psms easily fits in memory.
# Pass-1 sidecar merge. After `train_and_predict_pass1_oom!` has written each
# file's `.pass1_sidecar.arrow` (OOF `trace_prob_prepass`), fold the Pass-1 scores
# into each main file ONE FILE AT A TIME — the full experiment is never materialised.
# Reproduces the exact output columns the legacy in-memory path wrote:
# `:decoy`, `:trace_prob_prepass`, `:trace_prob` (= prepass),
# `:mbr_recovered` (= false).
function _merge_pass1_into_main!(
    file_paths::Vector{String},
    precursors::LibraryPrecursors;
    cleanup::Bool = true,
)
    for path in file_paths
        pass1_path = path * PASS1_SIDECAR_SUFFIX
        isfile(pass1_path) || continue
        # Materialise one file (Tables.columntable detaches from the mmap so the same
        # path can be safely rewritten below).
        main  = DataFrame(Tables.columntable(Arrow.Table(path)))
        pass1 = Arrow.Table(pass1_path)
        n = nrow(main)
        length(pass1.precursor_idx) == n ||
            error("_merge_pass1_into_main!: row-count mismatch at $path")
        @inbounds for i in 1:n
            (main.precursor_idx[i] == pass1.precursor_idx[i] &&
             main.scan_idx[i]      == pass1.scan_idx[i]) ||
                error("_merge_pass1_into_main!: sidecar misaligned at row $i of $path")
        end
        # :accession_numbers deliberately NOT added here. It was materialised per PSM row at this point
        # and then overwritten wholesale by process_final_psms! (IntegrateChromatogramsSearch/utils.jl),
        # so the early copy was carried through every intermediate read/write in between -- the MBR
        # feature pass, merge_recoveries, and the finalize materialise -- and then discarded. Protein
        # inference does not need it either: ProteinInferenceSearch reads getAccessionNumbers(precursors)
        # and indexes by precursor_idx (utils.jl:157/175/423) rather than using the table column.
        main[!, :decoy]              = main[!, :target] .== false
        main[!, :trace_prob_prepass] = collect(Float32.(pass1.trace_prob_prepass))
        if hasproperty(pass1, :trace_prob_infold)
            main[!, :trace_prob_infold] =
                collect(Float32.(pass1.trace_prob_infold))
        end
        main[!, :trace_prob]         = main[!, :trace_prob_prepass]
        main[!, :mbr_recovered]      = falses(n)
        pass1 = nothing; GC.gc(false)   # release sidecar mmap before rm + rewrite
        writeArrow(path, main)
        cleanup && safeRm(pass1_path, nothing; force=true)
    end
    return nothing
end

# MBR-off path. Streams Pass-1 LightGBM over the per-file Arrow tables via the same
# OOM trainer the MBR path uses, instead of materialising the whole experiment into
# one in-memory DataFrame (the legacy path's peak-RSS cost). ID-equivalent — not
# bit-identical — to the legacy path: the OOM trainer reservoir-samples for training,
# so scores differ in the last bits, but the full experiment is never held in RAM.
function _score_precursor_isotope_traces_no_mbr(
    second_pass_folder::String,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
)
    features = copy(ADVANCED_FEATURE_SET)
    pass1 = train_and_predict_pass1_oom!(
        file_paths;
        features        = features,
        compute_infold  = false,        # no MBR-FTR downstream consumes trace_prob_infold
        lgbm_hp         = SCORING_LGBM_HP,
        semisupervised  = true,
    )
    @debug_l1 "Pass-1 (no MBR, streamed) trained on $(length(pass1.available_features)) features"

    if pass1.last_classifier !== nothing
        lgbm_model = LightGBMModel(pass1.last_classifier, pass1.available_features, nothing)
        imp = importance(lgbm_model)
        if imp !== nothing
            sorted_imp = sort(imp, by = x -> -x[2])
            lines = ["ScoringSearch Pass-1 LGBM feature gains (all $(length(sorted_imp))):"]
            for (fname, gain) in sorted_imp
                push!(lines, "    $(rpad(string(fname), 40)) $(round(Int, gain))")
            end
            @debug_l1 join(lines, "\n")
        end
    end

    _merge_pass1_into_main!(file_paths, precursors)
    return nothing
end

"""
    get_psms_count(file_paths::Vector{String}) -> Int

Count total PSMs across the given Arrow files (used for diagnostic logging).
"""
function get_psms_count(file_paths::Vector{String})
    psms_count = 0
    for file_path in file_paths
        psms_count += length(Arrow.Table(file_path)[1])
    end
    return psms_count
end

"""
    load_psms_for_lightgbm(quant_psms_folder::String;
                           fold::Union{Nothing,UInt8}=nothing) -> DataFrame

Load PSMs from Arrow files in `quant_psms_folder` for experiment-wide
LightGBM training.

When `fold` is specified, only loads files ending with `_fold{fold}.arrow`
(memory-saving when only one fold is needed). When `fold` is `nothing`,
loads all `.arrow` files.
"""
function load_psms_for_lightgbm(quant_psms_folder::String;
                                fold::Union{Nothing,UInt8}=nothing)
    if fold !== nothing
        file_paths = [f for f in readdir(quant_psms_folder, join=true)
                      if endswith(f, "_fold$(fold).arrow")]
    else
        file_paths = [f for f in readdir(quant_psms_folder, join=true)
                      if endswith(f, ".arrow")]
    end
    return DataFrame(Tables.columntable(Arrow.Table(file_paths)))
end

"""
    write_scored_psms_to_files!(psms::DataFrame, file_paths::Vector{String})

Write scored PSMs back to Arrow files, grouped by ms_file_idx and cv_fold.
Supports both single-file-per-MS-run format (legacy) and fold-split format.
For fold-split files (containing "_fold" in path), groups by
(ms_file_idx, cv_fold). For legacy files, groups by ms_file_idx only.

# Arguments
- `psms`: DataFrame containing scored PSMs with `ms_file_idx` (and `cv_fold`
  for fold-split mode).
- `file_paths`: Vector of file paths for valid files only.
"""
function write_scored_psms_to_files!(psms::DataFrame, file_paths::Vector{String})
    dropVectorColumns!(psms)  # avoids writing issues

    # Detect if we're using fold-split files
    is_fold_split = any(p -> occursin("_fold", p), file_paths)

    if is_fold_split && hasproperty(psms, :cv_fold)
        # Fold-split mode: build (ms_file_idx, cv_fold) → path mapping by
        # peeking at each file's ms_file_idx.
        path_to_key = Dict{String, Tuple{UInt32, UInt8}}()
        for fpath in file_paths
            fold_match = match(r"_fold(\d)\.arrow$", fpath)
            if fold_match !== nothing
                fold_num = parse(UInt8, fold_match.captures[1])
                orig_df = DataFrame(Arrow.Table(fpath))
                if nrow(orig_df) > 0
                    ms_idx = first(orig_df.ms_file_idx)
                    path_to_key[fpath] = (ms_idx, fold_num)
                end
            end
        end

        key_to_path = Dict{Tuple{UInt32, UInt8}, String}()
        for (fpath, key) in path_to_key
            key_to_path[key] = fpath
        end

        for (key, gpsms) in pairs(groupby(psms, [:ms_file_idx, :cv_fold]))
            ms_idx = key[:ms_file_idx]
            cv_fold = key[:cv_fold]
            lookup_key = (ms_idx, cv_fold)
            if haskey(key_to_path, lookup_key)
                fpath = key_to_path[lookup_key]
                writeArrow(fpath, gpsms)
            else
                @warn "No output path found for ms_file_idx=$ms_idx, cv_fold=$cv_fold, skipping"
            end
        end
    else
        # Legacy mode: group by ms_file_idx only
        unique_file_indices = unique(psms[:, :ms_file_idx])
        sort!(unique_file_indices)

        if length(file_paths) != length(unique_file_indices)
            error("Mismatch: $(length(file_paths)) file paths provided but $(length(unique_file_indices)) unique file indices found in PSM data")
        end

        index_to_path = Dict(zip(unique_file_indices, file_paths))

        for (ms_file_idx, gpsms) in pairs(groupby(psms, :ms_file_idx))
            file_idx = ms_file_idx[:ms_file_idx]
            if haskey(index_to_path, file_idx)
                fpath = index_to_path[file_idx]
                writeArrow(fpath, gpsms)
            else
                @warn "No output path found for file index $file_idx, skipping"
            end
        end
    end
end
