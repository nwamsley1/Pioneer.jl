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
PSM scoring uses a fixed experiment-wide training pool with the existing
2-fold CV assignments. Both MBR modes fit and select LightGBM models on the
pool, then stream predictions over every source file.
==========================================================#

"""
    score_precursor_isotope_traces(file_paths; match_between_runs=true)

Fit semi-supervised LightGBM models on a fixed representative pool, then score
all input rows into aligned Pass-1 sidecars. MBR also requires in-fold scores;
without MBR only out-of-fold scores are computed. Predictions are attached during
fold merging, before experiment-wide FDR calibration.
"""
function score_precursor_isotope_traces(
    file_paths::Vector{String};
    match_between_runs::Bool = true,
)
    features = copy(ADVANCED_FEATURE_SET)
    pass1 = train_and_predict_pass1_oom!(
        file_paths;
        features        = features,
        compute_infold  = match_between_runs,
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

    return nothing
end

"""
    _merge_scored_folds!(fold_paths, merged_path)

Attach Pass-1 predictions while merging one run's folds in input order. Return
the row count and source paths for cleanup after the combined file is written.
"""
function _merge_scored_folds!(
    fold_paths::Vector{String},
    merged_path::String;
    sidecar_index = index_sidecar_paths(fold_paths),
)
    read_seconds = attach_seconds = 0.0
    fold_dfs = DataFrame[]
    cleanup_paths = String[]
    for path in fold_paths
        isfile(path) || continue
        pass1_path = path * PASS1_SIDECAR_SUFFIX
        isfile(pass1_path) || error("Missing Pass-1 predictions for $path")
        read_started = time()
        table = Arrow.Table(path)
        ref = PSMFileReference(path; table, sidecar_paths=sidecar_index[path])
        main = DataFrame(table; copycols=false)
        pass1 = Arrow.Table(pass1_path)
        read_seconds += time() - read_started
        attach_started = time()
        n = nrow(main)
        length(pass1.precursor_idx) == n ||
            error("Pass-1 row-count mismatch at $path")
        (main.precursor_idx == pass1.precursor_idx && main.scan_idx == pass1.scan_idx) ||
            error("Pass-1 sidecar misaligned at $path")
        main[!, :decoy]              = main[!, :target] .== false
        main[!, :trace_prob_prepass] = collect(Float32.(pass1.trace_prob_prepass))
        if hasproperty(pass1, :trace_prob_infold)
            main[!, :trace_prob_infold] =
                collect(Float32.(pass1.trace_prob_infold))
        end
        main[!, :trace_prob]         = main[!, :trace_prob_prepass]
        main[!, :mbr_recovered]      = falses(n)
        attach_seconds += time() - attach_started
        for sidecar in ref.sidecars
            read_started = time()
            table = Arrow.Table(sidecar.path)
            for name in sidecar.cols
                hasproperty(main, name) && continue
                main[!, name] = collect(Tables.getcolumn(table, name))
            end
            read_seconds += time() - read_started
        end
        push!(fold_dfs, main)
        append!(cleanup_paths, (path, pass1_path))
        append!(cleanup_paths, (sidecar.path for sidecar in ref.sidecars))
    end
    isempty(fold_dfs) && return nothing
    concatenate_started = time()
    combined = vcat(fold_dfs...)
    concatenate_seconds = time() - concatenate_started
    write_started = time()
    writeArrow(merged_path, combined; temp_dir=dirname(merged_path))
    write_seconds = time() - write_started
    return (; rows = nrow(combined), cleanup_paths,
            read_seconds, attach_seconds, concatenate_seconds, write_seconds)
end
