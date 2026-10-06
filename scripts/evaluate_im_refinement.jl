#!/usr/bin/env julia
# julia --project scripts/evaluate_im_refinement.jl anchors.arrow library.poin output_dir [tdfs_dir]
# Run on a baseline search (global.im_refinement="none"). Prefer its retained
# temp_data/main_search_psms directory (output.delete_temp=false): these anchors
# are selected before IM-dependent experiment-wide scoring. Final long tables
# are also supported, excluding MBR and withheld observations.
using Pioneer, Arrow, DataFrames, CSV, Statistics, LinearAlgebra

length(ARGS) in (3, 4) || error("Usage: anchors.arrow library.poin output_dir [tdfs_dir]")
anchors_path, library_path, out_dir = ARGS[1:3]
mkpath(out_dir)
BLAS.set_num_threads(2)
if isdir(anchors_path)
    parts = DataFrame[]
    cols = [:target, :precursor_idx, :ms_file_idx, :scan_idx, :charge, :im_obs, :lgbm_prob]
    files = sort(readdir(anchors_path; join = true))
    for file in files
        endswith(file, ".arrow") || continue
        endswith(file, "sidecar.arrow") && continue
        # ScoringSearch replaces fold tables with a merged table. Use one
        # representation per run, and ignore score-only sidecars.
        base = replace(file, r"_fold\d+\.arrow$" => ".arrow")
        base != file && isfile(base) && continue
        tbl = Arrow.Table(file)
        all(col -> hasproperty(tbl, col), cols) || continue
        part = DataFrame([col => collect(getproperty(tbl, col)) for col in cols])
        part[!, :file_name] = fill(replace(basename(file), r"(_fold\d+)?\.arrow$" => ""), nrow(part))
        push!(parts, part)
    end
    isempty(parts) && error("No MainSearch fold tables found in $anchors_path")
    df = reduce(vcat, parts)
    df[!, :qval] = zeros(Float32, nrow(df))
    for run in groupby(df, :ms_file_idx)
        q = zeros(Float32, nrow(run))
        Pioneer.get_qvalues!(run.lgbm_prob, collect(Bool, run.target), q; doSort = true)
        run.qval .= q
    end
    df = df[df.lgbm_prob .> 0.9f0, :]
    score_column = :lgbm_prob
else
    df = DataFrame(Arrow.Table(anchors_path))
    df = df[.!df.mbr_recovered .& .!df.quant_withheld, :]
    score_column = :prec_prob
end
lib = Arrow.Table(joinpath(library_path, "precursors_table.arrow"))
df = df[df.target .& (df.qval .<= 0.01) .& isfinite.(df.im_obs), :]
sort!(df, score_column; rev = true)
unique!(df, [:ms_file_idx, :precursor_idx])
metrics = DataFrame(file = String[], charge = Int[], mode = String[], n = Int[],
    median_abs_error = Float64[], p95_abs_error = Float64[], p99_abs_error = Float64[],
    rmse = Float64[], bias = Float64[])
selection = DataFrame(file = String[], fold = Int[], mode = String[])
heldout = DataFrame[]
for run in groupby(df, :ms_file_idx)
    ids = run.precursor_idx
    sequences = lib.sequence[ids]
    tokens = Pioneer.im_composition_tokens(sequences, lib.structural_mods[ids])
    pred = Float64.(lib.inv_ion_mobility[ids])
    obs = Float64.(run.im_obs)
    fname = String(first(run.file_name))
    if length(ARGS) == 4
        # Instrument calibration is independent of peptide identifications:
        # this avoids reusing a global regression fitted on the evaluation fold.
        tdfs = joinpath(ARGS[4], fname * ".tdfs")
        slices = Arrow.Table(joinpath(tdfs, "slices.arrow"))
        meta = Pioneer.JSON.parsefile(joinpath(tdfs, "meta.json"))
        obs = Float64(meta["im_scan0_1overK0"]) .+
              Float64(meta["im_slope_1overK0_per_scan"]) .* slices.im_scan[run.scan_idx]
    end
    detail = DataFrame(precursor_idx = ids, ms_file_idx = run.ms_file_idx,
                       sequence = sequences, charge = run.charge, im_obs = obs, im_pred = pred)
    for mode in Pioneer.IM_REFINEMENT_MODES
        result = Pioneer.crossfit_im_correction(pred, obs, run.charge, tokens, sequences; mode = mode)
        detail[!, Symbol("prediction_", mode)] = result.refined
        detail[!, :fold] = result.folds
        mode === :auto && foreach(p -> push!(selection, (fname, first(p), String(last(p)))), sort(collect(result.selected)))
        for z in [0; sort(unique(Int.(run.charge)))]
            rows = z == 0 ? collect(eachindex(pred)) : findall(==(z), run.charge)
            residual = result.refined[rows] .- obs[rows]
            push!(metrics, (fname, z, String(mode), length(rows), median(abs.(residual)),
                            quantile(abs.(residual), 0.95), quantile(abs.(residual), 0.99),
                            sqrt(mean(abs2, residual)), mean(residual)))
        end
    end
    push!(heldout, detail)
    println(fname, ": ", nrow(run), " anchors")
    CSV.write(joinpath(out_dir, "metrics.tsv"), metrics; delim = '\t')
    CSV.write(joinpath(out_dir, "selection.tsv"), selection; delim = '\t')
end
all_rows = reduce(vcat, heldout)
Arrow.write(joinpath(out_dir, "heldout_predictions.arrow"), all_rows)
summary = DataFrame(mode = String[], n = Int[], median_abs_error = Float64[],
                    p95_abs_error = Float64[], p99_abs_error = Float64[], rmse = Float64[])
for mode in Pioneer.IM_REFINEMENT_MODES
    residual = all_rows[!, Symbol("prediction_", mode)] .- all_rows.im_obs
    push!(summary, (String(mode), nrow(all_rows), median(abs.(residual)),
                    quantile(abs.(residual), 0.95), quantile(abs.(residual), 0.99), sqrt(mean(abs2, residual))))
    println(mode, " median=", median(abs.(residual)), " p95=", quantile(abs.(residual), 0.95))
end
CSV.write(joinpath(out_dir, "summary.tsv"), summary; delim = '\t')
