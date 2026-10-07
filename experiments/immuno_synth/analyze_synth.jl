# EXPERIMENT: summarize synthetic-index harness outputs. usage: julia analyze_synth.jl SYNTH_RUNS_DIR
using Arrow, DataFrames, JSON, CSV, Printf
d = ARGS[1]
ids = Arrow.Table(joinpath(d, "olsen_e40_01_ids.arrow")).precursor_idx
rows = DataFrame()
for o in sort(filter(x -> startswith(x, "out_"), readdir(d)))
    od = joinpath(d, o)
    sj = filter(x -> endswith(x, "_summary.json"), readdir(od)); isempty(sj) && continue
    s = JSON.parsefile(joinpath(od, sj[1]))
    c = Arrow.Table(joinpath(od, replace(sj[1], "_summary.json" => "_exact_candidates_synthLUT.arrow")))
    hits = Dict{UInt32, Int32}()
    for p in c.precursor_idx; hits[p] = get(hits, p, 0) + 1; end
    h = [get(hits, p, 0) for p in ids]
    nscan = s["n_ms2_scans"]
    push!(rows, (run = o[5:end], n_prec = s["n_precursors"], n_fake = s["n_fake"], lut_pass = s["lut_synth_pass"],
        emit_postLUT = s["emit_postLUT_total"], exact = s["exact_total"], exact_per_scan = s["exact_total"] / nscan,
        real_T = s["exact_real_target"], real_D = s["exact_real_decoy"], fake_T = s["exact_fake_target"], fake_D = s["exact_fake_decoy"],
        uniq_real_T = s["unique_real_target"], uniq_real_D = s["unique_real_decoy"],
        uniq_fake_T = s["unique_fake_target"], uniq_fake_D = s["unique_fake_decoy"],
        ids_ge1 = round(100 * count(>=(1), h) / length(ids), digits = 1), ids_ge3 = round(100 * count(>=(3), h) / length(ids), digits = 1),
        fi_s = round(s["frag_index_s"], digits = 1), load_s = round(s["load_s"], digits = 1), maxrss_gb = round(s["maxrss_gb"], digits = 1)))
end
CSV.write(joinpath(d, "synth_summary.csv"), rows)
show(stdout, rows; allcols = true, truncate = 0); println()
