# Reference for the streaming build: run the CURRENT (in-memory) BuildSpecLib stages up to the sorted precursor
# table (prepare_chronologer_input -> mock retention times -> parse_chronologer_output) and keep the intermediates.
# usage: julia --project=. old_stages.jl BUILD_CONFIG.json OUT_DIR
using Pioneer, JSON, Arrow, Tables
const P = Pioneer
cfg_path, out = ARGS[1], ARGS[2]
mkpath(out)
params = P.check_params_bsp(read(cfg_path, String))
lp = params["library_params"]
cin, cout = joinpath(out, "precursors_for_chronologer.arrow"), joinpath(out, "precursors_for_chronologer_rt.arrow")
P.with_koina_client(P.SyntheticKoinaClient(n_frags_per_prec = 20, realistic_ions = true)) do
    t = @timed P.prepare_chronologer_input(params, missing, Float32(lp["prec_mz_min"]), Float32(lp["prec_mz_max"]),
                                           cin, joinpath(out, "proteins_table.arrow"))
    println("prepare_chronologer_input: $(round(t.time, digits = 1)) s, peak RSS $(round(P.peak_rss() / 2^30, digits = 2)) GiB")
    P.predict_retention_times(cin, cout; rt_model = String(lp["rt_model"]))
    p = P.parse_chronologer_output(cout, out, Dict{String, Int8}(), Dict{String, Float32}(), params["isotope_mod_groups"], 3.0f0)
    println("sorted precursor table: $p")
end
t = Arrow.Table(joinpath(out, "precursors.arrow"))
println("rows: ", length(t.mz))
for (n, c) in zip(Tables.columnnames(t), Tables.columns(t)); println(rpad(string(n), 28), eltype(c)); end
