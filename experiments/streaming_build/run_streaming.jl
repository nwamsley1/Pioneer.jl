# Build only the sorted precursor table with the streaming builder (mock Koina), reporting time and memory.
# usage: julia --project=. run_streaming.jl BUILD.json OUT_DIR
using Pioneer, Arrow
const P = Pioneer
cfg_path, out = ARGS
mkpath(out)
params = P.check_params_bsp(read(cfg_path, String))
lp = params["library_params"]
t = @timed P.with_koina_client(P.SyntheticKoinaClient(n_frags_per_prec = 20, realistic_ions = true)) do
    P.build_precursor_table_streaming(params, Float32(lp["prec_mz_min"]), Float32(lp["prec_mz_max"]),
        joinpath(out, "precursors.arrow"), joinpath(out, "proteins_table.arrow"))
end
println("streaming: $(round(t.time / 60, digits = 1)) min, alloc $(round(t.bytes / 2^30, digits = 1)) GiB, peak RSS $(round(P.peak_rss() / 2^30, digits = 2)) GiB, rows $(length(Arrow.Table(joinpath(out, "precursors.arrow")).mz))")
