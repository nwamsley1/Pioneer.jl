# Build the sorted precursor table with the streaming builder and compare it, column by column, to the in-memory
# path's precursors.arrow (from old_stages.jl). usage: julia --project=. compare_streaming.jl BUILD.json REF_DIR OUT_DIR
using Pioneer, JSON, Arrow, Tables
const P = Pioneer
cfg_path, ref_dir, out = ARGS
mkpath(out)
params = P.check_params_bsp(read(cfg_path, String))
lp = params["library_params"]
t = @timed P.with_koina_client(P.SyntheticKoinaClient(n_frags_per_prec = 20, realistic_ions = true)) do
    P.build_precursor_table_streaming(params, Float32(lp["prec_mz_min"]), Float32(lp["prec_mz_max"]),
        joinpath(out, "precursors.arrow"), joinpath(out, "proteins_table.arrow"))
end
println("streaming: $(round(t.time, digits = 1)) s, alloc $(round(t.bytes / 2^30, digits = 2)) GiB, peak RSS $(round(P.peak_rss() / 2^30, digits = 2)) GiB")
a = Arrow.Table(joinpath(ref_dir, "precursors.arrow")); b = Arrow.Table(joinpath(out, "precursors.arrow"))
na, nb = Tables.columnnames(a), Tables.columnnames(b)
println("rows ref/new: $(length(a.mz)) / $(length(b.mz)); columns equal: $(collect(na) == collect(nb))")
ok = length(a.mz) == length(b.mz) && collect(na) == collect(nb)
for n in na
    ca, cb = Tables.getcolumn(a, n), Tables.getcolumn(b, n)
    same = length(ca) == length(cb) && all(isequal(ca[i], cb[i]) for i in eachindex(ca))
    tsame = eltype(ca) == eltype(cb)
    if !(same && tsame)
        global ok = false
        i = findfirst(i -> !isequal(ca[i], cb[i]), 1:min(length(ca), length(cb)))
        println("  MISMATCH $n  types $(eltype(ca)) / $(eltype(cb))  first diff row $i: $(i === nothing ? "-" : (ca[i], cb[i]))")
    end
end
pa = Arrow.Table(joinpath(ref_dir, "proteins_table.arrow")); pb = Arrow.Table(joinpath(out, "proteins_table.arrow"))
psame = all(n -> isequal(collect(Tables.getcolumn(pa, n)), collect(Tables.getcolumn(pb, n))), Tables.columnnames(pa))
println("proteins_table identical: $psame")
println(ok && psame ? "IDENTICAL" : "DIFFERENT")
