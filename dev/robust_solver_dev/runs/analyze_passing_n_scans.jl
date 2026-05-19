# For each file, fraction of PASSING precursors (07_n_scans_feature q<0.01)
# whose best-PSM n_scans is < 3 / < 2 / = 1.

using Arrow, DataFrames, Printf, Tables

const RUN07    = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/07_n_scans_feature"
const LIB_PATH = "/Users/nathanwamsley/Data/SPEC_LIBS/altimeter_3P_len7o40_ch2o3_mc1_OlsenAstral_mzsorted.poin/precursors_table.arrow"
const FILES = ["Sample-A_Rep1","Sample-A_Rep2","Sample-A_Rep3","Sample-B_Rep1","Sample-B_Rep2","Sample-B_Rep3"]

_norm(x) = ismissing(x) ? "" : String(x)
key3(seq,ch,sm,dc) = (String(seq), UInt8(ch), _norm(sm), Bool(dc))

println("[load] library + 07 precursors_long")
lib = DataFrame(Arrow.Table(LIB_PATH))
lib_key3 = Dict{Tuple{String,UInt8,String,Bool}, UInt32}()
@inbounds for i in 1:nrow(lib)
    lib_key3[key3(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i], lib.is_decoy[i])] = UInt32(i)
end

pl = DataFrame(Arrow.Table(joinpath(RUN07,"precursors_long.arrow")))
qcol = hasproperty(pl, :MBR_boosted_qval) ? :MBR_boosted_qval : :qval
pl = pl[(pl[!, qcol] .< 0.01) .& pl.target, :]
pass_by_file = Dict{String, Set{UInt32}}()
for sub in groupby(pl, :file_name)
    fn_full = String(first(sub.file_name))
    m = match(r"Sample-([AB]).*Rep(\d+)", fn_full)
    fn = m === nothing ? fn_full : "Sample-$(m.captures[1])_Rep$(m.captures[2])"
    s = Set{UInt32}()
    for r in eachrow(sub)
        pid = get(lib_key3, key3(r.sequence, r.charge, r.structural_mods, r.is_decoy), UInt32(0))
        pid != 0 && push!(s, pid)
    end
    pass_by_file[fn] = s
end

println("\n─── n_scans distribution among passing precursors (q<0.01) ───")
@printf "%-16s  %8s  %8s  %8s  %8s  %8s  %8s\n" "file" "n_pass" "=1" "%=1" "<2" "<3" "<5"

total_pass = 0; total_eq1 = 0; total_lt2 = 0; total_lt3 = 0; total_lt5 = 0
for fn in FILES
    P = pass_by_file[fn]
    f0 = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN07,"temp_data","main_search_psms","$(fn)_fold0.arrow"))))
    f1 = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN07,"temp_data","main_search_psms","$(fn)_fold1.arrow"))))
    main = vcat(f0, f1; cols=:union)
    sub = main[in.(main.precursor_idx, Ref(P)), :]
    ns = sub.n_scans
    n = length(ns)
    eq1 = count(==(UInt32(1)), ns)
    lt2 = count(<(UInt32(2)), ns)
    lt3 = count(<(UInt32(3)), ns)
    lt5 = count(<(UInt32(5)), ns)
    @printf "%-16s  %8d  %8d  %7.2f%%  %8d  %8d  %8d\n" fn n eq1 100*eq1/max(n,1) lt2 lt3 lt5
    global total_pass += n; global total_eq1 += eq1; global total_lt2 += lt2; global total_lt3 += lt3; global total_lt5 += lt5
end
println("─"^80)
@printf "%-16s  %8d  %8d  %7.2f%%  %8d  %8d  %8d\n" "TOTAL" total_pass total_eq1 100*total_eq1/max(total_pass,1) total_lt2 total_lt3 total_lt5
println()
@printf "Across-file: =1: %.2f%%, <2: %.2f%%, <3: %.2f%%, <5: %.2f%%\n" 100*total_eq1/total_pass 100*total_lt2/total_pass 100*total_lt3/total_pass 100*total_lt5/total_pass
