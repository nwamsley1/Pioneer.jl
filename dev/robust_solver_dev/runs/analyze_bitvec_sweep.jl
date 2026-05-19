using Arrow, DataFrames, Printf
_norm(x) = ismissing(x) ? "" : String(x)
key3(seq,ch,sm,dc) = (String(seq), UInt8(ch), _norm(sm), Bool(dc))
lib = DataFrame(Arrow.Table("/Users/nathanwamsley/Data/SPEC_LIBS/altimeter_3P_len7o40_ch2o3_mc1_OlsenAstral_mzsorted.poin/precursors_table.arrow"))
lib_key3 = Dict{Tuple{String,UInt8,String,Bool}, UInt32}()
@inbounds for i in 1:nrow(lib); lib_key3[key3(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i], lib.is_decoy[i])] = UInt32(i); end
function loadpass(p)
    df = DataFrame(Arrow.Table(p))
    qcol = hasproperty(df, :MBR_boosted_qval) ? :MBR_boosted_qval : :qval
    df[(df[!, qcol] .< 0.01) .& df.target, :]
end
hupo = loadpass("/Users/nathanwamsley/Data/HUPO2026/MtacThreeProteomeStandard_hupo2026_v0o6o6/precursors_long.arrow")
function pid_set(df)
    s = Dict{String, Set{UInt32}}()
    for sub in groupby(df, :file_name)
        fn_full = String(first(sub.file_name))
        m = match(r"Sample-([AB]).*Rep(\d+)", fn_full)
        fn = m === nothing ? fn_full : "Sample-$(m.captures[1])_Rep$(m.captures[2])"
        ss = Set{UInt32}()
        for r in eachrow(sub)
            pid = get(lib_key3, key3(r.sequence, r.charge, r.structural_mods, r.is_decoy), UInt32(0))
            pid != 0 && push!(ss, pid)
        end
        s[fn] = ss
    end; s
end
H = pid_set(hupo)
files = sort(collect(keys(H)))
hupo_n6 = reduce(intersect, [H[fn] for fn in files])

runs = [
    ("22 (0.20)",  "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/22_bitvec_20/precursors_long.arrow"),
    ("21 (0.10)",  "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/21_bitvec_10/precursors_long.arrow"),
    ("20 (0.05)",  "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/20_bitvec_05/precursors_long.arrow"),
    ("19 (0.03)",  "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/19_bitvec_03/precursors_long.arrow"),
    ("23 (0.02)",  "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/23_bitvec_02/precursors_long.arrow"),
    ("24 (0.01)",  "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/24_bitvec_01/precursors_long.arrow"),
    ("25 (0.005)", "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/25_bitvec_005/precursors_long.arrow"),
]

@printf "%-12s  %12s  %12s  %12s  %12s\n" "run (rate)" "tot_q<0.01" "tot_∩HUPO" "%N=6_cov" "Δ vs 0.03"
base_total = 0
for (label, path) in runs
    S = pid_set(loadpass(path))
    tot_pass = sum(length(s) for s in values(S))
    tot_in_hupo = 0
    n6_in = 0
    for fn in files
        s = get(S, fn, Set())
        h = get(H, fn, Set())
        tot_in_hupo += length(intersect(s, h))
        n6_in += length(intersect(s, hupo_n6))
    end
    cov = 100 * n6_in / (6 * length(hupo_n6))
    if label == "19 (0.03)"; global base_total = tot_pass; end
    delta = base_total > 0 ? tot_pass - base_total : 0
    @printf "%-12s  %12d  %12d  %11.2f%%  %+12d\n" label tot_pass tot_in_hupo cov delta
end
@printf "\nHUPO N=6 unique: %d\n" length(hupo_n6)
