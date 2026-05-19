using Arrow, DataFrames, Printf
_norm(x) = ismissing(x) ? "" : String(x)
key3(seq,ch,sm,dc) = (String(seq), UInt8(ch), _norm(sm), Bool(dc))

lib = DataFrame(Arrow.Table("/Users/nathanwamsley/Data/SPEC_LIBS/altimeter_3P_len7o40_ch2o3_mc1_OlsenAstral_mzsorted.poin/precursors_table.arrow"))
lib_key3 = Dict{Tuple{String,UInt8,String,Bool}, UInt32}()
@inbounds for i in 1:nrow(lib)
    lib_key3[key3(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i], lib.is_decoy[i])] = UInt32(i)
end

function loadpass(p)
    df = DataFrame(Arrow.Table(p))
    qcol = hasproperty(df, :MBR_boosted_qval) ? :MBR_boosted_qval : :qval
    df[(df[!, qcol] .< 0.01) .& df.target, :]
end

runs = [
    ("04_base",    "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/04_develop_baseline/precursors_long.arrow"),
    ("07_nscans",  "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/07_n_scans_feature/precursors_long.arrow"),
    ("16_combined","/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/16_combined/precursors_long.arrow"),
    ("hupo2026",   "/Users/nathanwamsley/Data/HUPO2026/MtacThreeProteomeStandard_hupo2026_v0o6o6/precursors_long.arrow"),
]

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
    end
    s
end

sets = Dict{String, Any}(name => pid_set(loadpass(p)) for (name, p) in runs)
files = sort(collect(keys(sets["hupo2026"])))
hupo_n6 = reduce(intersect, [sets["hupo2026"][fn] for fn in files])

@printf "%-16s  %10s  %10s  %10s  %10s  %10s  %10s\n" "file" "04_base" "07" "16_comb" "Δ16v07" "hupo" "16∩hupo"
for fn in files
    b04 = length(sets["04_base"][fn])
    s07 = length(sets["07_nscans"][fn])
    s16 = length(sets["16_combined"][fn])
    h   = length(sets["hupo2026"][fn])
    inter = length(intersect(sets["16_combined"][fn], sets["hupo2026"][fn]))
    @printf "%-16s  %10d  %10d  %10d  %+10d  %10d  %10d\n" fn b04 s07 s16 (s16-s07) h inter
end
@printf "\nHUPO N=6 unique: %d\n" length(hupo_n6)
println("Coverage by run (avg per-file):")
for (lbl, S) in [("04_base", sets["04_base"]), ("07_nscans", sets["07_nscans"]), ("16_combined", sets["16_combined"])]
    avg_cov = sum(length(intersect(S[fn], hupo_n6)) for fn in files) / (6 * length(hupo_n6))
    @printf "  %-12s  avg %% N=6 cov: %.2f%%\n" lbl 100*avg_cov
end
