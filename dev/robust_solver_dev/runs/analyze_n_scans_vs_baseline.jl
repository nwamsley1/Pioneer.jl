# Compare 07_n_scans_feature vs 04_develop_baseline vs HUPO 2026 N=6 targets.

using Arrow, DataFrames, Printf

const PATHS = [
    ("baseline_04", "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/04_develop_baseline/precursors_long.arrow"),
    ("n_scans_07",  "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/07_n_scans_feature/precursors_long.arrow"),
    ("hupo2026",    "/Users/nathanwamsley/Data/HUPO2026/MtacThreeProteomeStandard_hupo2026_v0o6o6/precursors_long.arrow"),
]
const LIB_PATH = "/Users/nathanwamsley/Data/SPEC_LIBS/altimeter_3P_len7o40_ch2o3_mc1_OlsenAstral_mzsorted.poin/precursors_table.arrow"

_norm(x) = ismissing(x) ? "" : String(x)
key3(seq,ch,sm,dc) = (String(seq), UInt8(ch), _norm(sm), Bool(dc))

# Library map (sequence, charge, struct_mods, is_decoy) -> our pid
println("[load] library")
lib = DataFrame(Arrow.Table(LIB_PATH))
lib_key3 = Dict{Tuple{String,UInt8,String,Bool}, UInt32}()
@inbounds for i in 1:nrow(lib)
    lib_key3[key3(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i], lib.is_decoy[i])] = UInt32(i)
end

function load_qpass(path, qcol_priority)
    df = DataFrame(Arrow.Table(path))
    qcol = first(filter(c -> hasproperty(df, c), qcol_priority))
    df = df[(df[!, qcol] .< 0.01) .& df.target, :]
    df
end

# Per-file (file_name, lib_pid) sets
println("[build] per-file id sets")
sets = Dict{String, Dict{String, Set{UInt32}}}()
for (name, p) in PATHS
    df = load_qpass(p, [:MBR_boosted_qval, :qval, :q_value])
    fmap = Dict{String, Set{UInt32}}()
    for sub in groupby(df, :file_name)
        fn = String(first(sub.file_name))
        # normalize: HUPO uses Sample-A_Rep1, baseline uses verbose name
        m = match(r"Sample-([AB]).*Rep(\d+)", fn)
        norm_fn = m === nothing ? fn : "Sample-$(m.captures[1])_Rep$(m.captures[2])"
        s = Set{UInt32}()
        for r in eachrow(sub)
            pid = get(lib_key3, key3(r.sequence, r.charge, r.structural_mods, r.is_decoy), UInt32(0))
            pid != 0 && push!(s, pid)
        end
        fmap[norm_fn] = s
    end
    sets[name] = fmap
end

println("\n─── Per-file q<0.01 target IDs ───")
@printf "%-16s  %12s  %12s  %12s  %12s  %12s\n" "file" "baseline_04" "n_scans_07" "Δ vs base" "hupo2026" "07∩hupo"
for fn in sort(unique(vcat([collect(keys(s)) for (_,s) in sets]...)))
    b04 = length(get(sets["baseline_04"], fn, Set()))
    s07 = length(get(sets["n_scans_07"], fn, Set()))
    h   = length(get(sets["hupo2026"], fn, Set()))
    inter = length(intersect(get(sets["n_scans_07"], fn, Set()), get(sets["hupo2026"], fn, Set())))
    @printf "%-16s  %12d  %12d  %+12d  %12d  %12d\n" fn b04 s07 (s07-b04) h inter
end

# HUPO N=6 per-file coverage
println("\n─── HUPO 2026 N=6 (in all 6 HUPO files) coverage ───")
file_keys = sort(collect(keys(sets["hupo2026"])))
hupo_n6 = reduce(intersect, [sets["hupo2026"][fn] for fn in file_keys])
println("HUPO N=6 unique precursors: ", length(hupo_n6))
for cell in ["baseline_04", "n_scans_07"]
    println("\n$cell ∩ HUPO N=6 per file:")
    @printf "  %-16s  %10s  %10s\n" "file" "in_run" "%cov"
    for fn in file_keys
        s = get(sets[cell], fn, Set())
        inter = length(intersect(s, hupo_n6))
        @printf "  %-16s  %10d  %9.2f%%\n" fn inter 100*inter/length(hupo_n6)
    end
end
