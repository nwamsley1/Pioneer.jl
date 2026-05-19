# For each file, profile feature distributions of MainSearch best-per-precursor
# PSMs for HUPO 2026 q<0.01 targets that were MISSED by 07_n_scans_feature at q<0.01.

using Arrow, DataFrames, Printf, Statistics, Tables

const RUN07     = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/07_n_scans_feature"
const HUPO_PATH = "/Users/nathanwamsley/Data/HUPO2026/MtacThreeProteomeStandard_hupo2026_v0o6o6/precursors_long.arrow"
const LIB_PATH  = "/Users/nathanwamsley/Data/SPEC_LIBS/altimeter_3P_len7o40_ch2o3_mc1_OlsenAstral_mzsorted.poin/precursors_table.arrow"
const FILES = ["Sample-A_Rep1","Sample-A_Rep2","Sample-A_Rep3","Sample-B_Rep1","Sample-B_Rep2","Sample-B_Rep3"]

const FILE_TO_VERBOSE = Dict(
    "Sample-A_Rep1" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-A_Standard_5min_Rep1",
    "Sample-A_Rep2" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-A_Standard_5min_Rep2",
    "Sample-A_Rep3" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-A_Standard_5min_Rep3",
    "Sample-B_Rep1" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-B_Standard_5min_Rep1",
    "Sample-B_Rep2" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-B_Standard_5min_Rep2",
    "Sample-B_Rep3" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-B_Standard_5min_Rep3",
)

const FEATURES = [:lgbm_prob, :n_scans, :weight, :gof, :poisson, :max_matched_residual,
                  :max_unmatched_residual, :irt_error, :total_ions, :longest_y, :y_count,
                  :log2_intensity_explained, :fitted_hellinger, :fitted_manhattan_distance,
                  :spectrum_peak_count]

_norm(x) = ismissing(x) ? "" : String(x)
key3(seq,ch,sm,dc) = (String(seq), UInt8(ch), _norm(sm), Bool(dc))
key_t(seq,ch,sm,im,dc) = (String(seq), UInt8(ch), _norm(sm), _norm(im), Bool(dc))
quant(x, p) = (xs = sort(x); xs[clamp(round(Int, p*length(xs)), 1, length(xs))])

println("[load] library + maps")
lib = DataFrame(Arrow.Table(LIB_PATH))
lib_key3 = Dict{Tuple{String,UInt8,String,Bool}, UInt32}()
lib_keyt = Dict{Tuple{String,UInt8,String,String,Bool}, UInt32}()
@inbounds for i in 1:nrow(lib)
    lib_key3[key3(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i], lib.is_decoy[i])] = UInt32(i)
    lib_keyt[key_t(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i], lib.isotopic_mods[i], lib.is_decoy[i])] = UInt32(i)
end

# HUPO q<0.01 targets per file (mapped to OUR pids)
hupo = DataFrame(Arrow.Table(HUPO_PATH))
hupo = hupo[(hupo.MBR_boosted_qval .< 0.01) .& hupo.target, :]
hupo_by_file = Dict{String, Set{UInt32}}()
for fn in FILES
    sub = hupo[hupo.file_name .== fn, :]
    s = Set{UInt32}()
    for r in eachrow(sub)
        pid = get(lib_keyt, key_t(r.sequence, r.charge, r.structural_mods, r.isotopic_mods, r.is_decoy), UInt32(0))
        pid != 0 && push!(s, pid)
    end
    hupo_by_file[fn] = s
end

# 07's q<0.01 targets per file
pl07 = DataFrame(Arrow.Table(joinpath(RUN07,"precursors_long.arrow")))
qcol_07 = hasproperty(pl07, :MBR_boosted_qval) ? :MBR_boosted_qval : :qval
pl07_pass = pl07[(pl07[!, qcol_07] .< 0.01) .& pl07.target, :]
pass07_by_file = Dict{String, Set{UInt32}}()
for sub in groupby(pl07_pass, :file_name)
    fn_full = String(first(sub.file_name))
    m = match(r"Sample-([AB]).*Rep(\d+)", fn_full)
    fn = m === nothing ? fn_full : "Sample-$(m.captures[1])_Rep$(m.captures[2])"
    s = Set{UInt32}()
    for r in eachrow(sub)
        pid = get(lib_key3, key3(r.sequence, r.charge, r.structural_mods, r.is_decoy), UInt32(0))
        pid != 0 && push!(s, pid)
    end
    pass07_by_file[fn] = s
end

# For each file, load MainSearch best-per-precursor (post-reduction; has lgbm_prob, n_scans)
# and look at distributions for the missing set vs the surviving set.
println("\n─── Feature distributions: missing (HUPO\\07) vs survived (in 07) ───")

function summarize(label, fn, df_subset)
    n = nrow(df_subset)
    n == 0 && return
    @printf "%-16s  %-9s  %6d" fn label n
    for f in FEATURES
        if hasproperty(df_subset, f)
            v = collect(skipmissing(df_subset[!, f]))
            isempty(v) && (print("   "*lpad("nan", 8)); continue)
            v = Float32.(v)
            @printf "  %8.3f" quant(v, 0.5)
        else
            print("   "*lpad("-", 8))
        end
    end
    println()
end

# Header
@printf "%-16s  %-9s  %6s" "file" "set" "n"
for f in FEATURES; @printf "  %8s" first(string(f), 8); end
println()

for fn in FILES
    f0 = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN07,"temp_data","main_search_psms","$(fn)_fold0.arrow"))))
    f1 = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN07,"temp_data","main_search_psms","$(fn)_fold1.arrow"))))
    main = vcat(f0, f1; cols=:union)
    H = hupo_by_file[fn]
    P = pass07_by_file[fn]
    miss = setdiff(H, P)
    surv = intersect(H, P)
    miss_df = main[in.(main.precursor_idx, Ref(miss)), :]
    surv_df = main[in.(main.precursor_idx, Ref(surv)), :]
    summarize("missing", fn, miss_df)
    summarize("survived", fn, surv_df)
end
