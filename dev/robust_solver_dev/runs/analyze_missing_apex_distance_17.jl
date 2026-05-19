# For each file: among HUPO 2026 q<0.01 targets that 16_combined misses at
# q<0.01, count how many of OUR best PSMs are at the same scan as HUPO's
# apex (=0 cycles) or within 1/2/3 cycles. Cycle = same-isolation-window
# spacing on the instrument.

using Arrow, DataFrames, Printf, Statistics, Tables

const RUN17     = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/17_max_weight"
const HUPO_PATH = "/Users/nathanwamsley/Data/HUPO2026/MtacThreeProteomeStandard_hupo2026_v0o6o6/precursors_long.arrow"
const LIB_PATH  = "/Users/nathanwamsley/Data/SPEC_LIBS/altimeter_3P_len7o40_ch2o3_mc1_OlsenAstral_mzsorted.poin/precursors_table.arrow"
const MS_DIR    = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard"
const FILES = ["Sample-A_Rep1","Sample-A_Rep2","Sample-A_Rep3","Sample-B_Rep1","Sample-B_Rep2","Sample-B_Rep3"]
const FILE_TO_VERBOSE = Dict(
    "Sample-A_Rep1" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-A_Standard_5min_Rep1.arrow",
    "Sample-A_Rep2" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-A_Standard_5min_Rep2.arrow",
    "Sample-A_Rep3" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-A_Standard_5min_Rep3.arrow",
    "Sample-B_Rep1" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-B_Standard_5min_Rep1.arrow",
    "Sample-B_Rep2" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-B_Standard_5min_Rep2.arrow",
    "Sample-B_Rep3" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-B_Standard_5min_Rep3.arrow",
)

_norm(x) = ismissing(x) ? "" : String(x)
key_t(seq,ch,sm,im,dc) = (String(seq), UInt8(ch), _norm(sm), _norm(im), Bool(dc))
key3(seq,ch,sm,dc)     = (String(seq), UInt8(ch), _norm(sm), Bool(dc))

println("[load] library + maps")
lib = DataFrame(Arrow.Table(LIB_PATH))
lib_keyt = Dict{Tuple{String,UInt8,String,String,Bool}, UInt32}()
lib_key3 = Dict{Tuple{String,UInt8,String,Bool}, UInt32}()
@inbounds for i in 1:nrow(lib)
    lib_keyt[key_t(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i], lib.isotopic_mods[i], lib.is_decoy[i])] = UInt32(i)
    lib_key3[key3(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i], lib.is_decoy[i])] = UInt32(i)
end

# HUPO q<0.01 target → per-file pid → (scan_idx, rt)
hupo = DataFrame(Arrow.Table(HUPO_PATH))
hupo = hupo[(hupo.MBR_boosted_qval .< 0.01) .& hupo.target, :]
hupo_by_file = Dict{String, Dict{UInt32, Tuple{Int, Float32}}}()
for fn in FILES
    sub = hupo[hupo.file_name .== fn, :]
    m = Dict{UInt32, Tuple{Int, Float32}}()
    for r in eachrow(sub)
        pid = get(lib_keyt, key_t(r.sequence, r.charge, r.structural_mods, r.isotopic_mods, r.is_decoy), UInt32(0))
        pid != 0 && (m[pid] = (Int(r.scan_idx), Float32(r.rt)))
    end
    hupo_by_file[fn] = m
end

# 16's q<0.01 set per file (mapped via 4-key — no isotopic_mods in our long-form)
pl17 = DataFrame(Arrow.Table(joinpath(RUN17,"precursors_long.arrow")))
qcol = hasproperty(pl17, :MBR_boosted_qval) ? :MBR_boosted_qval : :qval
pl17 = pl17[(pl17[!, qcol] .< 0.01) .& pl17.target, :]
pass17_by_file = Dict{String, Set{UInt32}}()
for sub in groupby(pl17, :file_name)
    fn_full = String(first(sub.file_name))
    m = match(r"Sample-([AB]).*Rep(\d+)", fn_full)
    fn = m === nothing ? fn_full : "Sample-$(m.captures[1])_Rep$(m.captures[2])"
    s = Set{UInt32}()
    for r in eachrow(sub)
        pid = get(lib_key3, key3(r.sequence, r.charge, r.structural_mods, r.is_decoy), UInt32(0))
        pid != 0 && push!(s, pid)
    end
    pass17_by_file[fn] = s
end

println("\n─── For HUPO q<0.01 targets that 16_combined missed: distance (in cycles) of our best PSM from HUPO apex ───")
@printf "%-16s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n" "file" "n_miss" "in_main" "=0" "%=0" "≤1" "≤2" "≤3" ">3" "no_psm"

for fn in FILES
    H = hupo_by_file[fn]                        # pid → (scan_idx, rt)
    P17 = pass17_by_file[fn]                    # passing in 16
    miss_pids = setdiff(Set(keys(H)), P17)

    # Load spectra to compute per-window cycle period
    spectra = DataFrame(Tables.columntable(Arrow.Table(joinpath(MS_DIR, FILE_TO_VERBOSE[fn]))))
    rt_arr = spectra.retentionTime
    is_ms2 = spectra.msOrder .== UInt8(2)
    ms2_idx = findall(is_ms2)
    # Group MS2 scans by isolation window key
    win_key(i) = (round(spectra.centerMz[i], digits=4), round(spectra.isolationWidthMz[i], digits=4))
    win_groups = Dict{Tuple{Float32,Float32}, Vector{Int}}()
    for i in ms2_idx
        k = win_key(i)
        push!(get!(win_groups, k, Int[]), i)
    end
    # Per-scan: window scan-position (cycle index within window)
    scan_to_winidx = Dict{Int, Tuple{Tuple{Float32,Float32}, Int}}()
    for (k, scans) in win_groups
        for (pos, s) in enumerate(scans)
            scan_to_winidx[s] = (k, pos)
        end
    end

    # Load our 16-combined main_search_psms: best PSM per precursor (post-reduction)
    f0 = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN17,"temp_data","main_search_psms","$(fn)_fold0.arrow"))))
    f1 = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN17,"temp_data","main_search_psms","$(fn)_fold1.arrow"))))
    main = vcat(f0, f1; cols=:union)
    pid_to_scan = Dict{UInt32, Int}()
    for r in eachrow(main)
        pid_to_scan[UInt32(r.precursor_idx)] = Int(r.scan_idx)
    end

    n_miss = length(miss_pids)
    in_main = 0
    eq0 = 0; le1 = 0; le2 = 0; le3 = 0; gt3 = 0; no_psm = 0
    for pid in miss_pids
        if !haskey(pid_to_scan, pid)
            no_psm += 1
            continue
        end
        in_main += 1
        ours = pid_to_scan[pid]
        hupo_scan = H[pid][1]
        if ours == hupo_scan
            eq0 += 1
            le1 += 1; le2 += 1; le3 += 1
            continue
        end
        # cycle distance: only valid if both scans are in same isolation window
        info_o = get(scan_to_winidx, ours, nothing)
        info_h = get(scan_to_winidx, hupo_scan, nothing)
        if info_o !== nothing && info_h !== nothing && info_o[1] == info_h[1]
            d = abs(info_o[2] - info_h[2])
            if d <= 1
                le1 += 1; le2 += 1; le3 += 1
            elseif d <= 2
                le2 += 1; le3 += 1
            elseif d <= 3
                le3 += 1
            else
                gt3 += 1
            end
        else
            gt3 += 1   # different window or scan not classifiable → counted as >3
        end
    end
    pct0 = n_miss == 0 ? 0.0 : 100*eq0/n_miss
    @printf "%-16s  %8d  %8d  %8d  %7.2f%%  %8d  %8d  %8d  %8d  %8d\n" fn n_miss in_main eq0 pct0 le1 le2 le3 gt3 no_psm
end
