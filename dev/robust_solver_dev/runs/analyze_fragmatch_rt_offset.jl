# For HUPO 2026 q<0.01 targets that ARE in our filtered fragment matches,
# compute |hupo_apex_rt − nearest_fragmatch_scan_rt| per (file, precursor).

using Arrow, DataFrames, Printf, Statistics

const RUN_DIR  = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/06_fragmatch_diag"
const MS_DIR   = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard"
const LIB_PATH = "/Users/nathanwamsley/Data/SPEC_LIBS/altimeter_3P_len7o40_ch2o3_mc1_OlsenAstral_mzsorted.poin/precursors_table.arrow"
const HUPO_PATH = "/Users/nathanwamsley/Data/HUPO2026/MtacThreeProteomeStandard_hupo2026_v0o6o6/precursors_long.arrow"

const FILE_MAP = Dict(
    "Sample-A_Rep1" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-A_Standard_5min_Rep1.arrow",
    "Sample-A_Rep2" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-A_Standard_5min_Rep2.arrow",
    "Sample-A_Rep3" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-A_Standard_5min_Rep3.arrow",
    "Sample-B_Rep1" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-B_Standard_5min_Rep1.arrow",
    "Sample-B_Rep2" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-B_Standard_5min_Rep2.arrow",
    "Sample-B_Rep3" => "20241211_bkc_25-0856_Goldfarb_Wamsley_Sample-B_Standard_5min_Rep3.arrow",
)

println("[load] library + HUPO + passing")
lib = DataFrame(Arrow.Table(LIB_PATH))
_norm(x) = ismissing(x) ? "" : String(x)
key_t(seq,ch,sm,im,dc) = (String(seq), UInt8(ch), _norm(sm), _norm(im), Bool(dc))

lib_key2pid = Dict{Tuple{String,UInt8,String,String,Bool}, UInt32}()
sizehint!(lib_key2pid, nrow(lib))
@inbounds for i in 1:nrow(lib)
    k = key_t(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i],
              lib.isotopic_mods[i], lib.is_decoy[i])
    lib_key2pid[k] = UInt32(i)
end

passing = Set{UInt32}(DataFrame(Arrow.Table(joinpath(RUN_DIR,"temp_data","passing_precursors.arrow"))).precursor_idx)

hupo = DataFrame(Arrow.Table(HUPO_PATH))
hupo = hupo[(hupo.MBR_boosted_qval .< 0.01) .& hupo.target, :]

println("\n─── RT offset (min) between HUPO apex and nearest filtered fragment-index scan ───")
@printf "%-16s  %8s  %10s  %10s  %8s  %8s  %8s  %8s  %8s  %8s\n" "file" "n" "p50" "p75" "p90" "p95" "p99" "max" "≤0.05" "≤0.10"

for fn in sort(collect(keys(FILE_MAP)))
    spectra_path = joinpath(MS_DIR, FILE_MAP[fn])
    spectra = DataFrame(Arrow.Table(spectra_path))
    rt_arr = spectra.retentionTime  # indexed by scan_idx

    # Filter fragment matches by passing
    fm = DataFrame(Arrow.Table(joinpath(RUN_DIR,"temp_data","fragment_index_matches","$fn.arrow")))
    keep = [pid in passing for pid in fm.precursor_idx]
    fm = fm[keep, :]

    # pid -> sorted vector of scan RTs
    pid2rts = Dict{UInt32, Vector{Float32}}()
    for r in eachrow(fm)
        pid = UInt32(r.precursor_idx)
        push!(get!(pid2rts, pid, Float32[]), Float32(rt_arr[r.scan_idx]))
    end
    for v in values(pid2rts); sort!(v); end

    # HUPO rows for this file
    sub = hupo[hupo.file_name .== fn, :]
    offsets = Float32[]
    sizehint!(offsets, nrow(sub))
    for r in eachrow(sub)
        k = key_t(r.sequence, r.charge, r.structural_mods, r.isotopic_mods, r.is_decoy)
        pid = get(lib_key2pid, k, UInt32(0))
        pid == 0 && continue
        rts = get(pid2rts, pid, nothing)
        rts === nothing && continue
        # Nearest via binary search
        target = Float32(r.rt)
        i = searchsortedfirst(rts, target)
        best = if i == 1
            abs(rts[1] - target)
        elseif i > length(rts)
            abs(rts[end] - target)
        else
            min(abs(rts[i] - target), abs(rts[i-1] - target))
        end
        push!(offsets, best)
    end
    if isempty(offsets)
        @printf "%-16s  no overlap\n" fn
        continue
    end
    sort!(offsets)
    n = length(offsets)
    qq(p) = offsets[clamp(round(Int, p*n), 1, n)]
    @printf "%-16s  %8d  %10.4f  %10.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n" fn n qq(0.50) qq(0.75) qq(0.90) qq(0.95) qq(0.99) maximum(offsets) (count(<=(0.05f0), offsets)/n) (count(<=(0.10f0), offsets)/n)
end
