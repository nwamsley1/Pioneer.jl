# For HUPO 2026 q<0.01 targets that ALSO scored at q<0.01 in 04_develop_baseline,
# compute |baseline_apex_rt − hupo_apex_rt| per (file, precursor).

using Arrow, DataFrames, Printf, Statistics

const BASE_PATH = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/04_develop_baseline/precursors_long.arrow"
const HUPO_PATH = "/Users/nathanwamsley/Data/HUPO2026/MtacThreeProteomeStandard_hupo2026_v0o6o6/precursors_long.arrow"

_norm(x) = ismissing(x) ? "" : String(x)
key_t(seq,ch,sm,dc) = (String(seq), UInt8(ch), _norm(sm), Bool(dc))

println("[load] baseline + HUPO")
base = DataFrame(Arrow.Table(BASE_PATH))
base = base[(base.qval .< 0.01) .& base.target, :]
hupo = DataFrame(Arrow.Table(HUPO_PATH))
hupo = hupo[(hupo.MBR_boosted_qval .< 0.01) .& hupo.target, :]

# (file, key) -> rt for baseline and hupo
function build_map(df, qcol_unused=nothing)
    m = Dict{Tuple{String,Tuple{String,UInt8,String,Bool}}, Float32}()
    for r in eachrow(df)
        k = key_t(r.sequence, r.charge, r.structural_mods, r.is_decoy)
        m[(String(r.file_name), k)] = Float32(r.rt)
    end
    m
end

base_map = build_map(base)
hupo_map = build_map(hupo)

# Per-file offset distribution over the intersection
files = sort(unique(hupo.file_name))
println("\n─── |baseline_apex_rt − HUPO_apex_rt| over the q<0.01 intersection (min) ───")
@printf "%-16s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n" "file" "n_∩" "p50" "p75" "p90" "p95" "p99" "max" "≤0.05" "≤0.10"
for fn in files
    keys_in_hupo = Set([k[2] for k in keys(hupo_map) if k[1] == fn])
    keys_in_base = Set([k[2] for k in keys(base_map) if k[1] == fn])
    inter = intersect(keys_in_hupo, keys_in_base)
    offsets = Float32[]
    for k in inter
        push!(offsets, abs(base_map[(fn, k)] - hupo_map[(fn, k)]))
    end
    isempty(offsets) && (println(fn, "  no overlap"); continue)
    sort!(offsets); n = length(offsets)
    qq(p) = offsets[clamp(round(Int, p*n), 1, n)]
    @printf "%-16s  %8d  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n" fn n qq(0.50) qq(0.75) qq(0.90) qq(0.95) qq(0.99) maximum(offsets) (count(<=(0.05f0), offsets)/n) (count(<=(0.10f0), offsets)/n)
end
