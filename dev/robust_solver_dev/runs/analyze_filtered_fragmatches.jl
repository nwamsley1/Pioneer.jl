# Diagnostic: per-file filtered fragment-index match composition + HUPO 2026 overlap.
#
# After running MainSearch (with the new dump hook + passing_precs.arrow), this:
#   - loads passing_precs (q≤0.10 globally, target ∪ decoy)
#   - for each file, loads temp_data/fragment_index_matches/<file>.arrow
#   - filters (scan_idx, precursor_idx) pairs to passing_precs
#   - reports total rows, unique target/decoy precursor counts
#   - computes overlap with HUPO 2026 q<0.01 IDs for the corresponding file
#
# Run:
#   julia --project=/private/tmp/pioneer-develop-latest analyze_filtered_fragmatches.jl

using Arrow, DataFrames, Printf, Statistics

const RUN_DIR = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/06_fragmatch_diag"
const LIB_PRECURSORS = "/Users/nathanwamsley/Data/SPEC_LIBS/altimeter_3P_len7o40_ch2o3_mc1_OlsenAstral_mzsorted.poin/precursors_table.arrow"
const HUPO_PATH = "/Users/nathanwamsley/Data/HUPO2026/MtacThreeProteomeStandard_hupo2026_v0o6o6/precursors_long.arrow"

# Map our verbose file_name -> HUPO short name
function short_name(s::AbstractString)
    m = match(r"Sample-([AB]).*Rep(\d+)", s)
    m === nothing && return s
    return "Sample-$(m.captures[1])_Rep$(m.captures[2])"
end

# ── Library decoy lookup ────────────────────────────────────────────────────
println("[load] library precursors_table")
lib = DataFrame(Arrow.Table(LIB_PRECURSORS))
n_prec = nrow(lib)
is_decoy = lib.is_decoy

# ── Passing set (q≤0.10) ────────────────────────────────────────────────────
passing_path = joinpath(RUN_DIR, "temp_data", "passing_precursors.arrow")
isfile(passing_path) || error("missing $passing_path — did MainSearch finish?")
passing = Set{UInt32}(DataFrame(Arrow.Table(passing_path)).precursor_idx)
n_pass_t = count(pid -> !is_decoy[pid], passing)
n_pass_d = count(pid -> is_decoy[pid], passing)
@printf "Global passing (q≤0.10): %d total = %d targets + %d decoys\n" length(passing) n_pass_t n_pass_d

# ── HUPO 2026 q<0.01 IDs per file (mapped to OUR library precursor_idx via
# (sequence, charge, structural_mods, isotopic_mods, is_decoy) key — the
# HUPO precursor_idx values come from a different library build) ────────────
println("\n[load] HUPO 2026 reference")
hupo = DataFrame(Arrow.Table(HUPO_PATH))
qcol = :MBR_boosted_qval
hupo_ok = hupo[(hupo[!, qcol] .< 0.01) .& hupo.target, :]

_norm(x) = ismissing(x) ? "" : String(x)
key_t(seq, ch, sm, im, dc) = (String(seq), UInt8(ch), _norm(sm), _norm(im), Bool(dc))

println("[build] library key → precursor_idx map")
lib_key2pid = Dict{Tuple{String,UInt8,String,String,Bool}, UInt32}()
sizehint!(lib_key2pid, nrow(lib))
@inbounds for i in 1:nrow(lib)
    k = key_t(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i],
              lib.isotopic_mods[i], lib.is_decoy[i])
    lib_key2pid[k] = UInt32(i)
end

hupo_by_file = Dict{String, Set{UInt32}}()
n_hupo_unmatched = 0
for sub in groupby(hupo_ok, :file_name)
    fn = String(first(sub.file_name))
    s = Set{UInt32}()
    for r in eachrow(sub)
        k = key_t(r.sequence, r.charge, r.structural_mods, r.isotopic_mods, r.is_decoy)
        pid = get(lib_key2pid, k, UInt32(0))
        if pid == 0
            global n_hupo_unmatched += 1
        else
            push!(s, pid)
        end
    end
    hupo_by_file[fn] = s
end
@printf "HUPO unmatched-to-library rows (sanity): %d\n" n_hupo_unmatched
avg_hupo = sum(length, values(hupo_by_file)) / max(length(hupo_by_file), 1)
@printf "HUPO ref files: %d (avg %.0f targets q<0.01 per file)\n" length(hupo_by_file) avg_hupo

# ── Per-file filtered fragment matches ──────────────────────────────────────
matches_dir = joinpath(RUN_DIR, "temp_data", "fragment_index_matches")
files = sort(filter(f -> endswith(f, ".arrow"), readdir(matches_dir)))
isempty(files) && error("no fragment_index_matches files found in $matches_dir")

println("\n─── Per-file filtered fragment matches ───")
@printf "%-20s  %12s  %12s  %12s  %10s  %10s  %10s  %10s  %10s  %10s\n"  "file" "rows_total" "rows_keep" "uniq_keep" "uniq_T" "uniq_D" "rows_T" "rows_D" "hupo_n" "hupo_∩keep"

summary_rows = NamedTuple[]
for f in files
    short_full = replace(f, ".arrow" => "")
    short = short_name(short_full)
    tbl = DataFrame(Arrow.Table(joinpath(matches_dir, f)))
    rows_total = nrow(tbl)
    keep_mask = [pid in passing for pid in tbl.precursor_idx]
    rows_keep = count(keep_mask)
    kept_pids = tbl.precursor_idx[keep_mask]
    uniq_keep = length(unique(kept_pids))
    target_mask_keep = [!is_decoy[pid] for pid in kept_pids]
    uniq_T = length(unique(kept_pids[target_mask_keep]))
    uniq_D = uniq_keep - uniq_T
    rows_T = count(target_mask_keep)
    rows_D = rows_keep - rows_T
    hupo_set = get(hupo_by_file, short, Set{UInt32}())
    hupo_n = length(hupo_set)
    hupo_in = count(pid -> pid in hupo_set, kept_pids)
    hupo_uniq = length(intersect(Set(kept_pids), hupo_set))
    push!(summary_rows, (file=short, rows_total=rows_total, rows_keep=rows_keep,
                         uniq_keep=uniq_keep, uniq_T=uniq_T, uniq_D=uniq_D,
                         rows_T=rows_T, rows_D=rows_D,
                         hupo_n=hupo_n, hupo_uniq=hupo_uniq))
    @printf "%-20s  %12d  %12d  %12d  %10d  %10d  %10d  %10d  %10d  %10d\n" short rows_total rows_keep uniq_keep uniq_T uniq_D rows_T rows_D hupo_n hupo_uniq
end

# ── Coverage summary ────────────────────────────────────────────────────────
println("\n─── HUPO 2026 coverage summary (per-file q<0.01 targets recovered in filtered list) ───")
@printf "%-20s  %10s  %10s  %8s\n" "file" "hupo_n" "in_keep" "%cov"
for r in summary_rows
    pct = r.hupo_n == 0 ? 0.0 : 100 * r.hupo_uniq / r.hupo_n
    @printf "%-20s  %10d  %10d  %7.2f%%\n" r.file r.hupo_n r.hupo_uniq pct
end

println()
