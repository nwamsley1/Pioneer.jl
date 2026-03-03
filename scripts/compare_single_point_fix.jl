##
# Compare precursor counts before/after single-point chromatogram fix
# Run after: SearchDIA("/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/olsen_eclipse_fix_single_point_chrom.json")
##
using Arrow, DataFrames, Tables

baseline_dir = "/Users/nathanwamsley/Data/For_Figures/OlsenEclipse/OlsenEclipse_LR_N200_d10_a01"
fix_dir      = "/Users/nathanwamsley/Data/For_Figures/OlsenEclipse/OlsenEclipse_FixSinglePointChrom"

# --- Precursors per file ---
println("="^80)
println("PRECURSOR COUNTS PER FILE (q < 0.01)")
println("="^80)

baseline_prec = DataFrame(Tables.columntable(Arrow.Table(joinpath(baseline_dir, "precursors_long.arrow"))))
fix_prec      = DataFrame(Tables.columntable(Arrow.Table(joinpath(fix_dir, "precursors_long.arrow"))))

# Get file names
baseline_files = sort(unique(baseline_prec.file_name))
fix_files      = sort(unique(fix_prec.file_name))

println("\n", rpad("File", 65), rpad("Baseline", 12), rpad("Fixed", 12), "Delta")
println("-"^100)

total_baseline = 0
total_fix = 0
for fname in union(baseline_files, fix_files)
    global total_baseline, total_fix
    n_base = nrow(filter(r -> r.file_name == fname, baseline_prec))
    n_fix  = nrow(filter(r -> r.file_name == fname, fix_prec))
    delta  = n_fix - n_base
    pct    = n_base > 0 ? round(100.0 * delta / n_base, digits=1) : 0.0
    short  = length(fname) > 60 ? fname[1:57] * "..." : fname
    println(rpad(short, 65), rpad(string(n_base), 12), rpad(string(n_fix), 12),
            delta >= 0 ? "+$delta (+$(pct)%)" : "$delta ($(pct)%)")
    total_baseline += n_base
    total_fix += n_fix
end

println("-"^100)
total_delta = total_fix - total_baseline
total_pct = total_baseline > 0 ? round(100.0 * total_delta / total_baseline, digits=1) : 0.0
println(rpad("TOTAL", 65), rpad(string(total_baseline), 12), rpad(string(total_fix), 12),
        total_delta >= 0 ? "+$total_delta (+$(total_pct)%)" : "$total_delta ($(total_pct)%)")

# --- Unique precursor IDs ---
println("\n\n", "="^80)
println("UNIQUE PRECURSOR IDs (experiment-wide)")
println("="^80)
baseline_unique = length(unique(baseline_prec.precursor_idx))
fix_unique      = length(unique(fix_prec.precursor_idx))
println("Baseline: $baseline_unique")
println("Fixed:    $fix_unique")
println("Delta:    +$(fix_unique - baseline_unique)")

# --- Protein groups ---
println("\n\n", "="^80)
println("PROTEIN GROUP COUNTS")
println("="^80)
baseline_prot = DataFrame(Tables.columntable(Arrow.Table(joinpath(baseline_dir, "protein_groups_long.arrow"))))
fix_prot      = DataFrame(Tables.columntable(Arrow.Table(joinpath(fix_dir, "protein_groups_long.arrow"))))

baseline_pg = length(unique(baseline_prot.protein))
fix_pg      = length(unique(fix_prot.protein))
println("Baseline protein groups: $baseline_pg")
println("Fixed protein groups:    $fix_pg")
println("Delta:                   +$(fix_pg - baseline_pg)")

# --- Check for single-point chromatograms in fix ---
if hasproperty(fix_prec, :points_integrated)
    println("\n\n", "="^80)
    println("POINTS INTEGRATED DISTRIBUTION (Fixed run)")
    println("="^80)
    pts = fix_prec.points_integrated
    for n in sort(unique(pts))
        cnt = count(==(n), pts)
        println("  $n points: $cnt precursors ($(round(100.0*cnt/length(pts), digits=1))%)")
    end
    single_point = count(==(1), pts)
    println("\nSingle-point chromatograms now retained: $single_point")
end
