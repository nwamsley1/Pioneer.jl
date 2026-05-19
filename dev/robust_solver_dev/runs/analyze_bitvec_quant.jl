# Quant precision (CV ECDF) + accuracy (log2 A/B) by species for runs
# 19 (bitvec 0.03) vs 24 (bitvec 0.01).

ENV["GKSwstype"] = "100"
using Arrow, DataFrames, Statistics, StatsBase, Plots, Printf

const OUT = joinpath(@__DIR__, "..", "results", "bitvec_03_vs_01")
mkpath(OUT); mkpath(joinpath(OUT, "plots"))

const RUNS = [
    ("0.03", "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/19_bitvec_03/precursors_long.arrow"),
    ("0.01", "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/24_bitvec_01/precursors_long.arrow"),
]
const COLORS = Dict("0.03" => :gray, "0.01" => :firebrick)

_sample(f) = let m = match(r"Sample-([AB])", f); m === nothing ? "?" : String(m.captures[1]); end

function load_cell(path)
    df = DataFrame(Arrow.Table(path))
    qcol = hasproperty(df, :MBR_boosted_qval) ? :MBR_boosted_qval : :qval
    df = df[(df[!, qcol] .< 0.01) .& df.target .& (.!occursin.(";", df.species)), :]
    df.sample = _sample.(df.file_name)
    df
end

# CV per precursor: 3 replicates per sample (A or B). Returns Dict{species, Vector{Float64}}
# of CVs computed within each sample (so HUMAN = 0 effective ratio, YEAST/ECOLI = nonzero).
function cv_by_species(df, sample_label)
    sub = df[df.sample .== sample_label, :]
    out = Dict{String, Vector{Float64}}()
    for sp_grp in groupby(sub, :species)
        sp = String(first(sp_grp.species))
        # collapse to one area per (precursor, file)
        collapsed = combine(groupby(sp_grp, [:precursor_idx, :file_name]),
                            :peak_area_normalized => (x -> begin
                                v = collect(skipmissing(x)); isempty(v) ? missing : Float64(first(v))
                            end) => :area)
        collapsed = collapsed[(.!ismissing.(collapsed.area)) .& (collapsed.area .> 0), :]
        cvs = Float64[]
        for prec_grp in groupby(collapsed, :precursor_idx)
            if nrow(prec_grp) == 3
                v = collect(skipmissing(prec_grp.area))
                push!(cvs, std(v) / mean(v))
            end
        end
        out[sp] = cvs
    end
    out
end

function compute_log2fc(df)
    g = combine(groupby(df, [:precursor_idx, :species, :sample]),
                :peak_area_normalized => (x -> median(skipmissing(x))) => :area)
    w = unstack(g, [:precursor_idx, :species], :sample, :area)
    w = w[(.!ismissing.(w.A)) .& (.!ismissing.(w.B)) .& (w.A .> 0) .& (w.B .> 0), :]
    w.log2fc = log2.(w.A ./ w.B); w
end

results = Dict{String, Any}()
for (n, p) in RUNS
    df = load_cell(p)
    cv_A = cv_by_species(df, "A")
    cv_B = cv_by_species(df, "B")
    fc   = compute_log2fc(df)
    results[n] = (df=df, cv_A=cv_A, cv_B=cv_B, fc=fc)
    println("[$n] rows=", nrow(df))
end

println("\n=== Per-species CV — Sample A 3-rep + Sample B 3-rep ===")
@printf "%-6s  %-7s  %-8s  %8s  %8s  %8s  %8s  %8s\n" "rate" "species" "sample" "n" "p50" "p75" "p90" "%≤0.20"
for (n, _) in RUNS
    r = results[n]
    for sp in ("HUMAN","YEAST","ECOLI")
        for (lbl, cvs) in (("A", get(r.cv_A, sp, Float64[])), ("B", get(r.cv_B, sp, Float64[])))
            isempty(cvs) && continue
            sort!(cvs)
            qq(p) = cvs[clamp(round(Int, p*length(cvs)), 1, length(cvs))]
            @printf "%-6s  %-7s  %-8s  %8d  %8.4f  %8.4f  %8.4f  %7.2f%%\n" n sp lbl length(cvs) qq(0.5) qq(0.75) qq(0.9) 100*count(<=(0.20), cvs)/length(cvs)
        end
    end
end

println("\n=== Per-species log2(A/B) accuracy ===")
@printf "%-6s  %-7s  %8s  %8s  %8s  %8s  %8s\n" "rate" "species" "n_pair" "log2_med" "log2_mad" "|>2|" "|>4|"
for (n, _) in RUNS
    r = results[n]
    for sp in ("HUMAN","YEAST","ECOLI","CONTAM")
        sub = r.fc[r.fc.species .== sp, :]; isempty(sub) && continue
        @printf "%-6s  %-7s  %8d  %8.4f  %8.4f  %7.4f  %7.4f\n" n sp nrow(sub) median(sub.log2fc) mad(sub.log2fc, normalize=false) mean(abs.(sub.log2fc) .> 2) mean(abs.(sub.log2fc) .> 4)
    end
end

# ── ECDF plots: CV per species, A and B combined into one curve per cell ──
for sp in ("HUMAN","YEAST","ECOLI")
    plt = plot(size=(900,500), title="CV ECDF — $sp (3-rep within sample)",
               xlabel="CV (std/mean)", ylabel="cumulative fraction",
               xlims=(0,0.5), legend=:bottomright, legendfontsize=10)
    for (n, _) in RUNS
        r = results[n]
        cvs = vcat(get(r.cv_A, sp, Float64[]), get(r.cv_B, sp, Float64[]))
        isempty(cvs) && continue
        sort!(cvs)
        ys = (1:length(cvs)) ./ length(cvs)
        plot!(plt, cvs, ys, label="bitvec=$n  n=$(length(cvs))", color=COLORS[n], linewidth=2)
    end
    savefig(plt, joinpath(OUT, "plots", "cv_ecdf_$(sp).png"))
end

# Combined log2FC density per species
for sp in ("HUMAN","YEAST","ECOLI")
    plt = plot(size=(900,500), title="MTAC_3P_Standard log2(A/B) — $sp",
               xlabel="log2 fold change", ylabel="density",
               xlims=(-6,6), legend=:topright, legendfontsize=10)
    for (n, _) in RUNS
        r = results[n]
        sub = r.fc[r.fc.species .== sp, :]
        isempty(sub) && continue
        stephist!(plt, sub.log2fc, bins=-6:0.1:6, label="bitvec=$n", color=COLORS[n], linewidth=2, normalize=:pdf)
    end
    savefig(plt, joinpath(OUT, "plots", "log2fc_$(sp).png"))
end

println("\nWrote $OUT")
