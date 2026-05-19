# 3-way comparison on MTAC_3P_Standard:
#   develop_baseline   — develop @ 932b8c77
#   develop_secondpass — develop + SecondPassConsensusSearch
#   hupo2026           — external reference: HUPO 0.6.6

ENV["GKSwstype"] = "100"
using Arrow, DataFrames, Statistics, StatsBase, Plots, Printf

const OUT = joinpath(@__DIR__, "..", "results", "mtac_3p_std_3way_secondpass")
mkpath(OUT); mkpath(joinpath(OUT, "plots"))

const CELLS = [
    ("develop_baseline",   "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/04_develop_baseline/precursors_long.arrow"),
    ("develop_secondpass", "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/05_develop_secondpass/precursors_long.arrow"),
    ("hupo2026",           "/Users/nathanwamsley/Data/HUPO2026/MtacThreeProteomeStandard_hupo2026_v0o6o6/precursors_long.arrow"),
]
const COLORS = Dict(
    "develop_baseline"   => :gray,
    "develop_secondpass" => :dodgerblue,
    "hupo2026"           => :firebrick,
)
const CV_BINS = [0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.0]

_sample(f) = let m = match(r"Sample-([AB])", f); m === nothing ? "?" : String(m.captures[1]); end

function load_cell(path)
    isfile(path) || return nothing
    df = DataFrame(Arrow.Table(path))
    qcol = hasproperty(df, :qval) ? :qval : (hasproperty(df, :MBR_boosted_qval) ? :MBR_boosted_qval : :q_value)
    df = df[(df[!, qcol] .< 0.01) .& df.target .& (.!occursin.(";", df.species)), :]
    df.sample = _sample.(df.file_name)
    df
end

function compute_log2fc(df)
    g = combine(groupby(df, [:precursor_idx, :species, :sample]),
                :peak_area_normalized => (x -> median(skipmissing(x))) => :area)
    w = unstack(g, [:precursor_idx, :species], :sample, :area)
    w = w[(.!ismissing.(w.A)) .& (.!ismissing.(w.B)) .& (w.A .> 0) .& (w.B .> 0), :]
    w.log2fc = log2.(w.A ./ w.B); w
end

function cv6_human(df)
    g = combine(groupby(df[df.species .== "HUMAN", :], [:precursor_idx, :file_name]),
                :peak_area_normalized => (x -> begin v = collect(skipmissing(x)); isempty(v) ? missing : Float64(first(v)); end) => :area)
    g = g[(.!ismissing.(g.area)) .& (g.area .> 0), :]
    cvs = Float64[]
    for sub in groupby(g, :precursor_idx)
        if nrow(sub) == 6
            v = collect(skipmissing(sub.area))
            push!(cvs, std(v) / mean(v))
        end
    end
    cvs
end

results = Dict{String, Any}()
for (n, p) in CELLS
    df = load_cell(p)
    df === nothing && (@warn "missing $n"; continue)
    fc  = compute_log2fc(df)
    cv6 = cv6_human(df)
    mp  = combine(groupby(df, :precursor_idx),
                  :file_name => (x -> length(unique(x))) => :n_files,
                  :species => first => :species)
    h = mp[mp.species .== "HUMAN", :]
    results[n] = (df=df, fc=fc, cv6=cv6, h=h, np=length(unique(df.precursor_idx)))
    println("[loaded] $n: rows=", nrow(df), " unique_prec=", results[n].np, " HUMAN-in-all-6=", length(cv6))
end

println("\n─── ID counts at q<0.01 (target, single-species) ───")
@printf "%-22s  %12s  %12s\n" "cell" "rows" "unique_prec"
for (n, _) in CELLS
    haskey(results, n) || continue
    @printf "%-22s  %12d  %12d\n" n nrow(results[n].df) results[n].np
end

println("\n─── HUMAN multiplicity ───")
@printf "%-22s  %7s  %7s  %7s  %7s  %7s  %7s  %8s\n" "cell" "N=6" "N=5" "N=4" "N=3" "N=2" "N=1" "total"
for (n, _) in CELLS
    haskey(results, n) || continue
    h = results[n].h
    cs = [count(==(k), h.n_files) for k in 1:6]
    @printf "%-22s  %7d  %7d  %7d  %7d  %7d  %7d  %8d\n" n cs[6] cs[5] cs[4] cs[3] cs[2] cs[1] nrow(h)
end

println("\n─── HUMAN-in-all-6 cumulative by CV ───")
@printf "%-22s  %6s  %5s  %5s  %5s  %5s  %5s  %5s  %8s\n" "cell" "CV≤.05" "≤.10" "≤.15" "≤.20" "≤.30" "≤.50" "≤1.0" "all-6"
for (n, _) in CELLS
    haskey(results, n) || continue
    cv6 = results[n].cv6
    cs = [count(<=(t), cv6) for t in CV_BINS]
    @printf "%-22s  %6d  %5d  %5d  %5d  %5d  %5d  %5d  %8d\n" n cs[1] cs[2] cs[3] cs[4] cs[5] cs[6] cs[7] length(cv6)
end

println("\n─── Per-species log2FC summary ───")
@printf "%-22s  %-7s  %8s  %8s  %8s  %8s  %8s\n" "cell" "species" "n_pair" "log2_med" "log2_mad" "|>2|" "|>4|"
for (n, _) in CELLS
    haskey(results, n) || continue
    f = results[n].fc
    for sp in ("HUMAN","YEAST","ECOLI","CONTAM")
        sub = f[f.species .== sp, :]; isempty(sub) && continue
        @printf "%-22s  %-7s  %8d  %8.4f  %8.4f  %8.4f  %8.4f\n" n sp nrow(sub) median(sub.log2fc) mad(sub.log2fc, normalize=false) mean(abs.(sub.log2fc) .> 2) mean(abs.(sub.log2fc) .> 4)
    end
end

# Plots
for sp in ("HUMAN","YEAST","ECOLI")
    plt = plot(size=(900,500), title="MTAC_3P_Standard log2(A/B) — $sp",
               xlabel="log2 fold change", ylabel="density",
               xlims=(-6,6), legend=:topright, legendfontsize=10)
    for (n,_) in CELLS
        haskey(results, n) || continue
        sub = results[n].fc; sub = sub[sub.species .== sp, :]
        isempty(sub) && continue
        stephist!(plt, sub.log2fc, bins=-6:0.1:6, label=n, color=COLORS[n], linewidth=2, normalize=:pdf)
    end
    savefig(plt, joinpath(OUT, "plots", "log2fc_$(sp).png"))
end
plt_cv = plot(size=(1000,500),
              title="HUMAN-in-all-6: cumulative count by CV",
              xlabel="CV threshold", ylabel="n unique HUMAN precursors",
              legend=:bottomright, legendfontsize=10)
for (n,_) in CELLS
    haskey(results, n) || continue
    cv6 = results[n].cv6
    ys = [count(<=(t), cv6) for t in CV_BINS]
    plot!(plt_cv, CV_BINS, ys, marker=:circle, linewidth=2, label=n, color=COLORS[n])
end
savefig(plt_cv, joinpath(OUT, "plots", "cv_in_6_human.png"))

println("\nWrote $OUT")
