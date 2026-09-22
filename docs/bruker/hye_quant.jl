# HYE benchmark: integration strategies compared on replicate CV (precision) and log2(A/B) accuracy per species.
# Inputs per (method tag): dump_<tag>_cond{A,B}/<file>_chrom_weights.arrow, out_<tag>_cond{A,B}/precursors_long.arrow
# (Pioneer's own quant) and temp_data/passing_psms/<file>.arrow.
# Usage: julia --project=<Pioneer> hye_quant.jl <tag: h50|h250> [lambda]
include(joinpath(@__DIR__, "chrom2d.jl"))
using Statistics, Printf, CSV

const H = expanduser("~/BrukerTims/pride_hye")
const LIBP = expanduser("~/BrukerTims/lib/HYE_canon_std_tims.poin/precursors_table.arrow")
const SUFFIX = Ref("")
const EXPECTED = Dict("HUMAN" => 0.0, "YEAST" => log2(30 / 15), "ECOLI" => log2(5 / 20))   # A/B

"Species per precursor index (only unambiguous single-species precursors)."
function species_map()
    t = Arrow.Table(LIBP)
    p = t.proteome_identifiers
    sp = Dict{UInt32, String}()
    for i in eachindex(p)
        s = String(p[i])
        s in ("HUMAN", "YEAST", "ECOLI") && (sp[UInt32(i)] = s)
    end
    sp
end

"Apex seeded at the best PSM's (cycle, IM slice) and hill-climbed over the 8-neighbourhood."
function seeded_apex(Zb, r0, c0)
    nr, nc = size(Zb); r = clamp(r0, 1, nr); c = clamp(c0, 1, nc)
    while true
        best = (r, c); bv = Zb[r, c]
        for dr in -1:1, dc in -1:1
            rr = r + dr; cc = c + dc
            (1 <= rr <= nr && 1 <= cc <= nc) || continue
            Zb[rr, cc] > bv && (bv = Zb[rr, cc]; best = (rr, cc))
        end
        best == (r, c) && return CartesianIndex(r, c)
        r, c = best
    end
end

"Region grown from the apex by descent only (a neighbour is included while it is <= 1.15x its parent and > 0)."
function grow_descent(Z, apex)
    nr, nc = size(Z); inreg = falses(nr, nc); stack = [apex]; inreg[apex] = true
    while !isempty(stack)
        p = pop!(stack)
        for d in (CartesianIndex(1, 0), CartesianIndex(-1, 0), CartesianIndex(0, 1), CartesianIndex(0, -1))
            q = p + d
            (1 <= q[1] <= nr && 1 <= q[2] <= nc) || continue
            inreg[q] && continue
            if Z[q] > 0 && Z[q] <= Z[p] * 1.15
                inreg[q] = true; push!(stack, q)
            end
        end
    end
    inreg
end

"Plane fitted on the cells just outside the region (its rim); falls back to the window border when the rim is small."
function rim_plane(Z, inreg)
    nr, nc = size(Z)
    rows = Float64[]; cols = Float64[]; vals = Float64[]
    for r in 1:nr, c in 1:nc
        inreg[r, c] && continue
        touches = any(((rr, cc),) -> 1 <= rr <= nr && 1 <= cc <= nc && inreg[rr, cc], ((r-1, c), (r+1, c), (r, c-1), (r, c+1)))
        touches || continue
        push!(rows, r); push!(cols, c); push!(vals, Z[r, c])
    end
    length(vals) < 4 && return border_plane(Z)
    β = hcat(ones(length(rows)), rows, cols) \ vals
    [β[1] + β[2] * r + β[3] * c for r in 1:nr, c in 1:nc]
end

"""
All strategies for one dump file. Returns a Dict precursor => NamedTuple of areas.
  raw     : sum of the raw weights in the window
  imsum1d : IM-summed trace, 1D WH (via the 2D solver with λ_im = 0), linear endpoint baseline, trapezoid
  v2      : 2D WH + window-border plane + 3%-threshold region + 2D trapezoid (first prototype)
  v3      : 2D WH + PSM-seeded apex + descent region + rim plane + 2D trapezoid
  apex    : v3's smoothed, baseline-subtracted apex height
"""
function integrate_file(dumpf, psmsf, λ)
    d = DataFrame(Arrow.Table(dumpf))
    p = DataFrame(Arrow.Table(psmsf)); p = p[(p.target) .& (p.qval .<= 0.01), :]
    best = Dict{UInt32, Int}(r.precursor_idx => Int(r.scan_idx) for r in eachrow(p))
    scan_cycle = Dict{Int, Int}(); scan_im = Dict{Int, Int}()
    for r in eachrow(d); scan_cycle[Int(r.scan_idx)] = Int(r.cycle_idx); scan_im[Int(r.scan_idx)] = Int(r.im_scan); end
    cyc_rt = Dict{Int, Float64}(); for r in eachrow(d); cyc_rt[Int(r.cycle_idx)] = Float64(r.rt); end
    out = Dict{UInt32, NamedTuple}()
    for sub in groupby(d, :precursor_idx)
        pid = sub.precursor_idx[1]; haskey(best, pid) || continue
        cyc = Int.(sub.cycle_idx); ims = Int.(sub.im_scan) .÷ STRIDE
        c0, c1 = extrema(cyc); i0, i1 = extrema(ims); nr = c1 - c0 + 1; nc = i1 - i0 + 1
        (nr >= 3 && nc >= 3) || continue
        Y = zeros(nr, nc); W = zeros(nr, nc)
        for k in eachindex(cyc); Y[cyc[k]-c0+1, ims[k]-i0+1] += Float64(sub.weight[k]); W[cyc[k]-c0+1, ims[k]-i0+1] = 1.0; end
        raw = sum(Y); raw > 0 || continue
        scale = maximum(Y); Yn = Y ./ scale
        rts = sort(unique(cyc)); drt = length(rts) > 1 ? median(diff([cyc_rt[c] for c in rts])) : 1.0
        # v2
        Z = wh2d(Yn, W, λ, λ); Zb2 = max.(Z .- border_plane(Z), 0.0); ap2 = argmax(Zb2); reg2 = grow_region(Zb2, ap2, 0.03)
        v2 = trapezoid2d(Zb2, reg2, drt, 1.0) * scale
        # v3
        bs = best[pid]; r0 = get(scan_cycle, bs, c0) - c0 + 1; cc0 = get(scan_im, bs, i0 * STRIDE) ÷ STRIDE - i0 + 1
        ap3 = seeded_apex(Z, r0, cc0); reg3 = grow_descent(Z, ap3)
        Zb3 = max.(Z .- rim_plane(Z, reg3), 0.0); reg3 = grow_descent(Zb3, seeded_apex(Zb3, ap3[1], ap3[2]))
        v3 = trapezoid2d(Zb3, reg3, drt, 1.0) * scale
        apex = Zb3[argmax(Zb3 .* reg3)] * scale
        # imsum1d: sum over IM, 1D smoothing (λ_im = 0 on a 1-column grid), endpoint baseline, trapezoid
        tr = vec(sum(Y; dims = 2)); trn = tr ./ maximum(tr)
        z1 = vec(wh2d(reshape(trn, nr, 1), ones(nr, 1), λ, 0.0))
        a1 = clamp(r0, 1, nr); while a1 < nr && z1[a1+1] >= z1[a1]; a1 += 1; end; while a1 > 1 && z1[a1-1] >= z1[a1]; a1 -= 1; end
        lo = a1; while lo > 1 && z1[lo-1] < z1[lo] && z1[lo-1] > 0; lo -= 1; end
        hi = a1; while hi < nr && z1[hi+1] < z1[hi] && z1[hi+1] > 0; hi += 1; end
        seg = z1[lo:hi] .- range(z1[lo], z1[hi], length = hi - lo + 1); seg = max.(seg, 0.0)
        imsum1d = (sum(seg) - (seg[1] + seg[end]) / 2) * drt * maximum(tr)
        # round 2: decompose the steps
        v4 = sum(Z) * scale                                               # smoothed, whole window, no baseline, no region
        core = [abs(r - ap3[1]) <= 2 && abs(c - ap3[2]) <= 2 for r in 1:nr, c in 1:nc]
        v5 = sum(Y[core])                                                 # raw sum in a 5x5 core around the seeded apex
        v6 = trapezoid2d(Z, grow_descent(Z, ap3), drt, 1.0) * scale       # smoothed + descent region, no baseline
        v7 = sum(max.(Z .- border_plane(Z), 0.0)) * scale                 # smoothed + border plane, whole window
        v8 = sum(Y .* core) - 25 * quantile(vec(Y[.!core]), 0.5)          # raw core minus the window's median cell (flat baseline)
        # round 3: adaptive footprint = cells where the smoothed surface exceeds a fraction of the seeded apex (connected),
        # raw weights summed inside, flat baseline = median of the raw cells outside the footprint
        function footprint_sum(frac)
            thr = frac * Z[ap3]; fp = falses(nr, nc); st = [ap3]; fp[ap3] = true
            while !isempty(st)
                q0 = pop!(st)
                for dd in (CartesianIndex(1, 0), CartesianIndex(-1, 0), CartesianIndex(0, 1), CartesianIndex(0, -1))
                    q = q0 + dd
                    (1 <= q[1] <= nr && 1 <= q[2] <= nc && !fp[q] && Z[q] >= thr) || continue
                    fp[q] = true; push!(st, q)
                end
            end
            outside = Y[.!fp]
            base = isempty(outside) ? 0.0 : median(outside)
            max(sum(Y[fp]) - count(fp) * base, 0.0), count(fp)
        end
        v9, n9 = footprint_sum(0.5); v10, n10 = footprint_sum(0.25); v11, _ = footprint_sum(0.1)
        # hybrid: v3's smoothed surface, the 5x5 rule's bounds and baseline (median SMOOTHED cell outside the block)
        smoothed = Z .* scale
        base_h = count(.!core) == 0 ? 0.0 : median(smoothed[.!core])
        v12 = max(sum(smoothed[core]) - count(core) * base_h, 0.0)
        # same, with a 7x7 block (does the leak outside 5x5 matter once the surface is smooth?)
        core7 = [abs(r - ap3[1]) <= 3 && abs(c - ap3[2]) <= 3 for r in 1:nr, c in 1:nc]
        base7 = count(.!core7) == 0 ? 0.0 : median(smoothed[.!core7])
        v13 = max(sum(smoothed[core7]) - count(core7) * base7, 0.0)
        out[pid] = (raw = raw, imsum1d = imsum1d, v2 = v2, v3 = v3, apex = apex, v4 = v4, v5 = v5, v6 = v6, v7 = v7, v8 = max(v8, 0.0), v9 = v9, v10 = v10, v11 = v11, v12 = v12, v13 = v13)
    end
    out
end

function main(tag, λ)
    sp = species_map()
    # joint run: all 12 files in one search (out_<tag>_joint, dump_<tag>_joint); condition from the file name
    run = "$(tag)_joint"
    dumpdir = joinpath(H, "dump_joint")      # one session, one dump dir; file names keep the gradients apart
    # precursors_long.arrow came out 0 bytes on the 12-file run (the .tsv of the same table is complete) --
    # read whichever is usable.
    pa = joinpath(H, "out_$run", "precursors_long.arrow"); pt = joinpath(H, "out_$run", "precursors_long.tsv")
    pl = (isfile(pa) && filesize(pa) > 0) ? DataFrame(Arrow.Table(pa)) : CSV.read(pt, DataFrame; select = [:file_name, :precursor_idx, :peak_area])
    files = Dict{String, Vector{String}}(); quant = Dict{String, Dict{UInt32, NamedTuple}}(); pioneer = Dict{String, Dict{UInt32, Float64}}()
    for C in ("A", "B"), R in 1:6
        f = "$(tag)_$(C)_REP$R"
        isfile(joinpath(dumpdir, f * "_chrom_weights.arrow")) || continue
        push!(get!(files, C, String[]), f)
        quant[f] = integrate_file(joinpath(dumpdir, f * "_chrom_weights.arrow"), joinpath(H, "out_$run", "temp_data", "passing_psms", f * ".arrow"), λ)
        sub = pl[pl.file_name .== f, :]
        pioneer[f] = Dict{UInt32, Float64}(r.precursor_idx => Float64(r.peak_area) for r in eachrow(sub) if !ismissing(r.peak_area) && r.peak_area > 0)
        @printf("%s: %d integrated, %d with Pioneer peak_area\n", f, length(quant[f]), length(pioneer[f]))
    end
    common = intersect([Set(keys(quant[f])) for C in ("A", "B") for f in files[C]]...)
    common_p = intersect([Set(keys(pioneer[f])) for C in ("A", "B") for f in files[C]]...)
    @printf("\nA: %d replicates, B: %d; precursors quantified in all: %d (dumps), %d (Pioneer peak_area)\n",
            length(files["A"]), length(files["B"]), length(common), length(common_p))
    strategies = [(:v8, "raw 5x5 core - flat baseline"), (:v12, "smoothed 5x5 core - flat base"),
                  (:v13, "smoothed 7x7 core - flat base"), (:v3, "2D v3 (seeded, descent, rim)"), (:pioneer, "Pioneer peak_area")]
    println("\n(log2 A/B centred on the human median of each strategy)\n", rpad("strategy", 36), " | ",
            join([rpad("$s: n  CV_A  CV_B  med log2A/B  MAD  |err|<0.5", 48) for s in ("HUMAN", "YEAST", "ECOLI")], " | "))
    for (key, name) in strategies
        row = rpad(name, 36)
        hcenter = 0.0
        for s in ("HUMAN", "YEAST", "ECOLI")
            pids = [pid for pid in (key == :pioneer ? common_p : common) if get(sp, pid, "") == s]
            get_v(f, pid) = key == :pioneer ? pioneer[f][pid] : getfield(quant[f][pid], key)
            cvs = Dict("A" => Float64[], "B" => Float64[]); ratios = Float64[]
            for pid in pids
                va = [get_v(f, pid) for f in files["A"]]; vb = [get_v(f, pid) for f in files["B"]]
                (all(>(0), va) && all(>(0), vb)) || continue
                length(va) > 1 && push!(cvs["A"], std(va) / mean(va)); length(vb) > 1 && push!(cvs["B"], std(vb) / mean(vb))
                push!(ratios, log2(mean(va) / mean(vb)))
            end
            isempty(ratios) && (row *= " | " * rpad("$s: none", 48); continue)
            s == "HUMAN" && (hcenter = median(ratios))
            ratios = ratios .- hcenter
            err = ratios .- EXPECTED[s]
            mcv(v) = isempty(v) ? NaN : 100median(v)
            row *= " | " * rpad(@sprintf("%s: %5d %5.1f%% %5.1f%% %6.2f (exp %5.2f) %5.2f %5.1f%%", s[1:1], length(ratios), mcv(cvs["A"]), mcv(cvs["B"]), median(ratios), EXPECTED[s], median(abs.(err .- median(err))) * 1.4826, 100mean(abs.(err) .< 0.5)), 48)
        end
        println(row)
    end
end
main(ARGS[1], length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1e-5)
