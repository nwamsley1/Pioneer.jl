# Empirical exploration: for each scan with two competing target localizations,
# what is the distribution of matched site-determining ions vs the number of
# possible ones? Two definitions of "possible":
#   - n_avail_theo : all theoretical b/y (charges 1,2) whose cumulative phospho
#                    count differs between the two localizations.
#   - n_avail_lib  : subset of the target's library-retained fragments that are
#                    distinguishing vs the competitor.
#
# Uses a "clean" analysis library (no positional-isomer decoys) for target-only
# competition, joined by (seq, phospho positions, CAM positions, charge) to the
# search-time library so we can reuse an existing search's precursors_long.arrow.
#
# Writes per-pair rows to `<out_prefix>.arrow` and prints stratified summaries.
#
# Usage:
#   julia --project=Pioneer.jl explore_matched_counts.jl \
#       <results_dir> <analysis_lib.poin> <search_lib.poin> <arrow_dir> \
#       [n_sample] [out_prefix]

using Pioneer  # for DetailedFrag / SplineCompactFrag types in the .jls sidecars
using Arrow, DataFrames, Statistics, Random, Serialization, CodecZlib

# ---------- constants + AA masses ----------
const PROTON = 1.0072764668
const H2O    = 18.0105646863
const PHOS   = 79.96633
const CAM    = 57.02146
const AA = Dict{Char,Float64}(
    'A'=>71.03711,'R'=>156.10111,'N'=>114.04293,'D'=>115.02694,'C'=>103.00919,
    'E'=>129.04259,'Q'=>128.05858,'G'=>57.02146,'H'=>137.05891,'I'=>113.08406,
    'L'=>113.08406,'K'=>128.09496,'M'=>131.04049,'F'=>147.06841,'P'=>97.05276,
    'S'=>87.03203,'T'=>101.04768,'W'=>186.07931,'Y'=>163.06333,'V'=>99.06841)

phpos(m)  = Int[parse(Int,x.captures[1]) for x in eachmatch(r"\((\d+),[A-Z],Unimod:21\)", String(m))]
campos(m) = Int[parse(Int,x.captures[1]) for x in eachmatch(r"\((\d+),[A-Z],Unimod:4\)",  String(m))]

@inline function frag_mz(res::Vector{Float64}, ph::Vector{Int}, cm::Vector{Int},
                         T::Symbol, i::Int, z::Int, L::Int)
    if T === :b
        m = 0.0
        @inbounds for k in 1:i; m += res[k]; end
        m += PHOS*count(<=(i), ph) + CAM*count(<=(i), cm)
    else
        m = H2O
        @inbounds for k in (L-i+1):L; m += res[k]; end
        m += PHOS*count(>=(L-i+1), ph) + CAM*count(>=(L-i+1), cm)
    end
    (m + z*PROTON) / z
end

@inline function has_match(mz::Float64, omz::Vector{Float64}, tol_ppm::Float64)
    d  = mz * tol_ppm * 1e-6
    lo = searchsortedfirst(omz, mz - d)
    hi = searchsortedlast(omz, mz + d)
    lo <= hi
end

# Deserialize the compressed .jls sidecars used by BuildSpecLib.
function load_jls(path)
    open(path, "r") do io
        return deserialize(GzipDecompressorStream(io))
    end
end

# ---------- CLI ----------
function main()
    length(ARGS) >= 4 || error("args: <results_dir> <analysis_lib.poin> <search_lib.poin> <arrow_dir> [n_sample] [out_prefix]")
    resdir     = ARGS[1]
    ana_lib    = ARGS[2]
    search_lib = ARGS[3]
    arrow_dir  = ARGS[4]
    nsample    = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 20_000
    out_prefix = length(ARGS) >= 6 ? ARGS[6] : joinpath(resdir, "explore_matched_counts")
    tol_ppm    = 20.0
    charges    = (1, 2)

    println("results_dir   = $resdir")
    println("analysis_lib  = $ana_lib")
    println("search_lib    = $search_lib")
    println("arrow_dir     = $arrow_dir")
    println("n_sample      = $nsample")

    # ---------- load analysis library (clean, target-only competition) ----------
    ap = DataFrame(Arrow.Table(joinpath(ana_lib, "precursors_table.arrow")))
    a_seq = String.(ap.sequence)
    a_mod = String.(coalesce.(ap.structural_mods, ""))
    a_chg = collect(ap.prec_charge)
    a_isd = collect(ap.is_decoy)
    a_ild = :is_loc_decoy in propertynames(ap) ? collect(ap.is_loc_decoy) : falses(nrow(ap))
    a_ph  = [phpos(m)  for m in a_mod]
    a_cm  = [campos(m) for m in a_mod]

    is_ana_target = (.!a_isd) .& (.!a_ild)
    println("analysis lib: n=$(nrow(ap))  targets=$(count(is_ana_target))  loc_decoys=$(count(a_ild))")

    # Per-precursor retained fragments from the analysis library (Prosit path
    # => Vector{DetailedFrag{Float32}}; explicit is_b/is_y/ion_position/charge).
    ana_frags = load_jls(joinpath(ana_lib, "detailed_fragments.jls"))
    ana_p2f   = load_jls(joinpath(ana_lib, "precursor_to_fragment_indices.jls"))
    println("analysis lib: $(length(ana_frags)) retained fragments (avg $(round(length(ana_frags)/nrow(ap), digits=1)) / precursor)")

    # (seq, chg) -> [analysis_lib pid, ...] for target-only competitor lookup.
    ana_grpT = Dict{Tuple{String,UInt8},Vector{Int}}()
    for i in 1:nrow(ap)
        is_ana_target[i] || continue
        push!(get!(ana_grpT, (a_seq[i], a_chg[i]), Int[]), i)
    end

    # Join key: (seq, sorted phospho positions, sorted CAM positions, charge)
    join_key(seq, ph, cm, chg) = (seq, Tuple(sort(ph)), Tuple(sort(cm)), UInt8(chg))
    ana_pid_by_key = Dict{Tuple{String,NTuple{N,Int} where N,NTuple{M,Int} where M,UInt8}, Int}()
    for i in 1:nrow(ap)
        is_ana_target[i] || continue
        ana_pid_by_key[join_key(a_seq[i], a_ph[i], a_cm[i], a_chg[i])] = i
    end

    # ---------- load search library (only for precursor_idx -> key lookup) ----------
    sp = DataFrame(Arrow.Table(joinpath(search_lib, "precursors_table.arrow")))
    s_seq = String.(sp.sequence)
    s_mod = String.(coalesce.(sp.structural_mods, ""))
    s_chg = collect(sp.prec_charge)
    s_ph  = [phpos(m)  for m in s_mod]
    s_cm  = [campos(m) for m in s_mod]

    # ---------- load PSMs ----------
    d = DataFrame(Arrow.Table(joinpath(resdir, "precursors_long.arrow")))
    println("PSMs total: $(nrow(d))")
    d = d[(d.target .== true) .& (d.qval .<= 0.01f0) .& (coalesce.(d.global_qval, 1f0) .<= 0.01f0), :]
    println("PSMs q<=0.01 targets: $(nrow(d))")

    # Restrict to PSMs whose target precursor exists in the clean analysis library
    # and has >=1 target sibling competitor there.
    keep = falses(nrow(d))
    ana_pid_of_psm = fill(0, nrow(d))
    for r in 1:nrow(d)
        sid = Int(d.precursor_idx[r])
        sid == 0 && continue
        k = join_key(s_seq[sid], s_ph[sid], s_cm[sid], s_chg[sid])
        aid = get(ana_pid_by_key, k, 0)
        aid == 0 && continue
        siblings = get(ana_grpT, (a_seq[aid], a_chg[aid]), Int[])
        length(siblings) >= 2 || continue
        keep[r] = true
        ana_pid_of_psm[r] = aid
    end
    d = d[keep, :]
    d.ana_pid = ana_pid_of_psm[keep]
    println("PSMs with >=1 target sibling in analysis lib: $(nrow(d))")

    # Best-scan-per-precursor is already what precursors_long emits; sample.
    Random.seed!(1)
    if nrow(d) > nsample
        d = d[shuffle(1:nrow(d))[1:nsample], :]
    end
    println("scoring $(nrow(d)) PSMs")

    # ---------- accumulate per-pair rows ----------
    out = DataFrame(
        precursor_idx  = Int[],
        ana_pid        = Int[],
        comp_pid       = Int[],
        file           = String[],
        scan_idx       = Int[],
        seq            = String[],
        charge         = UInt8[],
        C              = Int[],
        gsz            = Int[],
        n_peaks_scan   = Int[],
        n_peaks_iso    = Int[],
        n_avail_theo   = Int[],
        n_match_tgt_t  = Int[],
        n_match_cmp_t  = Int[],
        n_match_either_t = Int[],
        n_avail_lib    = Int[],
        n_match_tgt_l  = Int[],
        n_match_cmp_l  = Int[],
        n_match_either_l = Int[],
    )

    nSTY(s) = count(c -> c == 'S' || c == 'T' || c == 'Y', s)
    C_of(seq, nph) = (0 < nph <= nSTY(seq)) ? binomial(nSTY(seq), nph) : 1

    # Group by file so we only open each arrow once.
    for gi in groupby(d, :file_name)
        fname = gi.file_name[1]
        af = joinpath(arrow_dir, string(fname) * ".arrow")
        if !isfile(af)
            # short form in results, prefixed on disk -- pick any arrow matching *fname.arrow
            candidates = filter(f -> endswith(f, string(fname) * ".arrow"), readdir(arrow_dir))
            if !isempty(candidates)
                af = joinpath(arrow_dir, candidates[1])
            end
        end
        isfile(af) || (println("  missing $af  (skip $(nrow(gi)) rows)"); continue)
        t = Arrow.Table(af)
        MZ = t.mz_array; IN = t.intensity_array
        HAS_CENTER = :centerMz in propertynames(t) && :isolationWidthMz in propertynames(t)
        CMz = HAS_CENTER ? t.centerMz : nothing
        IW  = HAS_CENTER ? t.isolationWidthMz : nothing

        for row in eachrow(gi)
            aid = row.ana_pid
            seq = a_seq[aid]; L = length(seq)
            res = [AA[c] for c in seq]
            ph_t = a_ph[aid]; cm_t = a_cm[aid]
            chg  = a_chg[aid]
            sc = Int(row.scan_idx)

            omz = Float64.(MZ[sc]); oint = Float64.(IN[sc])
            p = sortperm(omz); omz = omz[p]; oint = oint[p]
            n_peaks = length(omz)

            n_iso = 0
            if HAS_CENTER
                cmz = Float64(CMz[sc]); iw = Float64(IW[sc])
                if !isnan(cmz) && !isnan(iw) && iw > 0
                    lo = cmz - iw/2; hi = cmz + iw/2
                    n_iso = searchsortedlast(omz, hi) - searchsortedfirst(omz, lo) + 1
                    n_iso = max(n_iso, 0)
                end
            end

            # Precompute library-retained (is_b/is_y, i, z) tuples for target.
            lo = Int(ana_p2f[aid]); hi = Int(ana_p2f[aid+1]) - 1
            lib_set = Set{Tuple{Symbol,Int,Int}}()
            @inbounds for k in lo:hi
                f = ana_frags[k]
                is_b_f = Pioneer.isB(f); is_y_f = Pioneer.isY(f)
                (is_b_f || is_y_f) || continue
                Pioneer.isP(f) && continue
                Pioneer.isIso(f) && continue
                z = Int(Pioneer.getFragCharge(f))
                z in charges || continue
                sym = is_b_f ? :b : :y
                push!(lib_set, (sym, Int(Pioneer.getIonPosition(f)), z))
            end

            C_here = C_of(seq, length(ph_t))
            gsz    = Int(coalesce(row.iso_group_size_at_scan, 1))
            siblings = ana_grpT[(seq, chg)]
            for cid in siblings
                cid == aid && continue
                ph_c = a_ph[cid]
                # Enumerate distinguishing ions.
                n_avail_theo = 0; n_mt_theo = 0; n_mc_theo = 0; n_me_theo = 0
                n_avail_lib  = 0; n_mt_lib  = 0; n_mc_lib  = 0; n_me_lib  = 0
                for T in (:b, :y), i in 1:(L-1)
                    if T === :b
                        cT = count(<=(i), ph_t); cC = count(<=(i), ph_c)
                    else
                        cT = count(>=(L-i+1), ph_t); cC = count(>=(L-i+1), ph_c)
                    end
                    cT == cC && continue
                    for z in charges
                        mz_t = frag_mz(res, ph_t, cm_t, T, i, z, L)
                        mz_c = frag_mz(res, ph_c, cm_t, T, i, z, L)
                        ht = has_match(mz_t, omz, tol_ppm)
                        hc = has_match(mz_c, omz, tol_ppm)
                        n_avail_theo += 1
                        ht && (n_mt_theo += 1)
                        hc && (n_mc_theo += 1)
                        (ht || hc) && (n_me_theo += 1)
                        if (T, i, z) in lib_set
                            n_avail_lib += 1
                            ht && (n_mt_lib += 1)
                            hc && (n_mc_lib += 1)
                            (ht || hc) && (n_me_lib += 1)
                        end
                    end
                end

                push!(out, (
                    Int(row.precursor_idx), aid, cid, String(fname), sc,
                    seq, UInt8(chg), C_here, gsz,
                    n_peaks, n_iso,
                    n_avail_theo, n_mt_theo, n_mc_theo, n_me_theo,
                    n_avail_lib,  n_mt_lib,  n_mc_lib,  n_me_lib,
                ))
            end
        end
    end

    outp = out_prefix * ".arrow"
    Arrow.write(outp, out)
    println("\nwrote $(nrow(out)) pair rows -> $outp")

    # ---------- summaries ----------
    println("\n=== distribution of n_avail_theo ===")
    for a in [2, 3, 4, 5, 6, 8, 10, 12, 16, 24]
        n = count(<=(a), out.n_avail_theo)
        println("  <=$a  : $n  ($(round(100*n/nrow(out),digits=1))%)")
    end

    stratify(col::Symbol, sub, num_tgt::Symbol, num_cmp::Symbol) = begin
        for b in [(2,2),(3,3),(4,4),(5,5),(6,6),(7,8),(9,10),(11,12),(13,16),(17,20),(21,24),(25,32),(33,999)]
            g = sub[(sub[!, col] .>= b[1]) .& (sub[!, col] .<= b[2]), :]
            isempty(g) && continue
            frac_t = mean(g[!, num_tgt] ./ g[!, col])
            frac_c = mean(g[!, num_cmp] ./ g[!, col])
            println("  n_avail=$(b[1])-$(b[2])  n=$(nrow(g))  <match_tgt/N>=$(round(frac_t,digits=3))  <match_cmp/N>=$(round(frac_c,digits=3))")
        end
    end
    println("\n=== match rate vs n_avail_theo (target side = 'signal'; competitor side = 'null') ===")
    stratify(:n_avail_theo, out[out.n_avail_theo .>= 1, :], :n_match_tgt_t, :n_match_cmp_t)
    println("\n=== match rate vs n_avail_lib ===")
    stratify(:n_avail_lib, out[out.n_avail_lib .>= 1, :], :n_match_tgt_l, :n_match_cmp_l)

    println("\n=== match rate stratified by scan peak count ===")
    peak_bins = [(0,500),(501,1000),(1001,1500),(1501,2000),(2001,3000),(3001,10^9)]
    o = out[out.n_avail_theo .>= 2, :]
    for pb in peak_bins
        g = o[(o.n_peaks_scan .>= pb[1]) .& (o.n_peaks_scan .<= pb[2]), :]
        isempty(g) && continue
        frac_c = mean(g.n_match_cmp_t ./ g.n_avail_theo)
        frac_t = mean(g.n_match_tgt_t ./ g.n_avail_theo)
        println("  peaks $(pb[1])-$(pb[2])  n=$(nrow(g))  <match_tgt/N>=$(round(frac_t,digits=3))  <match_cmp/N>=$(round(frac_c,digits=3))")
    end

    println("\n=== match rate stratified by isolation-window peak count ===")
    for pb in [(0,20),(21,50),(51,100),(101,200),(201,400),(401,10^9)]
        g = o[(o.n_peaks_iso .>= pb[1]) .& (o.n_peaks_iso .<= pb[2]), :]
        isempty(g) && continue
        frac_c = mean(g.n_match_cmp_t ./ g.n_avail_theo)
        frac_t = mean(g.n_match_tgt_t ./ g.n_avail_theo)
        println("  iso-peaks $(pb[1])-$(pb[2])  n=$(nrow(g))  <match_tgt/N>=$(round(frac_t,digits=3))  <match_cmp/N>=$(round(frac_c,digits=3))")
    end

    println("\n=== match rate stratified by C (number of possible localizations) ===")
    for cb in [(2,2),(3,3),(4,6),(7,10),(11,20),(21,10^6)]
        g = o[(o.C .>= cb[1]) .& (o.C .<= cb[2]), :]
        isempty(g) && continue
        frac_c = mean(g.n_match_cmp_t ./ g.n_avail_theo)
        frac_t = mean(g.n_match_tgt_t ./ g.n_avail_theo)
        println("  C=$(cb[1])-$(cb[2])  n=$(nrow(g))  <match_tgt/N>=$(round(frac_t,digits=3))  <match_cmp/N>=$(round(frac_c,digits=3))")
    end

    println("\n=== fraction distribution histogram (target-side match_tgt / n_avail_theo) ===")
    frac = out.n_match_tgt_t ./ max.(out.n_avail_theo, 1)
    for lo in 0.0:0.1:0.9
        hi = lo + 0.1
        n = count(x -> lo <= x < (hi == 1.0 ? Inf : hi), frac)
        println("  [$(lo)-$(hi)) : $n  ($(round(100*n/length(frac),digits=1))%)")
    end
end

main()
