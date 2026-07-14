# Optional precursor-level filters applied to the global-passing set BEFORE the
# Step-11 experiment-wide recompute. Gated by env vars (prototype / Phase-0 study):
#   IRT_FILTER_K  : cross-file observed-iRT dispersion filter at k*sigma (MAD from
#                   confident multi-file anchors). Grandfathers confident IDs
#                   (pre-global exp-wide qval <= 1%). Symmetric on targets+decoys.
#   FILEQ_PREC    : file-specific precursor q-value threshold (from the full per-file
#                   target+decoy distribution). Grandfathered, symmetric.
# When neither is set, returns passing_refs unchanged.
using Statistics

function apply_precursor_filters!(passing_refs::Vector{<:FileReference})
    ks  = get(ENV, "IRT_FILTER_K", "")
    fqs = get(ENV, "FILEQ_PREC", "")
    (isempty(ks) && isempty(fqs)) && return passing_refs
    gf_q = 1.0f-2

    paths = String[file_path(r) for r in passing_refs]
    # --- gather global-passing set in memory ---
    rc = Int[]; for p in paths; push!(rc, length(Arrow.Table(p).precursor_idx)); end
    off = cumsum([0; rc]); tot = off[end]
    pidx = Vector{UInt32}(undef, tot); tgt = Vector{Bool}(undef, tot)
    pp = Vector{Float32}(undef, tot); irt = Vector{Float32}(undef, tot)
    preq = Vector{Float32}(undef, tot); fidx = Vector{Int}(undef, tot)
    for (fi, p) in enumerate(paths)
        n = rc[fi]; n == 0 && continue; t = Arrow.Table(p); o = off[fi]
        @inbounds for i in 1:n
            pidx[o+i]=UInt32(t.precursor_idx[i]); tgt[o+i]=Bool(t.target[i])
            pp[o+i]=Float32(t.prec_prob[i]); irt[o+i]=Float32(t.irt_obs[i])
            preq[o+i]=Float32(t.qval[i]); fidx[o+i]=fi   # :qval here = pre-global exp-wide
        end
    end
    keep = trues(tot)
    grandf = preq .<= gf_q

    # --- iRT dispersion filter (bias-corrected) ---
    # (1) per-precursor consensus = median obs-iRT of its grandfathered (pre-global
    #     ew q <= 1%) instances; (2) per-file bias curve b_f(iRT) = adaptive ~1000-PSM
    #     quantile-bin median of (obs-iRT - consensus) over confident targets, applied
    #     to all instances -> corrected iRT (removes per-run offset + gradient warp);
    #     (3) reference = median of the top-3 scoring corrected-iRT instances;
    #     sigma = 1.4826*MAD of grandfathered-target corrected deltas. Grandfathered
    #     instances are never removed; non-grandfathered removed if |delta| > k*sigma.
    if !isempty(ks)
        k = parse(Float32, ks); N = 3; BIN = 1000
        cbuf = Dict{UInt32, Vector{Float32}}()
        for i in 1:tot; grandf[i] && push!(get!(cbuf, pidx[i], Float32[]), irt[i]); end
        cons = Dict{UInt32, Float32}()
        for (p, v) in cbuf; length(v) >= 2 && (cons[p] = median(v)); end
        cirt = copy(irt)
        fgi = Dict{Int, Vector{Int}}()
        for i in 1:tot; push!(get!(fgi, fidx[i], Int[]), i); end
        for (_, idxs) in fgi
            ci = [i for i in idxs if tgt[i] && grandf[i] && haskey(cons, pidx[i])]
            length(ci) < 50 && continue
            xs = Float32[irt[i] for i in ci]; ys = Float32[irt[i] - cons[pidx[i]] for i in ci]
            o = sortperm(xs); xs = xs[o]; ys = ys[o]; nb = max(1, div(length(xs), BIN))
            edges = Float32[xs[max(1, round(Int, q*length(xs)))] for q in range(0, 1, nb+1)]
            cx = Float32[]; cy = Float32[]
            for b in 1:nb
                m = (xs .>= edges[b]) .& (xs .<= edges[b+1])
                any(m) && (push!(cx, median(xs[m])); push!(cy, median(ys[m])))
            end
            isempty(cx) && continue
            for i in idxs
                x = irt[i]
                c = x <= cx[1] ? cy[1] : x >= cx[end] ? cy[end] : begin
                    j = searchsortedfirst(cx, x); t = (x - cx[j-1]) / (cx[j] - cx[j-1]); cy[j-1] + t*(cy[j] - cy[j-1])
                end
                cirt[i] = irt[i] - c
            end
        end
        pgi = Dict{UInt32, Vector{Int}}()
        for i in 1:tot; push!(get!(pgi, pidx[i], Int[]), i); end
        ref = Dict{UInt32, Float32}()
        for (p, idxs) in pgi
            ord = sortperm(view(pp, idxs); rev=true)
            topk = idxs[ord[1:min(N, length(ord))]]
            ref[p] = median(view(cirt, topk))
        end
        delta = Float32[cirt[i] - ref[pidx[i]] for i in 1:tot]
        tmask = tgt .& grandf
        med = median(@view delta[tmask]); sig = 1.4826f0*median(abs.(delta[tmask] .- med))
        irt_keep = .!((abs.(delta .- med) .> k*sig) .& .!grandf)
        keep .&= irt_keep
        @user_info "iRT filter (bias-corrected, median top-$N) k=$k sigma=$(round(sig,sigdigits=4)): removed $(sum(.!irt_keep)) / $tot instances"
    end

    # --- file-specific precursor q-value filter (from FULL per-file distribution) ---
    # Per-file q is computed over the full main_search_psms best-per-precursor set
    # (targets+decoys), whose prec_prob lives in a `.prec_prob.sidecar.arrow` sidecar.
    if !isempty(fqs)
        thr = parse(Float32, fqs)
        for (fi, p) in enumerate(paths)
            rc[fi] == 0 && continue
            msfile = replace(p, "passing_psms" => "main_search_psms")
            (isfile(msfile) && isfile(msfile * ".prec_prob.sidecar.arrow")) || continue
            mt = Arrow.Table(msfile)
            fprob = Float32.(Arrow.Table(msfile * ".prec_prob.sidecar.arrow").prec_prob)
            ftg = Bool.(mt.target)
            ord = sortperm(fprob; rev=true); q = Vector{Float32}(undef, length(fprob))
            get_qvalues!(fprob[ord], ftg[ord], q); qq = similar(q); qq[ord] = q
            lut = Dict{UInt32, Float32}()
            for j in eachindex(mt.precursor_idx); lut[UInt32(mt.precursor_idx[j])] = qq[j]; end
            o = off[fi]
            @inbounds for i in 1:rc[fi]
                fqv = get(lut, pidx[o+i], 0.0f0)
                if fqv > thr && !grandf[o+i]; keep[o+i] = false; end
            end
        end
        @user_info "file-specific precursor q filter thr=$thr applied"
    end

    # --- rewrite each passing file with kept rows ---
    GC.gc(false)
    for (fi, p) in enumerate(paths)
        rc[fi] == 0 && continue
        df = DataFrame(Arrow.Table(p); copycols=true)
        km = keep[off[fi]+1 : off[fi]+rc[fi]]
        df = df[km, :]
        tmp = p * ".pf"; writeArrow(tmp, df); safeRm(p, nothing; force=true); mv(tmp, p; force=true)
    end
    return [PSMFileReference(p) for p in paths]
end
