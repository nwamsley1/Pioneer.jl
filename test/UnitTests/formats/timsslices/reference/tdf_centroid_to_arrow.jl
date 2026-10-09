# Copyright (C) 2026 Nathan Wamsley
#
# This file is part of TimsSlices.jl
#
# TimsSlices.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# Bruker TDF (.d) -> Pioneer Arrow with IM smoothing + m/z centroiding per window, one row per kept IM slice.
#
#   For every frame and DIA window (MS1: the whole frame), output slices at scans s0, s0+stride, ... :
#     1. IM kernel: Gaussian (sigma im_sigma scans) over the window's raw entries around the slice scan,
#        accumulated per TOF bin (sum of equal bins);
#     2. m/z kernel: Gaussian (sigma mz_sigma bins) along the TOF axis within each cluster of nearby bins;
#     3. centroiding: local maxima with footprint walked down both sides (<= max_half bins); m/z = intensity-
#        weighted mean bin (wmean) or 3-point Gaussian apex (gauss); intensity = footprint sum;
#     4. cull: centroids with intensity below the raw file's cull_q intensity quantile are dropped.
#   Row = (frame, window, slice scan): imScan = slice scan, RT / eV at the slice scan, centerMz / width from the window.
#
# Usage: julia -t 12 --project=proto tdf_centroid_to_arrow.jl <run.d> <out_dir> [--im-sigma 5] [--mz-sigma 1] [--stride 8]
#        [--cull-q 0.01] [--centroid wmean|gauss] [--max-half 4] [--zstd] [--ms1-stride K] [--ms1-cull-q Q] [--split-cull]
#   MS1 frames use --ms1-stride / --ms1-cull-q (default: the MS2 values). --split-cull (implied by a different
#   --ms1-cull-q) derives each level's cull threshold from a 40-frame sample of that level only, instead of one
#   pooled sample (which MS1 entries dominate).
include(joinpath(@__DIR__, "tdf_to_arrow.jl"))
using Arrow, Statistics, Printf, Base.Threads

struct CParams
    im_sigma::Float64; mz_sigma::Float64; stride::Int; cull_q::Float64; centroid::Symbol; max_half::Int
    # sum_scale: scale the IM kernel to sum to `stride` so a centroid's intensity is the ion's summed
    # intensity over the scans the slice stands for (comparable to a bin-N merge); false = kernel sums to 1
    # (per-scan average), which makes a raw-quantile cull far too strict on low-input data.
    sum_scale::Bool
    # min_scans: keep a centroid only if at least this many raw entries (from distinct scans, since a
    # scan holds one entry per ion) fall inside its footprint -- a persistence cull that removes
    # single-scan events regardless of intensity. 1 = off.
    min_scans::Int
end
CParams(a, b, c, d, e, f) = CParams(a, b, c, d, e, f, false, 1)
CParams(a, b, c, d, e, f, g) = CParams(a, b, c, d, e, f, g, 1)
gauss_kernel(sigma) = sigma <= 0 ? [1.0] : (h = ceil(Int, 3sigma); k = [exp(-(x^2) / (2sigma^2)) for x in -h:h]; k ./ sum(k))

# Per-thread scratch
mutable struct Scratch
    tof::Vector{Int}; val::Vector{Float64}; ord::Vector{Int}; dense::Vector{Float64}; sm::Vector{Float64}
    cnt::Vector{Int}                                   # raw entries per bin of the current dense run
    out_mz::Vector{Float32}; out_int::Vector{Float32}
end
Scratch() = Scratch(Int[], Float64[], Int[], Float64[], Float64[], Int[], Float32[], Float32[])

# Centroid one dense run `sm` (bins t0 .. t0+L-1) into (tof position, intensity) pairs appended to out.
# `cnt` holds the number of raw entries per bin (before smoothing) for the persistence cull.
function centroid_dense!(out_mz, out_int, sm, cnt, L, t0, thr, p::CParams, cal)
    @inbounds for j in 2:L-1
        # local maximum (any positive apex); the cull threshold applies to the footprint SUM below
        (sm[j] > 0 && sm[j] > sm[j-1] && sm[j] >= sm[j+1]) || continue
        lo = j; while lo > 1 && lo > j - p.max_half && sm[lo-1] < sm[lo] && sm[lo-1] > 0; lo -= 1; end
        hi = j; while hi < L && hi < j + p.max_half && sm[hi+1] < sm[hi] && sm[hi+1] > 0; hi += 1; end
        s = 0.0; ws = 0.0; n_entries = 0
        for k in lo:hi; s += sm[k]; ws += sm[k] * (t0 + k - 1); n_entries += cnt[k]; end
        s >= thr || continue
        n_entries >= p.min_scans || continue
        if p.centroid == :gauss
            l0 = log(sm[j]); lm = log(max(sm[j-1], 1e-9)); lp = log(max(sm[j+1], 1e-9))
            den = lm - 2l0 + lp
            off = den < 0 ? clamp(0.5 * (lm - lp) / den, -1.0, 1.0) : 0.0
            tpos = t0 + j - 1 + off
        else
            tpos = ws / s
        end
        push!(out_mz, Float32(tof_to_mz(cal, tpos))); push!(out_int, Float32(s))
    end
end

# One window of one frame -> slices. scan_start: 1-based peak start per 0-based scan (length nscans+1).
function centroid_window!(slices, sc::Scratch, scan_start, tofs, ints, s0, s1, p::CParams, kim, kmz, thr, cal)
    h_im = length(kim) ÷ 2; h_mz = length(kmz) ÷ 2
    gap_max = 2h_mz + 2
    for slice in s0:p.stride:(s1 - 1)
        empty!(sc.tof); empty!(sc.val)
        for s in max(s0, slice - h_im):min(s1 - 1, slice + h_im)
            w = kim[s - slice + h_im + 1]
            @inbounds for k in scan_start[s + 1]:(scan_start[s + 2] - 1)
                push!(sc.tof, Int(tofs[k])); push!(sc.val, w * ints[k])
            end
        end
        isempty(sc.tof) && continue
        resize!(sc.ord, length(sc.tof)); sortperm!(sc.ord, sc.tof)
        empty!(sc.out_mz); empty!(sc.out_int)
        # walk runs of nearby bins
        i = 1; n = length(sc.ord)
        while i <= n
            t_first = sc.tof[sc.ord[i]]; jj = i; t_last = t_first
            while jj < n && sc.tof[sc.ord[jj + 1]] - t_last <= gap_max
                jj += 1; t_last = sc.tof[sc.ord[jj]]
            end
            t0 = t_first - h_mz; L = t_last - t_first + 1 + 2h_mz
            length(sc.dense) < L && (resize!(sc.dense, L); resize!(sc.sm, L); resize!(sc.cnt, L))
            fill!(view(sc.dense, 1:L), 0.0); fill!(view(sc.sm, 1:L), 0.0); fill!(view(sc.cnt, 1:L), 0)
            @inbounds for q in i:jj
                k = sc.ord[q]; sc.dense[sc.tof[k] - t0 + 1] += sc.val[k]; sc.cnt[sc.tof[k] - t0 + 1] += 1
            end
            if h_mz == 0
                @inbounds for a in 1:L; sc.sm[a] = sc.dense[a]; end
            else
                @inbounds for a in 1:L
                    acc = 0.0
                    for (q, w) in enumerate(kmz)
                        b = a + q - h_mz - 1
                        1 <= b <= L && (acc += w * sc.dense[b])
                    end
                    sc.sm[a] = acc
                end
            end
            centroid_dense!(sc.out_mz, sc.out_int, sc.sm, sc.cnt, L, t0, thr, p, cal)
            i = jj + 1
        end
        isempty(sc.out_mz) || push!(slices, (slice, copy(sc.out_mz), copy(sc.out_int)))
    end
end

function convert_centroid(dir::AbstractString, out_dir::AbstractString, p::CParams; zstd::Bool = false, batch::Int = 48,
                          ms1_stride::Int = p.stride, ms1_cull_q::Float64 = p.cull_q, split_cull::Bool = false)
    b = open_tdf(dir); fr = b.frames; nfr = nrow(fr)
    groups, frame_group = read_windows(b.db)
    cal, resid = regressed_mz_cal(b.db, Int(fr.MzCalibration[1]))
    println("mz cal residuals on reference peaks (ppm): ", round.(resid, digits = 2))
    imcal = boundary_im_cal(b.meta, fr)
    mz_lo = Float32(parse(Float64, b.meta["MzAcqRangeLower"])); mz_hi = Float32(parse(Float64, b.meta["MzAcqRangeUpper"]))
    ev_line = let xs = Float64[], ys = Float64[]
        for ws in values(groups), w in ws; push!(xs, (w.scan_begin + w.scan_end) / 2); push!(ys, w.ce); end
        xm = mean(xs); ym = mean(ys); sl = sum((xs .- xm) .* (ys .- ym)) / sum((xs .- xm) .^ 2); (ym - sl * xm, sl)
    end
    ev_at(s) = Float32(ev_line[1] + ev_line[2] * s)
    kim = gauss_kernel(p.im_sigma) .* (p.sum_scale ? p.stride : 1); kmz = gauss_kernel(p.mz_sigma)
    # MS1 frames: own stride (and kernel scale) and own cull quantile
    p1 = CParams(p.im_sigma, p.mz_sigma, ms1_stride, ms1_cull_q, p.centroid, p.max_half, p.sum_scale, p.min_scans)
    kim1 = gauss_kernel(p.im_sigma) .* (p.sum_scale ? ms1_stride : 1)
    split_cull = split_cull || ms1_cull_q != p.cull_q

    # cull threshold(s) from the raw intensity quantile of a 40-frame sample: one pooled sample of both levels,
    # or (split_cull) one sample per level -> thr1 (MS1), thr2 (MS2)
    valid = [i for i in 1:nfr if fr.MsMsType[i] == 0 || fr.MsMsType[i] == 9]
    sample_of(idxs) = idxs[round.(Int, range(1, length(idxs); length = min(40, length(idxs))))]
    function raw_of(idxs)
        v = Float64[]
        for i in idxs
            _, _, ints = read_frame(b, Int(fr.Id[i])); append!(v, Float64.(ints))
        end
        v
    end
    if split_cull
        raw1 = raw_of(sample_of(filter(i -> fr.MsMsType[i] == 0, valid)))
        raw2 = raw_of(sample_of(filter(i -> fr.MsMsType[i] == 9, valid)))
        thr1 = ms1_cull_q > 0 ? quantile(raw1, ms1_cull_q) : 0.0
        thr2 = p.cull_q > 0 ? quantile(raw2, p.cull_q) : 0.0
        println("params: ", p, "; MS1 stride $ms1_stride; split cull: MS1 q$ms1_cull_q = $(round(thr1, digits = 1)) counts ($(length(raw1)) MS1 entries), MS2 q$(p.cull_q) = $(round(thr2, digits = 1)) counts ($(length(raw2)) MS2 entries)")
    else
        raw_sample = raw_of(sample_of(valid))
        thr1 = thr2 = p.cull_q > 0 ? quantile(raw_sample, p.cull_q) : 0.0
        println("params: ", p, "; MS1 stride $ms1_stride; raw intensity q$(p.cull_q) = $(round(thr2, digits = 1)) counts from $(length(raw_sample)) sampled entries -> cull threshold")
    end

    # outputs per frame
    frame_slices = Vector{Vector{Tuple{Int, Vector{Float32}, Vector{Float32}}}}(undef, nfr)
    for i in 1:nfr; frame_slices[i] = Tuple{Int, Vector{Float32}, Vector{Float32}}[]; end
    scratch = [Scratch() for _ in 1:Threads.maxthreadid()]
    t0 = time(); done = 0
    idx = 1
    while idx <= length(valid)
        ids = valid[idx:min(idx + batch - 1, length(valid))]
        decoded = [read_frame(b, Int(fr.Id[i])) for i in ids]      # serial decode (shared file handle)
        @threads for q in eachindex(ids)
            i = ids[q]; scans, tofs, ints = decoded[q]
            nscans = Int(fr.NumScans[i])
            scan_start = zeros(Int, nscans + 1)
            @inbounds for k in eachindex(scans); scan_start[scans[k] + 2] += 1; end
            scan_start[1] = 1
            @inbounds for s in 2:nscans+1; scan_start[s] += scan_start[s-1]; end
            sc = scratch[threadid()]
            if fr.MsMsType[i] == 0
                centroid_window!(frame_slices[i], sc, scan_start, tofs, ints, 0, nscans, p1, kim1, kmz, thr1, cal)
            else
                for w in groups[frame_group[Int(fr.Id[i])]]
                    centroid_window!(frame_slices[i], sc, scan_start, tofs, ints, w.scan_begin, w.scan_end, p, kim, kmz, thr2, cal)
                end
            end
        end
        done += length(ids); idx += batch
        done % (batch * 20) == 0 && println("  $done / $(length(valid)) frames, $(round(time() - t0, digits = 1)) s")
    end
    close(b.bin)
    println("processed $(length(valid)) frames in $(round(time() - t0, digits = 1)) s")

    # flatten (slices are in frame order; within a frame in window/scan order -> RT ascending within the frame)
    nrows = sum(length(v) for v in frame_slices); npk = sum(length(s[2]) for v in frame_slices for s in v)
    T = Union{Missing, Float32}
    mz_flat = Vector{T}(undef, npk); int_flat = Vector{T}(undef, npk); starts = Vector{Int}(undef, nrows + 1)
    retentionTime = Vector{Float32}(undef, nrows); TIC = Vector{Float32}(undef, nrows)
    centerMz = Vector{T}(undef, nrows); isolationWidthMz = Vector{T}(undef, nrows); collisionEnergyField = Vector{T}(undef, nrows)
    collisionEnergyEv = Vector{Float32}(undef, nrows); msOrder = Vector{UInt8}(undef, nrows); cycle_idx = Vector{Int32}(undef, nrows)
    frameId = Vector{Int32}(undef, nrows); imScan = Vector{UInt16}(undef, nrows); windowGroup = Vector{UInt8}(undef, nrows)
    row = 0; pk = 0; cycle = Int32(0)
    for i in 1:nfr
        mst = fr.MsMsType[i]; (mst == 0 || mst == 9) || continue
        mst == 0 && (cycle += Int32(1))
        fid = Int(fr.Id[i]); nscans = Int(fr.NumScans[i])
        rt0 = Float64(fr.Time[i]); dt = Float64(fr.RampTime[i]) / 1000 / nscans
        wins = mst == 0 ? (DiaWindow(0, nscans, 0f0, 0f0, 0f0),) : groups[frame_group[fid]]
        for (s, mzs, its) in frame_slices[i]
            w = mst == 0 ? wins[1] : wins[findfirst(w -> w.scan_begin <= s < w.scan_end, wins)]
            row += 1; starts[row] = pk + 1
            @inbounds for k in eachindex(mzs); pk += 1; mz_flat[pk] = mzs[k]; int_flat[pk] = its[k]; end
            retentionTime[row] = Float32((rt0 + s * dt) / 60); TIC[row] = Float32(sum(its))
            msOrder[row] = mst == 0 ? 0x01 : 0x02
            centerMz[row] = mst == 0 ? missing : w.center; isolationWidthMz[row] = mst == 0 ? missing : w.width
            collisionEnergyField[row] = mst == 0 ? missing : w.ce; collisionEnergyEv[row] = mst == 0 ? 0f0 : ev_at(s)
            cycle_idx[row] = cycle; frameId[row] = Int32(fid); imScan[row] = UInt16(s)
            windowGroup[row] = mst == 0 ? 0x00 : UInt8(frame_group[fid])
        end
    end
    starts[row + 1] = pk + 1
    mz_views = [view(mz_flat, starts[r]:starts[r+1]-1) for r in 1:nrows]
    int_views = [view(int_flat, starts[r]:starts[r+1]-1) for r in 1:nrows]
    tbl = (mz_array = mz_views, intensity_array = int_views,
           scanHeader = fill("", nrows), scanNumber = Int32.(1:nrows), packetType = zeros(Int32, nrows),
           retentionTime = retentionTime, lowMz = fill(mz_lo, nrows), highMz = fill(mz_hi, nrows), TIC = TIC,
           centerMz = centerMz, isolationWidthMz = isolationWidthMz,
           collisionEnergyField = collisionEnergyField, collisionEnergyEvField = collisionEnergyEv,
           msOrder = msOrder, cycle_idx = cycle_idx, frameId = frameId, imScan = imScan, windowGroup = windowGroup)
    meta = Dict("source" => basename(rstrip(dir, '/')), "instrument" => get(b.meta, "InstrumentName", ""),
                "mz_cal_sqrt_intercept" => string(cal.intercept), "mz_cal_sqrt_slope" => string(cal.slope),
                "im_scan0_1overK0" => string(imcal.intercept), "im_slope_1overK0_per_scan" => string(imcal.slope),
                "ce_ev_intercept" => string(ev_line[1]), "ce_ev_slope_per_scan" => string(ev_line[2]),
                "NumScans" => string(maximum(fr.NumScans)), "OneOverK0AcqRangeLower" => b.meta["OneOverK0AcqRangeLower"],
                "OneOverK0AcqRangeUpper" => b.meta["OneOverK0AcqRangeUpper"],
                "centroid_im_sigma" => string(p.im_sigma), "centroid_mz_sigma" => string(p.mz_sigma), "centroid_stride" => string(p.stride),
                "centroid_cull_q" => string(p.cull_q), "centroid_cull_thr" => string(thr2), "centroid_method" => string(p.centroid),
                "centroid_sum_scale" => string(p.sum_scale), "centroid_min_scans" => string(p.min_scans),
                "centroid_ms1_stride" => string(ms1_stride), "centroid_ms1_cull_q" => string(ms1_cull_q),
                "centroid_split_cull" => string(split_cull), "centroid_cull_thr_ms1" => string(thr1), "centroid_cull_thr_ms2" => string(thr2))
    mkpath(out_dir)
    name = replace(basename(rstrip(dir, '/')), r"\.d$" => "") *
           @sprintf("_cen_s%g_m%g_k%d_q%g_%s", p.im_sigma, p.mz_sigma, p.stride, p.cull_q, p.centroid) * (p.sum_scale ? "_sum" : "") *
           (p.min_scans > 1 ? "_n$(p.min_scans)" : "") *
           (ms1_stride != p.stride ? "_ms1k$(ms1_stride)" : "") *
           (ms1_cull_q != p.cull_q ? @sprintf("_ms1q%g", ms1_cull_q) : (split_cull ? "_splitq" : ""))
    out = joinpath(out_dir, name * (zstd ? ".zstd.arrow" : ".arrow"))
    t1 = time()
    # Arrow list columns use Int32 offsets per record batch (max ~2.1 G peaks); split the rows into several
    # record batches when needed. Pioneer reads multi-batch files (ChainedVector columns).
    max_peaks_per_batch = 1_000_000_000
    if npk > max_peaks_per_batch
        bounds = Int[1]
        for r in 1:nrows
            (starts[r + 1] - starts[bounds[end]]) > max_peaks_per_batch && push!(bounds, r)
        end
        push!(bounds, nrows + 1)
        chunks = [map(c -> view(c, bounds[i]:bounds[i+1]-1), tbl) for i in 1:length(bounds)-1]
        src = Arrow.Tables.partitioner(identity, chunks)
        println("peaks $npk > $max_peaks_per_batch: writing $(length(chunks)) record batches")
    else
        src = tbl
    end
    zstd ? Arrow.write(out, src; metadata = meta, compress = :zstd) : Arrow.write(out, src; metadata = meta)
    println("rows $nrows, peaks $npk; wrote $out  $(round(filesize(out) / 1e9, digits = 3)) GB in $(round(time() - t1, digits = 1)) s")
    return out
end

if abspath(PROGRAM_FILE) == @__FILE__
    getarg(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i + 1])
    p = CParams(parse(Float64, getarg("--im-sigma", "5")), parse(Float64, getarg("--mz-sigma", "1")), parse(Int, getarg("--stride", "8")),
                parse(Float64, getarg("--cull-q", "0.01")), Symbol(getarg("--centroid", "wmean")), parse(Int, getarg("--max-half", "4")),
                "--sum-scale" in ARGS, parse(Int, getarg("--min-scans", "1")))
    convert_centroid(ARGS[1], ARGS[2], p; zstd = "--zstd" in ARGS,
                     ms1_stride = parse(Int, getarg("--ms1-stride", getarg("--stride", "8"))),
                     ms1_cull_q = parse(Float64, getarg("--ms1-cull-q", getarg("--cull-q", "0.01"))),
                     split_cull = "--split-cull" in ARGS)
end
