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

# Prototype: Bruker TDF (.d) -> Pioneer Arrow, one row per (frame, IM scan) packet. No summing, no merging.
#
#   row = (frame f, scan s) for s inside some quad window w of frame f (MS1: every scan)
#     msOrder          1 (MS1 frame) / 2 (diaPASEF frame)
#     retentionTime    (Frames.Time + s * RampTime/NumScans) in MINUTES (Pioneer convention)
#     centerMz, width  window w's IsolationMz / IsolationWidth (MS1: missing)
#     imScan           s (UInt16); 1/K0 = OneOverK0AcqRangeUpper + (Lower-Upper)/NumScans * s, stored in file metadata
#     mz_array         that scan's peaks, TOF bin -> m/z via the line fitted to CalibrationInfo reference peaks
#     frameId, windowGroup  kept as columns
#
# Usage: julia --project=. tdf_to_arrow.jl <run.d> <out_dir> [--no-zstd]

include(joinpath(@__DIR__, "tdf_reader.jl"))
using Arrow

# --- m/z calibration fitted to the tune-mix reference peaks stored in the file ---------------------
function regressed_mz_cal(db::SQLite.DB, mzcal_id::Integer)
    function blob(key)
        r = first(DBInterface.execute(db, "SELECT Value FROM CalibrationInfo WHERE KeyPolarity='+' AND KeyName='$key'"))
        v = r[:Value]
        bytes = v isa Vector{UInt8} ? v : Vector{UInt8}(v)
        return collect(reinterpret(Float64, bytes))
    end
    masses = blob("ReferencePeakMasses"); t = blob("MeasuredTimesOfFlight")
    r = first(DBInterface.execute(db, "SELECT DigitizerTimebase, DigitizerDelay FROM MzCalibration WHERE Id=$mzcal_id"))
    x = (t .- r[:DigitizerDelay]) ./ r[:DigitizerTimebase]
    y = sqrt.(masses)
    xm = sum(x) / length(x); ym = sum(y) / length(y)
    b = sum((x .- xm) .* (y .- ym)) / sum((x .- xm) .^ 2)
    a = ym - b * xm
    resid_ppm = [((a + b * xi)^2 - m) / m * 1e6 for (xi, m) in zip(x, masses)]
    return LinearMzCal(a, b), resid_ppm
end

struct DiaWindow
    scan_begin::Int   # inclusive
    scan_end::Int     # exclusive
    center::Float32
    width::Float32
    ce::Float32
end

function read_windows(db::SQLite.DB)
    groups = Dict{Int,Vector{DiaWindow}}()
    for r in DBInterface.execute(db, "SELECT WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth, CollisionEnergy FROM DiaFrameMsMsWindows ORDER BY WindowGroup, ScanNumBegin")
        push!(get!(groups, Int(r[:WindowGroup]), DiaWindow[]),
              DiaWindow(Int(r[:ScanNumBegin]), Int(r[:ScanNumEnd]), Float32(r[:IsolationMz]), Float32(r[:IsolationWidth]), Float32(r[:CollisionEnergy])))
    end
    frame_group = Dict{Int,Int}()
    for r in DBInterface.execute(db, "SELECT Frame, WindowGroup FROM DiaFrameMsMsInfo")
        frame_group[Int(r[:Frame])] = Int(r[:WindowGroup])
    end
    return groups, frame_group
end

# im_bin > 1: sum blocks of im_bin adjacent IM scans inside each window (MS1: inside the frame) into one
# row. TOF bins are one global axis per run, so merging = sort by TOF bin and add intensities of equal
# bins. The row's imScan is the block's centre scan; RT / eV are evaluated at the centre scan.
function convert(dir::AbstractString, out_dir::AbstractString; zstd::Bool=true, im_bin::Int=1, plain::Bool=true, merge_tol::Int=0)
    im_bin >= 1 || error("im_bin must be >= 1")
    merge_tol >= 0 || error("merge_tol must be >= 0")
    (zstd || plain) || error("nothing to write: zstd and plain both false")
    b = open_tdf(dir)
    fr = b.frames
    nfr = nrow(fr)
    groups, frame_group = read_windows(b.db)
    mzcal_id = Int(fr.MzCalibration[1])
    cal, resid = regressed_mz_cal(b.db, mzcal_id)
    println("mz cal: sqrt(mz) = $(cal.intercept) + $(cal.slope) * tof   (residual on reference peaks, ppm: $(round.(resid, digits=2)))")
    imcal = boundary_im_cal(b.meta, fr)
    mz_lo = Float32(parse(Float64, b.meta["MzAcqRangeLower"])); mz_hi = Float32(parse(Float64, b.meta["MzAcqRangeUpper"]))

    # ---- collision energy per scan: CE is a linear ramp in the IM scan index; fit it from the
    #      window table (each window stores the ramp value at its mid-scan) --------------------
    ev_line = let xs = Float64[], ys = Float64[]
        for ws in values(groups), w in ws
            push!(xs, (w.scan_begin + w.scan_end) / 2); push!(ys, w.ce)
        end
        xm = sum(xs) / length(xs); ym = sum(ys) / length(ys)
        slope_ev = sum((xs .- xm) .* (ys .- ym)) / sum((xs .- xm) .^ 2); icpt_ev = ym - slope_ev * xm
        println("collision energy ramp: eV = $icpt_ev + $slope_ev * scan  (max residual $(maximum(abs.(ys .- (icpt_ev .+ slope_ev .* xs)))) eV over $(length(xs)) windows)")
        (icpt_ev, slope_ev)
    end
    ev_at(s) = Float32(ev_line[1] + ev_line[2] * s)

    # ---- row count ------------------------------------------------------------------------------
    nrows = 0
    for i in 1:nfr
        if fr.MsMsType[i] == 0
            nrows += cld(Int(fr.NumScans[i]), im_bin)
        elseif fr.MsMsType[i] == 9
            nrows += sum(w -> cld(w.scan_end - w.scan_begin, im_bin), groups[frame_group[Int(fr.Id[i])]])
        end
    end
    npeaks_total = Int(sum(fr.NumPeaks))
    println("frames $nfr -> rows $nrows (im_bin = $im_bin), peaks <= $npeaks_total")
    blk_tof = Int[]; blk_int = Float64[]      # scratch for merging a block of IM scans

    # ---- flat peak buffers + per-row columns ----------------------------------------------------
    T = Union{Missing,Float32}
    mz_flat = Vector{T}(undef, npeaks_total); int_flat = Vector{T}(undef, npeaks_total)
    starts = Vector{Int}(undef, nrows + 1)
    retentionTime = Vector{Float32}(undef, nrows); TIC = Vector{Float32}(undef, nrows)
    centerMz = Vector{T}(undef, nrows); isolationWidthMz = Vector{T}(undef, nrows); collisionEnergyField = Vector{T}(undef, nrows)
    collisionEnergyEv = Vector{Float32}(undef, nrows)   # per-scan ramp value (eV); 0 for MS1
    msOrder = Vector{UInt8}(undef, nrows); cycle_idx = Vector{Int32}(undef, nrows)
    frameId = Vector{Int32}(undef, nrows); imScan = Vector{UInt16}(undef, nrows); windowGroup = Vector{UInt8}(undef, nrows)

    row = 0; pk = 0; cycle = Int32(0)
    t0 = time()
    for i in 1:nfr
        mst = fr.MsMsType[i]
        (mst == 0 || mst == 9) || continue
        fid = Int(fr.Id[i]); nscans = Int(fr.NumScans[i])
        scans, tofs, ints = read_frame(b, fid)
        # per-scan ranges within this frame (peaks are emitted in ascending scan order)
        scan_start = zeros(Int, nscans + 1)          # 1-based start index of scan s = scan_start[s+1]
        @inbounds for k in eachindex(scans); scan_start[scans[k] + 2] += 1; end
        scan_start[1] = 1
        @inbounds for s in 2:nscans+1; scan_start[s] += scan_start[s-1]; end
        rt0 = Float64(fr.Time[i]); dt = Float64(fr.RampTime[i]) / 1000 / nscans   # seconds per scan
        if mst == 0
            cycle += Int32(1)
            wins = (DiaWindow(0, nscans, 0f0, 0f0, 0f0),)
        else
            wins = groups[frame_group[fid]]
        end
        for w in wins
            s0 = w.scan_begin
            while s0 < w.scan_end
                s1 = min(s0 + im_bin, w.scan_end)            # block = scans s0 .. s1-1
                row += 1
                starts[row] = pk + 1
                tic = 0.0
                if s1 - s0 == 1
                    @inbounds for k in scan_start[s0+1]:(scan_start[s0+2] - 1)
                        pk += 1
                        mz_flat[pk] = Float32(tof_to_mz(cal, tofs[k]))
                        int_flat[pk] = Float32(ints[k])
                        tic += ints[k]
                    end
                else
                    # merge the block's scans on the shared TOF-bin axis: sort by bin, then cluster
                    # entries within merge_tol bins of the running cluster's last bin (0 = exact bin
                    # only). The centroid across scans jitters by ~1-2 bins, so exact merging splits
                    # one ion into adjacent-bin doublets; a cluster is written at its intensity-
                    # weighted mean bin. Peaks within one scan are >= 5 bins apart, so tol <= 2 never
                    # merges two ions of the same scan.
                    empty!(blk_tof); empty!(blk_int)
                    @inbounds for k in scan_start[s0+1]:(scan_start[s1+1] - 1)
                        push!(blk_tof, Int(tofs[k])); push!(blk_int, Float64(ints[k]))
                    end
                    prev_t = -1; c_wsum = 0.0; c_isum = 0.0
                    @inbounds for j in sortperm(blk_tof)
                        t = blk_tof[j]; it = blk_int[j]
                        if prev_t >= 0 && t - prev_t <= merge_tol
                            c_wsum += t * it; c_isum += it
                        else
                            if prev_t >= 0
                                pk += 1
                                mz_flat[pk] = Float32(tof_to_mz(cal, c_wsum / c_isum)); int_flat[pk] = Float32(c_isum)
                            end
                            c_wsum = t * it; c_isum = it
                        end
                        prev_t = t
                        tic += it
                    end
                    if prev_t >= 0
                        pk += 1
                        mz_flat[pk] = Float32(tof_to_mz(cal, c_wsum / c_isum)); int_flat[pk] = Float32(c_isum)
                    end
                end
                sc = (s0 + s1 - 1) / 2                        # block centre scan
                retentionTime[row] = Float32((rt0 + sc * dt) / 60)
                TIC[row] = Float32(tic)
                msOrder[row] = mst == 0 ? 0x01 : 0x02
                centerMz[row] = mst == 0 ? missing : w.center
                isolationWidthMz[row] = mst == 0 ? missing : w.width
                collisionEnergyField[row] = mst == 0 ? missing : w.ce      # window value (ramp at window mid-scan)
                collisionEnergyEv[row] = mst == 0 ? 0f0 : ev_at(sc)        # this packet's actual eV on the ramp
                cycle_idx[row] = cycle
                frameId[row] = Int32(fid); imScan[row] = UInt16(round(Int, sc))
                windowGroup[row] = mst == 0 ? 0x00 : UInt8(frame_group[fid])
                s0 = s1
            end
        end
    end
    starts[row + 1] = pk + 1
    @assert row == nrows
    println("decoded + calibrated $row rows, $pk peaks in $(round(time() - t0, digits=1)) s")

    mz_views = [view(mz_flat, starts[r]:starts[r+1]-1) for r in 1:nrows]
    int_views = [view(int_flat, starts[r]:starts[r+1]-1) for r in 1:nrows]
    tbl = (
        mz_array = mz_views, intensity_array = int_views,
        scanHeader = fill("", nrows), scanNumber = Int32.(1:nrows), packetType = zeros(Int32, nrows),
        retentionTime = retentionTime, lowMz = fill(mz_lo, nrows), highMz = fill(mz_hi, nrows), TIC = TIC,
        centerMz = centerMz, isolationWidthMz = isolationWidthMz,
        collisionEnergyField = collisionEnergyField, collisionEnergyEvField = collisionEnergyEv,
        msOrder = msOrder, cycle_idx = cycle_idx,
        frameId = frameId, imScan = imScan, windowGroup = windowGroup,
    )
    meta = Dict(
        "source" => basename(rstrip(dir, '/')), "instrument" => get(b.meta, "InstrumentName", ""),
        "mz_cal_sqrt_intercept" => string(cal.intercept), "mz_cal_sqrt_slope" => string(cal.slope),
        "im_scan0_1overK0" => string(imcal.intercept), "im_slope_1overK0_per_scan" => string(imcal.slope),
        "ce_ev_intercept" => string(ev_line[1]), "ce_ev_slope_per_scan" => string(ev_line[2]),
        "NumScans" => string(maximum(fr.NumScans)), "OneOverK0AcqRangeLower" => b.meta["OneOverK0AcqRangeLower"],
        "OneOverK0AcqRangeUpper" => b.meta["OneOverK0AcqRangeUpper"], "im_bin" => string(im_bin),
        "merge_tol_bins" => string(merge_tol),
    )
    mkpath(out_dir)
    name = replace(basename(rstrip(dir, '/')), r"\.d$" => "") * (im_bin > 1 ? "_imbin$(im_bin)" : "") *
           (im_bin > 1 && merge_tol > 0 ? "_tol$(merge_tol)" : "")
    out = joinpath(out_dir, name * ".arrow")
    if plain
        t1 = time(); Arrow.write(out, tbl; metadata = meta);
        println("wrote $out  $(round(filesize(out)/1e9, digits=3)) GB in $(round(time() - t1, digits=1)) s")
    end
    if zstd
        outz = joinpath(out_dir, name * ".zstd.arrow")
        t2 = time(); Arrow.write(outz, tbl; metadata = meta, compress = :zstd)
        println("wrote $outz  $(round(filesize(outz)/1e9, digits=3)) GB in $(round(time() - t2, digits=1)) s")
        plain || (out = outz)
    end
    close(b.bin)
    return out
end

# CLI: tdf_to_arrow.jl <run.d> <out_dir> [--no-zstd | --zstd-only] [--im-bin N] [--merge-tol K]
if abspath(PROGRAM_FILE) == @__FILE__
    ib = let i = findfirst(==("--im-bin"), ARGS); i === nothing ? 1 : parse(Int, ARGS[i + 1]) end
    mt = let i = findfirst(==("--merge-tol"), ARGS); i === nothing ? 0 : parse(Int, ARGS[i + 1]) end
    zo = "--zstd-only" in ARGS
    convert(ARGS[1], ARGS[2]; zstd = zo || !("--no-zstd" in ARGS), im_bin = ib, plain = !zo, merge_tol = mt)
end
