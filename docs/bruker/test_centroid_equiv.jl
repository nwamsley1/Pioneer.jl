# Check that the sparse converter (tdf_centroid_to_arrow.jl) reproduces the dense-matrix centroids of
# im_mz_centroid.jl / centroid_lib.jl on one window (frame 4092, window 1 of the E. coli file, m/z 905-930).
# Usage: julia --project=proto test_centroid_equiv.jl <run.d> <map.arrow>
include(joinpath(@__DIR__, "tdf_centroid_to_arrow.jl"))
include(joinpath(@__DIR__, "centroid_lib.jl"))
using Arrow, Printf
dpath, mappath = ARGS[1], ARGS[2]
p = CParams(5.0, 1.0, 8, 0.0, :wmean, 4)
kim = gauss_kernel(p.im_sigma); kmz = gauss_kernel(p.mz_sigma)
b = open_tdf(dpath); fr = b.frames
groups, frame_group = read_windows(b.db)
cal, _ = regressed_mz_cal(b.db, Int(fr.MzCalibration[4092]))
fid = 4092; w = groups[frame_group[fid]][1]
scans, tofs, ints = read_frame(b, fid); nscans = Int(fr.NumScans[fid])
scan_start = zeros(Int, nscans + 1)
for k in eachindex(scans); scan_start[scans[k] + 2] += 1; end
scan_start[1] = 1
for s in 2:nscans+1; scan_start[s] += scan_start[s-1]; end
slices = Tuple{Int, Vector{Float32}, Vector{Float32}}[]
centroid_window!(slices, Scratch(), scan_start, tofs, ints, w.scan_begin, w.scan_end, p, kim, kmz, 0.1, cal)
close(b.bin)
# dense reference on the map (905-930): same slices (window scans start at w.scan_begin, stride 8)
t = Arrow.Table(mappath)
M = reshape(Float32.(t.matrix[1]), Int(t.nscan[1]), Int(t.ntof[1])); scan0 = Int(t.scan0[1]); mz_all = Float64.(t.mz_axis[1])
S = smooth_dim(smooth_dim(M, gauss_kernel(5.0), 1), gauss_kernel(1.0), 2)
lo_edge, hi_edge = mz_all[1] + 0.05, mz_all[end] - 0.05      # ignore the map's edges (dense map has no data beyond them)
n_cmp = 0; max_dmz_ppm = 0.0; max_dint_rel = 0.0; n_missing = 0
for (s, mzs, its) in slices
    i = s - scan0 + 1
    # converter rule: any positive apex, footprint SUM >= thr -> same rule for the dense reference
    dense = [c for c in centroid_row(@view(S[i, :]), mz_all, 0.0; min_width = 1) if c.intensity >= 0.1 && lo_edge <= c.mz_wmean <= hi_edge]
    sparse = [(mzs[k], its[k]) for k in eachindex(mzs) if lo_edge <= mzs[k] <= hi_edge]
    length(dense) == length(sparse) || (global n_missing += abs(length(dense) - length(sparse)))
    for (c, (mz, it)) in zip(sort(dense; by = c -> c.mz_wmean), sort(sparse))
        global n_cmp += 1
        global max_dmz_ppm = max(max_dmz_ppm, 1e6 * abs(c.mz_wmean - mz) / mz)
        global max_dint_rel = max(max_dint_rel, abs(c.intensity - it) / max(c.intensity, 1e-9))
    end
end
println(@sprintf("slices %d; compared %d centroids; count mismatches %d; max |dm/z| %.3f ppm; max relative intensity diff %.2e",
                 length(slices), n_cmp, n_missing, max_dmz_ppm, max_dint_rel))
