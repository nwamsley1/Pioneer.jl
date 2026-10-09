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

# tdfs -> Pioneer slice Arrow (today's schema). Validates the container end to end: expand(convert_run(x)) must equal
# the Arrow that convert writes directly.

"""
    expand(tdfs_dir, arrow_path) -> arrow_path
"""
function expand(tdfs_dir::AbstractString, arrow_path::AbstractString; log::IO = stdout)
    t = open_tdfs(tdfs_dir)
    meta = t.meta
    pd = meta["params"]
    # Parameters a later TimsSlices no longer has (e.g. the quantile culls, removed in 0.1) are ignored.
    p = ConvertParams(; (Symbol(k) => (k in ("centroid", "format") ? Symbol(v) : k == "frames" ? nothing : v)
                         for (k, v) in pd if Symbol(k) in fieldnames(ConvertParams))...)
    aw = SliceArrowWriter(arrow_path, arrow_metadata(meta, p),
                          t.mz_cal, t.bin_scale, t.int_scale, Float64(meta["mz_lo"]), Float64(meta["mz_hi"]))
    blk = SliceBlock(); codec = BlockCodec(); rows = SliceRows()
    fr = t.frames; sl = t.slices
    t0 = time()
    for i in 1:n_frames(t)
        read_frame_block!(blk, codec, t, i)
        fm = FrameMeta(fr.frame_id[i], fr.ms_order[i], fr.window_group[i], fr.rt_s[i], 0.0, fr.n_scans[i])
        empty!(rows)
        r = fr.first_slice[i]:fr.first_slice[i] + fr.n_slices[i] - 1
        append!(rows.im_scan, view(sl.im_scan, r)); append!(rows.window, view(sl.window, r))
        append!(rows.retention_time, view(sl.retention_time, r)); append!(rows.center_mz, view(sl.center_mz, r))
        append!(rows.isolation_width, view(sl.isolation_width, r)); append!(rows.collision_energy_ev, view(sl.collision_energy_ev, r))
        append!(rows.window_ce, view(sl.window_ce, r))
        append!(rows.tic, view(sl.tic, r))
        # format 2 has no base peak column: recompute it from the decoded peaks, as the converter does
        hasproperty(sl, :base_peak_intensity) ? append!(rows.base_peak_intensity, view(sl.base_peak_intensity, r)) :
            append!(rows.base_peak_intensity, (Float32(maximum(view(blk.intensity, slice_range(blk, j)); init = zero(eltype(blk.intensity))) / t.int_scale)
                                               for j in 1:blk.n_slices))
        append!(rows.n_peaks, view(sl.n_peaks, r)); append!(rows.peak_offset, view(sl.peak_offset, r))
        append!(rows.block_size, view(sl.block_size, r))
        write_frame!(aw, fm, rows, blk)
    end
    close(aw)
    @printf(log, "expanded %d frames / %d slices to %s (%.3f GB) in %.1f s\n", n_frames(t), n_slices(t), arrow_path, filesize(arrow_path) / 1e9, time() - t0)
    arrow_path
end
