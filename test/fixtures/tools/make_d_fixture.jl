# Build the truncated timsTOF diaPASEF `.d` test fixture (published on Zenodo, record 23023765) from a full run.
#
#   julia --project=<Pioneer> make_d_fixture.jl <source.d> <out.d> [rt_lo_s rt_hi_s prec_lo prec_hi]
#
# The published fixture (ecoli_tims_fixture.d) is the defaults applied to the E. coli 50 ng 5 min diaPASEF run
# (PXD070049, LFQ_Ultra2_diaPASEF_5min_50ng_Ecoli_01): 420-510 s around the densest identifications, and the two
# isolation windows per cycle inside the committed E. coli test library's precursor range (test/integration/
# ecoli_lib.poin, 495.5-601.4): 537.5 and 562.5 m/z, i.e. precursors 526-574. 90 s is needed: with 45 s (even with
# every in-range window) Precursor Scoring has too few PSMs to train and ends with ~0 IDs; at 90 s the fixture
# searches to ~750 precursors / ~530 protein groups at 1% FDR with every stage running. The fixture keeps:
#   - MS2 frames in the RT range, each reduced to the isolation windows overlapping [prec_lo, prec_hi] (frames with
#     none are dropped); the other windows' scans are emptied and the frame re-encoded (TimsSlices.encode_codec2);
#   - MS1 frames in the RT range, reduced to the kept windows' IM scans and to m/z [prec_lo - 25, prec_hi + 15];
#   - `analysis.tdf` with every table that has a `Frame` column, the window tables, Segments and the rewritten
#     Frames fields (TimsId, NumPeaks, MaxIntensity, SummedIntensities) trimmed to match. NumScans is unchanged.
# The script then reopens the fixture and checks every frame decodes to exactly the filtered source peaks.

using Pioneer, SQLite, Printf
const TS = Pioneer.TimsSlices
const DBI = SQLite.DBInterface

src = ARGS[1]; out = ARGS[2]
rt_lo, rt_hi, prec_lo, prec_hi = length(ARGS) >= 6 ? parse.(Float64, ARGS[3:6]) : (420.0, 510.0, 526.0, 574.0)
isdir(out) && error("$out exists")
mkpath(out)

f = TS.open_tdf(src)
fr = f.frames; groups = f.dia.groups; fgroup = f.dia.frame_group
overlaps(w) = w.center - w.width / 2 <= prec_hi && w.center + w.width / 2 >= prec_lo
kept_windows = [filter(overlaps, g) for g in groups]           # per window group
in_rt = [rt_lo <= fr.time[i] <= rt_hi for i in eachindex(fr.id)]
keep_row = [in_rt[i] && (TS.is_ms1(f, i) || (TS.is_dia(f, i) && fgroup[i] > 0 && !isempty(kept_windows[fgroup[i]])))
            for i in eachindex(fr.id)]
rows = findall(keep_row)
all_kept = reduce(vcat, [kept_windows[fgroup[i]] for i in rows if TS.is_dia(f, i)])
ms1_scans = minimum(w.scan_begin for w in all_kept):(maximum(w.scan_end for w in all_kept) - 1)
ms1_bins = floor(UInt32, TS.mz_to_bin(f.mz_cal, prec_lo - 25)):ceil(UInt32, TS.mz_to_bin(f.mz_cal, prec_hi + 15))
@printf("source: %d frames; keeping %d (%d MS1, %d MS2) at %.0f-%.0f s; windows kept per group: %s\n",
        length(fr.id), length(rows), count(i -> TS.is_ms1(f, i), rows), count(i -> TS.is_dia(f, i), rows), rt_lo, rt_hi,
        join([string(length(k)) for k in kept_windows], ","))

"The filtered peaks of frame row i: (tof, intensity, scan_start)."
function filtered_frame(buf, i)
    TS.read_frame!(buf, f, i)
    ms1 = TS.is_ms1(f, i)
    scan_ok = falses(buf.n_scans)
    if ms1
        for s in ms1_scans; s < buf.n_scans && (scan_ok[s + 1] = true); end
    else
        for w in kept_windows[fgroup[i]], s in w.scan_begin:(w.scan_end - 1); s < buf.n_scans && (scan_ok[s + 1] = true); end
    end
    tof = UInt32[]; it = UInt32[]; ss = Int[1]
    for s in 0:buf.n_scans-1
        if scan_ok[s + 1]
            for k in TS.scan_range(buf, s)
                (!ms1 || buf.tof[k] in ms1_bins) || continue
                push!(tof, buf.tof[k]); push!(it, buf.intensity[k])
            end
        end
        push!(ss, length(tof) + 1)
    end
    tof, it, ss, buf.n_scans
end

# --- analysis.tdf_bin -------------------------------------------------------------------------------
new_offset = Dict{Int32, Int64}(); new_stats = Dict{Int32, NTuple{3, Int64}}()   # frame id => (peaks, max, sum)
buf = TS.FrameBuffer()
open(joinpath(out, "analysis.tdf_bin"), "w") do io
    open(joinpath(src, "analysis.tdf_bin")) do s; write(io, read(s, fr.tims_id[1])); end   # keep the file preamble
    for i in rows
        tof, it, ss, ns = filtered_frame(buf, i)
        new_offset[fr.id[i]] = position(io)
        write(io, TS.encode_codec2(tof, it, ss, ns))
        new_stats[fr.id[i]] = (length(tof), isempty(it) ? 0 : Int(maximum(it)), Int(sum(Int, it; init = 0)))
    end
end

# --- analysis.tdf ------------------------------------------------------------------------------------
cp(joinpath(src, "analysis.tdf"), joinpath(out, "analysis.tdf"))
db = SQLite.DB(joinpath(out, "analysis.tdf"))
kept_ids = [fr.id[i] for i in rows]
DBI.execute(db, "BEGIN")
DBI.execute(db, "CREATE TEMP TABLE kept(Id INTEGER PRIMARY KEY)")
for id in kept_ids; DBI.execute(db, "INSERT INTO kept VALUES ($id)"); end
tables = [string(r[:name]) for r in DBI.execute(db, "SELECT name FROM sqlite_master WHERE type = 'table'")]
for t in tables
    cols = [string(r[:name]) for r in DBI.execute(db, "PRAGMA table_info(\"$t\")")]
    "Frame" in cols && DBI.execute(db, "DELETE FROM \"$t\" WHERE Frame NOT IN (SELECT Id FROM kept)")
end
DBI.execute(db, "DELETE FROM Frames WHERE Id NOT IN (SELECT Id FROM kept)")
for (g, ws) in enumerate(groups), w in ws
    w in kept_windows[g] && continue
    DBI.execute(db, "DELETE FROM DiaFrameMsMsWindows WHERE WindowGroup = $g AND ScanNumBegin = $(w.scan_begin)")
end
DBI.execute(db, "DELETE FROM DiaFrameMsMsWindowGroups WHERE Id NOT IN (SELECT DISTINCT WindowGroup FROM DiaFrameMsMsWindows)")
DBI.execute(db, "UPDATE Segments SET FirstFrame = $(minimum(kept_ids)), LastFrame = $(maximum(kept_ids))")
for id in kept_ids
    np, mx, sm = new_stats[id]
    DBI.execute(db, "UPDATE Frames SET TimsId = $(new_offset[id]), NumPeaks = $np, MaxIntensity = $mx, SummedIntensities = $sm WHERE Id = $id")
end
DBI.execute(db, "COMMIT")
DBI.execute(db, "VACUUM")
close(db)

# --- verify ------------------------------------------------------------------------------------------
g = TS.open_tdf(out)
@assert length(g.frames.id) == length(rows) "frame count"
bufg = TS.FrameBuffer(); n_bad = 0
for (j, i) in enumerate(rows)
    tof, it, ss, ns = filtered_frame(buf, i)
    TS.read_frame!(bufg, g, j)
    ok = bufg.n_scans == ns && bufg.n_peaks == length(tof) && bufg.scan_start[1:ns + 1] == ss &&
         bufg.tof[1:length(tof)] == tof && bufg.intensity[1:length(it)] == it
    global n_bad += !ok
end
@assert n_bad == 0 "$n_bad frames differ from the filtered source"
@assert all(length(TS.windows(g, j)) == length(kept_windows[fgroup[rows[j]]]) for j in eachindex(rows) if TS.is_dia(g, j))
@assert g.mz_cal == f.mz_cal && g.im_cal == f.im_cal "calibration changed"
sz(p) = round(filesize(joinpath(out, p)) / 1e6; digits = 2)
@printf("wrote %s: analysis.tdf %.2f MB, analysis.tdf_bin %.2f MB; %d frames verified against the filtered source\n",
        out, sz("analysis.tdf"), sz("analysis.tdf_bin"), length(rows))

# --- checksums (committed: test/UnitTests/formats/timsslices/fixtures/) -----------------------------------
# Raw frames, same columns as hela_first50_checksums.csv, and SHA-256 of the .tdfs a default conversion writes.
using SHA
open(out * "_frames_checksums.csv", "w") do io
    println(io, "row,frame_id,msms_type,n_scans,n_peaks,sum_tof,sum_intensity,sum_scan,max_tof")
    for j in eachindex(g.frames.id)
        TS.read_frame!(bufg, g, j)
        np = bufg.n_peaks
        sum_scan = sum(Int(s) * (bufg.scan_start[s + 2] - bufg.scan_start[s + 1]) for s in 0:bufg.n_scans-1)
        @printf(io, "%d,%d,%d,%d,%d,%d,%d,%d,%d\n", j, g.frames.id[j], g.frames.msms_type[j], bufg.n_scans, np,
                sum(Int, view(bufg.tof, 1:np); init = 0), sum(Int, view(bufg.intensity, 1:np); init = 0), sum_scan,
                np == 0 ? 0 : Int(maximum(view(bufg.tof, 1:np))))
    end
end
tmp = mktempdir()
conv = TS.convert_run(out, tmp; name = "fixture", log = devnull)
open(out * "_tdfs_sha256.csv", "w") do io
    println(io, "file,sha256")
    for fn in ("blocks.bin", "slices.arrow", "frames.arrow")
        println(io, fn, ",", bytes2hex(open(sha256, joinpath(conv.tdfs, fn))))
    end
end
println("wrote ", out, "_frames_checksums.csv and ", out, "_tdfs_sha256.csv")
