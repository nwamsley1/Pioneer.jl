# Build the truncated SCIEX SWATH `.wiff` + `.wiff.scan` test fixture (published on Zenodo; see
# dev_docs/formats/FOLD_TIMSSLICES_SCIEXWIFF.md) from a full run.
#
#   julia --project=<Pioneer> make_wiff_fixture.jl <source.wiff> <out_dir> [rt_lo_min rt_hi_min prec_lo prec_hi]
#
# Source: JPST002949 BenchSample_B_nswath4_25ng (three-proteome, ZenoTOF 7600; 173 windows of 2.9 Da per cycle,
# 5.2 s cycles, 80 min). The fixture keeps, inside the RT range, every cycle's MS1 scan and the MS2 scans whose
# isolation window overlaps [prec_lo, prec_hi] (defaults: the committed E. coli test library's precursor range,
# test/integration/ecoli_lib.poin, so that library identifies peptides in it). Everything else becomes an empty
# record.
#
# How, without rewriting the compound file:
#   .wiff.scan  the 44-byte header, then each kept record's bytes copied in order: its protobuf metadata (which the
#               reader locates from the end of the previous non-empty block) and its ffffffff-marked block.
#   .wiff       the `SampleSubtree/Sample1/Idx` stream (one 54-byte record per cycle x experiment slot: u32 block
#               offset relative to the scan header, u32 block size, RT, TIC, ...) is patched IN PLACE, over the
#               stream's own sectors: kept records get their new offsets, dropped records size 0 (the reader treats
#               those as empty scans) at a non-decreasing offset. The stream keeps its length, so no other byte of
#               the compound file changes and the file stays the size of the original.
# The script then reopens the fixture and checks every kept record decodes identically, every dropped one is empty,
# and every other stream is byte-identical.

using Pioneer, Printf, SHA
const S = Pioneer.SciexWiff
const C = S.CFB

src_wiff = ARGS[1]; out_dir = ARGS[2]
rt_lo, rt_hi, prec_lo, prec_hi = length(ARGS) >= 6 ? parse.(Float64, ARGS[3:6]) : (30.0, 50.0, 495.5, 601.4)
isdir(out_dir) && error("$out_dir exists")
mkpath(out_dir)
name = splitext(basename(src_wiff))[1]
out_wiff = joinpath(out_dir, name * ".wiff"); out_scan = out_wiff * ".scan"

run = S.WiffRun(src_wiff)
ix = run.index; n = length(ix)
function keep(k)
    ix.block_size[k] > 0 || return false
    rt_lo <= S.retention_time_min(run, k) <= rt_hi || return false
    w = S.window(run, k)
    w === nothing || (w.lo <= prec_hi && w.hi >= prec_lo)
end
kept = [keep(k) for k in 1:n]
nwin = count(w -> w.lo <= prec_hi && w.hi >= prec_lo, run.windows)
@printf("source: %d records (%d cycles x %d experiments); keeping %d at %.1f-%.1f min: MS1 + %d of %d windows per cycle\n",
        n, n ÷ (length(run.windows) + 1), length(run.windows) + 1, count(kept), rt_lo, rt_hi, nwin, length(run.windows))

# --- .wiff.scan -------------------------------------------------------------------------------------
s = run.scan
new_rel = Vector{UInt32}(undef, n)                 # new u32 offset field per record
open(out_scan, "w") do io
    write(io, view(s, 1:S.SCAN_FILE_HEADER))
    for k in 1:n
        if kept[k]
            ms = ix.meta_start[k]; ff = ix.block_offset[k]; sz = Int64(ix.block_size[k])
            newff = position(io) + (ff - ms)
            write(io, view(s, ms + 1:ff + sz))
            new_rel[k] = UInt32(newff - S.SCAN_FILE_HEADER)
        else
            new_rel[k] = UInt32(position(io) - S.SCAN_FILE_HEADER)   # empty: size 0 at the current end
        end
    end
end

# --- .wiff: patch the Idx stream in place --------------------------------------------------------------
data = read(src_wiff)
cf = C.CompoundFile(copy(data))
idx_path = "SampleSubtree/Sample$(run.sample)/Idx"
idx = C.read_stream(cf, idx_path)
for k in 1:n
    o = S.IDX_HEADER + S.IDX_RECORD * (k - 1)
    idx[o + 1:o + 4] .= reinterpret(UInt8, [htol(new_rel[k])])
    kept[k] || (idx[o + 5:o + 8] .= 0x00)
end
e = cf.entries[cf.paths[idx_path]]
length(idx) >= cf.mini_cutoff || error("Idx stream is in the mini-stream; in-place patch not implemented")
chain = C.follow_chain(cf.fat, e.start_sector, length(cf.fat))
for (j, sid) in enumerate(chain)
    lo = (j - 1) * cf.sector_size + 1; hi = min(j * cf.sector_size, length(idx))
    lo > hi && break
    fo = (Int(sid) + 1) * cf.sector_size
    data[fo + 1:fo + (hi - lo + 1)] .= view(idx, lo:hi)
end
write(out_wiff, data)

# --- verify -------------------------------------------------------------------------------------------
cf2 = C.CompoundFile(read(out_wiff))
@assert C.read_stream(cf2, idx_path) == idx "Idx stream not written as intended"
for p in C.stream_paths(cf)
    p == idx_path && continue
    @assert C.read_stream(cf2, p) == C.read_stream(cf, p) "stream $p changed"
end
fx = S.WiffRun(out_wiff)
@assert length(fx.index) == n
b1 = S.ScanBuffer(); b2 = S.ScanBuffer(); n_bad = 0
for k in 1:n
    if kept[k]
        S.read_scan!(b1, run, k); S.read_scan!(b2, fx, k)
        ok = b1.n == b2.n && b1.bin[1:b1.n] == b2.bin[1:b2.n] && b1.intensity[1:b1.n] == b2.intensity[1:b2.n] &&
             (b1.cal_a, b1.cal_b) == (b2.cal_a, b2.cal_b) && fx.index.rt_ms[k] == ix.rt_ms[k]
        global n_bad += !ok
    else
        global n_bad += !S.isempty_scan(fx, k)
    end
end
@assert n_bad == 0 "$n_bad records differ"
@printf("wrote %s (%.2f MB) and %s (%.2f MB); %d kept records verified, %d empty\n", out_wiff, filesize(out_wiff) / 1e6,
        out_scan, filesize(out_scan) / 1e6, count(kept), n - count(kept))

# --- checksums (committed: test/UnitTests/formats/sciexwiff/fixtures/) ----------------------------------
tmp = mktempdir()
conv = S.convert_run(out_wiff, tmp; params = S.ConvertParams(format = :scxs), name = "fixture", log = devnull)
scxs = conv isa NamedTuple ? conv.scxs : joinpath(tmp, "fixture.scxs")
open(joinpath(out_dir, "fixture_scxs_sha256.csv"), "w") do io
    println(io, "file,sha256")
    for fn in ("blocks.bin", "scans.arrow")
        println(io, fn, ",", bytes2hex(open(sha256, joinpath(scxs, fn))))
    end
end
println("wrote ", joinpath(out_dir, "fixture_scxs_sha256.csv"))
