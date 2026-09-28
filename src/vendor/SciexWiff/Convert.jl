# Conversion driver: `.wiff` -> Pioneer Arrow and/or `.scxs`, scans processed in batches on all threads,
# written in acquisition order.

using Printf

"""
    ConvertParams(; centroid = CentroidParams(), bin_scale = 4, int_scale = 10.0, zstd_level = 3,
                    format = :both, batch_scans = 4096)

- `bin_scale`: stored bin = round(bin_scale · centroid position in raw TOF units). One raw unit is ≈ 1.2 ppm
  at m/z 400, so the default 4 rounds positions to ≈ 0.3 ppm.
- `int_scale`: stored intensity = round(int_scale · intensity); centroids that round to 0 are dropped.
- `format`: `:arrow`, `:scxs` or `:both`. Both outputs hold exactly the same (quantised) peaks.
"""
Base.@kwdef struct ConvertParams
    centroid::CentroidParams = CentroidParams()
    bin_scale::Int = 4
    int_scale::Float64 = 10.0
    zstd_level::Int = 3
    format::Symbol = :both
    batch_scans::Int = 4096
end

"Everything one worker task owns."
mutable struct Worker
    sb::ScanBuffer
    cb::CentroidBuffer
    codec::BlockCodec
    qbin::Vector{UInt32}; qint::Vector{UInt32}
    mz::Vector{Float32}; intensity::Vector{Float32}
    t_decode::Float64; t_centroid::Float64; t_encode::Float64
end
Worker() = Worker(ScanBuffer(), CentroidBuffer(), BlockCodec(), UInt32[], UInt32[], Float32[], Float32[], 0.0, 0.0, 0.0)

"What one scan yields for the writers."
struct ScanResult
    row::NamedTuple
    n::Int
    mz::Vector{Float32}
    intensity::Vector{Float32}
    zbytes::Vector{UInt8}
end

"Round centroids to stored integers (merging equal bins, dropping zero intensity); returns the count."
function quantize!(wk::Worker, p::ConvertParams)
    cb = wk.cb; n = 0
    length(wk.qbin) < cb.n && (resize!(wk.qbin, cb.n); resize!(wk.qint, cb.n))
    @inbounds for k in 1:cb.n
        q = round(UInt32, p.bin_scale * cb.pos[k])
        v = round(UInt64, p.int_scale * cb.intensity[k])
        if n > 0 && wk.qbin[n] == q
            wk.qint[n] = UInt32(min(UInt64(wk.qint[n]) + v, typemax(UInt32)))
        elseif v > 0
            n += 1; wk.qbin[n] = q; wk.qint[n] = UInt32(min(v, typemax(UInt32)))
        end
    end
    n
end

function process_scan!(wk::Worker, run::WiffRun, r::Int, p::ConvertParams, want_arrow::Bool, want_scxs::Bool)
    t0 = time_ns()
    read_scan!(wk.sb, run, r)
    t1 = time_ns()
    centroid!(wk.cb, wk.sb, p.centroid)
    t2 = time_ns()
    n = quantize!(wk, p)
    a, b = wk.sb.cal_a, wk.sb.cal_b
    mz = Float32[]; it = Float32[]
    if want_arrow
        mz = Vector{Float32}(undef, n); it = Vector{Float32}(undef, n)
        @inbounds for k in 1:n
            mz[k] = Float32(bin_to_mz(Float64(wk.qbin[k]) / p.bin_scale, a, b))
            it[k] = Float32(wk.qint[k] / p.int_scale)
        end
    end
    z = want_scxs ? wk.codec.zbuf[1:encode_block!(wk.codec, wk.qbin, wk.qint, n, p.zstd_level)] : UInt8[]
    t3 = time_ns()
    wk.t_decode += (t1 - t0) / 1e9; wk.t_centroid += (t2 - t1) / 1e9; wk.t_encode += (t3 - t2) / 1e9
    cyc, ex = cycle_experiment(run, r)
    w = window(run, r); lo, hi = scan_range(run, r)
    row = (record = Int32(r), cycle = Int32(cyc), experiment = Int16(ex), ms_order = UInt8(ex == 1 ? 1 : 2),
           retention_time = Float32(retention_time_min(run, r)), low_mz = Float32(lo), high_mz = Float32(hi),
           tic = Float32(run.index.tic[r]),
           center_mz = w === nothing ? NaN32 : Float32(center(w)), isolation_width = w === nothing ? NaN32 : Float32(width(w)),
           cal_a = a, cal_b = b, n_peaks = Int32(n))
    ScanResult(row, n, mz, it, z)
end

"""
    convert_run(wiff, out_dir; params = ConvertParams(), name = <file stem>, sample = 1, log = stdout)
        -> (arrow = path | nothing, scxs = path | nothing, timings)
"""
function convert_run(wiff::AbstractString, out_dir::AbstractString; params::ConvertParams = ConvertParams(),
                 name::Union{Nothing, AbstractString} = nothing, sample::Integer = 1, log::IO = stdout)
    p = params
    p.format in (:arrow, :scxs, :both) || throw(ArgumentError("format must be :arrow, :scxs or :both"))
    t_start = time()
    run = WiffRun(wiff; sample = sample)
    t_open = time() - t_start
    name === nothing && (name = splitext(basename(wiff))[1])
    records = [r for r in 1:length(run) if run.index.block_size[r] > 0]
    mkpath(out_dir)
    want_arrow = p.format in (:arrow, :both); want_scxs = p.format in (:scxs, :both)
    arrow_path = want_arrow ? joinpath(out_dir, name * ".arrow") : nothing
    scxs_path = want_scxs ? joinpath(out_dir, name * ".scxs") : nothing
    cp = p.centroid
    meta = Dict{String, Any}(
        "source" => basename(wiff), "source_scan" => basename(run.scan_path), "sample" => sample,
        "n_records" => length(run), "n_windows" => length(run.windows),
        "windows" => [[w.lo, w.hi] for w in run.windows], "scan_ranges" => [[lo, hi] for (lo, hi) in run.scan_ranges],
        "centroid" => Dict(string(f) => string(getfield(cp, f)) for f in fieldnames(CentroidParams)),
        "bin_scale" => p.bin_scale, "int_scale" => p.int_scale, "zstd_level" => p.zstd_level,
        "mz_formula" => "(cal_a * (stored_bin / bin_scale / 40 - cal_b))^2",
        "converter" => "Pioneer.jl $(pkgversion(parentmodule(@__MODULE__))) (SciexWiff, from e4d9097)", "converted_at" => Libc.strftime("%Y-%m-%dT%H:%M:%SZ", time()))
    aw = want_arrow ? PioneerArrowWriter(arrow_path; metadata = Dict("source" => basename(wiff), "converter" => meta["converter"],
        "centroid" => string(cp), "bin_scale" => string(p.bin_scale), "int_scale" => string(p.int_scale))) : nothing
    sw = want_scxs ? ScxsWriter(scxs_path, meta) : nothing

    nt = Threads.nthreads()
    workers = [Worker() for _ in 1:nt]
    # double-buffered: the main task writes batch i while the workers compute batch i+1
    bufs = (Vector{ScanResult}(undef, p.batch_scans), Vector{ScanResult}(undef, p.batch_scans))
    t_write_arrow = 0.0; t_write_scxs = 0.0; t_compute = 0.0; npk = 0
    function write_batch(results, m)
        ta = 0.0; ts = 0.0; np = 0
        for j in 1:m
            res = results[j]; np += res.n
            if want_arrow
                t0 = time(); write_scan!(aw, res.row, res.mz, res.intensity, res.n); ta += time() - t0
            end
            if want_scxs
                t0 = time(); write_scan!(sw, res.row, res.zbytes); ts += time() - t0
            end
        end
        (ta, ts, np)
    end
    writer = nothing
    for (bi, bstart) in enumerate(1:p.batch_scans:length(records))
        batch = bstart:min(bstart + p.batch_scans - 1, length(records))
        results = bufs[isodd(bi) ? 1 : 2]
        tc = time()
        next = Threads.Atomic{Int}(1)
        @sync for t in 1:nt
            Threads.@spawn begin
                wk = workers[t]
                while true
                    j = Threads.atomic_add!(next, 1)
                    j > length(batch) && break
                    results[j] = process_scan!(wk, run, records[batch[j]], p, want_arrow, want_scxs)
                end
            end
        end
        t_compute += time() - tc
        if writer !== nothing                      # the previous batch must be written before this one
            ta, ts, np = fetch(writer); t_write_arrow += ta; t_write_scxs += ts; npk += np
        end
        writer = Threads.@spawn write_batch(results, length(batch))
    end
    if writer !== nothing
        ta, ts, np = fetch(writer); t_write_arrow += ta; t_write_scxs += ts; npk += np
    end
    tc = time(); want_arrow && close(aw); t_write_arrow += time() - tc
    tc = time(); want_scxs && close(sw); t_write_scxs += time() - tc
    total = time() - t_start
    timings = (open = t_open, compute_wall = t_compute, write_arrow = t_write_arrow, write_scxs = t_write_scxs, total = total,
               cpu_decode = sum(w.t_decode for w in workers), cpu_centroid = sum(w.t_centroid for w in workers),
               cpu_encode = sum(w.t_encode for w in workers), threads = nt, n_scans = length(records), n_peaks = npk)
    @printf(log, "%s: %d scans, %d centroids | open %.2f s, compute %.2f s (CPU decode %.2f, centroid %.2f, quantize+encode %.2f), write arrow %.2f s, scxs %.2f s | total %.2f s on %d threads\n",
        basename(wiff), length(records), npk, t_open, t_compute, timings.cpu_decode, timings.cpu_centroid, timings.cpu_encode,
        t_write_arrow, t_write_scxs, total, nt)
    want_arrow && @printf(log, "  wrote %s (%.1f MB)\n", arrow_path, filesize(arrow_path) / 1e6)
    want_scxs && @printf(log, "  wrote %s (%.1f MB)\n", scxs_path, sum(filesize(joinpath(scxs_path, f)) for f in readdir(scxs_path)) / 1e6)
    (arrow = arrow_path, scxs = scxs_path, timings = timings)
end
