# A `.wiff` + `.wiff.scan` pair opened for reading.

using Mmap

struct WiffRun
    wiff_path::String
    scan_path::String
    sample::Int
    index::ScanIndex
    windows::Vector{SwathWindow}
    scan_ranges::Vector{Tuple{Float64,Float64}}   # acquired m/z range per experiment
    method_name::String
    scan::Vector{UInt8}          # memory-mapped `.wiff.scan`
    calibration::Vector{NTuple{2,Float64}}   # per-record (a, b) from TOFCalibrationData; used when a block has none
    bin_step::Int                            # bins per delta unit (8 ZenoTOF; 4 TripleTOF, detected)
    instrument_model::String                 # from the run log ("TripleTOF 5600"/"TripleTOF 6600"), or ""
end

# TripleTOF 5600/5600+ TOFCalibrationData reads every m/z low by a constant factor 0.99068 (y1 K/R ions on three
# instruments: PXD028735, PXD001587, PXD010912); m/z ∝ a², so scale a by 1/sqrt(0.99068). The 6600 needs none.
const TTOF5600_A_SCALE = 1 / sqrt(0.99068)

"""
    find_scan_file(wiff_path) -> path

`X.wiff.scan` next to `X.wiff`; also accepts `X.scan` (some repositories rename it).
"""
function find_scan_file(wiff::AbstractString)
    for p in (wiff * ".scan", splitext(wiff)[1] * ".scan")
        isfile(p) && return p
    end
    throw(ArgumentError("no .wiff.scan (or .scan) next to $wiff"))
end

"Sample numbers present in the container (`SampleSubtree/Sample{N}/Idx`), sorted numerically."
function sample_numbers(cf::CFB.CompoundFile)
    ns = Int[]
    for p in CFB.stream_paths(cf)
        m = match(r"^SampleSubtree/Sample(\d+)/Idx$", p)
        m === nothing || push!(ns, parse(Int, m[1]))
    end
    sort!(ns)
end

function WiffRun(wiff::AbstractString; sample::Integer = 1, scan_path::AbstractString = find_scan_file(wiff))
    endswith(lowercase(wiff), ".wiff2") &&
        throw(ArgumentError(".wiff2 files are encrypted and unsupported; convert with msConvert on Windows"))
    cf = CFB.CompoundFile(wiff)
    sample in sample_numbers(cf) || throw(ArgumentError("sample $sample not in $wiff (has $(sample_numbers(cf)))"))
    index = ScanIndex(CFB.read_stream(cf, "SampleSubtree/Sample$sample/Idx"))
    windows = read_swath_windows(cf)
    # A run stopped mid-cycle (seen on TripleTOF 5600/6600) ends with a partial cycle; records are
    # cycle-aligned from the first, so drop the incomplete last cycle rather than mislabel its scans.
    full, partial = divrem(length(index), length(windows) + 1)
    full > 0 || throw(ArgumentError("$(length(index)) Idx records is less than one $(length(windows) + 1)-experiment cycle"))
    calibration = read_tof_calibration(cf, sample)
    model = read_instrument_model(cf, sample)
    if !isempty(calibration)
        is5600 = model == "TripleTOF 5600"
        # 5600 tables keep b near 0, 6600 tables near -14: a cross-check on the logged model
        is5600 == (_median(last.(calibration)) > -5) || @warn "$(basename(wiff)): run log names \"$model\" but its " *
            "TOFCalibrationData b values look like the other TripleTOF generation; check the m/z calibration"
        is5600 && (calibration = [(a * TTOF5600_A_SCALE, b) for (a, b) in calibration])
    end
    if partial > 0
        @warn "$(basename(wiff)): dropping $partial records of an incomplete final cycle ($(length(windows) + 1) experiments per cycle)"
        keep = 1:(full * (length(windows) + 1))
        index = ScanIndex(index.block_offset[keep], index.block_size[keep], index.meta_start[keep], index.rt_ms[keep],
                          index.tic[keep], index.base_peak_intensity[keep], index.base_peak_bin[keep])
        length(calibration) >= length(keep) && (calibration = calibration[keep])
    end
    ranges = read_scan_ranges(cf, length(windows) + 1)
    scan = Mmap.mmap(scan_path)
    step = detect_bin_step(scan, index)
    WiffRun(String(wiff), String(scan_path), sample, index, windows, ranges, read_method_name(cf), scan, calibration, step, model)
end

Base.length(run::WiffRun) = length(run.index)
experiments_per_cycle(run::WiffRun) = length(run.windows) + 1
"1-based (cycle, experiment) of Idx record `r` (1-based)."
cycle_experiment(run::WiffRun, r::Integer) = divrem(r - 1, experiments_per_cycle(run)) .+ 1
ms_order(run::WiffRun, r::Integer) = cycle_experiment(run, r)[2] == 1 ? 1 : 2
"The SWATH window of record `r`, or `nothing` for MS1."
function window(run::WiffRun, r::Integer)
    e = cycle_experiment(run, r)[2]
    e == 1 ? nothing : run.windows[e-1]
end
scan_range(run::WiffRun, r::Integer) = run.scan_ranges[cycle_experiment(run, r)[2]]
isempty_scan(run::WiffRun, r::Integer) = run.index.block_size[r] == 0
retention_time_min(run::WiffRun, r::Integer) = run.index.rt_ms[r] / 60_000

"Decode record `r` (1-based) into `buf`; an empty record gives `buf.n == 0`."
function read_scan!(buf::ScanBuffer, run::WiffRun, r::Integer)
    ix = run.index
    if ix.block_size[r] == 0
        buf.n = 0
        return buf
    end
    if ix.meta_start[r] == ix.block_offset[r]
        # no per-block metadata (TripleTOF): calibration comes from the .wiff's TOFCalibrationData table
        r <= length(run.calibration) || throw(ScanFormatError("record $r has no block metadata and no TOFCalibrationData entry"))
        buf.cal_a, buf.cal_b = run.calibration[r]
    end
    decode_block!(buf, run.scan, ix.meta_start[r], ix.block_offset[r], ix.block_size[r], run.bin_step)
end

"""
    detect_bin_step(scan, index) -> Int

Bins per delta unit. ZenoTOF blocks carry metadata and use 8. Blocks without metadata (TripleTOF) are checked against the
Idx: the step is the one (4 or 8) whose most intense decoded peak lands within 64 bins of the record's base-peak bin
in clearly more of 50 sampled records.
"""
function detect_bin_step(scan::AbstractVector{UInt8}, index::ScanIndex)
    recs = filter(r -> index.block_size[r] > 0, 1:min(length(index), 2000))
    isempty(recs) && return 8
    index.meta_start[recs[1]] < index.block_offset[recs[1]] && return 8
    sample = recs[round.(Int, range(1, length(recs), length = min(50, length(recs))))]
    buf = ScanBuffer(); hits = Dict(4 => 0, 8 => 0)
    for step in (4, 8), r in sample
        decode_block!(buf, scan, index.meta_start[r], index.block_offset[r], index.block_size[r], step)
        buf.n > 0 && abs(Float64(buf.bin[argmax(view(buf.intensity, 1:buf.n))]) - index.base_peak_bin[r]) <= 64 && (hits[step] += 1)
    end
    best = hits[4] >= hits[8] ? 4 : 8
    # the Idx base peak is the vendor's centroided value, so the raw maximum often misses it (16/50 on a 5600+);
    # the wrong step essentially never matches, so require a clear winner rather than a majority
    hits[best] >= max(5, 4 * hits[12 - best]) ||
        throw(ScanFormatError("cannot determine the bin step: base-peak matches $(hits) of $(length(sample)) records"))
    best
end

"""
    zt_candidate_windows(windows) -> Bool

The window table ZT Scan DIA writes: at least 100 contiguous windows (its Q1 reporting bins). A sanity check on a
run the user declares ZT, not a detector: stepped SWATH with contiguous narrow windows passes it too (nSWATH:
173 × 2.9 Da). The `.wiff` does not record the scan mode (notes/zt_scan_format.md), so the user says which it is.
"""
function zt_candidate_windows(w::AbstractVector{SwathWindow})
    length(w) >= 100 || return false
    all(k -> abs(w[k+1].lo - w[k].hi) < 1e-3, 1:length(w)-1)
end

_median(v) = (s = sort(v); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)

"""
    acquisition_metadata(run; zt_scan = false) -> Dict{String,String}

Facts about the acquisition written into the converter's outputs (Arrow schema metadata, `.scxs` meta.json):

- `acquisition_type`: `zt_scan_dia` or `swath`, as the user declared it (`zt_scan`); `acquisition_type_source`:
  `user`. The `.wiff` does not record the scan mode, so it is never guessed.
- `acquisition_method`: the method name.
- ZT only: `q1_bin_width_mz`, `q1_bin_step_mz` (mean over the method's Q1 bins), `q1_bins_per_cycle`,
  `q1_first_bin_center_mz`, `q1_last_bin_center_mz`, `q1_bin_dwell_ms` (median MS2 spacing within a cycle),
  `q1_scan_rate_mz_per_s`. The quad window width is not recorded: the `.wiff` does not store it, and method names
  are not a reliable source (`…_5-3DaQ1`).

Throws `ArgumentError` when `zt_scan` is set on a run whose window table cannot be ZT (`zt_candidate_windows`).
"""
function acquisition_metadata(run::WiffRun; zt_scan::Bool = false)
    m = Dict{String,String}("acquisition_method" => run.method_name, "acquisition_type_source" => "user")
    if !zt_scan
        m["acquisition_type"] = "swath"
        return m
    end
    zt_candidate_windows(run.windows) || throw(ArgumentError(
        "$(basename(run.wiff_path)) was declared ZT Scan DIA, but its method has $(length(run.windows)) windows " *
        "that are not a contiguous Q1 bin table (ZT has hundreds). Convert it without ZT."))
    w = run.windows; E = experiments_per_cycle(run)
    step = (w[end].lo - w[1].lo) / (length(w) - 1)      # bins alternate slightly in width: average over the ramp
    dwell = Float64[]
    for c in 1:min(20, length(run) ÷ E), e in 3:E
        r = (c - 1) * E + e
        push!(dwell, run.index.rt_ms[r] - run.index.rt_ms[r-1])
    end
    dw = _median(dwell)
    m["acquisition_type"] = "zt_scan_dia"
    m["q1_bin_width_mz"] = string(sum(x -> x.hi - x.lo, w) / length(w))
    m["q1_bin_step_mz"] = string(step)
    m["q1_bins_per_cycle"] = string(length(w))
    m["q1_first_bin_center_mz"] = string(center(w[1])); m["q1_last_bin_center_mz"] = string(center(w[end]))
    m["q1_bin_dwell_ms"] = string(dw)
    m["q1_scan_rate_mz_per_s"] = string(step / (dw / 1000))
    m
end
