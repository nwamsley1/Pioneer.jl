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
end

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
    length(index) % (length(windows) + 1) == 0 ||
        throw(ArgumentError("$(length(index)) Idx records is not a whole number of $(length(windows) + 1)-experiment cycles"))
    ranges = read_scan_ranges(cf, length(windows) + 1)
    scan = Mmap.mmap(scan_path)
    WiffRun(String(wiff), String(scan_path), sample, index, windows, ranges, read_method_name(cf), scan)
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
    decode_block!(buf, run.scan, ix.meta_start[r], ix.block_offset[r], ix.block_size[r])
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
