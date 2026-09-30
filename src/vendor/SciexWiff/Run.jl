# A `.wiff` + `.wiff.scan` pair opened for reading.

using Mmap

struct WiffRun
    wiff_path::String
    scan_path::String
    sample::Int
    index::ScanIndex
    windows::Vector{SwathWindow}
    scan_ranges::Vector{Tuple{Float64,Float64}}   # acquired m/z range per experiment
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
    WiffRun(String(wiff), String(scan_path), sample, index, windows, ranges, scan)
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
    base_peak(bpi, bin, sb) -> (mz, intensity)

A scan's base peak: the vendor's, from its Idx record (intensity `bpi` and TOF bin `bin`, converted to m/z
with the scan's calibration), or, when the Idx lacks either (the bin is 0 on every ZT Scan DIA MS2 record),
the most intense point of the decoded profile scan in `sb`, which is how the vendor defines it. Both are in
the profile's intensity units, not the centroids' (area) units. `(NaN, 0.0)` for an empty scan.
"""
function base_peak(bpi::Real, bin::Real, sb::ScanBuffer)
    bpi > 0 && bin > 0 && return bin_to_mz(Float64(bin), sb.cal_a, sb.cal_b), Float64(bpi)
    sb.n == 0 && return NaN, 0.0
    k = argmax(view(sb.intensity, 1:sb.n))
    return mz(sb, k), Float64(sb.intensity[k])
end
base_peak(run::WiffRun, r::Integer, sb::ScanBuffer) =
    base_peak(run.index.base_peak_intensity[r], run.index.base_peak_bin[r], sb)
