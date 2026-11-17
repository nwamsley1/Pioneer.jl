# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

"""
    Ms1Coverage

The precursor m/z ranges each MS1 scan of a file measures. Normally one range per scan, [lowMz, highMz]. Two cases
need more: segmented MS1 (consecutive MS1 scans covering different parts of the precursor range) and multiplexed
MS1 (BoxCar / msx), where one scan holds several narrow boxes listed in its filter string,
`FTMS + p NSI Full msx ms [400.0000-416.3000, 441.2000-454.2000, ...]`. Ranges of scan `s` are
`starts[s]:starts[s + 1] - 1` in `lo` / `hi`; MS2 scans have none. `uniform` is true when every MS1 scan covers the
same ranges, i.e. no scan selection by m/z is needed.
"""
struct Ms1Coverage
    starts::Vector{Int32}
    lo::Vector{Float32}
    hi::Vector{Float32}
    uniform::Bool
end

"""
    parse_filter_ranges!(lo, hi, header) -> Int

Appends the m/z ranges of a Thermo filter string's last bracket (`[a-b, c-d, ...]`) and returns how many were found
(0 for a header without one).
"""
function parse_filter_ranges!(lo::Vector{Float32}, hi::Vector{Float32}, header::AbstractString)::Int
    a = findlast('[', header); b = findlast(']', header)
    (a === nothing || b === nothing || b < a) && return 0
    n = 0
    for part in eachsplit(SubString(header, nextind(header, a), prevind(header, b)), ',')
        d = findnext('-', part, nextind(part, firstindex(part)))   # skip a leading sign
        d === nothing && continue
        l = tryparse(Float32, strip(SubString(part, 1, prevind(part, d))))
        h = tryparse(Float32, strip(SubString(part, nextind(part, d))))
        (l === nothing || h === nothing || h <= l) && continue
        push!(lo, l); push!(hi, h); n += 1
    end
    return n
end

"""
    Ms1Coverage(spectra)

Per-scan MS1 coverage. A multiplexed scan (its filter string lists more than one range) covers its boxes; any other
MS1 scan covers [lowMz, highMz], or everything when that range is not valid (high <= low).
"""
function Ms1Coverage(spectra)
    n_scans = length(spectra)
    starts = Vector{Int32}(undef, n_scans + 1)
    lo = Float32[]; hi = Float32[]
    blo = Float32[]; bhi = Float32[]
    for s in 1:n_scans
        starts[s] = Int32(length(lo) + 1)
        getMsOrder(spectra, s) == 1 || continue
        empty!(blo); empty!(bhi)
        if parse_filter_ranges!(blo, bhi, getScanHeader(spectra, s)) > 1
            append!(lo, blo); append!(hi, bhi)
            continue
        end
        l = Float32(getLowMz(spectra, s)); h = Float32(getHighMz(spectra, s))
        h > l || (l = -Inf32; h = Inf32)
        push!(lo, l); push!(hi, h)
    end
    starts[n_scans + 1] = Int32(length(lo) + 1)
    uniform = true
    first_ms1 = findfirst(s -> starts[s + 1] > starts[s], 1:n_scans)
    if first_ms1 !== nothing
        r0 = starts[first_ms1]:(starts[first_ms1 + 1] - 1)
        for s in (first_ms1 + 1):n_scans
            r = starts[s]:(starts[s + 1] - 1)
            isempty(r) && continue
            if length(r) != length(r0) || any(k -> lo[r[k]] != lo[r0[k]] || hi[r[k]] != hi[r0[k]], eachindex(r))
                uniform = false
                break
            end
        end
    end
    return Ms1Coverage(starts, lo, hi, uniform)
end

"""Whether MS1 scan `s` measures m/z `mz`."""
@inline function covers(c::Ms1Coverage, s::Integer, mz::Float32)::Bool
    @inbounds for k in c.starts[s]:(c.starts[s + 1] - 1)
        c.lo[k] <= mz <= c.hi[k] && return true
    end
    return false
end
