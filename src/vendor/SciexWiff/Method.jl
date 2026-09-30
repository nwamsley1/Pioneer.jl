# SWATH windows from `MethodSubtree/Method1/DeviceMethod0/SWATHMethod` (notes/format.md §5):
# a 40-byte preamble (u16 window count at byte 38) then 20-byte records {f64 lo, f64 hi, u32}.
# Experiment k+1 of each cycle is window k; experiment 1 is the MS1 survey.

const SWATH_STREAM = "MethodSubtree/Method1/DeviceMethod0/SWATHMethod"

struct SwathWindow
    lo::Float64
    hi::Float64
end
center(w::SwathWindow) = (w.lo + w.hi) / 2
width(w::SwathWindow) = w.hi - w.lo

function read_swath_windows(cf::CFB.CompoundFile)
    CFB.has_stream(cf, SWATH_STREAM) || throw(ArgumentError("no SWATHMethod stream: not a SWATH acquisition"))
    s = CFB.read_stream(cf, SWATH_STREAM)
    n, rem = divrem(length(s) - 40, 20)
    (length(s) >= 40 && rem == 0) || throw(ArgumentError("SWATHMethod has $(length(s)) bytes, not 40 + 20·n"))
    Int(_rd(UInt16, s, 38)) == n || throw(ArgumentError("SWATHMethod count $(_rd(UInt16, s, 38)) != $n records"))
    wins = SwathWindow[]
    for k in 0:n-1
        lo, hi = _rd(Float64, s, 40 + 20k), _rd(Float64, s, 48 + 20k)
        (isfinite(lo) && isfinite(hi) && 0 < lo < hi) || throw(ArgumentError("SWATHMethod record $k is ($lo, $hi)"))
        push!(wins, SwathWindow(lo, hi))
    end
    wins
end

"""
    read_scan_ranges(cf, n_experiments) -> Vector{Tuple{Float64,Float64}}

The acquired m/z range of each experiment (`Experiment{e}/ExperimentHeaderEx`, f64 at bytes 66 and 74).
Matches the msConvert scan window (paper run: MS1 400–1250, MS2 100–1500).
"""
function read_scan_ranges(cf::CFB.CompoundFile, n::Integer)
    out = Tuple{Float64,Float64}[]
    for e in 0:n-1
        p = "MethodSubtree/Method1/DeviceMethod0/Period0/Experiment$e/ExperimentHeaderEx"
        CFB.has_stream(cf, p) || throw(ArgumentError("no $p"))
        s = CFB.read_stream(cf, p)
        length(s) >= 82 || throw(ArgumentError("$p has $(length(s)) bytes"))
        lo, hi = _rd(Float64, s, 66), _rd(Float64, s, 74)
        (isfinite(lo) && isfinite(hi) && 0 < lo < hi < 1e5) || throw(ArgumentError("$p scan range ($lo, $hi)"))
        push!(out, (lo, hi))
    end
    out
end

"UTF-16LE runs of printable ASCII (length ≥ 4) in a stream."
function _utf16_strings(s::Vector{UInt8})
    out = String[]; cur = Char[]
    for o in 1:2:length(s)-1
        u = UInt16(s[o]) | UInt16(s[o+1]) << 8
        if 0x20 <= u < 0x7f
            push!(cur, Char(u))
        else
            length(cur) >= 4 && push!(out, String(cur)); empty!(cur)
        end
    end
    length(cur) >= 4 && push!(out, String(cur))
    out
end

"The acquisition method's name (the file name of the method, without its folder), or \"\"."
function read_method_name(cf::CFB.CompoundFile)
    p = "MethodSubtree/Method1/AcqMethodFileInfoStm"
    CFB.has_stream(cf, p) || return ""
    strs = _utf16_strings(CFB.read_stream(cf, p))
    isempty(strs) && return ""
    name = last(split(first(strs), '\\'))
    # the stream runs the save date on after the name
    name = split(name, ',')[1]
    replace(name, r"\d?(Monday|Tuesday|Wednesday|Thursday|Friday|Saturday|Sunday).*$" => "")
end
