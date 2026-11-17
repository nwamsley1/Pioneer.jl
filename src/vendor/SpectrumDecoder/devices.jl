# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Optional auxiliary data: LC devices (traces, status logs) and per-scan extras. Nothing here is needed for spectra,
# and nothing here throws: every failure becomes a warning string.

"""Counted UTF-16 strings (u32 length + code units) found in [o, o + len)."""
function counted_strings(b, o::Integer, len::Integer)
    out = String[]; k = o; stop = min(o + len, length(b)) - 8
    while k < stop
        n = Int(ru32(b, k))
        if 1 <= n <= 200 && k + 4 + 2n <= length(b) && all(j -> 0x20 <= b[k+4+2j+1] < 0x7f && b[k+4+2j+2] == 0, 0:n-1)
            push!(out, String(Char.(b[k+5:2:k+4+2n]))); k += 4 + 2n
        else
            k += 2
        end
    end
    out
end

devoff(f::Acquisition, dev::Integer, k::Integer) =
    f.layout.dev_offset_width == 8 ? Int64(ru64(f.bytes, dev + k)) : Int64(ru32(f.bytes, dev + k))

"""Kind-2 traces: f64 value per sample; 72-byte index per sample with f64 time (min) at +32, value copy at +56,
byte offset at +64."""
function traces(f::Acquisition; warnings = String[])
    b = f.bytes; L = f.layout; out = []
    _, dirs = device_directory(b, L.count_back)
    for (kind, idx, off) in dirs
        kind in (2, 4) || continue
        try
            # kind 2: f64 value per sample, time in the index (+32); kind 4 (LC link device): (f64 value, f64 time)
            # per sample. Both: index record +64 = byte offset of the sample.
            n = Int(ru32(b, off + 0xc)); st = devoff(f, off, L.dev_scan_table); d = devoff(f, off, L.dev_packet_start)
            ss = kind == 4 ? 16 : 8
            t = Vector{Float64}(undef, n); v = Vector{Float64}(undef, n)
            for j in 0:n-1
                r = st + 72j
                ru32(b, r + 64) == ss * j || error("index record $j offset disagrees with the data")
                kind == 2 && rf64(b, r + 56) != rf64(b, d + 8j) && error("index record $j value disagrees with the data")
                v[j+1] = rf64(b, d + ss * j); t[j+1] = kind == 4 ? rf64(b, d + ss * j + 8) : rf64(b, r + 32)
            end
            s = counted_strings(b, off, 24576)
            j = isempty(s) ? nothing : findfirst(x -> startswith(x, s[1] * " "), s)
            push!(out, (device = Int(idx), name = isempty(s) ? "device $idx" : j === nothing ? s[1] : s[j], times = t, values = v))
        catch e
            push!(warnings, "trace device $idx: " * sprint(showerror, e))
        end
    end
    out
end

"""
    status_logs(f; warnings) -> Vector of (device, labels, times, values::Matrix{Any})

Kind-5 device logs. Each record = 4-byte prefix + the fields of a typed list found in the device's block, sized by
`extras_field_size` (e.g. an LC event log of mostly-text fields with a record-number prefix, or an LC status log with
an f32-time prefix). The list is chosen so that 4 + Σ sizes equals the record size. If no typed list fits, f32 rows
[time, field…] labelled by the strings after the "None" placeholders are tried. `times` holds the prefix as f32 (time in minutes for numeric logs, record number for the
event log).
"""
function status_logs(f::Acquisition; warnings = String[])
    b = f.bytes; L = f.layout; out = []
    _, dirs = device_directory(b, L.count_back)
    for (kind, idx, off) in dirs
        kind in (4, 5) || continue   # kind-4 LC link devices can carry a log next to their trace
        try
            n = Int(ru32(b, off + 0x10)); n == 0 && continue
            a = devoff(f, off, L.dev_region_a); z = devoff(f, off, L.dev_region_b)
            ((z - a) % n == 0) || (push!(warnings, "log device $idx: region not a whole number of records"); continue)
            rs = Int((z - a) ÷ n)
            # single-record placeholder devices: nothing to log
            n == 1 && (push!(warnings, "log device $idx: single placeholder record (no log data)"); continue)
            fits = filter(l -> any(x -> length(x[3]) > 1, l) && 4 + sum(extras_field_size(t, fm) for (t, fm, _) in l) == rs,
                          typed_lists(b, off + 0x100, min(off + 0x8000, length(b))))
            if !isempty(fits)
                sch = fits[1]; sz = [extras_field_size(t, fm) for (t, fm, _) in sch]; offs = cumsum([0; sz])
                keep = [c for c in eachindex(sch) if sz[c] > 0]
                t = [Float64(rf32(b, a + rs * j)) for j in 0:n-1]
                V = Any[extras_value(b, a + rs * j + 4 + offs[c], sch[c][1], sch[c][2]) for j in 0:n-1, c in keep]
                push!(out, (device = Int(idx), labels = [sch[c][3] for c in keep], times = t, values = V))
                continue
            end
            (rs % 4 == 0 && rs <= 400) || (push!(warnings, "log device $idx: no typed list fits $rs-byte records (skipped)"); continue)
            k = rs ÷ 4 - 1
            s = counted_strings(b, off, 24576); i = findlast(==("None"), s)
            labels = i === nothing ? String[] : s[i+1:min(end, i + k)]
            length(labels) == k || (push!(warnings, "log device $idx: $k fields but $(length(labels)) labels (skipped)"); continue)
            t = [Float64(rf32(b, a + rs * j)) for j in 0:n-1]
            issorted(t) || (push!(warnings, "log device $idx: times not increasing (skipped)"); continue)
            V = Any[Float64(rf32(b, a + rs * j + 4c)) for j in 0:n-1, c in 1:k]
            push!(out, (device = Int(idx), labels = labels, times = t, values = V))
        catch e
            push!(warnings, "log device $idx: " * sprint(showerror, e))
        end
    end
    out
end

# ---------------- per-scan extras (region D, schema in region B) ----------------
# Byte size of a schema field by type code:
# 0 section header, 1-5 one byte (5 = charge), 6-7 two bytes, 8-10 four bytes (10 = f32), 11 f64,
# 12 text of fmt bytes, 13 UTF-16 text of fmt characters.
extras_field_size(t, fm) = t == 0 ? 0 : t <= 5 ? 1 : t <= 7 ? 2 : t <= 10 ? 4 : t == 11 ? 8 : t == 12 ? fm : t == 13 ? 2fm : -1

function typed_lists(b, lo, hi)
    out = []; k = lo
    while k < hi - 16
        n = Int(ru32(b, k)); ok = 3 <= n <= 2000; x = k + 4; L = Tuple{Int,Int,String}[]
        if ok
            for _ in 1:n
                (x + 12 <= hi) || (ok = false; break)
                t = Int(ru32(b, x)); fm = Int(ru32(b, x + 4)); m = Int(ru32(b, x + 8))
                (t < 64 && 0 <= m <= 300 && x + 12 + 2m <= hi) || (ok = false; break)
                # labels are UTF-16 and may hold non-ASCII characters
                u = reinterpret(UInt16, b[x+13:x+12+2m])
                all(c -> c != 0x0000 && !(0xd800 <= c <= 0xdfff), u) || (ok = false; break)
                push!(L, (t, fm, transcode(String, collect(u)))); x += 12 + 2m
            end
        end
        ok && !any(e -> length(e[3]) > 1, L) && (ok = false)   # runs of zeros look like lists of empty labels
        if ok; push!(out, L); k = x; else k += 2; end
    end
    out
end

function extras_value(b, x, t, fm)
    t == 11 ? rf64(b, x) : t == 10 ? Float64(rf32(b, x)) : t in (8, 9) ? Float64(ru32(b, x)) :
    t <= 5 ? Float64(b[x+1]) : t <= 7 ? Float64(reinterpret(Int16, b[x+1:x+2])[1]) :
    t == 12 ? strip(String(filter(c -> 0x20 <= c < 0x7f, b[x+1:x+fm]))) :
    t == 13 ? strip(transcode(String, filter(c -> c != 0x0000 && !(0xd800 <= c <= 0xdfff), reinterpret(UInt16, b[x+1:x+2fm])))) : missing
end

"""
    scan_extras(f; warnings) -> Union{Nothing, NamedTuple{(:labels, :columns)}}

Per-scan labelled values (e.g. "Ion Injection Time (ms):", "Charge State:", "Monoisotopic M/Z:", "FAIMS CV:",
"Master Scan Number:"). The schema is the typed list in region B that contains "Charge State:"; records start at
region D, one per scan. Returns `nothing` (with a warning) if the schema or the size check fails.
"""
function scan_extras(f::Acquisition; warnings = String[])
    try
        b = f.bytes; L = f.layout; dev = f.device_block; n = length(f.scans)
        lo = devoff(f, dev, L.dev_region_b); rd = devoff(f, dev, L.dev_region_d)
        sch = filter(l -> any(x -> x[3] == "Charge State:", l), typed_lists(b, lo, f.scan_table))
        isempty(sch) && (push!(warnings, "scan extras: no schema list found"); return nothing)
        sch = sch[1]; sz = [extras_field_size(t, fm) for (t, fm, _) in sch]
        any(<(0), sz) && (push!(warnings, "scan extras: unknown field type $(unique(t for (t, _, _) in sch if extras_field_size(t, 0) < 0))"); return nothing)
        rs = sum(sz); offs = cumsum([0; sz])
        rd + n * rs <= length(b) || (push!(warnings, "scan extras: $n × $rs bytes run past the end of the file"); return nothing)
        keep = [c for c in eachindex(sch) if sch[c][1] != 0]
        cols = Dict{String,Vector}(strip(sch[c][3], [':', ' ']) => [extras_value(b, rd + (j - 1) * rs + offs[c], sch[c][1], sch[c][2]) for j in 1:n] for c in keep)
        (labels = [strip(sch[c][3], [':', ' ']) for c in keep], columns = cols)
    catch e
        push!(warnings, "scan extras: " * sprint(showerror, e)); nothing
    end
end

"""
    instrument_status(f; warnings) -> Union{Nothing, NamedTuple{(:labels, :times, :columns)}}

Region A of the mass-spec device: periodic instrument status records. Schema = the longest typed list in the device
block (after its offset group); record = f32 time (min) + fields sized by `extras_field_size`; record count = u32 at
device block +0x10 and must tile region A exactly.
"""
function instrument_status(f::Acquisition; warnings = String[])
    try
        b = f.bytes; L = f.layout; dev = f.device_block
        ra = devoff(f, dev, L.dev_region_a); rb = devoff(f, dev, L.dev_region_b); n = Int(ru32(b, dev + 0x10))
        lists = typed_lists(b, dev + 0x100, ra)
        isempty(lists) && (push!(warnings, "instrument status: no schema list"); return nothing)
        sch = first(sort(lists, by = length, rev = true)); sz = [extras_field_size(t, fm) for (t, fm, _) in sch]
        any(<(0), sz) && (push!(warnings, "instrument status: unknown field type"); return nothing)
        rs = 4 + sum(sz); offs = cumsum([0; sz])
        rb - ra == n * rs || (push!(warnings, "instrument status: region A $(rb - ra) bytes != $n × $rs"); return nothing)
        keep = [c for c in eachindex(sch) if sch[c][1] != 0]
        t = [Float64(rf32(b, ra + j * rs)) for j in 0:n-1]
        cols = Dict{String,Vector}(strip(sch[c][3], [':', ' ']) => [extras_value(b, ra + j * rs + 4 + offs[c], sch[c][1], sch[c][2]) for j in 0:n-1] for c in keep)
        (labels = [strip(sch[c][3], [':', ' ']) for c in keep], times = t, columns = cols)
    catch e
        push!(warnings, "instrument status: " * sprint(showerror, e)); nothing
    end
end

"""
    read_auxiliary(f; lc = true, extras = true, status = true) -> (traces, logs, extras, status, warnings)

Everything optional in one call; never throws.
"""
function read_auxiliary(f::Acquisition; lc::Bool = true, extras::Bool = true, status::Bool = true)
    w = String[]
    tr = lc ? (try traces(f; warnings = w) catch e; push!(w, sprint(showerror, e)); [] end) : []
    lg = lc ? (try status_logs(f; warnings = w) catch e; push!(w, sprint(showerror, e)); [] end) : []
    ex = extras ? scan_extras(f; warnings = w) : nothing
    st = status ? instrument_status(f; warnings = w) : nothing
    (traces = tr, logs = lg, extras = ex, status = st, warnings = w)
end
