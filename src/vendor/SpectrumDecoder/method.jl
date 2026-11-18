# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# The acquisition method embedded in the file head: a generic compound-document container (public, general-purpose
# container format, read through its own sector table and directory).
# Its "Text" streams are readable method summaries; the "Data" streams (binary settings) are not decoded.

"""
    method_streams(f) -> Vector{(name, bytes)}; empty if no container is found. Never throws.
"""
function method_streams(f::Acquisition)
    try
        b = f.bytes; hi = min(f.packet_start, length(b))
        sig = UInt8[0xd0, 0xcf, 0x11, 0xe0, 0xa1, 0xb1, 0x1a, 0xe1]
        base = findfirst(i -> view(b, i+1:i+8) == sig, 0:hi-8); base === nothing && return Tuple{String,Vector{UInt8}}[]
        base -= 1
        u16(k) = Int(ru16(b, k)); u32(k) = Int(ru32(b, k))
        ss = 1 << u16(base + 0x1e); ms = 1 << u16(base + 0x20); nfat = u32(base + 0x2c); dir0 = u32(base + 0x30)
        cutoff = u32(base + 0x38); mfat0 = u32(base + 0x3c)
        sec(n) = base + (n + 1) * ss
        fat = vcat([[u32(sec(s) + 4i) for i in 0:ss÷4-1] for s in [u32(base + 0x4c + 4i) for i in 0:min(nfat, 109)-1]]...)
        chain(s) = (out = Int[]; while 0 <= s < 0xfffffffa && length(out) < 1_000_000; push!(out, s); s = fat[s+1]; end; out)
        readchain(s, n) = (d = vcat([b[sec(x)+1:sec(x)+ss] for x in chain(s)]...); d[1:min(n, length(d))])
        dirb = readchain(dir0, typemax(Int)); ents = []
        for e in 0:length(dirb)÷128-1
            r = e * 128; nl = Int(reinterpret(UInt16, dirb[r+0x41:r+0x42])[1]); t = dirb[r+0x43]; t == 0 && continue
            name = transcode(String, reinterpret(UInt16, dirb[r+1:r+max(nl - 2, 0)]))
            push!(ents, (name, t, Int(reinterpret(UInt32, dirb[r+0x75:r+0x78])[1]), Int(reinterpret(UInt64, dirb[r+0x79:r+0x80])[1])))
        end
        root = ents[findfirst(e -> e[2] == 5, ents)]; mini = readchain(root[3], root[4])
        mfat = mfat0 < 0xfffffffa ? reinterpret(UInt32, readchain(mfat0, typemax(Int))) : UInt32[]
        function stream(e)
            e[4] >= cutoff && return readchain(e[3], e[4])
            out = UInt8[]; s = e[3]
            while s < 0xfffffffa && length(out) < e[4]; append!(out, mini[s*ms+1:min(end, (s+1)*ms)]); s = Int(mfat[s+1]); end
            out[1:min(e[4], length(out))]
        end
        [(e[1], stream(e)) for e in ents if e[2] == 2]
    catch
        Tuple{String,Vector{UInt8}}[]
    end
end

"""
    method_text(f) -> Vector{String}

The readable method summaries ("Text" streams, UTF-16) of the embedded method; empty if unavailable.
"""
function method_text(f::Acquisition)
    out = String[]
    for (name, d) in method_streams(f)
        name == "Text" || continue
        u = reinterpret(UInt16, d[1:2*(length(d)÷2)])
        push!(out, transcode(String, filter(c -> c != 0x0000 && !(0xd800 <= c <= 0xdfff), u)))
    end
    out
end

"""
    method_windows(f) -> Vector{(centre, width)}

Isolation windows listed in the method text: "|m/z|lo-hi|" mass-list rows (Exploris/Astral style), or a table with an
"m/z" column and a width column headed "Isolation Window (m/z)" (Fusion style) or "Q1 Resolution (m/z)" (custom-window
tMS2 mass lists). Empty without a table.
"""
function method_windows(f::Acquisition)
    out = Tuple{Float64,Float64}[]
    for t in method_text(f)
        # edges may be whole numbers ("415.3-430", "429-443.7")
        for m in eachmatch(r"\|\s*([0-9]+(?:\.[0-9]+)?)-([0-9]+(?:\.[0-9]+)?)\|", t)
            lo = parse(Float64, m.captures[1]); hi = parse(Float64, m.captures[2]); push!(out, ((lo + hi) / 2, hi - lo))
        end
        lines = split(t, '\n')
        h = findfirst(l -> occursin("|", l) && (occursin("Isolation Window (m/z)", l) || occursin("Q1 Resolution (m/z)", l)), lines)
        if h !== nothing
            cols = strip.(split(lines[h], '|')); cm = findfirst(==("m/z"), cols)
            cw = something(findfirst(==("Isolation Window (m/z)"), cols), findfirst(==("Q1 Resolution (m/z)"), cols), Some(nothing))
            (cm === nothing || cw === nothing) && continue
            for l in lines[h+1:end]
                c = strip.(split(l, '|')); length(c) < max(cm, cw) && break
                x = tryparse(Float64, c[cm]); w = tryparse(Float64, c[cw]); (x === nothing || w === nothing) && break
                push!(out, (x, w))
            end
        end
    end
    unique(out)
end

"""Isolation widths stated per DIA experiment in Q Exactive-style method text ("Isolation window  45.0 m/z"), in order."""
method_experiment_widths(f::Acquisition) =
    [parse(Float64, m.captures[1]) for t in method_text(f) for m in eachmatch(r"Isolation window\s+([0-9.]+)\s*m/z", t)]

"""
Single isolation-window statement in Fusion/Exploris/Stellar-style method text ("Isolation Window = 10",
"Isolation Window (m/z) = 2"), plus a stated "Window Overlap (m/z) = o" if any; returns the expected width w + o,
or nothing if no single statement is found or the windows are custom (table-defined).
"""
function method_single_width(f::Acquisition)
    t = join(method_text(f), "\n")
    # "Isolation Window = Custom": windows come from the mass-list table; the stated number is not used
    occursin(r"Isolation Window\s*=\s*Custom", t) && return nothing
    ws = unique([parse(Float64, m.captures[1]) for m in eachmatch(r"Isolation Window(?: \(m/z\))?\s*=\s*([0-9]+\.?[0-9]*)", t)])
    length(ws) == 1 || return nothing
    ov = [parse(Float64, m.captures[1]) for m in eachmatch(r"Window Overlap \(m/z\)\s*=\s*([0-9]+\.?[0-9]*)", t)]
    ws[1] + (length(unique(ov)) == 1 ? ov[1] : 0.0)
end
