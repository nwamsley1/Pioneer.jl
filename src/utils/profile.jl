# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

### ───────────── POSIX (Linux, macOS, …) ─────────────
if Sys.isunix()

    """
        peak_rss  →  Int64

    Maximum resident‑set size (bytes) of the current process so far.
    """
    function peak_rss()::Int64
        return Sys.maxrss()
    end

### ───────────── Windows ─────────────
elseif Sys.iswindows()
    # bring in WinAPI function GetProcessMemoryInfo
    const PROCESS_MEMORY_COUNTERS = NTuple{11, UInt64}

    function peak_rss()::Int64
        pmc   = Ref(ntuple(_ -> UInt64(0), 11))
        hProc = ccall(("GetCurrentProcess", "kernel32"), Ptr{Cvoid}, ())
        ok = ccall(("GetProcessMemoryInfo", "psapi"),
                   Int32,
                   (Ptr{Cvoid}, Ref{PROCESS_MEMORY_COUNTERS}, UInt32),
                   hProc, pmc, sizeof(pmc))
        ok == 0 && error("GetProcessMemoryInfo failed")
        return Int64(pmc[][3])    # WorkingSetSize (bytes)
    end
end

"""
    process_private_bytes() -> Int

Memory this process holds that the OS cannot drop for free: private commit (`PrivateUsage`) on
Windows, anonymous resident memory (`RssAnon`) on Linux, Julia's live heap elsewhere or if those
cannot be read. File-backed pages, such as memory-mapped Arrow files, are not counted: the OS can
drop and reread them.
"""
function process_private_bytes()::Int
    if Sys.iswindows()
        # PROCESS_MEMORY_COUNTERS_EX: two UInt32 then nine size_t; PrivateUsage is the last
        pmc = Ref(ntuple(_ -> UInt64(0), 10))
        h = ccall((:GetCurrentProcess, "kernel32"), Ptr{Cvoid}, ())
        ok = ccall((:GetProcessMemoryInfo, "psapi"), Int32,
                   (Ptr{Cvoid}, Ref{NTuple{10, UInt64}}, UInt32), h, pmc, sizeof(pmc))
        ok != 0 && return Int(pmc[][10])
    elseif Sys.islinux()
        try
            for line in eachline("/proc/self/status")
                startswith(line, "RssAnon:") && return parse(Int, split(line)[2]) * 1024
            end
        catch
        end
    end
    return Int(Base.gc_live_bytes())
end

"""
    commit_headroom_bytes() -> Int

How much more memory this process can commit before Windows refuses allocations (which crashes the
search): the commit charge still available now (`ullAvailPageFile` of `GlobalMemoryStatusEx`) plus what
the page file may still grow by (`pagefile_growth_bytes`). Physical RAM is not the limit on Windows:
memory beyond it is paged out (slower, not fatal), but committing past the commit limit fails.
`typemax(Int)` elsewhere (on Linux the limit is `Sys.total_memory()`, which respects cgroup limits)
or if it cannot be read.
"""
function commit_headroom_bytes()::Int
    Sys.iswindows() || return typemax(Int)
    # MEMORYSTATUSEX: dwLength and dwMemoryLoad (two UInt32, packed in the first UInt64), then
    # ullTotalPhys, ullAvailPhys, ullTotalPageFile, ullAvailPageFile, ullTotalVirtual, ullAvailVirtual,
    # ullAvailExtendedVirtual
    st = Ref(ntuple(i -> i == 1 ? UInt64(64) : UInt64(0), 8))
    ok = ccall((:GlobalMemoryStatusEx, "kernel32"), Int32, (Ref{NTuple{8, UInt64}},), st)
    ok == 0 && return typemax(Int)
    total_phys, commit_limit, commit_avail = Int(st[][2]), Int(st[][4]), Int(st[][5])
    return commit_avail + pagefile_growth_bytes(commit_limit, total_phys)
end

"""Page files from the registry (`PagingFiles`): (path, initial MB, maximum MB); read once per process."""
const _PAGEFILES = Ref{Union{Nothing, Vector{Tuple{String, Int, Int}}}}(nothing)

function windows_pagefiles()
    _PAGEFILES[] === nothing || return _PAGEFILES[]
    entries = Tuple{String, Int, Int}[]
    try
        out = read(`reg query "HKLM\SYSTEM\CurrentControlSet\Control\Session Manager\Memory Management" /v PagingFiles`, String)
        m = match(r"PagingFiles\s+REG_MULTI_SZ\s*(.*)", out)
        if m !== nothing
            for e in split(strip(m.captures[1]), "\0"; keepempty = false)
                f = split(strip(e))
                isempty(f) && continue
                push!(entries, (String(f[1]), length(f) >= 3 ? parse(Int, f[2]) : 0,
                                length(f) >= 3 ? parse(Int, f[3]) : 0))
            end
        end
    catch
        entries = Tuple{String, Int, Int}[]    # unreadable: assume the page file cannot grow
    end
    return _PAGEFILES[] = entries
end

"""
    pagefile_growth_bytes(commit_limit, total_phys) -> Int

How much the page files may still grow: up to their configured maximum, or for a system-managed page
file (no sizes, `0 0`, or a `?:` drive) up to Windows' maximum of 3 x RAM (at least 4 GB), at most
1/8 of its volume; never more than the volume's free space less 1 GB. 0 when there is no page file or
the configuration cannot be read.
"""
function pagefile_growth_bytes(commit_limit::Integer, total_phys::Integer)
    files = windows_pagefiles()
    isempty(files) && return 0
    current = max(commit_limit - total_phys, 0)      # current page-file size, about
    max_total = 0; free = 0; seen = Set{String}()
    for (path, _, max_mb) in files
        drive = startswith(path, "?") ? get(ENV, "SystemDrive", "C:") * "\\" : first(path, 3)
        vol = try diskstat(drive) catch; nothing end
        if max_mb > 0
            max_total += max_mb * 2^20
        else
            cap = max(3 * Int(total_phys), 4 * 2^30)
            vol === nothing || (cap = min(cap, fld(Int(vol.total), 8)))
            max_total += cap
        end
        vol === nothing || drive in seen || (free += Int(vol.available); push!(seen, drive))
    end
    return max(0, min(max_total - current, free - 2^30))
end