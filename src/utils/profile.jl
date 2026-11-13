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
Windows, anonymous resident memory (`RssAnon`) on Linux, the physical footprint (`phys_footprint`, what
Activity Monitor shows as "Memory", including compressed pages) on macOS, Julia's live heap elsewhere or if
those cannot be read. File-backed pages, such as memory-mapped Arrow files, are not counted: the OS can
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
    elseif Sys.isapple()
        # task_info(mach_task_self(), TASK_VM_INFO = 22) fills task_vm_info_data_t; phys_footprint is at byte
        # offset 144 (the 19th UInt64). The count is in 4-byte units; 38 of them reach phys_footprint.
        try
            buf = zeros(UInt64, 128)
            cnt = Ref{UInt32}(UInt32(sizeof(buf) ÷ 4))
            task = unsafe_load(cglobal(:mach_task_self_, UInt32))
            kr = ccall(:task_info, Cint, (UInt32, Cint, Ptr{UInt64}, Ref{UInt32}), task, 22, buf, cnt)
            kr == 0 && cnt[] >= 38 && return Int(buf[19])
        catch
        end
    end
    return Int(Base.gc_live_bytes())
end

"""
    memory_limit_bytes() -> Int

Memory this process may use: the smallest of `Sys.total_memory()`, every cgroup memory limit from this
process's cgroup up to the root (`cgroup_memory_limit_bytes`), and the Slurm job's allocation
(`slurm_memory_limit_bytes`). `Sys.total_memory()` alone is not enough on a cluster: it reads only the
process's own cgroup, but Slurm sets the job's limit on a parent cgroup (`job_<id>`), and inside a container
with its own cgroup namespace (enroot, Apptainer, Docker) those parents are not visible at all.
"""
memory_limit_bytes()::Int = memory_limit_info().bytes

"""
    memory_limit_info() -> (bytes, source, details)

`memory_limit_bytes()` with where it came from, for the log: `source` names the smallest limit, `details` lists
every limit found (physical RAM, Julia's `Sys.total_memory()` when a container limit lowers it, cgroup limit,
Slurm allocation).
"""
function memory_limit_info()
    phys = Int(Sys.total_physical_memory())
    found = Tuple{Int, String}[(phys, "physical RAM")]
    jl = Int(Sys.total_memory())
    jl < phys && push!(found, (jl, "container limit"))
    if Sys.islinux()
        cg = cgroup_memory_limit_bytes()
        cg < typemax(Int) && push!(found, (cg, "cgroup limit"))
        sl = slurm_memory_limit_bytes()
        sl < typemax(Int) && push!(found, (sl, "Slurm allocation"))
    end
    bytes, source = found[argmin(first.(found))]
    details = join(["$name $(round(b / 1e9; digits = 1)) GB" for (b, name) in found], ", ")
    return (; bytes, source, details)
end

"""
    set_gc_memory_target!(bytes) -> (before, after)

Tells Julia's garbage collector to keep its heap under `bytes` (what `--heap-size-hint` sets), collecting more
often as the heap nears it. Julia sets no target by default, so the GC lets freed memory pile up regardless of the
machine's or the job's limit; with many threads allocating at once this added up to 1.4 x the PSM table during
Main Search scoring. Only lowers the target (a smaller `--heap-size-hint` is kept). Returns the targets before and
after, or `nothing` if this Julia has no such setting.
"""
function set_gc_memory_target!(bytes::Integer)
    before = try
        ccall(:jl_gc_get_max_memory, UInt64, ())
    catch
        return nothing
    end
    bytes > 0 && UInt64(bytes) < before && ccall(:jl_gc_set_max_memory, Cvoid, (UInt64,), UInt64(bytes))
    return (Int(before), Int(ccall(:jl_gc_get_max_memory, UInt64, ())))
end

"""
    in_cgroup_namespace() -> Bool

Whether this process sees itself at the root of the cgroup tree (`/proc/self/cgroup` is `0::/`), as inside a
container with its own cgroup namespace (enroot, Docker, Podman). Limits set on parents outside the container, such
as a batch job's, are then invisible.
"""
function in_cgroup_namespace()
    Sys.islinux() || return false
    try
        lines = filter(!isempty, readlines("/proc/self/cgroup"))
        return !isempty(lines) && all(l -> endswith(l, ":/"), lines)
    catch
        return false
    end
end

"""How `process_private_bytes()` measures memory on this platform, for the log."""
process_memory_source() = Sys.iswindows() ? "private commit (PrivateUsage)" :
                          Sys.islinux() ? "anonymous resident memory (RssAnon)" :
                          Sys.isapple() ? "physical footprint (phys_footprint)" :
                          "Julia's live heap (no OS reading on this platform; may undercount)"

"""
    batch_job() -> String or nothing

The batch scheduler job this process runs in (Slurm, LSF, PBS/Torque, Grid Engine), from its environment
variables, or `nothing`.
"""
function batch_job()
    haskey(ENV, "SLURM_JOB_ID") && return "Slurm job $(ENV["SLURM_JOB_ID"])"
    haskey(ENV, "LSB_JOBID") && return "LSF job $(ENV["LSB_JOBID"])"
    haskey(ENV, "PBS_JOBID") && return "PBS job $(ENV["PBS_JOBID"])"
    haskey(ENV, "JOB_ID") && haskey(ENV, "SGE_ROOT") && return "Grid Engine job $(ENV["JOB_ID"])"
    return nothing
end

"""
    cgroup_memory_limit_bytes() -> Int

Smallest memory limit on this process's cgroup or any parent of it: `memory.max` (cgroup v2) or
`memory.limit_in_bytes` (cgroup v1). `typemax(Int)` if there is none or it cannot be read.
"""
function cgroup_memory_limit_bytes()::Int
    lim = typemax(Int)
    try
        for line in eachline("/proc/self/cgroup")
            parts = split(line, ':'; limit = 3)
            length(parts) == 3 || continue
            hier, ctrls, path = parts
            if hier == "0" && isempty(ctrls)                 # cgroup v2 (unified)
                root, file = "/sys/fs/cgroup", "memory.max"
            elseif "memory" in split(ctrls, ',')             # cgroup v1 memory controller
                root, file = "/sys/fs/cgroup/memory", "memory.limit_in_bytes"
            else
                continue
            end
            d = rstrip(path, '/')                            # "" is the (namespace) root
            while true
                f = root * d * "/" * file
                if isfile(f)
                    v = tryparse(Int, strip(read(f, String)))   # "max" (no limit) does not parse
                    v !== nothing && v > 0 && (lim = min(lim, v))
                end
                isempty(d) && break
                d = dirname(d); d == "/" && (d = "")
            end
        end
    catch
    end
    return lim
end

"""
    slurm_memory_limit_bytes() -> Int

The Slurm job's memory allocation on this node, from `SLURM_MEM_PER_NODE`, or `SLURM_MEM_PER_CPU` times
`SLURM_CPUS_ON_NODE` (both in MiB). This is the only place a containerized job can see its limit.
`typemax(Int)` outside Slurm or when the job has the whole node's memory (`--mem=0`).
"""
function slurm_memory_limit_bytes()::Int
    mb = tryparse(Int, get(ENV, "SLURM_MEM_PER_NODE", ""))
    mb !== nothing && mb > 0 && return mb * 2^20
    per = tryparse(Int, get(ENV, "SLURM_MEM_PER_CPU", ""))
    n = tryparse(Int, get(ENV, "SLURM_CPUS_ON_NODE", ""))
    per !== nothing && n !== nothing && per > 0 && n > 0 && return per * n * 2^20
    return typemax(Int)
end

"""
    commit_headroom_bytes() -> Int

How much more memory this process can commit before Windows refuses allocations (which crashes the
search): the commit charge still available now (`ullAvailPageFile` of `GlobalMemoryStatusEx`) plus what
the page file may still grow by (`pagefile_growth_bytes`). Physical RAM is not the limit on Windows:
memory beyond it is paged out (slower, not fatal), but committing past the commit limit fails.
`typemax(Int)` elsewhere (on Linux the limit is `memory_limit_bytes()`: cgroup and Slurm limits)
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