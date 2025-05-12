
### ───────────── POSIX (Linux, macOS, …) ─────────────
if Sys.isunix()
    const RUSAGE_SELF = 0    # from sys/resource.h

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
        pmc   = zero(PROCESS_MEMORY_COUNTERS)
        hProc = ccall(("GetCurrentProcess", "kernel32"), Ptr{Cvoid}, ())
        ok = ccall(("GetProcessMemoryInfo", "psapi"),
                   Int32,
                   (Ptr{Cvoid}, Ref{PROCESS_MEMORY_COUNTERS}, UInt32),
                   hProc, pmc, sizeof(pmc))
        ok == 0 && error("GetProcessMemoryInfo failed")
        return Int64(pmc[3])    # WorkingSetSize (bytes)
    end
end