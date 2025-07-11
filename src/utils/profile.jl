# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.


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