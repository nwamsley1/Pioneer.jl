# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Safe file operations for cross-platform compatibility

const WINDOWS_DELETE_MAX_ATTEMPTS = 3
const WINDOWS_DELETE_SETTLE_SECONDS = 2.0
const _WINDOWS_DELETE_GC_LOCK = ReentrantLock()

const _WINDOWS_FILE_ATTRIBUTE_NORMAL = UInt32(0x80)
const _WINDOWS_ERROR_FILE_NOT_FOUND = UInt32(2)
const _WINDOWS_ERROR_PATH_NOT_FOUND = UInt32(3)

"""
    _windows_path(fpath)

Absolute backslash form of `fpath` for Win32 calls. Paths of 260 characters or more
get the `\\\\?\\` prefix (`\\\\?\\UNC\\` for network shares), which lifts the MAX_PATH limit.
"""
function _windows_path(fpath::AbstractString)
    win_path = replace(abspath(normpath(String(fpath))), "/" => "\\")
    (length(win_path) < 260 || startswith(win_path, "\\\\?\\")) && return win_path
    startswith(win_path, "\\\\") && return "\\\\?\\UNC\\" * win_path[3:end]
    return "\\\\?\\" * win_path
end

"""
    _windows_delete_file(win_path) -> UInt32

Delete one file with kernel32 `DeleteFileW`, after clearing a read-only attribute
(what `del /f` does), and return 0 or the Windows error code. Like `cmd.exe del`, it
removes a file this process still has memory-mapped, which Julia's `rm` refuses with
EACCES, but it needs no child process (about 1 ms per file instead of about 50 ms).
"""
# Defined on Windows only. Its callers run only there, but a system image compiles
# every concrete method it holds, and LLVM has no `stdcall` on arm64: building the
# app on Apple Silicon aborted with "LLVM ERROR: Unsupported calling convention".
@static if Sys.iswindows()
    function _windows_delete_file(win_path::String)
        ccall((:SetFileAttributesW, "kernel32"), stdcall, Cint, (Cwstring, UInt32),
              win_path, _WINDOWS_FILE_ATTRIBUTE_NORMAL)
        ccall((:DeleteFileW, "kernel32"), stdcall, Cint, (Cwstring,), win_path) != 0 &&
            return UInt32(0)
        return UInt32(Libc.GetLastError())
    end
end

"""
    safeRm(fpath; force=false)

Safely remove a file with Windows-specific handling for file locks and permissions.

# Arguments
- `fpath`: Path to file to remove
- `force`: Force removal on Unix. Windows preserves its historical forced-delete behavior.

# Implementation
- On Windows: Normalize once and delete with `DeleteFileW`, waiting for a delete-pending
  name to clear. If that fails (a still-mapped file on a network share), release dropped
  mappings (`_release_dropped_mappings`) and retry. If the file is still held it warns and
  returns rather than throwing, as the earlier `cmd.exe del` path effectively did.
- On Unix: Standard rm() call

This function handles common Windows file locking issues that occur with Arrow files
and other binary formats that may have lingering memory mappings. Callers must let
their Arrow or IO references leave scope before invoking this function; rebinding an
argument inside `safeRm` cannot release a reference held by the caller.
"""
function safeRm(fpath::AbstractString; force::Bool=false)
    path = abspath(normpath(String(fpath)))
    (isfile(path) || islink(path)) || return nothing

    if !Sys.iswindows()
        rm(path; force=force)
        return nothing
    end

    win_path = _windows_path(path)
    code = _windows_delete(path, win_path)
    code == 0 && return nothing
    # On local NTFS a file this process has mapped can still be deleted; on a network
    # share it cannot. A dropped Arrow table's mapping is released by a collection once
    # no worker thread still holds the finished task that read it, so release those
    # tasks, collect, and try again. Local deletes almost never get here.
    @debug_l1 "Windows deletion failed for $path ($(_windows_error_message(code))); releasing dropped mappings"
    _release_dropped_mappings()
    for attempt in 1:WINDOWS_DELETE_MAX_ATTEMPTS
        code = _windows_delete(path, win_path)
        code == 0 && return nothing
        attempt < WINDOWS_DELETE_MAX_ATTEMPTS && sleep(0.1 * attempt)   # e.g. a virus scanner
    end
    # Still held: by a table that is still referenced, or by another process. The
    # cmd.exe del this replaced reported success even then, so callers never saw a
    # failure here; keep it non-fatal and leave a warning.
    @user_warn "Could not delete $path: $(_windows_error_message(code))"
    return nothing
end

"""
    _release_dropped_mappings()

Give every default-pool thread a fresh task, so none still holds a finished task (and the
data it captured, such as an Arrow record batch over a mapped file), then run a full
collection, which unmaps files whose tables are no longer referenced. Collections are
serialized: parallel MaxLFQ writes have crashed Julia on Windows when several threads
entered GC.gc() at once.
"""
function _release_dropped_mappings()
    try
        Threads.@threads :static for _ in 1:Threads.nthreads(:default)
            nothing
        end
    catch
        # `:static` cannot run inside another threaded loop; queue one task per thread
        # instead (best effort: the scheduler may not place one on every thread).
        @sync for _ in 1:Threads.nthreads(:default)
            Threads.@spawn nothing
        end
    end
    lock(_WINDOWS_DELETE_GC_LOCK) do
        GC.gc(true)
    end
    return nothing
end

"""
    _windows_delete(path, win_path) -> UInt32

Delete with `DeleteFileW`; return 0 once `path` is gone (a file that was already missing
counts) or the Windows error code. After a successful call a network share can keep the
name briefly in "delete pending" state, listed but unopenable, so a rename onto the same
path would fail; wait up to `WINDOWS_DELETE_SETTLE_SECONDS` for the name to clear.
"""
function _windows_delete(path::String, win_path::String)
    code = _windows_delete_file(win_path)
    code in (_WINDOWS_ERROR_FILE_NOT_FOUND, _WINDOWS_ERROR_PATH_NOT_FOUND) && return UInt32(0)
    code == 0 || return code
    deadline = time() + WINDOWS_DELETE_SETTLE_SECONDS
    while !_windows_name_gone(path)
        if time() > deadline
            @debug_l1 "Windows deletion of $path succeeded but the name was still listed after $(WINDOWS_DELETE_SETTLE_SECONDS) s"
            break
        end
        sleep(0.01)
    end
    return UInt32(0)
end

# A delete-pending name makes stat fail (EACCES); that is not gone yet.
_windows_name_gone(path::String) =
    try
        !(isfile(path) || islink(path))
    catch
        false
    end

_windows_error_message(code::UInt32) = "DeleteFileW error $code: $(rstrip(strip(Libc.FormatMessage(code)), '.'))"

const TEMP_DATA_RETRY_SECONDS = 1.0

"""
    remove_temp_data(temp_path; rm_entry) -> Vector{String}

Delete a search's `temp_data` folder and return the paths that could not be removed
(empty when it is gone).

A single `rm(temp_path; recursive=true)` stops at the first entry it cannot delete.
On NFS a file that is still open is renamed to `.nfsXXXX` instead of deleted, so its
directory stays non-empty; the score calibrations' workspace (`jl_*`) sorts before
`main_search_psms`, so every PSM table behind it was left on disk. Open calibrations
are therefore closed first, each entry is removed on its own, and entries that still
fail are retried once after garbage collection. `rm_entry` exists for tests.
"""
function remove_temp_data(temp_path::AbstractString;
                          rm_entry = p -> rm(p; recursive=true, force=true))
    isdir(temp_path) || return String[]
    close_open_score_calibrations!()
    GC.gc()
    failed = String[]
    for attempt in 1:2
        empty!(failed)
        for entry in readdir(temp_path; join=true)
            try
                rm_entry(entry)
            catch e
                @debug_l1 "temp_data: could not remove $entry (attempt $attempt): $e"
                push!(failed, entry)
            end
        end
        isempty(failed) && break
        if attempt == 1
            GC.gc(true)
            sleep(TEMP_DATA_RETRY_SECONDS)
        end
    end
    if isempty(failed)
        try
            rm_entry(temp_path)
        catch e
            @debug_l1 "temp_data: could not remove $temp_path: $e"
            push!(failed, temp_path)
        end
    end
    return failed
end
