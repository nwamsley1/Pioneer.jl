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

# Safe file operations for cross-platform compatibility

const WINDOWS_DELETE_MAX_ATTEMPTS = 3
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
function _windows_delete_file(win_path::String)
    ccall((:SetFileAttributesW, "kernel32"), stdcall, Cint, (Cwstring, UInt32),
          win_path, _WINDOWS_FILE_ATTRIBUTE_NORMAL)
    ccall((:DeleteFileW, "kernel32"), stdcall, Cint, (Cwstring,), win_path) != 0 &&
        return UInt32(0)
    return UInt32(Libc.GetLastError())
end

"""
    safeRm(fpath; force=false)

Safely remove a file with Windows-specific handling for file locks and permissions.

# Arguments
- `fpath`: Path to file to remove
- `force`: Force removal on Unix. Windows preserves its historical forced-delete behavior.

# Implementation
- On Windows: Normalize once, retry `DeleteFileW`, fall back to rename, then
  use serialized garbage collection and forced removal as a last resort
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
    delete_error = nothing
    for attempt in 1:WINDOWS_DELETE_MAX_ATTEMPTS
        code = _windows_delete_file(win_path)
        # Not found means another caller removed it first.
        code in (0, _WINDOWS_ERROR_FILE_NOT_FOUND, _WINDOWS_ERROR_PATH_NOT_FOUND) && return nothing
        delete_error = ErrorException("DeleteFileW error $code: $(strip(Libc.FormatMessage(code)))")
        @debug_l1 "Windows deletion failed on attempt $attempt for $path: $(sprint(showerror, delete_error))"
        attempt < WINDOWS_DELETE_MAX_ATTEMPTS && sleep(0.1 * attempt)
    end

    backup_path = path * ".backup_" * string(time_ns())
    try
        mv(path, backup_path; force=true)
        @user_warn "Could not delete $path, renamed to $backup_path"
        return nothing
    catch rename_error
        @debug_l1 "Windows rename fallback failed for $path: $(sprint(showerror, rename_error))"

        # Parallel MaxLFQ writes have previously crashed Julia on Windows when
        # several threads entered GC.gc() concurrently. Keep this last-resort
        # collection serialized even though the deletion logic is centralized.
        try
            lock(_WINDOWS_DELETE_GC_LOCK) do
                GC.gc(true)
            end
            rm(path; force=true)
            return nothing
        catch final_error
            error(
                "Unable to remove or rename file $path. " *
                "Delete error: $(sprint(showerror, delete_error)). " *
                "Rename error: $(sprint(showerror, rename_error)). " *
                "Final removal error: $(sprint(showerror, final_error))."
            )
        end
    end
end
