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

"""
    _windows_delete_command(fpath)

Build a deterministic `cmd.exe del` command for `fpath`. Windows command
built-ins interpret `/` as the start of a switch, so paths must be made
absolute and converted to backslashes before they reach `del`. The `/f` flag
preserves the historical Windows behavior of removing read-only files.
"""
function _windows_delete_command(fpath::AbstractString)
    win_path = replace(abspath(normpath(String(fpath))), "/" => "\\")
    args = String["cmd.exe", "/d", "/c", "del", "/f", "/q"]
    push!(args, win_path)
    return Cmd(args)
end

"""
    safeRm(fpath; force=false)

Safely remove a file with Windows-specific handling for file locks and permissions.

# Arguments
- `fpath`: Path to file to remove
- `force`: Force removal on Unix. Windows preserves its historical forced-delete behavior.

# Implementation
- On Windows: Normalize once, retry `cmd.exe del`, fall back to rename, then
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

    delete_cmd = _windows_delete_command(path)
    delete_error = nothing
    for attempt in 1:WINDOWS_DELETE_MAX_ATTEMPTS
        try
            run(delete_cmd)
            return nothing
        catch delete_attempt_error
            delete_error = delete_attempt_error
            @debug_l1 "Windows deletion failed on attempt $attempt for $path: $(sprint(showerror, delete_attempt_error))"
            attempt < WINDOWS_DELETE_MAX_ATTEMPTS && sleep(0.1 * attempt)
        end
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
