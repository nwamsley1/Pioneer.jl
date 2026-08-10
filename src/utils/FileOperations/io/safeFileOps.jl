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

"""
    _windows_delete_command(fpath; force=false)

Build a deterministic `cmd.exe del` command for `fpath`. Windows command
built-ins interpret `/` as the start of a switch, so paths must be made
absolute and converted to backslashes before they reach `del`.
"""
function _windows_delete_command(fpath::AbstractString; force::Bool=false)
    win_path = replace(abspath(normpath(String(fpath))), "/" => "\\")
    args = String["cmd.exe", "/d", "/c", "del", "/q"]
    force && push!(args, "/f")
    push!(args, win_path)
    return Cmd(args)
end

"""
    safeRm(fpath; force=false)

Safely remove a file with Windows-specific handling for file locks and permissions.

# Arguments
- `fpath`: Path to file to remove
- `force`: Force removal even if file is read-only

# Implementation
- On Windows: Normalize once, retry `cmd.exe del`, then fall back to rename
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

    delete_cmd = _windows_delete_command(path; force=force)
    for attempt in 1:WINDOWS_DELETE_MAX_ATTEMPTS
        try
            run(delete_cmd)
            return nothing
        catch delete_error
            @debug_l1 "Windows deletion failed on attempt $attempt for $path: $(sprint(showerror, delete_error))"
            if attempt == WINDOWS_DELETE_MAX_ATTEMPTS
                backup_path = path * ".backup_" * string(time_ns())
                try
                    mv(path, backup_path; force=force)
                    @user_warn "Could not delete $path, renamed to $backup_path"
                    return nothing
                catch rename_error
                    error(
                        "Unable to remove or rename file $path. " *
                        "Delete error: $(sprint(showerror, delete_error)). " *
                        "Rename error: $(sprint(showerror, rename_error))."
                    )
                end
            end

            # Arrow memory maps can keep a file locked on Windows until their
            # finalizers run. Collect only after a real deletion failure.
            attempt == 1 && GC.gc(true)
            sleep(0.05 * 2^(attempt - 1))
        end
    end

    return nothing
end
