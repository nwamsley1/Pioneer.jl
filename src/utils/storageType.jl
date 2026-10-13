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

# ── What kind of drive a file is on ──────────────────────────────────────────
#
# SearchDIA reads each MS data file out of order (the fragment-index search walks the library partition by
# partition, and BitVec calibration samples scans at random). On a mechanical drive every such read is a seek:
# on one lab PC a 119-file search took 32 h with the data on a hard drive, against an expected 5-6 h on an SSD.
# `storage_kind` asks the operating system what a path is stored on so the search can warn about it. It never
# throws: anything it cannot determine is `:unknown`, and an unknown drive is never warned about.

const STORAGE_KINDS = (:ssd, :hdd, :optical, :network, :unknown)
const STORAGE_PROBE_TIMEOUT_S = 10.0

"""
    storage_kind(path) -> Symbol

The kind of storage `path` is on: `:ssd`, `:hdd` (a spinning disk), `:optical`, `:network` (a network share) or
`:unknown`. Asks the operating system: the kernel's `rotational` flag on Linux, `diskutil`'s `SolidState` on
macOS, `Get-PhysicalDisk`'s `MediaType` on Windows. Never throws.
"""
function storage_kind(path::AbstractString)::Symbol
    p = try
        realpath(expanduser(String(path)))
    catch
        return :unknown
    end
    try
        Sys.islinux() && return _storage_kind_linux(p)
        Sys.isapple() && return _storage_kind_macos(p)
        Sys.iswindows() && return _storage_kind_windows(p)
    catch
    end
    return :unknown
end

"Human-readable name of a storage kind, for messages."
storage_kind_label(kind::Symbol) =
    kind === :hdd ? "a mechanical hard drive" :
    kind === :optical ? "an optical drive" :
    kind === :ssd ? "an SSD" :
    kind === :network ? "a network share" : "an unknown drive"

"Whether searching data on `kind` storage should be warned about (mechanical and optical drives)."
storage_is_slow(kind::Symbol) = kind === :hdd || kind === :optical

"Run `cmd` and return its stdout, or `nothing` if it fails or takes longer than `timeout` seconds."
function _read_with_timeout(cmd::Cmd; timeout::Real = STORAGE_PROBE_TIMEOUT_S)
    out = Pipe()
    proc = try
        run(pipeline(ignorestatus(cmd); stdout = out, stderr = devnull); wait = false)
    catch
        return nothing
    end
    close(out.in)
    reader = @async read(out, String)
    if timedwait(() -> process_exited(proc), Float64(timeout)) !== :ok
        kill(proc)
        return nothing
    end
    text = fetch(reader)
    return success(proc) ? text : nothing
end

# ── Linux ────────────────────────────────────────────────────────────────────

const _LINUX_NETWORK_FS = ("nfs", "nfs4", "cifs", "smb3", "smbfs", "sshfs", "fuse.sshfs", "9p", "afs",
                           "glusterfs", "fuse.glusterfs", "ceph", "lustre", "gpfs", "beegfs", "fuse.s3fs")

"Unescape a mountinfo field (`\\040` for a space and the like)."
_unescape_mountinfo(s::AbstractString) = replace(s, r"\\([0-7]{3})" => m -> string(Char(parse(Int, m[2:end]; base = 8))))

"""
    _linux_mount_for(path, mountinfo) -> (fstype, source) or nothing

The filesystem type and source of the mount containing `path`, from the text of `/proc/self/mountinfo` (the
mount point that is the longest prefix of `path`).
"""
function _linux_mount_for(path::AbstractString, mountinfo::AbstractString)
    best, best_len = nothing, -1
    for line in eachsplit(mountinfo, '\n')
        isempty(line) && continue
        halves = split(line, " - "; limit = 2)
        length(halves) == 2 || continue
        pre, post = split(halves[1]), split(halves[2])
        (length(pre) >= 5 && length(post) >= 2) || continue
        mp = _unescape_mountinfo(pre[5])
        inside = path == mp || mp == "/" || startswith(path, rstrip(mp, '/') * "/")
        if inside && length(mp) > best_len
            best, best_len = (String(post[1]), _unescape_mountinfo(post[2])), length(mp)
        end
    end
    return best
end

"""
    _linux_block_kind(sysdev) -> Symbol

Kind of the block device at `/sys/dev/block/MAJ:MIN` (or its `/sys/class/block/NAME` directory): partitions are
resolved to their disk; device-mapper, RAID and similar stacked devices are `:hdd` if any underlying disk is
rotational, `:ssd` if all are not.
"""
function _linux_block_kind(sysdev::AbstractString; depth::Int = 0)::Symbol
    depth > 8 && return :unknown
    dev = realpath(sysdev)
    disk = isfile(joinpath(dev, "partition")) ? dirname(dev) : dev
    name = basename(disk)
    startswith(name, "sr") && return :optical
    typefile = joinpath(disk, "device", "type")
    isfile(typefile) && strip(read(typefile, String)) == "5" && return :optical   # SCSI type 5: CD/DVD
    slaves = joinpath(disk, "slaves")
    if isdir(slaves) && !isempty(readdir(slaves))
        kinds = [_linux_block_kind(joinpath(slaves, s); depth = depth + 1) for s in readdir(slaves)]
        any(==(:hdd), kinds) && return :hdd
        any(==(:optical), kinds) && return :optical
        all(==(:ssd), kinds) && return :ssd
        return :unknown
    end
    (startswith(name, "loop") || startswith(name, "zram") || startswith(name, "ram")) && return :unknown
    rot = joinpath(disk, "queue", "rotational")
    isfile(rot) || return :unknown
    flag = strip(read(rot, String))
    return flag == "1" ? :hdd : flag == "0" ? :ssd : :unknown
end

_linux_major(dev::UInt64) = Int(((dev >> 8) & 0xfff) | ((dev >> 32) & ~UInt64(0xfff)))
_linux_minor(dev::UInt64) = Int((dev & 0xff) | ((dev >> 12) & ~UInt64(0xff)))

function _storage_kind_linux(p::String)::Symbol
    mount = isfile("/proc/self/mountinfo") ? _linux_mount_for(p, read("/proc/self/mountinfo", String)) : nothing
    if mount !== nothing
        fstype, source = mount
        (fstype in _LINUX_NETWORK_FS || startswith(source, "//")) && return :network
    end
    dev = UInt64(stat(p).device)
    major, minor = _linux_major(dev), _linux_minor(dev)
    # btrfs, overlayfs and similar report an anonymous device (major 0): fall back to the mount's source device
    if major == 0 && mount !== nothing && startswith(mount[2], "/dev/") && ispath(mount[2])
        rdev = UInt64(stat(mount[2]).rdev)
        major, minor = _linux_major(rdev), _linux_minor(rdev)
    end
    major == 0 && return :unknown
    sysdev = "/sys/dev/block/$major:$minor"
    ispath(sysdev) || return :unknown
    return _linux_block_kind(sysdev)
end

# ── macOS ────────────────────────────────────────────────────────────────────

"""
    _macos_plist_kind(plist) -> Symbol

Kind of a disk from the text of `diskutil info -plist`: optical media first, then `SolidState`.
"""
function _macos_plist_kind(plist::AbstractString)::Symbol
    occursin(r"<key>OpticalMediaType</key>|<key>OpticalDeviceType</key>", plist) && return :optical
    m = match(r"<key>SolidState</key>\s*<(true|false)\s*/>", plist)
    m === nothing && return :unknown
    return m[1] == "true" ? :ssd : :hdd
end

function _storage_kind_macos(p::String)::Symbol
    df = _read_with_timeout(`df -P $p`)
    df === nothing && return :unknown
    lines = split(strip(df), '\n')
    length(lines) >= 2 || return :unknown
    source = String(first(split(lines[end])))
    (startswith(source, "//") || occursin(r"^[^/].*:/", source)) && return :network   # smb/afp "//host/x", nfs "host:/x"
    startswith(source, "/dev/") || return :unknown
    plist = _read_with_timeout(`diskutil info -plist $source`)
    plist === nothing && return :unknown
    return _macos_plist_kind(plist)
end

# ── Windows ──────────────────────────────────────────────────────────────────

"""
    _windows_probe_script(drive_letter) -> String

PowerShell that prints one word for a drive letter: NETWORK, OPTICAL, or the physical disk's MediaType (SSD,
HDD, SCM, Unspecified).
"""
_windows_probe_script(letter::Char) = """
\$ErrorActionPreference = 'Stop'
\$ld = Get-CimInstance Win32_LogicalDisk -Filter "DeviceID='$(letter):'"
if (\$ld.DriveType -eq 4) { 'NETWORK'; exit }
if (\$ld.DriveType -eq 5) { 'OPTICAL'; exit }
\$disk = Get-Partition -DriveLetter $(letter) | Get-Disk
\$pd = Get-PhysicalDisk | Where-Object { \$_.DeviceId -eq [string]\$disk.Number }
if (\$pd) { [string]\$pd.MediaType } else { 'Unspecified' }
"""

"Kind from the word `_windows_probe_script` prints."
function _windows_probe_kind(word::AbstractString)::Symbol
    w = uppercase(strip(word))
    w == "NETWORK" && return :network
    w == "OPTICAL" && return :optical
    w == "HDD" && return :hdd
    (w == "SSD" || w == "SCM") && return :ssd
    return :unknown
end

function _storage_kind_windows(p::String)::Symbol
    (startswith(p, "\\\\") || startswith(p, "//")) && return :network        # UNC path
    m = match(r"^([A-Za-z]):", p)
    m === nothing && return :unknown
    script = _windows_probe_script(uppercase(m[1][1]))
    out = _read_with_timeout(`powershell -NoProfile -NonInteractive -ExecutionPolicy Bypass -Command $script`)
    out === nothing && return :unknown
    return _windows_probe_kind(out)
end

# ── Search-time check ────────────────────────────────────────────────────────

"A key that is the same for paths on the same volume, so each volume is probed once."
function _storage_volume_key(path::AbstractString)
    p = realpath(path)
    if Sys.iswindows()
        m = match(r"^([A-Za-z]:|\\\\[^\\]+\\[^\\]+)", p)     # drive letter, or \\server\share
        return uppercase(m === nothing ? p : m.match)
    end
    return stat(p).device
end

"""
    slow_storage_paths(paths) -> Vector{Tuple{String, Symbol}}

`(path, kind)` for one representative path per volume among `paths` whose volume is a mechanical or optical
drive. Symlinks are resolved first, so a folder of links is judged by where the files really are.
"""
function slow_storage_paths(paths::AbstractVector{<:AbstractString})
    seen = Dict{Any, Symbol}()
    slow = Tuple{String, Symbol}[]
    for path in paths
        key = try
            _storage_volume_key(path)
        catch
            continue
        end
        haskey(seen, key) && continue
        kind = storage_kind(path)
        seen[key] = kind
        storage_is_slow(kind) && push!(slow, (realpath(path), kind))
    end
    return slow
end

"""
    warn_slow_storage(ms_paths, library_dir, results_dir)

Warn (never fail) when the MS data files, the spectral library or the output folder are on a mechanical or
optical drive. The output folder matters as much as the inputs: a search writes and re-reads large temporary
files there.
"""
function warn_slow_storage(ms_paths::AbstractVector{<:AbstractString}, library_dir::AbstractString,
                           results_dir::AbstractString)
    try
        for (path, kind) in slow_storage_paths(ms_paths)
            @user_warn "The MS data ($(dirname(path))) are on $(storage_kind_label(kind)). " *
                       "This may result in substantially longer search times."
        end
        for (path, kind) in slow_storage_paths([library_dir])
            @user_warn "The spectral library ($path) is on $(storage_kind_label(kind)). " *
                       "This may result in substantially longer search times."
        end
        for (path, kind) in slow_storage_paths([results_dir])
            @user_warn "The output folder ($path) is on $(storage_kind_label(kind)). " *
                       "This may result in substantially longer search times."
        end
    catch
    end
    return nothing
end
