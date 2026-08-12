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

# =============================================================================
# Windows LightGBM performance fix (MSVC DLL override)
# =============================================================================
#
# The `LightGBM_jll` binary that Julia ships on Windows is built with mingw-w64,
# and its OpenMP runtime (GNU libgomp on winpthreads) does NOT parallelize the
# LightGBM tree-learner FIT: on an 18-core box it runs at ~1 core. Microsoft's
# MSVC build of the SAME LightGBM version (3.3.5) saturates all cores with
# bit-identical results (LightGBM is seeded/deterministic in Pioneer). See
# `src/build/windows/LIGHTGBM_MSVC.md` for the full write-up, benchmarks, and
# the packaged-app equivalent of this fix.
#
# `setup_windows_lightgbm()` performs the developer / `Pkg.add` install-time fix:
# it fetches Microsoft's verified MSVC `lib_lightgbm.dll` v3.3.5 and points
# `LightGBM_jll` at it through Julia's supported artifact-override mechanism
# (`~/.julia/artifacts/Overrides.toml`). The packaged Windows app applies the
# equivalent swap at build time (see the CI workflow), so END USERS never run
# this — it is for source installs.
#
# Windows-only. No new Julia dependencies: download / SHA-256 / unzip are done
# through PowerShell (always present on Windows), TOML through `Base.TOML`.

# LightGBM_jll package UUID (the override is keyed by this).
const _LIGHTGBM_JLL_UUID = "0e4427ef-1ff7-5cd7-8faa-8ff0877bb2ec"

# Pinned Microsoft MSVC build of LightGBM 3.3.5, distributed in the official
# PyPI wheel (a zip). Both hashes are verified before anything is installed.
# NOTE: this pin is coupled to `LightGBM = "2.2"` in Project.toml, whose C API
# expects LightGBM 3.3.5. If the LightGBM.jl compat is bumped, update BOTH this
# pin and the CI swap in .github/workflows/build_app_windows.yml.
const _MSVC_WHEEL_URL =
    "https://files.pythonhosted.org/packages/6a/31/13f80e5abac627c53375c1dc49404d622d929272a231c05b2f4f7bb98b9e/lightgbm-3.3.5-py3-none-win_amd64.whl"
const _MSVC_WHEEL_SHA256 = "02a40745c1972cf4a2cde764c7739228f45178c2237af2df40fde7063a58ac6a"
const _MSVC_DLL_SHA256   = "a52a9494c894f9987f9be24c97e0dec7c37fd23b1184aa0ac39d3200ffa51699"

# Small PowerShell helpers (avoid adding Downloads/SHA/zip deps).
function _pwsh(cmd::String)
    return read(`powershell -NoProfile -NonInteractive -Command $cmd`, String)
end
_sha256_hex(path::String) =
    lowercase(strip(_pwsh("(Get-FileHash -Algorithm SHA256 -LiteralPath '$path').Hash")))

# Locate the currently-installed (MinGW) LightGBM artifact so we can copy its
# `lightgbm.exe`: the JLL wrapper declares two products, so the override dir
# must contain BOTH lib_lightgbm.dll and lightgbm.exe or __init__ fails. Pioneer
# only calls the C API, so a MinGW CLI exe alongside the MSVC library is fine.
function _find_current_lightgbm_exe()
    for depot in DEPOT_PATH
        adir = joinpath(depot, "artifacts")
        isdir(adir) || continue
        for h in readdir(adir)
            exe = joinpath(adir, h, "bin", "lightgbm.exe")
            isfile(exe) && return exe
        end
    end
    return nothing
end

# Merge (never clobber) the override into Overrides.toml. Uses Base.TOML to read
# any existing entries; writes a minimal, valid TOML back.
function _register_override!(override_root::String; force::Bool)
    ovr = joinpath(first(DEPOT_PATH), "artifacts", "Overrides.toml")
    data = Dict{String,Any}()
    if isfile(ovr)
        try
            data = Base.TOML.parsefile(ovr)
        catch e
            @user_warn "Could not parse existing $ovr ($e); it will be rewritten with only the LightGBM entry."
            data = Dict{String,Any}()
        end
    end
    entry = get!(() -> Dict{String,Any}(), data, _LIGHTGBM_JLL_UUID)
    entry isa AbstractDict || (entry = Dict{String,Any}())
    if !force && get(entry, "LightGBM", nothing) == override_root
        return (ovr, false)   # already registered at this path
    end
    entry["LightGBM"] = override_root
    data[_LIGHTGBM_JLL_UUID] = entry

    mkpath(dirname(ovr))
    io = IOBuffer()
    for (uuid, val) in data
        if val isa AbstractDict
            println(io, "[", uuid, "]")
            for (k, v) in val
                println(io, k, " = \"", replace(String(v), "\\" => "\\\\"), "\"")
            end
            println(io)
        else
            println(io, uuid, " = \"", replace(String(val), "\\" => "\\\\"), "\"")
        end
    end
    write(ovr, take!(io))
    return (ovr, true)
end

"""
    setup_windows_lightgbm(; force::Bool=false)

Install the Microsoft MSVC build of `lib_lightgbm` (v3.3.5) and register it as a
`LightGBM_jll` artifact override, so Pioneer's LightGBM stages actually use all
CPU cores on Windows. No-op on Linux/macOS. Idempotent. **Restart Julia** after
running for the override to take effect.

Why: the shipped mingw-w64 LightGBM binary barely threads the tree-learner FIT
(~1 of N cores) because of libgomp-on-winpthreads barrier cost; the MSVC build
saturates all cores with bit-identical (seeded, deterministic) results.

The pinned wheel and its DLL are SHA-256 verified before installation. Requires
the Visual C++ Redistributable (for `VCOMP140.DLL` etc.); if it is missing this
prints where to get it and continues laying out the files. For the packaged
Windows app the equivalent swap is done at build time, so end users do not run
this — see `src/build/windows/LIGHTGBM_MSVC.md`.

Pass `force=true` to reinstall even if an override is already registered.
"""
function setup_windows_lightgbm(; force::Bool=false)
    if !Sys.iswindows()
        @user_info "setup_windows_lightgbm: not on Windows; nothing to do (the shipped LightGBM binary is fine off Windows)."
        return nothing
    end

    # (1) Visual C++ Redistributable check (the MSVC DLL imports these).
    sys32 = joinpath(get(ENV, "WINDIR", "C:\\Windows"), "System32")
    missing_vc = filter(f -> !isfile(joinpath(sys32, f)),
                        ["vcomp140.dll", "vcruntime140.dll", "msvcp140.dll"])
    if !isempty(missing_vc)
        @user_warn """Visual C++ Redistributable appears to be missing: $(join(missing_vc, ", ")).
        The MSVC LightGBM DLL will fail to load without it. Install:
            https://aka.ms/vs/17/release/vc_redist.x64.exe
        then re-run setup_windows_lightgbm()."""
    end

    override_root = joinpath(first(DEPOT_PATH), "pioneer_lgbm_msvc")
    bin = joinpath(override_root, "bin")
    mkpath(bin)
    dll_dest = joinpath(bin, "lib_lightgbm.dll")

    # (2) Download the pinned MSVC wheel to a temp dir and verify its hash.
    work = mktempdir()
    whl = joinpath(work, "lightgbm-3.3.5-win_amd64.whl")
    @user_info "Downloading MSVC LightGBM 3.3.5 (verified by SHA-256)..."
    _pwsh("Invoke-WebRequest -Uri '$_MSVC_WHEEL_URL' -OutFile '$whl' -UseBasicParsing")
    got = _sha256_hex(whl)
    got == _MSVC_WHEEL_SHA256 ||
        error("wheel SHA-256 mismatch: expected $_MSVC_WHEEL_SHA256, got $got. Aborting (nothing installed).")

    # (3) Extract lib_lightgbm.dll (the wheel is a zip) and verify the DLL hash.
    exdir = joinpath(work, "unz")
    _pwsh("Expand-Archive -LiteralPath '$whl' -DestinationPath '$exdir' -Force")
    src_dll = joinpath(exdir, "lightgbm", "lib_lightgbm.dll")
    isfile(src_dll) || error("lib_lightgbm.dll not found in wheel at $src_dll")
    got_dll = _sha256_hex(src_dll)
    got_dll == _MSVC_DLL_SHA256 ||
        error("DLL SHA-256 mismatch: expected $_MSVC_DLL_SHA256, got $got_dll. Aborting.")
    cp(src_dll, dll_dest; force=true)

    # (4) Satisfy the JLL's second product (lightgbm.exe) from the current artifact.
    exe = _find_current_lightgbm_exe()
    if exe === nothing
        @user_warn "Could not find an existing lightgbm.exe to copy; load LightGBM once (`using LightGBM`) then re-run so the artifact is materialized."
    else
        cp(exe, joinpath(bin, "lightgbm.exe"); force=true)
    end

    # (5) Register the override.
    ovr, changed = _register_override!(override_root; force=force)
    if !changed
        @user_info "LightGBM MSVC override already registered at $override_root (use force=true to reinstall)."
        return nothing
    end

    @user_info """LightGBM MSVC override installed.
        DLL:      $dll_dest
        Override: $ovr
    RESTART Julia for it to take effect. Verify with:
        using Libdl; filter(h -> occursin("lightgbm", lowercase(h)), Libdl.dllist())
    which should now point at $bin."""
    return nothing
end
