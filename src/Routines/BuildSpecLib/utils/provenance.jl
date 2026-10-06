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

# Making a library carry its own inputs and identify them.
#
# `fasta_paths` names the inputs only as paths on the build machine, which tell a
# downloader nothing. UniProt serves per-proteome FASTAs only for the current
# release, so after a release rolls over the exact file is recoverable only from
# a ~100 GB archive. The library therefore keeps a copy of every input FASTA in
# `<lib>.poin/fasta/`, plus a content-addressed record of each in config.json.

"""
Suffix of the optional sidecar describing where a FASTA came from: next to
`UP000005640_9606.fasta.gz`, a `UP000005640_9606.fasta.gz.provenance.json` with
any fields worth keeping (UniProt release and date, source URL, download time).
FASTA headers carry no release information, so this is the only way it reaches
the library.
"""
const FASTA_PROVENANCE_SUFFIX = ".provenance.json"

"""
    fasta_provenance_record(path, name) -> Dict{String,Any}

Identify one input FASTA by content: its name, base name, size and the SHA-256
of its bytes, merged with the fields of its provenance sidecar if one exists.

Throws if the sidecar records a `sha256` that does not match the file: the
sidecar would then describe some other download.
"""
function fasta_provenance_record(path::AbstractString, name::AbstractString)
    isfile(path) || throw(ArgumentError("FASTA file not found: $path"))
    sha = open(io -> bytes2hex(SHA.sha256(io)), path)
    record = Dict{String, Any}()
    sidecar = path * FASTA_PROVENANCE_SUFFIX
    if isfile(sidecar)
        merge!(record, JSON.parse(read(sidecar, String)))
        recorded = get(record, "sha256", sha)
        recorded == sha || throw(ArgumentError(
            "$sidecar records sha256 $recorded, but $path has sha256 $sha"))
    end
    record["name"] = String(name)
    record["file"] = basename(String(path))
    record["bytes"] = filesize(path)
    record["sha256"] = sha
    return record
end

"""
    bundle_fastas!(params, lib_dir) -> params

Copy every input FASTA (and its provenance sidecar, if any) into
`lib_dir/fasta/`, and record each in `params["fasta_provenance"]` with
`bundled_file` giving its path relative to the library. Two inputs with the same
base name are kept apart by prefixing the later one with its position.
"""
function bundle_fastas!(params::AbstractDict, lib_dir::AbstractString)
    paths = params["fasta_paths"]
    names = params["fasta_names"]
    fasta_dir = joinpath(lib_dir, "fasta")
    mkpath(fasta_dir)
    used = Set{String}()
    records = Dict{String, Any}[]
    for (i, path) in enumerate(paths)
        record = fasta_provenance_record(path, i <= length(names) ? string(names[i]) : "")
        dest = basename(path)
        dest in used && (dest = "$(i)_$(dest)")
        push!(used, dest)
        cp(path, joinpath(fasta_dir, dest); force = true)
        sidecar = path * FASTA_PROVENANCE_SUFFIX
        isfile(sidecar) && cp(sidecar, joinpath(fasta_dir, dest * FASTA_PROVENANCE_SUFFIX); force = true)
        record["bundled_file"] = "fasta/" * dest
        push!(records, record)
    end
    params["fasta_provenance"] = records
    return params
end

"""
    stamp_build_provenance!(params, lib_dir) -> params

Record what built the library and which inputs went into it, in the parameters
that become config.json: `pioneer_version`, `build_date` (UTC), the resolved
`library_params["prediction_model"]` (BuildSpecLib defaults it to "altimeter"
after config.json is written, so an omitted key would otherwise leave the
library unable to say what predicted its fragments), and the bundled FASTAs.
"""
function stamp_build_provenance!(params::AbstractDict, lib_dir::AbstractString)
    params["pioneer_version"] = get_pioneer_version()
    # "Z" appended rather than placed in the format string, where Dates reads it
    # as a timezone token.
    params["build_date"] = Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SS") * "Z"
    library_params = params["library_params"]
    library_params["prediction_model"] = String(get(library_params, "prediction_model", "altimeter"))
    bundle_fastas!(params, lib_dir)
    return params
end
