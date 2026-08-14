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

# Reading the Hugging Face dataset repository that hosts prebuilt Pioneer
# libraries. Read-only by construction: nothing here writes to the remote.
#
# The network is confined to `hf_get` / `hf_stream`; everything else is pure and
# is unit-tested against recorded JSON.

"""
Dataset repository consulted when `--repo` is not given.
"""
const HF_DEFAULT_REPO = "nwamsley/altimeter_spectral_libraries"

"""
Optional curation manifest at the repo root. Absent is normal: metadata then
comes from each library's own `config.json`.
"""
const HF_MANIFEST_FILE = "libraries.json"

"""
Per-library checksum file, as written by the publisher.
"""
const HF_CHECKSUM_FILE = "SHA256SUMS"

"""
One file inside the repository.

`path` is repo-relative (`"Foo.poin/config.json"`); `size` is the byte count the
tree API reports, which for LFS-backed files is the real size rather than the
size of the pointer.
"""
struct HFFile
    path::String
    size::Int
end

"""
One spectral library: a repository directory whose name ends in `.poin` or
`.pion`, plus every file beneath it.
"""
struct HFLibrary
    name::String
    files::Vector{HFFile}
    total_bytes::Int
end

hf_tree_url(repo::AbstractString) =
    "https://huggingface.co/api/datasets/$(repo)/tree/main?recursive=true"

hf_resolve_url(repo::AbstractString, path::AbstractString) =
    "https://huggingface.co/datasets/$(repo)/resolve/main/$(path)"

"""
    is_library_dir(name) -> Bool

A directory is a library iff its name ends in `.poin` or `.pion`. Both are in
the wild: `BuildSpecLib` writes `.poin`, but `.pion`-named libraries exist and
load fine, so both are accepted.
"""
is_library_dir(name::AbstractString) = endswith(name, ".poin") || endswith(name, ".pion")

"""
    _entry_size(entry) -> Int

Byte size of a tree entry. LFS-backed files carry the real size under `lfs`,
while the top-level `size` is the pointer; small files have only `size`.
"""
function _entry_size(entry::AbstractDict)
    lfs = get(entry, "lfs", nothing)
    if lfs isa AbstractDict
        sz = get(lfs, "size", nothing)
        sz isa Number && return Int(sz)
    end
    sz = get(entry, "size", nothing)
    return sz isa Number ? Int(sz) : 0
end

"""
    parse_tree(tree_json) -> Vector{HFLibrary}

Group a recursive tree listing into libraries. Files outside any library
directory (a repo-root README, `.gitattributes`) are ignored.

Pure: the caller supplies the JSON text.
"""
function parse_tree(tree_json::AbstractString)
    entries = JSON.parse(tree_json)
    by_lib = Dict{String, Vector{HFFile}}()
    for entry in entries
        entry isa AbstractDict || continue
        get(entry, "type", "") == "file" || continue
        path = String(get(entry, "path", ""))
        isempty(path) && continue
        lib_name = String(first(split(path, '/')))
        is_library_dir(lib_name) || continue
        push!(get!(by_lib, lib_name, HFFile[]), HFFile(path, _entry_size(entry)))
    end

    libraries = HFLibrary[]
    for name in sort!(collect(keys(by_lib)))
        files = sort!(by_lib[name]; by = f -> f.path)
        push!(libraries, HFLibrary(name, files, sum(f -> f.size, files; init = 0)))
    end
    return libraries
end

"""
    hf_get(url) -> String

Fetch a small text resource. This and `hf_stream` are the only functions in the
feature that touch the network; everything above them is pure, and everything
below takes one of them as an injectable keyword so tests can serve fixtures
from a local socket.
"""
function hf_get(url::AbstractString)
    response = HTTP.get(url; status_exception = true)
    return String(response.body)
end

"""
    list_libraries(repo; fetch=hf_get) -> Vector{HFLibrary}

Every library in `repo`, newest listing each call. A library that has just been
uploaded appears here with no metadata file to update.
"""
list_libraries(repo::AbstractString = HF_DEFAULT_REPO; fetch = hf_get) =
    parse_tree(fetch(hf_tree_url(repo)))

"""
    fetch_json(repo, path; fetch=hf_get) -> Union{Dict, Nothing}

Fetch and parse one JSON file, returning `nothing` when it is absent or
unparseable. Callers treat a missing file as "no curation", never as an error:
the manifest is optional and a library without one still lists.
"""
function fetch_json(repo::AbstractString, path::AbstractString; fetch = hf_get)
    text = try
        fetch(hf_resolve_url(repo, path))
    catch
        return nothing
    end
    parsed = try
        JSON.parse(text)
    catch
        return nothing
    end
    return parsed isa AbstractDict ? Dict{String, Any}(parsed) : nothing
end
