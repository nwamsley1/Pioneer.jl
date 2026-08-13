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

# Transferring a library out of the repository.
#
# Two invariants the rest of Pioneer depends on:
#
#   * A library directory exists only once it is complete. Everything lands in
#     `<name>.partial` and is renamed at the end, because `library_info` (and
#     the search itself) accept a directory by its marker files -- a half-copied
#     one would pass that check and then fail deep inside a search.
#   * A resumed transfer never silently concatenates. A server that ignores our
#     Range header answers 200 with the whole body, and appending that to what
#     we already had would produce a corrupt file whose size looks plausible.

"""
Most hops any one file may take. huggingface.co uses one (to its CDN); the cap
exists so a redirect loop fails rather than spinning.
"""
const HF_MAX_REDIRECTS = 5

"""
    hf_download_file(url, path, expected_size; on_chunk) -> Int

Fetch one file to `path`, resuming when a partial copy is already there. Returns
the byte count on disk.

Resume is a `Range` request; a 206 continues the file, anything else restarts it
from scratch rather than appending to a body that starts at byte zero. A local
file longer than `expected_size` is treated as junk and discarded -- it cannot
be a prefix of the real thing.
"""
function hf_download_file(url::AbstractString, path::AbstractString, expected_size::Integer;
                          on_chunk = nothing)
    existing = isfile(path) ? filesize(path) : 0

    if expected_size > 0
        existing == expected_size && return existing          # already complete
        existing > expected_size && (rm(path; force = true); existing = 0)
    end

    current = String(url)
    for _ in 1:HF_MAX_REDIRECTS
        # Identity encoding is required, not merely preferred. HTTP.jl asks for
        # gzip by default and does not decompress in streaming mode, so the
        # bytes reaching the file would be the compressed ones -- a short file
        # that starts with 1f 8b. It also keeps the byte offsets below honest:
        # a Range is relative to the representation the server sends, so
        # resuming a compressed stream by identity offset is meaningless.
        headers = Pair{String, String}["Accept-Encoding" => "identity"]
        existing > 0 && push!(headers, "Range" => "bytes=$(existing)-")

        next_url = ""
        restart = false
        total = existing

        # `redirect = false` is load-bearing. Left to itself HTTP.jl runs this
        # callback once per hop -- including the 307 that huggingface.co answers
        # with -- and the body of that hop is the string "Temporary Redirect.
        # Redirecting to ...". Streaming every hop into the same file therefore
        # prepends ~1 KB of prose to a multi-gigabyte binary, producing a file
        # of plausible size and wrong contents. Following redirects here means
        # only the final response is ever written.
        HTTP.open("GET", current, headers; redirect = false, status_exception = false) do stream
            response = HTTP.startread(stream)

            if response.status in (301, 302, 303, 307, 308)
                location = HTTP.header(response, "Location", "")
                isempty(location) && error("redirect with no Location from $current")
                # Relative, in HF's case: /api/resolve-cache/datasets/...
                next_url = string(HTTP.URIs.resolvereference(HTTP.URI(current), location))
                while !eof(stream)                       # discard, never write
                    readavailable(stream)
                end
                return
            end

            if response.status >= 400
                error("HTTP $(response.status) fetching $current")
            end

            # 206 = our range was honoured. Anything else restarts the file,
            # since a 200 body begins at byte zero and appending it would
            # corrupt what is already there.
            if existing > 0 && response.status != 206
                restart = true
                while !eof(stream)
                    readavailable(stream)
                end
                return
            end

            open(path, existing > 0 ? "a" : "w") do out
                while !eof(stream)
                    chunk = readavailable(stream)
                    total += write(out, chunk)
                    on_chunk === nothing || on_chunk(length(chunk))
                end
            end
        end

        if !isempty(next_url)
            current = next_url
            continue
        end
        if restart
            rm(path; force = true)
            existing = 0
            continue
        end
        return total
    end
    error("too many redirects fetching $url")
end

"""
    parse_checksums(text) -> Dict{String,String}

Parse `SHA256SUMS` (`<hex>  <filename>` per line) into filename => hex. Pure.
"""
function parse_checksums(text::AbstractString)
    sums = Dict{String, String}()
    for line in eachline(IOBuffer(text))
        stripped = strip(line)
        isempty(stripped) && continue
        fields = split(stripped)
        length(fields) >= 2 || continue
        sums[String(last(fields))] = String(lowercase(first(fields)))
    end
    return sums
end

"""
    verify_checksums(dir, sums) -> Vector{String}

Names under `dir` whose contents do not match `sums`, plus any that are missing.
Empty means everything verified.

Hashing streams the file: `detailed_fragments.jls` is comfortably over a
gigabyte and must never be slurped into memory to be checked.
"""
function verify_checksums(dir::AbstractString, sums::AbstractDict)
    bad = String[]
    for name in sort!(collect(keys(sums)))
        path = joinpath(dir, name)
        if !isfile(path)
            push!(bad, name)
            continue
        end
        actual = open(io -> bytes2hex(SHA.sha256(io)), path)
        actual == sums[name] || push!(bad, name)
    end
    return bad
end

"""
    download_library(library, dest, repo; force, on_progress, on_file) -> String

Download every file of `library` into `dest`, verify it, and return the final
directory path.

`on_progress(bytes_done, bytes_total)` fires as bytes land; `on_file(name, i, n)`
fires as each file starts. The CLI turns these into progress output, which the
GUI's log drawer renders.
"""
function download_library(library::HFLibrary, dest::AbstractString,
                          repo::AbstractString = HF_DEFAULT_REPO;
                          force::Bool = false,
                          on_progress = (done, total) -> nothing,
                          on_file = (name, i, n) -> nothing,
                          url_for = path -> hf_resolve_url(repo, path))
    isdir(dest) || mkpath(dest)
    final = joinpath(dest, library.name)

    if ispath(final)
        force || throw(ArgumentError(
            "$(final) already exists. Choose another destination or pass --force to replace it."))
    end

    staging = final * ".partial"
    isdir(staging) || mkpath(staging)

    done = 0
    total = library.total_bytes
    for (i, file) in enumerate(library.files)
        # Path within the library, so nested files keep their layout.
        relative = String(chopprefix(file.path, library.name * "/"))
        target = joinpath(staging, relative)
        mkpath(dirname(target))
        on_file(relative, i, length(library.files))
        hf_download_file(url_for(file.path), target, file.size;
                         on_chunk = n -> (done += n; on_progress(done, total)))
    end

    sums_path = joinpath(staging, HF_CHECKSUM_FILE)
    if isfile(sums_path)
        bad = verify_checksums(staging, parse_checksums(read(sums_path, String)))
        isempty(bad) || throw(ErrorException(
            "checksum verification failed for: " * join(bad, ", ") *
            ". The partial download is kept at $(staging) for inspection."))
    end

    # Only now does anything named like a library exist.
    force && ispath(final) && rm(final; recursive = true)
    mv(staging, final)
    return final
end
