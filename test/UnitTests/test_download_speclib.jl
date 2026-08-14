# Unit tests for DownloadSpecLib (src/Routines/DownloadSpecLib/).
#
# No test here reaches the internet. The parsing layer runs against JSON
# recorded from the real repository (test/data/download_speclib/), and the
# transfer layer runs against an HTTP server on 127.0.0.1 — so resume,
# truncation and checksum failure are exercised for real rather than mocked.
#
# Run standalone: julia --project=. test/UnitTests/test_download_speclib.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
end

# Outside the guard: the suite loads Pioneer but does not import these into
# Main, and this file uses both directly.
using HTTP
using SHA

const DL_FIXTURES = joinpath(@__DIR__, "..", "data", "download_speclib")

@testset "DownloadSpecLib" begin

# ═══════════════════════════════════════════════════════════════════════════════
# Tree parsing — against JSON recorded from the live repository
# ═══════════════════════════════════════════════════════════════════════════════

@testset "parse_tree" begin
    tree = read(joinpath(DL_FIXTURES, "tree.json"), String)
    libs = Pioneer.parse_tree(tree)

    @testset "finds the library, ignores repo-root files" begin
        # README.md and .gitattributes sit at the root and are not library files.
        @test length(libs) == 1
        @test libs[1].name == "Pioneer_Human_canon_std.poin"
        @test all(f -> startswith(f.path, "Pioneer_Human_canon_std.poin/"), libs[1].files)
    end

    @testset "reports real LFS sizes" begin
        # The tree's top-level `size` for an LFS file is the pointer; the real
        # size is what a download has to plan for.
        @test libs[1].total_bytes == 3456781421
        @test length(libs[1].files) == 15
        frag = libs[1].files[findfirst(f -> endswith(f.path, "detailed_fragments.jls"), libs[1].files)]
        @test frag.size == 1702026297
    end

    @testset "library directories are recognised by extension" begin
        @test Pioneer.is_library_dir("foo.poin")
        @test Pioneer.is_library_dir("foo.pion")   # exists in the wild
        @test !Pioneer.is_library_dir("foo")
        @test !Pioneer.is_library_dir("README.md")
    end

    @testset "empty and malformed trees do not throw" begin
        @test isempty(Pioneer.parse_tree("[]"))
        @test isempty(Pioneer.parse_tree("""[{"type":"file","path":"README.md","size":10}]"""))
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Catalog — derived metadata and manifest precedence
# ═══════════════════════════════════════════════════════════════════════════════

@testset "catalog" begin
    tree = read(joinpath(DL_FIXTURES, "tree.json"), String)
    libs = Pioneer.parse_tree(tree)
    config = Pioneer.JSON.parse(read(joinpath(DL_FIXTURES, "config.json"), String))

    @testset "describe_config derives the build parameters" begin
        details = Dict(Pioneer.describe_config(config))
        @test details["FASTA"] == "HUMAN, CONTAM"
        @test details["Missed cleavages"] == "1"
        @test details["Peptide length"] == "7-30"
        @test details["Fixed mods"] == "Unimod:4 (C)"
        @test details["Variable mods"] == "Unimod:35 (M)"
        @test details["NCE"] == "26.0"
        @test details["Precursor m/z"] == "359.0-1034.0"
    end

    @testset "a library with no manifest still lists" begin
        catalog = Pioneer.build_catalog(libs, nothing; config_for = _ -> config)
        @test length(catalog) == 1
        @test catalog[1].title == "Pioneer_Human_canon_std.poin"   # falls back to the name
        @test !isempty(catalog[1].details)
        @test catalog[1].total_bytes == 3456781421
    end

    @testset "manifest supplies what a config cannot know" begin
        manifest = Dict("libraries" => Dict(
            "Pioneer_Human_canon_std.poin" => Dict(
                "title" => "Human canonical",
                "model" => "Altimeter",
                "description" => "UniProt human canonical + contaminants.",
                "recommended_for" => "General human DIA.")))
        catalog = Pioneer.build_catalog(libs, manifest; config_for = _ -> config)
        @test catalog[1].title == "Human canonical"
        @test catalog[1].model == "Altimeter"
        @test catalog[1].recommended_for == "General human DIA."
        # No "details" in the manifest, so the config still supplied them.
        @test !isempty(catalog[1].details)
    end

    @testset "manifest details skip the config fetch entirely" begin
        manifest = Dict("libraries" => Dict(
            "Pioneer_Human_canon_std.poin" => Dict(
                "title" => "Human canonical",
                "details" => Dict("NCE" => "26", "FASTA" => "HUMAN"))))
        fetched = Ref(false)
        catalog = Pioneer.build_catalog(libs, manifest;
                                        config_for = _ -> (fetched[] = true; config))
        @test !fetched[]                       # the single-request fast path
        @test Dict(catalog[1].details)["NCE"] == "26"
    end

    @testset "an unreadable config degrades to no details" begin
        catalog = Pioneer.build_catalog(libs, nothing; config_for = _ -> nothing)
        @test length(catalog) == 1
        @test isempty(catalog[1].details)
    end

    @testset "json payload carries the contract fields" begin
        catalog = Pioneer.build_catalog(libs, nothing; config_for = _ -> config)
        payload = Pioneer.JSON.parse(Pioneer.catalog_to_json(catalog))
        entry = payload["libraries"][1]
        for key in ("name", "title", "model", "description", "recommended_for",
                    "total_bytes", "size_human", "n_files", "details")
            @test haskey(entry, key)
        end
        @test entry["total_bytes"] == 3456781421
    end

    @testset "format_bytes" begin
        @test Pioneer.format_bytes(512) == "512 B"
        @test Pioneer.format_bytes(3456781421) == "3.22 GiB"
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Checksums
# ═══════════════════════════════════════════════════════════════════════════════

@testset "checksums" begin
    @testset "parse_checksums" begin
        sums = Pioneer.parse_checksums("abc123  config.json\ndef456  data.jls\n\n")
        @test sums["config.json"] == "abc123"
        @test sums["data.jls"] == "def456"
        @test length(sums) == 2
    end

    @testset "verify_checksums detects corruption and absence" begin
        mktempdir() do dir
            write(joinpath(dir, "good.txt"), "hello")
            write(joinpath(dir, "bad.txt"), "hello")
            good = bytes2hex(SHA.sha256("hello"))
            sums = Dict("good.txt" => good,
                        "bad.txt" => bytes2hex(SHA.sha256("goodbye")),
                        "missing.txt" => good)
            bad = Pioneer.verify_checksums(dir, sums)
            @test "bad.txt" in bad
            @test "missing.txt" in bad
            @test !("good.txt" in bad)
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Transfer — against a real HTTP server on the loopback interface
# ═══════════════════════════════════════════════════════════════════════════════

# Ports are handed out in sequence rather than drawn at random: two testsets
# drawing the same port let a pooled client connection outlive the server that
# opened it, and the next testset's request was answered with the previous
# testset's data.
_DL_PORT = Ref(rand(20000:39000))

"""Serve `files` (name => bytes) with Range support, honouring `range_mode`.

`redirect = true` puts every path behind a 307 to /r/<path>, the way
huggingface.co bounces a resolve URL to its CDN. That hop is why real downloads
arrived with ~1 KB of "Temporary Redirect. Redirecting to ..." written into
them while every loopback test passed: HTTP.jl runs the response callback once
per hop, so the redirect's own body was being streamed into the file.
"""
function _serve(files::Dict{String, Vector{UInt8}}, range_mode::Symbol;
                redirect::Bool = false)
    port = (_DL_PORT[] += 1)
    server = HTTP.serve!("127.0.0.1", port; verbose = false) do request
        path = HTTP.URI(request.target).path
        if redirect && !startswith(path, "/r/")
            return HTTP.Response(307, ["Location" => "/r" * path];
                                 body = "Temporary Redirect. Redirecting to /r$path")
        end
        name = lstrip(startswith(path, "/r/") ? path[3:end] : path, '/')
        haskey(files, name) || return HTTP.Response(404, "no")
        body = files[name]
        range_header = HTTP.header(request, "Range", "")
        if !isempty(range_header) && range_mode === :honour
            first_byte = parse(Int, split(split(range_header, "=")[2], "-")[1])
            return HTTP.Response(206, ["Content-Range" => "bytes $first_byte-$(length(body)-1)/$(length(body))"];
                                 body = body[(first_byte + 1):end])
        end
        # range_mode === :ignore -> full 200 even though a Range was asked for
        return HTTP.Response(200, body)
    end
    return server, port
end

@testset "hf_download_file" begin
    payload = Vector{UInt8}(rand(UInt8, 5000))
    files = Dict("f.bin" => payload)

    @testset "downloads a whole file" begin
        server, port = _serve(files, :honour)
        try
            mktempdir() do dir
                path = joinpath(dir, "f.bin")
                n = Pioneer.hf_download_file("http://127.0.0.1:$port/f.bin", path, length(payload))
                @test n == length(payload)
                @test read(path) == payload
            end
        finally
            close(server)
        end
    end

    @testset "an existing file is overwritten, never continued" begin
        # Resume was removed deliberately: a partial of unknown provenance can
        # be a wrong prefix, and continuing it yields a plausible size with
        # corrupt contents.
        server, port = _serve(files, :honour)
        try
            mktempdir() do dir
                path = joinpath(dir, "f.bin")
                write(path, payload[1:2000])                  # short partial
                Pioneer.hf_download_file("http://127.0.0.1:$port/f.bin", path, length(payload))
                @test read(path) == payload
                @test filesize(path) == length(payload)       # not 7000

                write(path, vcat(payload, rand(UInt8, 100)))  # over-long junk
                Pioneer.hf_download_file("http://127.0.0.1:$port/f.bin", path, length(payload))
                @test read(path) == payload
            end
        finally
            close(server)
        end
    end

    @testset "a redirect body is never written into the file" begin
        # The bug this exists for: huggingface.co answers a resolve URL with a
        # 307 whose body is prose. Streaming every hop produced a file of
        # plausible size whose first ~1 KB was that prose.
        server, port = _serve(files, :honour; redirect = true)
        try
            mktempdir() do dir
                path = joinpath(dir, "f.bin")
                n = Pioneer.hf_download_file("http://127.0.0.1:$port/f.bin", path, length(payload))
                @test n == length(payload)
                @test filesize(path) == length(payload)
                @test read(path) == payload
            end
        finally
            close(server)
        end
    end

    @testset "a stale partial is replaced across a redirect too" begin
        server, port = _serve(files, :honour; redirect = true)
        try
            mktempdir() do dir
                path = joinpath(dir, "f.bin")
                write(path, payload[1:2000])
                Pioneer.hf_download_file("http://127.0.0.1:$port/f.bin", path, length(payload))
                @test read(path) == payload
                @test filesize(path) == length(payload)
            end
        finally
            close(server)
        end
    end

    @testset "a short body is rejected rather than accepted" begin
        # Nothing else catches truncation for a file absent from SHA256SUMS.
        truncated = Dict("f.bin" => payload[1:1000])
        server, port = _serve(truncated, :honour)
        try
            mktempdir() do dir
                path = joinpath(dir, "f.bin")
                @test_throws ErrorException Pioneer.hf_download_file(
                    "http://127.0.0.1:$port/f.bin", path, length(payload))
            end
        finally
            close(server)
        end
    end
end

@testset "download_library" begin
    config = read(joinpath(DL_FIXTURES, "config.json"))
    fragments = Vector{UInt8}(rand(UInt8, 3000))
    sums = "$(bytes2hex(SHA.sha256(config)))  config.json\n" *
           "$(bytes2hex(SHA.sha256(fragments)))  detailed_fragments.jls\n"
    files = Dict(
        "lib.poin/config.json" => config,
        "lib.poin/detailed_fragments.jls" => fragments,
        "lib.poin/SHA256SUMS" => Vector{UInt8}(sums),
    )
    library = Pioneer.HFLibrary("lib.poin",
        [Pioneer.HFFile("lib.poin/config.json", length(config)),
         Pioneer.HFFile("lib.poin/detailed_fragments.jls", length(fragments)),
         Pioneer.HFFile("lib.poin/SHA256SUMS", length(sums))],
        length(config) + length(fragments) + length(sums))

    @testset "downloads, verifies and renames atomically" begin
        server, port = _serve(files, :honour)
        try
            mktempdir() do dir
                progress = Ref(0)
                path = Pioneer.download_library(library, dir, "";
                    url_for = p -> "http://127.0.0.1:$port/$p",
                    on_progress = (done, total) -> progress[] = done)
                @test path == joinpath(dir, "lib.poin")
                @test isdir(path)
                @test read(joinpath(path, "config.json")) == config
                @test !ispath(joinpath(dir, "lib.poin.partial"))   # staging is gone
                @test progress[] == library.total_bytes
            end
        finally
            close(server)
        end
    end

    @testset "a bad checksum leaves no library directory" begin
        corrupt = copy(files)
        # Flip one bit rather than redrawing at random: @testset reseeds the RNG
        # on entry, so a fresh rand() here reproduces `fragments` exactly and the
        # "corrupt" file would match its own checksum.
        corrupt["lib.poin/detailed_fragments.jls"] =
            vcat(fragments[1:(end - 1)], UInt8[fragments[end] ⊻ 0xff])
        server, port = _serve(corrupt, :honour)
        try
            mktempdir() do dir
                @test_throws ErrorException Pioneer.download_library(library, dir, "";
                    url_for = p -> "http://127.0.0.1:$port/$p")
                # The invariant that matters: nothing that looks like a library.
                @test !ispath(joinpath(dir, "lib.poin"))
                @test isdir(joinpath(dir, "lib.poin.partial"))     # kept for inspection
            end
        finally
            close(server)
        end
    end

    @testset "an existing library is not silently replaced" begin
        server, port = _serve(files, :honour)
        try
            mktempdir() do dir
                mkpath(joinpath(dir, "lib.poin"))
                @test_throws ArgumentError Pioneer.download_library(library, dir, "";
                    url_for = p -> "http://127.0.0.1:$port/$p")
            end
        finally
            close(server)
        end
    end
end

end # top-level testset
