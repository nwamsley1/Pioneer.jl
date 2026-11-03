# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# A built library carries copies of its input FASTAs and a content-addressed
# record of each (see src/Routines/BuildSpecLib/utils/provenance.jl).

using SHA, JSON

@testset "Library FASTA provenance" begin
    mktempdir() do dir
        a = joinpath(dir, "a"); b = joinpath(dir, "b"); mkpath(a); mkpath(b)
        fa = joinpath(a, "proteome.fasta"); write(fa, ">sp|P1|X_HUMAN x OS=Homo sapiens GN=X\nPEPTIDEK\n")
        fb = joinpath(b, "proteome.fasta"); write(fb, ">sp|P2|Y_HUMAN y OS=Homo sapiens GN=Y\nKPEPTIDE\n")
        sha_a = bytes2hex(SHA.sha256(read(fa)))

        @testset "record without a sidecar" begin
            r = Pioneer.fasta_provenance_record(fa, "HUMAN")
            @test r["name"] == "HUMAN"
            @test r["file"] == "proteome.fasta"
            @test r["bytes"] == filesize(fa)
            @test r["sha256"] == sha_a
            @test !haskey(r, "uniprot_release")
        end

        @testset "sidecar fields are merged" begin
            write(fa * ".provenance.json", JSON.json(Dict(
                "uniprot_release" => "2026_03", "source_url" => "https://example.org/p.fasta",
                "sha256" => sha_a)))
            r = Pioneer.fasta_provenance_record(fa, "HUMAN")
            @test r["uniprot_release"] == "2026_03"
            @test r["source_url"] == "https://example.org/p.fasta"
            @test r["sha256"] == sha_a
        end

        @testset "a sidecar for different bytes is refused" begin
            write(fb * ".provenance.json", JSON.json(Dict("sha256" => sha_a)))
            @test_throws ArgumentError Pioneer.fasta_provenance_record(fb, "OTHER")
            rm(fb * ".provenance.json")
        end

        @testset "missing FASTA" begin
            @test_throws ArgumentError Pioneer.fasta_provenance_record(joinpath(dir, "nope.fasta"), "X")
        end

        @testset "bundle copies FASTAs and sidecars, keeping same-named inputs apart" begin
            lib = joinpath(dir, "lib.poin")
            params = Dict{String, Any}("fasta_paths" => [fa, fb], "fasta_names" => ["HUMAN", "OTHER"],
                                       "library_params" => Dict{String, Any}())
            Pioneer.stamp_build_provenance!(params, lib)
            recs = params["fasta_provenance"]
            @test [r["bundled_file"] for r in recs] == ["fasta/proteome.fasta", "fasta/2_proteome.fasta"]
            @test read(joinpath(lib, "fasta", "proteome.fasta")) == read(fa)
            @test read(joinpath(lib, "fasta", "2_proteome.fasta")) == read(fb)
            @test isfile(joinpath(lib, "fasta", "proteome.fasta.provenance.json"))
            @test !isfile(joinpath(lib, "fasta", "2_proteome.fasta.provenance.json"))
            @test recs[1]["uniprot_release"] == "2026_03"
            @test params["library_params"]["prediction_model"] == "altimeter"
            @test params["pioneer_version"] == Pioneer.get_pioneer_version()
            @test occursin(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$", params["build_date"])
            # JSON round trip, as config.json is written.
            @test JSON.parse(JSON.json(params))["fasta_provenance"][2]["sha256"] == bytes2hex(SHA.sha256(read(fb)))
        end
    end
end
