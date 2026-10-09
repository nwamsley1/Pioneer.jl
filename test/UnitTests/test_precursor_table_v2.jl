# Compact precursor table (schema 2): packed sequences, mod entries, side tables, schema-1 getters.

using Test, Random, Arrow
using Pioneer

@testset "precursor table schema 2" begin
    @testset "packed sequences" begin
        rng = Xoshiro(3)
        seqs = vcat(["A", "Y", "PEP", "PEPA", "PEPTIDEK", "UOXBZJ", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"],
                    [String(rand(rng, collect("ACDEFGHIKLMNPQRSTVWY"), rand(rng, 1:40))) for _ in 1:2000])
        packed = Pioneer.pack_sequence.(seqs)
        @test all(Pioneer.unpack_sequence(p) == s for (p, s) in zip(packed, seqs))
        @test all(length(p) == cld(5 * length(s), 8) for (p, s) in zip(packed, seqs))
        @test sortperm(packed) == sortperm(seqs)                      # byte order is String order
        @test_throws ErrorException Pioneer.pack_sequence("PEPtIDE")
    end
    @testset "mod entries" begin
        names = Dict{String, UInt8}()
        seq = "MPEPCTIDEK"
        for mods in ("", "(1,M,Unimod:35)", "(1,n,Unimod:1)(1,M,Unimod:35)(5,C,Unimod:4)", "(10,c,Unimod:2)")
            e = Pioneer.encode_mods(mods, seq, names)
            name_list = Vector{String}(undef, length(names)); for (k, v) in names; name_list[v] = k; end
            @test Pioneer.decode_mods(e, Pioneer.pack_sequence(seq), name_list) == mods
        end
        @test ismissing(Pioneer.encode_mods(missing, seq, names))
        @test_throws ErrorException Pioneer.encode_mods("(2,M,Unimod:35)", seq, names)     # residue 2 is P
        @test_throws ErrorException Pioneer.encode_mods("(64,M,Unimod:35)", repeat("M", 70), names)
    end
    @testset "conversion and getters" begin
        dir = mktempdir()
        tbl = (accession_numbers = ["P1;P2", "P1;P2", "P3", "P3", "P1;P2"], sequence = ["PEPTIDEK", "PEPTIDEK", "MAGIC", "CMAGI", "KEDITPEP"],
               start_idx = [UInt32[3, 9], UInt32[3, 9], UInt32[1], UInt32[1], UInt32[3, 9]],
               structural_mods = Union{Missing, String}["", "", "(1,M,Unimod:35)", "(1,C,Unimod:4)(2,M,Unimod:35)", missing],
               proteome_identifiers = ["HUMAN;HUMAN", "HUMAN;HUMAN", "YEAST", "YEAST", "HUMAN;HUMAN"],
               prec_charge = UInt8[2, 3, 2, 2, 2], is_decoy = [false, false, false, true, true],
               mz = Float32[1, 2, 3, 4, 5], irt = Float32[5, 4, 3, 2, 1], base_pep_id = UInt32[1, 1, 3, 3, 1],
               num_variable_modifications = UInt8[0, 0, 1, 1, 0])
        Arrow.write(joinpath(dir, "precursors_table.arrow"), tbl)
        v1 = Pioneer.SetPrecursors(dir)
        ref = (collect(Pioneer.getSequence(v1)), collect(Pioneer.getStructuralMods(v1)), collect(Pioneer.getAccessionNumbers(v1)),
               collect(Pioneer.getProteomeIdentifiers(v1)), map(collect, Pioneer.getStartIdx(v1)), copy(v1.pid_to_cv_fold))
        Pioneer.convert_precursor_table_v2(dir)
        t = Arrow.Table(joinpath(dir, "precursors_table.arrow"))
        @test Pioneer.precursor_schema(t) == 2
        @test !hasproperty(t, :sequence) && !hasproperty(t, :accession_numbers) && hasproperty(t, :sequence_packed)
        @test_throws ArgumentError Pioneer.SetPrecursors(t)                 # side tables required
        v2 = Pioneer.SetPrecursors(dir)
        @test isequal((collect(Pioneer.getSequence(v2)), collect(Pioneer.getStructuralMods(v2)),
                       collect(Pioneer.getAccessionNumbers(v2)), collect(Pioneer.getProteomeIdentifiers(v2)),
                       map(collect, Pioneer.getStartIdx(v2)), v2.pid_to_cv_fold), ref)
        @test_throws ErrorException Pioneer.convert_precursor_table_v2(dir)  # already schema 2
        bad = mktempdir()                                                     # one base_pep_id, two accession sets
        Arrow.write(joinpath(bad, "precursors_table.arrow"), merge(tbl, (base_pep_id = UInt32[1, 1, 1, 1, 1],)))
        @test_throws ErrorException Pioneer.convert_precursor_table_v2(bad)
        @test Pioneer.precursor_schema(Arrow.Table(joinpath(bad, "precursors_table.arrow"))) == 1   # left unchanged
    end
end
