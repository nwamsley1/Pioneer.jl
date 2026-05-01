# Tests for src/Routines/BuildSpecLib/utils/parse_mods.jl. The functions here
# parse the (pos,aa,mod_name) modification strings that drive precursor and
# fragment m/z calculation; bugs surface as wrong masses everywhere
# downstream. Targets the difficult / edge-case paths:
# - empty / one-mod / multi-mod strings
# - n- and c-terminal mods (regex must accept lowercase n/c)
# - unknown mod names (silently skipped, not crash)
# - shared structural-mod names across iso channels
# - reverseSequence preserves first/last positions and relocates internal mods
# - sulfur-count modifications (positive and negative deltas)

using Pioneer: get_aa_masses!, get_structural_mod_masses!, getIsoModMasses!,
               get_sulfur_counts!, get_fragment_mz, get_precursor_mz,
               reverseSequence
using Pioneer: AA_to_mass, H2O, PROTON

@testset "parse_mods" begin

    @testset "get_aa_masses!" begin
        @testset "all 20 standard amino acids" begin
            seq = "ACDEFGHIKLMNPQRSTVWY"
            buf = zeros(Float32, length(seq))
            get_aa_masses!(buf, seq)
            for (i, aa) in enumerate(seq)
                @test buf[i] == Float32(AA_to_mass[aa])
            end
        end

        @testset "buffer is reset to zero before fill" begin
            buf = ones(Float32, 5)
            get_aa_masses!(buf, "AAAAA")
            @test all(buf .== Float32(AA_to_mass['A']))
            # Now fill with shorter sequence (only 3 chars touched in iter)
            buf2 = ones(Float32, 5)
            get_aa_masses!(buf2, "GGG")
            @test buf2[1] == buf2[2] == buf2[3] == Float32(AA_to_mass['G'])
            @test buf2[4] == 0.0f0
            @test buf2[5] == 0.0f0
        end

        @testset "accepts SubString" begin
            full = "PEPTIDE"
            sub = SubString(full, 1, 4)   # "PEPT"
            buf = zeros(Float32, 4)
            get_aa_masses!(buf, sub)
            @test buf[1] == Float32(AA_to_mass['P'])
            @test buf[4] == Float32(AA_to_mass['T'])
        end
    end

    @testset "get_structural_mod_masses!" begin
        mod_to_mass = Dict{String, Float32}(
            "Phospho" => 79.96633f0,
            "Acetyl"  => 42.01057f0,
            "Carbamidomethyl" => 57.021464f0,
            "Oxidation" => 15.99491f0,
        )

        @testset "empty mod string is a no-op" begin
            buf = ones(Float32, 7)   # initialized to 1 to verify reset
            get_structural_mod_masses!(buf, "", mod_to_mass)
            @test all(buf .== 0.0f0)
        end

        @testset "single mod placed at correct index" begin
            buf = zeros(Float32, 7)
            get_structural_mod_masses!(buf, "(4,T,Phospho)", mod_to_mass)
            @test buf[4] == 79.96633f0
            @test count(==(0.0f0), buf) == 6
        end

        @testset "multiple mods including n-terminal" begin
            buf = zeros(Float32, 7)
            get_structural_mod_masses!(
                buf,
                "(1,n,Acetyl)(3,C,Carbamidomethyl)(7,M,Oxidation)",
                mod_to_mass)
            @test buf[1] == 42.01057f0
            @test buf[3] == 57.021464f0
            @test buf[7] == 15.99491f0
            @test buf[2] == buf[4] == buf[5] == buf[6] == 0.0f0
        end

        @testset "unknown mod names are silently skipped" begin
            buf = zeros(Float32, 5)
            get_structural_mod_masses!(buf,
                "(2,P,NotAMod)(4,T,Phospho)", mod_to_mass)
            @test buf[2] == 0.0f0   # NotAMod skipped
            @test buf[4] == 79.96633f0
        end
    end

    @testset "getIsoModMasses!" begin
        iso_mods_dict = Dict{String, Dict{String, Float32}}(
            "exTag1" => Dict("d0" => 0.0f0, "d4" => 4.0f0,  "d8" => 8.0f0),
            "exTag3" => Dict("d0" => 0.0f0, "d4" => 4.0f0),
        )

        @testset "empty structural OR isotopic string is a no-op" begin
            buf = ones(Float32, 8)   # verify reset to zero
            getIsoModMasses!(buf, "", "(1, d4)", iso_mods_dict)
            @test all(buf .== 0.0f0)

            buf = ones(Float32, 8)
            getIsoModMasses!(buf, "(3,K,exTag1)", "", iso_mods_dict)
            @test all(buf .== 0.0f0)
        end

        @testset "iso index points into structural mod sequence" begin
            buf = zeros(Float32, 7)
            getIsoModMasses!(buf,
                "(1,I,exTag1)(5,L,exTag3)(7,K,exTag1)",
                "(3, d4)",   # third mod (position 7) gets d4
                iso_mods_dict)
            @test buf[7] == 4.0f0
            @test all(buf[1:6] .== 0.0f0)
        end

        @testset "n-terminal and c-terminal structural mods are indexable" begin
            buf = zeros(Float32, 16)
            getIsoModMasses!(buf,
                "(1,n,exTag1)(5,L,exTag3)(16,c,exTag1)",
                "(1, d0)(2, d4)(3, d8)",
                iso_mods_dict)
            @test buf[1]  == 0.0f0   # exTag1 d0 = 0
            @test buf[5]  == 4.0f0   # exTag3 d4 = 4
            @test buf[16] == 8.0f0   # exTag1 d8 = 8
        end

        @testset "out-of-range iso index is silently skipped" begin
            buf = zeros(Float32, 5)
            # Only 1 structural mod, but iso refers to index 3
            getIsoModMasses!(buf,
                "(2,K,exTag1)",
                "(3, d4)",
                iso_mods_dict)
            @test all(buf .== 0.0f0)

            # Index 0 also out of range
            buf2 = zeros(Float32, 5)
            getIsoModMasses!(buf2,
                "(2,K,exTag1)",
                "(0, d4)",
                iso_mods_dict)
            @test all(buf2 .== 0.0f0)
        end

        @testset "unknown channel for known mod silently skipped" begin
            buf = zeros(Float32, 5)
            getIsoModMasses!(buf,
                "(2,K,exTag1)",
                "(1, d99)",   # d99 not in iso_mods_dict["exTag1"]
                iso_mods_dict)
            @test all(buf .== 0.0f0)
        end

        @testset "structural mod with name not in iso dict is excluded from indexing" begin
            # If a structural mod's name has no entry in iso_mods_dict, it should
            # not consume an iso-index slot.
            buf = zeros(Float32, 7)
            # Two structural mods: first not in iso dict, second is.
            # Iso index 1 should map to the FIRST mod that *is* in the iso dict
            # (since the unknown one was filtered out).
            getIsoModMasses!(buf,
                "(2,K,Carbamidomethyl)(5,K,exTag1)",
                "(1, d4)",
                iso_mods_dict)
            @test buf[5] == 4.0f0
            @test buf[2] == 0.0f0
        end
    end

    @testset "get_sulfur_counts!" begin
        mods_to_sulfur = Dict{String, Int8}(
            "Carbamidomethyl" => Int8(0),  # adds C to side chain but on existing S, net 0
            "DiSulfide"       => Int8(2),
            "Loss"            => Int8(-1),
        )

        @testset "sequence-only sulfurs (C and M)" begin
            buf = zeros(Int8, 7)
            get_sulfur_counts!(buf, "ACMSTCM", "", mods_to_sulfur)
            @test buf == Int8[0, 1, 1, 0, 0, 1, 1]
        end

        @testset "no C/M, no mods → all zero" begin
            buf = zeros(Int8, 5)
            get_sulfur_counts!(buf, "GAGAG", "", mods_to_sulfur)
            @test all(buf .== Int8(0))
        end

        @testset "modification adds sulfur" begin
            buf = zeros(Int8, 5)
            get_sulfur_counts!(buf, "GAGAG", "(2,A,DiSulfide)", mods_to_sulfur)
            @test buf == Int8[0, 2, 0, 0, 0]
        end

        @testset "modification removes sulfur (negative delta)" begin
            buf = zeros(Int8, 4)
            # M at position 2 contributes 1; Loss subtracts 1 → net 0.
            get_sulfur_counts!(buf, "AMAA", "(2,M,Loss)", mods_to_sulfur)
            @test buf == Int8[0, 0, 0, 0]
        end

        @testset "buffer is reset to zero before fill" begin
            buf = fill(Int8(99), 5)
            get_sulfur_counts!(buf, "GAGAG", "", mods_to_sulfur)
            @test all(buf .== Int8(0))
        end
    end

    @testset "get_fragment_mz / get_precursor_mz" begin
        # Use a known peptide: PEPTIDE, no mods.
        # Sum of monoisotopic AA residue masses for PEPTIDE:
        seq = "PEPTIDE"
        aa_masses = zeros(Float32, length(seq))
        get_aa_masses!(aa_masses, seq)
        struct_mods = zeros(Float32, length(seq))
        iso_mods = zeros(Float32, length(seq))
        residue_sum = sum(aa_masses)

        @testset "precursor mz at +1 charge" begin
            mz = get_precursor_mz(length(seq), UInt8(1),
                                  aa_masses, struct_mods, iso_mods)
            expected = residue_sum + Float32(H2O) + Float32(PROTON) * 1.0f0
            @test mz ≈ expected
        end

        @testset "precursor mz scales with charge as (M + zH) / z" begin
            mz1 = get_precursor_mz(length(seq), UInt8(1),
                                   aa_masses, struct_mods, iso_mods)
            mz2 = get_precursor_mz(length(seq), UInt8(2),
                                   aa_masses, struct_mods, iso_mods)
            mz3 = get_precursor_mz(length(seq), UInt8(3),
                                   aa_masses, struct_mods, iso_mods)
            # M+H = mz1; M+2H = 2*mz2; M+3H = 3*mz3 — all should differ by H mass each
            M = mz1 - Float32(PROTON)
            @test 2 * mz2 - 2 * Float32(PROTON) ≈ M atol=1f-3
            @test 3 * mz3 - 3 * Float32(PROTON) ≈ M atol=1f-3
        end

        @testset "y-ion includes H2O, b-ion does not" begin
            # y3 covers residues 5-7 (IDE), b3 covers residues 1-3 (PEP)
            y3 = get_fragment_mz(5, 7, 'y', UInt8(1),
                                 aa_masses, struct_mods, iso_mods)
            b3 = get_fragment_mz(1, 3, 'b', UInt8(1),
                                 aa_masses, struct_mods, iso_mods)
            ide_residue = sum(aa_masses[5:7])
            pep_residue = sum(aa_masses[1:3])
            @test y3 ≈ ide_residue + Float32(H2O) + Float32(PROTON)
            @test b3 ≈ pep_residue + Float32(PROTON)
        end

        @testset "fragment charge halves the m/z difference" begin
            y3_z1 = get_fragment_mz(5, 7, 'y', UInt8(1),
                                    aa_masses, struct_mods, iso_mods)
            y3_z2 = get_fragment_mz(5, 7, 'y', UInt8(2),
                                    aa_masses, struct_mods, iso_mods)
            # M = y3_z1 - PROTON; y3_z2 = (M + 2*PROTON) / 2
            M = y3_z1 - Float32(PROTON)
            @test y3_z2 ≈ (M + 2 * Float32(PROTON)) / 2
        end
    end

    @testset "reverseSequence" begin
        @testset "no mods: keeps last residue, reverses prefix" begin
            rev_seq, rev_mods = reverseSequence("PEPTIDE", "")
            @test rev_seq == "DITPEPE"
            @test rev_mods == ""
        end

        @testset "internal mod relocates by seq_length - pos" begin
            rev_seq, rev_mods = reverseSequence("PEPTIDE", "(4,T,Phospho)")
            @test rev_seq == "DITPEPE"
            # pos 4 → new pos 7-4=3; rev_seq[3] = 'T'
            @test rev_mods == "(3,T,Phospho)"
        end

        @testset "c-terminal mod (pos == seq_length) stays put" begin
            rev_seq, rev_mods = reverseSequence("PEPTIDE", "(7,E,Acetyl)")
            @test rev_seq == "DITPEPE"
            @test rev_mods == "(7,E,Acetyl)"
        end

        @testset "n-terminal mod stays at pos 1" begin
            rev_seq, rev_mods = reverseSequence("PEPTIDE", "(1,n,Acetyl)")
            @test rev_seq == "DITPEPE"
            @test rev_mods == "(1,n,Acetyl)"
        end

        @testset "mix of n-terminal, internal, c-terminal mods" begin
            rev_seq, rev_mods = reverseSequence(
                "PEPTIDE",
                "(1,n,nterm)(1,P,internal)(4,T,Phospho)(7,E,Acetyl)(7,c,cterm)")
            @test rev_seq == "DITPEPE"
            # n-terminal stays at 1; (1,P,internal) → new pos 7-1=6, rev_seq[6]='P';
            # (4,T,Phospho) → new pos 3; (7,E,Acetyl) stays at 7; (7,c,cterm) stays at 7
            # Sorted by new position: 1, 3, 6, 7, 7
            @test rev_mods == "(1,n,nterm)(3,T,Phospho)(6,P,internal)(7,E,Acetyl)(7,c,cterm)"
        end

        @testset "round-trip: reverse twice over rest gives original prefix" begin
            # Reversing the prefix twice (last AA fixed both times) returns to the
            # original sequence.
            seq = "ACDEFGHIK"
            once_seq, _ = reverseSequence(seq, "")
            twice_seq, _ = reverseSequence(once_seq, "")
            @test twice_seq == seq
        end
    end

end
