using Test
using Pioneer

# Access internal functions
const compute_diann_mutation_shifts = Pioneer.compute_diann_mutation_shifts
const create_diann_decoy_fragments = Pioneer.create_diann_decoy_fragments
const extract_mod_positions = Pioneer.extract_mod_positions
const DIANN_MUTATION_TABLE = Pioneer.DIANN_MUTATION_TABLE
const AA_to_mass = Pioneer.AA_to_mass
const DetailedFrag = Pioneer.DetailedFrag

@testset "DIA-NN Style Decoy Generation" begin

    @testset "Mutation table completeness" begin
        # All 20 standard amino acids should be in the table
        standard_aas = "GAVLIFMPWSCTYHKRQEND"
        for aa in standard_aas
            @test haskey(DIANN_MUTATION_TABLE, aa)
            # Mutated AA should also be a valid amino acid
            @test haskey(AA_to_mass, DIANN_MUTATION_TABLE[aa])
            # Mutation should produce a DIFFERENT amino acid
            @test DIANN_MUTATION_TABLE[aa] != aa
        end
    end

    @testset "Mutation shifts for known peptides" begin
        # PEPTIDEK: P(1) E(2) P(3) T(4) I(5) D(6) E(7) K(8)
        # N-terminal mutation: position 2 = 'E', E→Q, shift = mass(Q) - mass(E)
        # C-terminal mutation: position 7 = 'E' (n-1=7), E→Q, shift = mass(Q) - mass(E)
        n_shift, c_shift = compute_diann_mutation_shifts("PEPTIDEK")

        expected_n_shift = AA_to_mass['Q'] - AA_to_mass['E']  # E→Q at position 2
        expected_c_shift = AA_to_mass['Q'] - AA_to_mass['E']  # E→Q at position 7

        @test n_shift ≈ expected_n_shift atol=1e-8
        @test c_shift ≈ expected_c_shift atol=1e-8

        # Manually verify: E=129.04259, Q=128.05858 → shift = -0.98401
        @test n_shift ≈ -0.98401 atol=1e-4
        @test c_shift ≈ -0.98401 atol=1e-4

        println("PEPTIDEK:")
        println("  N-terminal: pos 2, E→Q, shift = $(round(n_shift, digits=5)) Da")
        println("  C-terminal: pos 7, E→Q, shift = $(round(c_shift, digits=5)) Da")
    end

    @testset "Mutation shifts for ACDEFGHIK" begin
        # ACDEFGHIK: 9 AAs
        # N-terminal: position 2 = 'C', C→S, shift = mass(S) - mass(C)
        # C-terminal: position 8 = 'I' (n-1=8), I→V, shift = mass(V) - mass(I)
        n_shift, c_shift = compute_diann_mutation_shifts("ACDEFGHIK")

        expected_n_shift = AA_to_mass['S'] - AA_to_mass['C']  # C→S
        expected_c_shift = AA_to_mass['V'] - AA_to_mass['I']  # I→V

        @test n_shift ≈ expected_n_shift atol=1e-8
        @test c_shift ≈ expected_c_shift atol=1e-8

        println("ACDEFGHIK:")
        println("  N-terminal: pos 2, C→S, shift = $(round(n_shift, digits=5)) Da")
        println("  C-terminal: pos 8, I→V, shift = $(round(c_shift, digits=5)) Da")
    end

    @testset "Position selection with modifications" begin
        # If position 2 is modified, should fall to min(3, n-1)
        # Sequence ABCDEFK (7 AAs), mod at position 2
        n_shift, c_shift = compute_diann_mutation_shifts("GADEFGK", Set{Int}([2]))

        # Position 2 ('A') is modified → skip to min(3, 6) = 3 → 'D'
        # D→Q
        expected_n_shift = AA_to_mass['Q'] - AA_to_mass['D']
        @test n_shift ≈ expected_n_shift atol=1e-8

        println("GADEFGK (mod at pos 2):")
        println("  N-terminal: pos 3 (skipped 2), D→Q, shift = $(round(n_shift, digits=5)) Da")
    end

    @testset "Short peptide (3 AAs)" begin
        # Minimum: 3 AAs (n=3, m=n-2=1 in DIA-NN 0-indexed)
        # N-terminal: pos 2 (only interior position)
        # C-terminal: pos 2 (n-1=2, same position!)
        n_shift, c_shift = compute_diann_mutation_shifts("GAK")

        # N-terminal: pos 2 = 'A', A→L
        expected_n = AA_to_mass['L'] - AA_to_mass['A']
        # C-terminal: pos 2 = 'A' (n-1=2), A→L
        expected_c = AA_to_mass['L'] - AA_to_mass['A']

        @test n_shift ≈ expected_n atol=1e-8
        @test c_shift ≈ expected_c atol=1e-8

        println("GAK (3 AAs):")
        println("  N-terminal: pos 2, A→L, shift = $(round(n_shift, digits=5)) Da")
        println("  C-terminal: pos 2, A→L, shift = $(round(c_shift, digits=5)) Da")
    end

    @testset "Too-short peptide" begin
        # < 3 AAs: returns (0, 0)
        n_shift, c_shift = compute_diann_mutation_shifts("AK")
        @test n_shift == 0.0
        @test c_shift == 0.0
    end

    @testset "Fragment m/z shifting" begin
        # Create synthetic target fragments for a simple peptide
        # Precursor 1: 3 y-ions and 2 b-ions
        target_frags = DetailedFrag{Float32}[
            # Y-ions (should shift by c_shift/charge)
            DetailedFrag{Float32}(UInt32(1), 500.0f0, Float16(1.0), UInt16(0),
                true, false, false, false, UInt8(1), UInt8(3), UInt8(2), UInt8(1), UInt8(0)),
            DetailedFrag{Float32}(UInt32(1), 400.0f0, Float16(0.8), UInt16(0),
                true, false, false, false, UInt8(1), UInt8(2), UInt8(2), UInt8(2), UInt8(0)),
            DetailedFrag{Float32}(UInt32(1), 250.0f0, Float16(0.5), UInt16(0),
                true, false, false, false, UInt8(2), UInt8(3), UInt8(2), UInt8(3), UInt8(0)),
            # B-ions (should shift by n_shift/charge)
            DetailedFrag{Float32}(UInt32(1), 300.0f0, Float16(0.6), UInt16(0),
                false, true, false, false, UInt8(1), UInt8(2), UInt8(2), UInt8(4), UInt8(0)),
            DetailedFrag{Float32}(UInt32(1), 200.0f0, Float16(0.3), UInt16(0),
                false, true, false, false, UInt8(1), UInt8(3), UInt8(2), UInt8(5), UInt8(0)),
        ]

        n_shift = -0.98401  # E→Q
        c_shift = 13.03164  # D→Q

        decoy_frags = create_diann_decoy_fragments(
            target_frags, 1:5, n_shift, c_shift, UInt32(2)
        )

        @test length(decoy_frags) == 5

        # Check Y-ions shifted by c_shift/charge
        @test decoy_frags[1].mz ≈ 500.0f0 + Float32(c_shift / 1) atol=0.001  # charge 1
        @test decoy_frags[2].mz ≈ 400.0f0 + Float32(c_shift / 1) atol=0.001  # charge 1
        @test decoy_frags[3].mz ≈ 250.0f0 + Float32(c_shift / 2) atol=0.001  # charge 2

        # Check B-ions shifted by n_shift/charge
        @test decoy_frags[4].mz ≈ 300.0f0 + Float32(n_shift / 1) atol=0.001  # charge 1
        @test decoy_frags[5].mz ≈ 200.0f0 + Float32(n_shift / 1) atol=0.001  # charge 1

        # Intensities should be UNCHANGED
        @test decoy_frags[1].intensity == target_frags[1].intensity
        @test decoy_frags[2].intensity == target_frags[2].intensity
        @test decoy_frags[3].intensity == target_frags[3].intensity
        @test decoy_frags[4].intensity == target_frags[4].intensity
        @test decoy_frags[5].intensity == target_frags[5].intensity

        # Ion types should be preserved
        @test decoy_frags[1].is_y == true
        @test decoy_frags[4].is_b == true

        # Precursor ID should be the decoy's
        @test all(f.prec_id == UInt32(2) for f in decoy_frags)

        println("\nFragment shifting (n_shift=$(round(n_shift, digits=4)), c_shift=$(round(c_shift, digits=4))):")
        for (i, (t, d)) in enumerate(zip(target_frags, decoy_frags))
            ion = t.is_y ? "y" : "b"
            println("  $ion-ion charge=$(t.frag_charge): $(round(t.mz, digits=3)) → $(round(d.mz, digits=3))  " *
                    "(Δ=$(round(d.mz - t.mz, digits=4))), intensity: $(t.intensity) → $(d.intensity)")
        end
    end

    @testset "extract_mod_positions" begin
        @test extract_mod_positions(missing) == Set{Int}()
        @test extract_mod_positions("") == Set{Int}()
        @test extract_mod_positions("2:Carbamidomethyl") == Set{Int}([2])
        @test extract_mod_positions("1:Acetyl,5:Oxidation") == Set{Int}([1, 5])
        @test extract_mod_positions("3:Phospho;7:Oxidation") == Set{Int}([3, 7])
    end

    @testset "Verify against DIA-NN mass values" begin
        # Cross-check our AA masses match DIA-NN's (diann.cpp lines 1375-1395)
        diann_masses = Dict(
            'G' => 57.021464, 'A' => 71.037114, 'V' => 99.068414,
            'L' => 113.084064, 'I' => 113.084064, 'F' => 147.068414,
            'M' => 131.040485, 'P' => 97.052764, 'W' => 186.079313,
            'S' => 87.032028, 'C' => 103.009185, 'T' => 101.047679,
            'Y' => 163.063329, 'H' => 137.058912, 'K' => 128.094963,
            'R' => 156.101111, 'Q' => 128.058578, 'E' => 129.042593,
            'N' => 114.042927, 'D' => 115.026943
        )

        println("\nAA mass comparison (Pioneer vs DIA-NN):")
        max_diff = 0.0
        for (aa, diann_mass) in sort(collect(diann_masses))
            pioneer_mass = AA_to_mass[aa]
            diff = abs(pioneer_mass - diann_mass)
            max_diff = max(max_diff, diff)
            if diff > 0.001
                println("  WARNING: $aa: Pioneer=$(pioneer_mass), DIA-NN=$(diann_mass), diff=$(diff)")
            end
        end
        println("  Max mass difference: $(round(max_diff, digits=6)) Da")
        @test max_diff < 0.01  # Should be very close
    end
end
