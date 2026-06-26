using Test
using Arrow
using DataFrames

import Pioneer
using Pioneer: FragBoundModel, SplineCoefficientModel, buildPionLib, serialize_to_jls

@testset "buildPionLib separate detailed/index fragment filters" begin
    function write_index_filter_fixture!(test_dir::String)
        # Spline coefficient layout: 4-tuple of Float32 per fragment.
        coef_zero = (0.0f0, 0.0f0, 0.0f0, 0.0f0)
        fragments_table = (
            mz = Float32[120.0, 140.0, 210.0, 230.0, 500.0],
            coefficients = NTuple{4, Float32}[coef_zero, coef_zero, coef_zero, coef_zero, coef_zero],
            intensity = Float16[0.9, 0.8, 0.7, 0.6, 0.5],
            is_y = Bool[true, true, false, false, false],
            is_b = Bool[false, false, true, true, false],
            is_p = Bool[false, false, false, false, true],
            fragment_index = UInt8[2, 4, 1, 3, 0],
            charge = UInt8[1, 1, 1, 1, 1],
            sulfur_count = UInt8[0, 0, 0, 0, 0],
            ion_type = UInt16[1, 1, 2, 2, 3],
            isotope = UInt8[0, 0, 0, 0, 0],
            is_internal = Bool[false, false, false, false, false],
            is_immonium = Bool[false, false, false, false, false],
            has_neutral_diff = Bool[false, false, false, false, false],
        )

        precursors_table = (
            sequence = ["PEPTIDE"],
            accession_numbers = [["P1"]],
            mz = Float32[600.0],
            irt = Float32[10.0],
            prec_charge = UInt8[2],
            length = UInt8[7],
            structural_mods = [""],
            isotopic_mods = [""],
            collision_energy = Float32[0.0],
            is_decoy = Bool[false],
            entrapment_group_id = UInt8[0],
            missed_cleavages = UInt8[0],
            sulfur_count = UInt8[0],
            base_pep_id = UInt32[1],
            pair_id = UInt32[1],
        )

        proteins_table = DataFrame(
            accession = ["P1"],
            gene_name = ["GENE1"],
            protein_name = ["Protein 1"],
        )

        Arrow.write(joinpath(test_dir, "fragments_table.arrow"), fragments_table)
        Arrow.write(joinpath(test_dir, "precursors_table.arrow"), precursors_table)
        Arrow.write(joinpath(test_dir, "proteins_table.arrow"), proteins_table)
        Arrow.write(joinpath(test_dir, "prec_to_frag.arrow"), DataFrame(start_idx = UInt64[1, 6]))
        # SplineCoefficientModel buildPionLib loads spline_knots.jls; emit a minimal one.
        serialize_to_jls(joinpath(test_dir, "spline_knots.jls"), Float32[0.0, 0.5, 1.0])
    end

    test_dir = mktempdir()
    try
        write_index_filter_fixture!(test_dir)

        frag_bounds = FragBoundModel(
            Pioneer.ImmutablePolynomial(Float32[0.0, 0.0]),
            Pioneer.ImmutablePolynomial(Float32[1000.0, 0.0]),
        )

        @test buildPionLib(
            test_dir,
            UInt8(4),       # y_start_index: index keeps y4, drops y2
            UInt8(2),       # y_start: detailed fragments keep y2 and y4
            UInt8(3),       # b_start_index: index keeps b3, drops b1
            UInt8(1),       # b_start: detailed fragments keep b1 and b3
            false,          # include_p_index: index drops precursor ion
            true,           # include_p: detailed fragments keep precursor ion
            false,
            false,
            false,
            true,
            UInt8(2),
            UInt8(10),
            Float32(1000),
            0.0f0,
            frag_bounds,
            0.0f0,
            3.0f0,
            SplineCoefficientModel("altimeter"),
        ) === nothing

        detailed_frags = Pioneer.deserialize_from_jls(joinpath(test_dir, "detailed_fragments.jls"))
        @test length(detailed_frags) == 5
        @test sort(Float32[Pioneer.getMz(f) for f in detailed_frags]) == Float32[120.0, 140.0, 210.0, 230.0, 500.0]

        partitioned_index = Pioneer.deserialize_from_jls(joinpath(test_dir, "partitioned_fragment_index.jls"))
        partition = Pioneer.getPartition(partitioned_index, 1)
        frag_bins = Pioneer.getFragBins(partition)
        indexed_fragments = Pioneer.getFragments(partition)

        @test length(indexed_fragments) == 2
        @test frag_bins.lows == Float32[140.0, 230.0]
        @test frag_bins.highs[1:2] == Float32[140.0, 230.0]
        @test UInt8[Pioneer.getScore(f) for f in indexed_fragments] == UInt8[1, 2]
    finally
        safe_rmdir(test_dir)
    end
end

@testset "fragment index filter validation" begin
    @test_throws ArgumentError Pioneer.validate_fragment_index_filters(
        UInt8(1), UInt8(2), UInt8(1), UInt8(1), false, false)
    @test_throws ArgumentError Pioneer.validate_fragment_index_filters(
        UInt8(2), UInt8(2), UInt8(1), UInt8(2), false, false)
    @test_throws ArgumentError Pioneer.validate_fragment_index_filters(
        UInt8(2), UInt8(2), UInt8(2), UInt8(2), true, false)
    @test Pioneer.validate_fragment_index_filters(
        UInt8(4), UInt8(2), UInt8(3), UInt8(1), false, true) === nothing
end
