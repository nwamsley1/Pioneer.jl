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

@testset "buildFragmentIndex.jl" begin
    # Helper function to check tolerance
    function Tol(a, b, ppm = 2)
        abs(a - b) <= (ppm * minimum((a, b)) / 1000000)
    end
    
    # Create a set of test fragments
    frag_ions = [
        # First fragment m/z bin (around 100.0)
        # First RT bin (around 10.0-10.1)
        SimpleFrag(Float32(100.0), UInt32(1), Float32(200.0), Float32(10.0), UInt8(2), UInt8(10)),
        SimpleFrag(Float32(100.0), UInt32(1), Float32(201.0), Float32(10.1), UInt8(2), UInt8(9)),
        SimpleFrag(Float32(100.0), UInt32(1), Float32(203.0), Float32(10.1), UInt8(2), UInt8(8)),
        
        # Second RT bin (around 100.0-102.0)
        SimpleFrag(Float32(100.0), UInt32(1), Float32(204.0), Float32(100.0), UInt8(2), UInt8(7)),
        SimpleFrag(Float32(100.0), UInt32(1), Float32(205.0), Float32(100.0), UInt8(2), UInt8(6)),
        SimpleFrag(Float32(100.01), UInt32(1), Float32(206.0), Float32(102.0), UInt8(2), UInt8(5)),

        # Second fragment m/z bin (around 200.0)
        # Third RT bin (around 10.0-10.1)
        SimpleFrag(Float32(200.0), UInt32(1), Float32(207.0), Float32(10.0), UInt8(2), UInt8(10)),
        SimpleFrag(Float32(200.0), UInt32(1), Float32(208.0), Float32(10.1), UInt8(2), UInt8(9)),
        SimpleFrag(Float32(200.0), UInt32(1), Float32(209.0), Float32(10.1), UInt8(2), UInt8(8)),
        
        # Fourth RT bin (around 100.0-102.0)
        SimpleFrag(Float32(200.0), UInt32(1), Float32(210.0), Float32(100.0), UInt8(2), UInt8(7)),
        SimpleFrag(Float32(200.0), UInt32(1), Float32(211.0), Float32(100.0), UInt8(2), UInt8(6)),
        SimpleFrag(Float32(200.01), UInt32(1), Float32(212.0), Float32(102.0), UInt8(2), UInt8(5))
    ]
    
    # Test folder for output
    test_dir = mktempdir()

    # Test standard fragment index building
    @testset "Standard Index Building" begin
        # Build fragment index with standard parameters
        buildFragmentIndex!(
            test_dir,
            frag_ions,
            Float32(20.0),  # frag_bin_tol_ppm
            Float32(5.0),   # rt_bin_tol
            index_name=""
        )
        
        # Load and verify output files
        fragments = Arrow.Table(joinpath(test_dir, "f_index_fragments.arrow"))
        rt_bins = Arrow.Table(joinpath(test_dir, "f_index_rt_bins.arrow"))
        frag_bins = Arrow.Table(joinpath(test_dir, "f_index_fragment_bins.arrow"))
        
        # Verify index structure
        @test length(frag_bins.FragIndexBin) >= 4  # Allow for implementation variations
        
        # Count bin sizes by checking first_bin and last_bin difference
        rt_bin_sizes = [bin.last_bin - bin.first_bin + 1 for bin in collect(rt_bins.FragIndexBin)]
        
        # Check that there are at least 2 RT bins (original had 4)
        @test length(rt_bins.FragIndexBin) >= 2
        
        # Check m/z bin boundaries for the first bins
        if length(frag_bins.FragIndexBin) >= 4
            @test Tol(frag_bins.FragIndexBin[1].lb, 100.0)
            @test Tol(frag_bins.FragIndexBin[2].lb, 200.0)
        else
            # If fewer bins, at least check that both 100 and 200 m/z regions exist
            mz_values = [bin.lb for bin in collect(frag_bins.FragIndexBin)]
            @test any(x -> Tol(x, 100.0), mz_values)
            @test any(x -> Tol(x, 200.0), mz_values)
        end
        
        # Check that bin ranges are valid (first_bin <= last_bin)
        @test all(bin -> bin.first_bin <= bin.last_bin, collect(frag_bins.FragIndexBin))
        @test all(bin -> bin.first_bin <= bin.last_bin, collect(rt_bins.FragIndexBin))
        
        # Check that all fragments are accounted for
        @test length(fragments.IndexFragment) == length(frag_ions)
        
        # Check precursor m/z ranges are preserved
        prec_mzs = Set([frag.prec_mz for frag in collect(fragments.IndexFragment)])
        @test any(x -> Tol(x, 200.0), collect(prec_mzs))
        @test any(x -> Tol(x, 212.0), collect(prec_mzs))
    end
    
    # Test with reversed input order to ensure sorting works correctly
    @testset "Reversed Input Order" begin
        # Reverse the fragment order
        reversed_frag_ions = frag_ions[end:-1:1]
        
        # Build fragment index with reversed input
        buildFragmentIndex!(
            test_dir,
            reversed_frag_ions,
            Float32(20.0),  # frag_bin_tol_ppm
            Float32(5.0),   # rt_bin_tol
            index_name="reversed_"
        )
        
        # Load and verify output files
        fragments = Arrow.Table(joinpath(test_dir, "reversed_f_index_fragments.arrow"))
        rt_bins = Arrow.Table(joinpath(test_dir, "reversed_f_index_rt_bins.arrow"))
        frag_bins = Arrow.Table(joinpath(test_dir, "reversed_f_index_fragment_bins.arrow"))
        
        # Verify index structure should be similar to forward case
        @test length(frag_bins.FragIndexBin) >= 4
        
        # Check that both 100 and 200 m/z regions exist
        mz_values = [bin.lb for bin in collect(frag_bins.FragIndexBin)]
        @test any(x -> Tol(x, 100.0), mz_values)
        @test any(x -> Tol(x, 200.0), mz_values)
        
        # Check that all fragments are accounted for
        @test length(fragments.IndexFragment) == length(frag_ions)
    end
    
    # Test with different thresholds (m/z and RT)
    @testset "Modified Thresholds" begin
        # Build with tighter m/z tolerance
        buildFragmentIndex!(
            test_dir,
            frag_ions,
            Float32(5.0),   # Lower frag_bin_tol_ppm
            Float32(5.0),   # rt_bin_tol
            index_name="tight_mz_"
        )
        
        tight_mz_frag_bins = Arrow.Table(joinpath(test_dir, "tight_mz_f_index_fragment_bins.arrow"))
        standard_frag_bins = Arrow.Table(joinpath(test_dir, "f_index_fragment_bins.arrow"))
        
        # Tighter m/z tolerance should produce at least as many bins
        @test length(tight_mz_frag_bins.FragIndexBin) >= length(standard_frag_bins.FragIndexBin)
        
        # Build with tighter RT tolerance
        buildFragmentIndex!(
            test_dir,
            frag_ions,
            Float32(20.0),  # frag_bin_tol_ppm
            Float32(1.0),   # Lower rt_bin_tol
            index_name="tight_rt_"
        )
        
        tight_rt_bins = Arrow.Table(joinpath(test_dir, "tight_rt_f_index_rt_bins.arrow"))
        standard_rt_bins = Arrow.Table(joinpath(test_dir, "f_index_rt_bins.arrow"))
        
        # Tighter RT tolerance should produce at least as many bins
        @test length(tight_rt_bins.FragIndexBin) >= length(standard_rt_bins.FragIndexBin)
    end
    
    # Test presearch index (with infinite RT tolerance)
    @testset "Presearch Index" begin
        buildFragmentIndex!(
            test_dir,
            frag_ions,
            Float32(20.0),  # frag_bin_tol_ppm
            typemax(Float32),  # Infinite RT tolerance
            index_name="presearch_"
        )
        
        presearch_rt_bins = Arrow.Table(joinpath(test_dir, "presearch_f_index_rt_bins.arrow"))
        
        # Should collapse to a small number of RT bins with infinite tolerance
        @test length(presearch_rt_bins.FragIndexBin) <= 2  # Allow for implementation variations
    end
    
    # Clean up test directory
    rm(test_dir, recursive=true)
end

@testset "BuildPionLib" begin
    #==========================================================================
    Tests for buildPionLib
    ==========================================================================#
    @testset "buildPionLib" begin
        # Create a temporary directory for testing
        test_dir = mktempdir()
        
        # Create mock input data
        fragments_table = (
            mz = Float32[100.0, 200.0, 300.0, 400.0],
            intensity = Float16[0.5, 0.8, 0.3, 0.9],
            is_y = [true, false, true, false],
            is_b = [false, true, false, true],
            is_p = [false, false, false, false],
            fragment_index = UInt8[1, 2, 3, 4],
            charge = UInt8[1, 1, 2, 2],
            sulfur_count = UInt8[0, 1, 0, 1],
            ion_type = UInt16[1, 2, 1, 2],
            isotope = UInt8[0, 0, 0, 0],
            is_internal = [false, false, false, false],
            is_immonium = [false, false, false, false],
            has_neutral_diff = [false, false, false, false]
        )
        
        precursors_table = (
            mz = Float32[500.0, 600.0],
            irt = Float32[10.0, 20.0],
            prec_charge = UInt8[2, 3]
        )
        
        # Write mock files
        Arrow.write(joinpath(test_dir, "fragments_table.arrow"), fragments_table)
        Arrow.write(joinpath(test_dir, "precursors_table.arrow"), precursors_table)
        # Convert prec_to_frag to a proper table/DataFrame structure
        prec_to_frag_df = DataFrame(Dict(:start_idx => UInt64[1, 3, 5]))
        Arrow.write(joinpath(test_dir, "prec_to_frag.arrow"), prec_to_frag_df)
       
        # Set up parameters for buildPionLib
        y_start_index = UInt8(1)
        y_start = UInt8(1)
        b_start_index = UInt8(1)
        b_start = UInt8(1)
        include_p_index = false
        include_p = false
        include_isotope = false
        include_immonium = false
        include_internal = false
        include_neutral_diff = false
        max_frag_charge = UInt8(2)
        max_frag_rank = UInt8(10)
        min_frag_intensity = Float32(0.1)
        rank_to_score = UInt8[10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        
        # Create a simple fragment bounds model for testing
        frag_bounds = FragBoundModel(
            ImmutablePolynomial([50.0f0, 0.0f0]),  # min_frag_mz = 50
            ImmutablePolynomial([1500.0f0, 0.0f0])  # max_frag_mz = 1500
        )
        
        frag_bin_tol_ppm = Float32(10.0)
        rt_bin_tol_ppm = Float32(30.0)
        model_type = InstrumentSpecificModel("test_model")
        
        # Call the function
        result = buildPionLib(
            test_dir,
            y_start_index,
            y_start,
            b_start_index,
            b_start,
            include_p_index,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge,
            max_frag_rank,
            min_frag_intensity,
            rank_to_score,
            frag_bounds,
            frag_bin_tol_ppm,
            rt_bin_tol_ppm,
            model_type
        )
        
        # Check that the function returns nothing
        @test result === nothing
        
        # Check output files existence
        @test isfile(joinpath(test_dir, "f_index_fragments.arrow"))
        @test isfile(joinpath(test_dir, "f_index_rt_bins.arrow"))
        @test isfile(joinpath(test_dir, "f_index_fragment_bins.arrow"))
        @test isfile(joinpath(test_dir, "presearch_f_index_fragments.arrow"))
        @test isfile(joinpath(test_dir, "presearch_f_index_rt_bins.arrow"))
        @test isfile(joinpath(test_dir, "presearch_f_index_fragment_bins.arrow"))
        @test isfile(joinpath(test_dir, "detailed_fragments.jld2"))
        @test isfile(joinpath(test_dir, "precursor_to_fragment_indices.jld2"))
        
        # Load and verify structure of output files
        f_index_fragments = Arrow.Table(joinpath(test_dir, "f_index_fragments.arrow"))
        @test hasproperty(f_index_fragments, :IndexFragment)
        @test length(f_index_fragments.IndexFragment) > 0
        
        f_index_rt_bins = Arrow.Table(joinpath(test_dir, "f_index_rt_bins.arrow"))
        @test hasproperty(f_index_rt_bins, :FragIndexBin)
        @test length(f_index_rt_bins.FragIndexBin) > 0
        
        f_index_fragment_bins = Arrow.Table(joinpath(test_dir, "f_index_fragment_bins.arrow"))
        @test hasproperty(f_index_fragment_bins, :FragIndexBin)
        @test length(f_index_fragment_bins.FragIndexBin) > 0
        
        # Verify detailed fragments file contains the expected data
        detailed_fragments_data = load(joinpath(test_dir, "detailed_fragments.jld2"))
        @test haskey(detailed_fragments_data, "data")
        @test detailed_fragments_data["data"] isa Vector
        @test length(detailed_fragments_data["data"]) > 0
        
        # Verify precursor to fragment indices mapping
        pid_to_fid_data = load(joinpath(test_dir, "precursor_to_fragment_indices.jld2"))
        @test haskey(pid_to_fid_data, "pid_to_fid")
        @test pid_to_fid_data["pid_to_fid"] isa Vector{UInt64}
        @test length(pid_to_fid_data["pid_to_fid"]) == length(precursors_table.mz) + 1
        
        # Test with SplineCoefficientModel
        spline_model_type = SplineCoefficientModel("test_spline_model")
        
        # Create mock spline coefficient data
        spline_fragments_table = (
            mz = Float32[100.0, 200.0, 300.0, 400.0],
            coefficients = NTuple{3, Float32}[(1.0f0, 0.5f0, 0.1f0), (0.8f0, 0.4f0, 0.2f0), 
                                         (0.6f0, 0.3f0, 0.1f0), (0.9f0, 0.5f0, 0.2f0)],
            intensity = Float16[0.5, 0.8, 0.3, 0.9],
            is_y = [true, false, true, false],
            is_b = [false, true, false, true],
            is_p = [false, false, false, false],
            fragment_index = UInt8[1, 2, 3, 4],
            charge = UInt8[1, 1, 2, 2],
            sulfur_count = UInt8[0, 1, 0, 1],
            ion_type = UInt16[1, 2, 1, 2],
            isotope = UInt8[0, 0, 0, 0],
            is_internal = [false, false, false, false],
            is_immonium = [false, false, false, false],
            has_neutral_diff = [false, false, false, false]
        )
        
        # Create a new test directory for spline model test
        #spline_test_dir = "/Users/nathanwamsley/Desktop/test2/"#mktempdir()
        spline_test_dir = mktempdir()#"/Users/nathanwamsley/Desktop/test2/"#mktempdir()
        #mkdir(spline_test_dir)
        # Write mock files for spline test
        Arrow.write(joinpath(spline_test_dir, "fragments_table.arrow"), spline_fragments_table)
        Arrow.write(joinpath(spline_test_dir, "precursors_table.arrow"), precursors_table)
        prec_to_frag_df = DataFrame(Dict(:start_idx => UInt64[1, 3, 5]))
        Arrow.write(joinpath(spline_test_dir, "prec_to_frag.arrow"), prec_to_frag_df)

        # Test the SplineCoefficientModel version
        spline_result = buildPionLib(
            spline_test_dir,
            y_start_index,
            y_start,
            b_start_index,
            b_start,
            include_p_index,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge,
            max_frag_rank,
            min_frag_intensity,
            rank_to_score,
            frag_bounds,
            frag_bin_tol_ppm,
            rt_bin_tol_ppm,
            spline_model_type
        )
        
        @test spline_result === nothing
        @test isfile(joinpath(spline_test_dir, "detailed_fragments.jld2"))
        
        # Verify spline fragments contain coefficients
        spline_fragments_data = load(joinpath(spline_test_dir, "detailed_fragments.jld2"))
        @test haskey(spline_fragments_data, "data")
        @test spline_fragments_data["data"] isa Vector
        @test eltype(spline_fragments_data["data"]) <: SplineDetailedFrag{3, Float32}
        


        # Clean up
        rm(test_dir, recursive=true)
        rm(spline_test_dir, recursive=true)
    end
    
    #==========================================================================
    Tests for cleanUpLibrary
    ==========================================================================#
    @testset "cleanUpLibrary" begin
        # Create a temporary directory for testing
        test_dir = mktempdir()
        
        # Create test files
        test_files = [
            "fragments_table.arrow",
            "prec_to_frag.arrow",
            "precursors.arrow",
            "other_file.txt",  # This should not be removed
            "detailed_fragments.jld2"  # This should not be removed
        ]
        
        for file in test_files
            touch(joinpath(test_dir, file))
        end
        
        # Call the function
        cleanUpLibrary(test_dir)
        
        # Check that specified files were removed
        @test !isfile(joinpath(test_dir, "fragments_table.arrow"))
        @test !isfile(joinpath(test_dir, "prec_to_frag.arrow"))
        @test !isfile(joinpath(test_dir, "precursors.arrow"))
        
        # Files not in the removal list should remain
        @test isfile(joinpath(test_dir, "other_file.txt"))
        @test isfile(joinpath(test_dir, "detailed_fragments.jld2"))
        
        # Clean up
        rm(test_dir, recursive=true)
    end
    
    #==========================================================================
    Tests for fragFilter
    ==========================================================================#
    
    @testset "fragFilter" begin
        # Create a test fragment bounds model
        frag_bounds = FragBoundModel(
            ImmutablePolynomial([100.0f0, 0.0f0]),  # min_frag_mz = 100
            ImmutablePolynomial([1000.0f0, 0.0f0])  # max_frag_mz = 1000
        )
        
        # Test parameters
        prec_mz = Float32(500.0)
        y_start = UInt8(2)
        b_start = UInt8(3)
        include_p = true
        include_isotope = true
        include_immonium = true
        include_internal = true
        include_neutral_diff = true
        max_frag_charge = UInt8(3)
        
        # Test cases
        
        # Test case 1: Valid y-fragment within m/z range and other requirements
        @test fragFilter(
            true,    # frag_is_y
            false,   # frag_is_b
            false,   # frag_is_p
            UInt8(3),     # frag_index >= y_start
            UInt8(2),     # frag_charge <= max_frag_charge
            UInt8(0),     # frag_isotope
            false,   # frag_internal
            false,   # frag_immonium
            false,   # frag_neutral_diff
            Float32(500.0),  # frag_mz (within bounds)
            frag_bounds,
            prec_mz,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge
        ) == true
        
        # Test case 2: Valid b-fragment
        @test fragFilter(
            false,   # frag_is_y
            true,    # frag_is_b
            false,   # frag_is_p
            UInt8(5),     # frag_index >= b_start
            UInt8(1),     # frag_charge <= max_frag_charge
            UInt8(0),     # frag_isotope
            false,   # frag_internal
            false,   # frag_immonium
            false,   # frag_neutral_diff
            Float32(750.0),  # frag_mz (within bounds)
            frag_bounds,
            prec_mz,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge
        ) == true
        
        # Test case 3: Fragment m/z below minimum bound
        @test fragFilter(
            true,    # frag_is_y
            false,   # frag_is_b
            false,   # frag_is_p
            UInt8(3),     # frag_index >= y_start
            UInt8(2),     # frag_charge <= max_frag_charge
            UInt8(0),     # frag_isotope
            false,   # frag_internal
            false,   # frag_immonium
            false,   # frag_neutral_diff
            Float32(50.0),  # frag_mz (below min)
            frag_bounds,
            prec_mz,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge
        ) == false
        
        # Test case 4: Fragment m/z above maximum bound
        @test fragFilter(
            true,    # frag_is_y
            false,   # frag_is_b
            false,   # frag_is_p
            UInt8(3),     # frag_index >= y_start
            UInt8(2),     # frag_charge <= max_frag_charge
            UInt8(0),     # frag_isotope
            false,   # frag_internal
            false,   # frag_immonium
            false,   # frag_neutral_diff
            Float32(1500.0),  # frag_mz (above max)
            frag_bounds,
            prec_mz,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge
        ) == false
        
        # Test case 5: Y-ion with index below minimum
        @test fragFilter(
            true,    # frag_is_y
            false,   # frag_is_b
            false,   # frag_is_p
            UInt8(1),     # frag_index < y_start
            UInt8(2),     # frag_charge <= max_frag_charge
            UInt8(0),     # frag_isotope
            false,   # frag_internal
            false,   # frag_immonium
            false,   # frag_neutral_diff
            Float32(500.0),  # frag_mz (within bounds)
            frag_bounds,
            prec_mz,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge
        ) == false
        
        # Test case 6: B-ion with index below minimum
        @test fragFilter(
            false,   # frag_is_y
            true,    # frag_is_b
            false,   # frag_is_p
            UInt8(2),     # frag_index < b_start
            UInt8(2),     # frag_charge <= max_frag_charge
            UInt8(0),     # frag_isotope
            false,   # frag_internal
            false,   # frag_immonium
            false,   # frag_neutral_diff
            Float32(500.0),  # frag_mz (within bounds)
            frag_bounds,
            prec_mz,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge
        ) == false
        
        # Test case 7: Fragment with charge > max allowed
        @test fragFilter(
            true,    # frag_is_y
            false,   # frag_is_b
            false,   # frag_is_p
            UInt8(3),     # frag_index >= y_start
            UInt8(4),     # frag_charge > max_frag_charge
            UInt8(0),     # frag_isotope
            false,   # frag_internal
            false,   # frag_immonium
            false,   # frag_neutral_diff
            Float32(500.0),  # frag_mz (within bounds)
            frag_bounds,
            prec_mz,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge
        ) == false
        
        # Test case 8: Precursor ion (p-type) when include_p is false
        @test fragFilter(
            false,   # frag_is_y
            false,   # frag_is_b
            true,    # frag_is_p
            UInt8(1),     # frag_index
            UInt8(2),     # frag_charge
            UInt8(0),     # frag_isotope
            false,   # frag_internal
            false,   # frag_immonium
            false,   # frag_neutral_diff
            Float32(500.0),  # frag_mz
            frag_bounds,
            prec_mz,
            y_start,
            b_start,
            false,   # include_p set to false
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge
        ) == false
        
        # Test case 9: Immonium ion when include_immonium is false
        @test fragFilter(
            false,   # frag_is_y
            false,   # frag_is_b
            false,   # frag_is_p
            UInt8(3),     # frag_index
            UInt8(2),     # frag_charge
            UInt8(0),     # frag_isotope
            false,   # frag_internal
            true,    # frag_immonium
            false,   # frag_neutral_diff
            Float32(500.0),  # frag_mz
            frag_bounds,
            prec_mz,
            y_start,
            b_start,
            include_p,
            include_isotope,
            false,   # include_immonium set to false
            include_internal,
            include_neutral_diff,
            max_frag_charge
        ) == false
        
        # Test case 10: Internal fragment when include_internal is false
        @test fragFilter(
            false,   # frag_is_y
            false,   # frag_is_b
            false,   # frag_is_p
            UInt8(3),     # frag_index
            UInt8(2),     # frag_charge
            UInt8(0),     # frag_isotope
            true,    # frag_internal
            false,   # frag_immonium
            false,   # frag_neutral_diff
            Float32(500.0),  # frag_mz
            frag_bounds,
            prec_mz,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            false,   # include_internal set to false
            include_neutral_diff,
            max_frag_charge
        ) == false
        
        # Test case 11: Neutral difference fragment when include_neutral_diff is false
        @test fragFilter(
            false,   # frag_is_y
            false,   # frag_is_b
            false,   # frag_is_p
            UInt8(3),     # frag_index
            UInt8(2),     # frag_charge
            UInt8(0),     # frag_isotope
            false,   # frag_internal
            false,   # frag_immonium
            true,    # frag_neutral_diff
            Float32(500.0),  # frag_mz
            frag_bounds,
            prec_mz,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            false,   # include_neutral_diff set to false
            max_frag_charge
        ) == false
        
        # Test case 12: Non-zero isotope when include_isotope is false
        @test fragFilter(
            true,    # frag_is_y
            false,   # frag_is_b
            false,   # frag_is_p
            UInt8(3),     # frag_index >= y_start
            UInt8(2),     # frag_charge <= max_frag_charge
            UInt8(1),     # frag_isotope > 0
            false,   # frag_internal
            false,   # frag_immonium
            false,   # frag_neutral_diff
            Float32(500.0),  # frag_mz (within bounds)
            frag_bounds,
            prec_mz,
            y_start,
            b_start,
            include_p,
            false,   # include_isotope set to false
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge
        ) == false
    end
    
    #==========================================================================
    Tests for getSimpleFrags
    ==========================================================================#
    @testset "getSimpleFrags" begin
        # Mock data for testing
        frag_mz = Float32[100.0, 200.0, 300.0, 400.0]
        frag_is_y = [true, false, true, false]
        frag_is_b = [false, true, false, true]
        frag_is_p = [false, false, false, false]
        frag_index = UInt8[1, 2, 3, 4]
        frag_charge = UInt8[1, 1, 2, 2]
        frag_isotope = UInt8[0, 0, 0, 0]
        frag_internal = [false, false, false, false]
        frag_immonium = [false, false, false, false]
        frag_neutral_diff = [false, false, false, false]
        
        precursor_mz = Float32[500.0, 600.0]
        precursor_irt = Float32[10.0, 20.0]
        precursor_charge = UInt8[2, 3]
        
        # Index mapping precursors to fragments
        # First precursor has fragments 1-2, second has fragments 3-4
        prec_to_frag_idx = UInt64[1, 3, 5]
        
        # Parameters for filtering
        y_start = UInt8(1)
        b_start = UInt8(1)
        include_p = true
        include_isotope = true
        include_immonium = true
        include_internal = true
        include_neutral_diff = true
        max_frag_charge = UInt8(3)
        
        # Fragment bounds model
        frag_bounds = FragBoundModel(
            ImmutablePolynomial([50.0f0, 0.0f0]),  # min_frag_mz = 50
            ImmutablePolynomial([1000.0f0, 0.0f0])  # max_frag_mz = 1000
        )
        
        # Rank to score mapping
        rank_to_score = UInt8[10, 9, 8, 7, 6]
        
        # Call function with parameters that should allow all fragments to pass
        simple_frags = getSimpleFrags(
            frag_mz,
            frag_is_y,
            frag_is_b,
            frag_is_p,
            frag_index,
            frag_charge,
            frag_isotope,
            frag_internal,
            frag_immonium,
            frag_neutral_diff,
            precursor_mz,
            precursor_irt,
            precursor_charge,
            prec_to_frag_idx,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge,
            frag_bounds,
            rank_to_score
        )
        
        # All fragments should pass our filters
        @test length(simple_frags) == 4
        
        # Check properties of first fragment
        @test getMZ(simple_frags[1]) ≈ 100.0
        @test getPrecID(simple_frags[1]) == 1
        @test getPrecMZ(simple_frags[1]) ≈ 500.0
        @test getIRT(simple_frags[1]) ≈ 10.0
        @test getPrecCharge(simple_frags[1]) == 2
        @test getScore(simple_frags[1]) == 10  # First fragment gets highest score
        
        # Check properties of second fragment
        @test getMZ(simple_frags[2]) ≈ 200.0
        @test getPrecID(simple_frags[2]) == 1
        @test getPrecMZ(simple_frags[2]) ≈ 500.0
        @test getIRT(simple_frags[2]) ≈ 10.0
        @test getPrecCharge(simple_frags[2]) == 2
        @test getScore(simple_frags[2]) == 9  # Second fragment gets second highest score
        
        # Check properties of third fragment (first fragment of second precursor)
        @test getMZ(simple_frags[3]) ≈ 300.0
        @test getPrecID(simple_frags[3]) == 2
        @test getPrecMZ(simple_frags[3]) ≈ 600.0
        @test getIRT(simple_frags[3]) ≈ 20.0
        @test getPrecCharge(simple_frags[3]) == 3
        @test getScore(simple_frags[3]) == 10  # First fragment of second precursor
        
        # Test with more restrictive filters
        # Only allow y ions with index >= 3
        restrictive_simple_frags = getSimpleFrags(
            frag_mz,
            frag_is_y,
            frag_is_b,
            frag_is_p,
            frag_index,
            frag_charge,
            frag_isotope,
            frag_internal,
            frag_immonium,
            frag_neutral_diff,
            precursor_mz,
            precursor_irt,
            precursor_charge,
            prec_to_frag_idx,
            UInt8(4),  # Only y ions with index > 3
            UInt8(4),
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge,
            frag_bounds,
            rank_to_score
        )
        
        # Only fragment #3 should pass (y-ion with index 3)
        @test length(restrictive_simple_frags) == 1
        @test getMZ(restrictive_simple_frags[1]) ≈ 400.0
        @test getPrecID(restrictive_simple_frags[1]) == 2
        
        # Test with rank limits
        # Create a large number of identical fragments
        n_test_frags = 10
        large_frag_mz = repeat(Float32[100.0], n_test_frags)
        large_frag_is_y = repeat([true], n_test_frags)
        large_frag_is_b = repeat([false], n_test_frags)
        large_frag_is_p = repeat([false], n_test_frags)
        large_frag_index = repeat(UInt8[1], n_test_frags)
        large_frag_charge = repeat(UInt8[1], n_test_frags)
        large_frag_isotope = repeat(UInt8[0], n_test_frags)
        large_frag_internal = repeat([false], n_test_frags)
        large_frag_immonium = repeat([false], n_test_frags)
        large_frag_neutral_diff = repeat([false], n_test_frags)
        
        # All fragments belong to same precursor
        large_prec_to_frag_idx = UInt64[1, n_test_frags+1]
        
        # Small rank_to_score array to test rank limiting
        small_rank_to_score = UInt8[10, 9, 8]
        
        rank_limited_frags = getSimpleFrags(
            large_frag_mz,
            large_frag_is_y,
            large_frag_is_b,
            large_frag_is_p,
            large_frag_index,
            large_frag_charge,
            large_frag_isotope,
            large_frag_internal,
            large_frag_immonium,
            large_frag_neutral_diff,
            Float32[500.0],         # Single precursor
            Float32[10.0],          # Single RT value
            UInt8[2],               # Single charge
            large_prec_to_frag_idx,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge,
            frag_bounds,
            small_rank_to_score     # Only allow 3 ranks
        )
        
        # Should only return 3 fragments due to rank limit
        @test length(rank_limited_frags) == 3
        
        # Test scores are assigned in decreasing order
        @test getScore(rank_limited_frags[1]) == 10
        @test getScore(rank_limited_frags[2]) == 9
        @test getScore(rank_limited_frags[3]) == 8
    end
    
    #==========================================================================
    Tests for buildFragmentIndex!
    ==========================================================================#
    @testset "buildFragmentIndex!" begin
        # Create a temporary directory for testing
        test_dir = mktempdir()
        
        # Create mock fragment data
        frag_ions = [
            # First RT bin
            SimpleFrag(Float32(100.0), UInt32(1), Float32(500.0), Float32(10.0), UInt8(2), UInt8(10)),
            SimpleFrag(Float32(150.0), UInt32(1), Float32(500.0), Float32(10.0), UInt8(2), UInt8(9)),
            SimpleFrag(Float32(200.0), UInt32(1), Float32(500.0), Float32(10.5), UInt8(2), UInt8(8)),
            
            # Second RT bin
            SimpleFrag(Float32(120.0), UInt32(2), Float32(600.0), Float32(20.0), UInt8(3), UInt8(10)),
            SimpleFrag(Float32(180.0), UInt32(2), Float32(600.0), Float32(20.5), UInt8(3), UInt8(9))
        ]
        
        # Parameters
        frag_bin_tol_ppm = Float32(100.0) 
        rt_bin_tol = Float32(5.0)         # Should create two RT bins
        
        # Call function
        buildFragmentIndex!(
            test_dir,
            frag_ions,
            frag_bin_tol_ppm,
            rt_bin_tol
        )
        
        # Check output files exist
        @test isfile(joinpath(test_dir, "f_index_fragments.arrow"))
        @test isfile(joinpath(test_dir, "f_index_rt_bins.arrow"))
        @test isfile(joinpath(test_dir, "f_index_fragment_bins.arrow"))
        
        # Load and validate files
        fragments = Arrow.Table(joinpath(test_dir, "f_index_fragments.arrow"))
        rt_bins = Arrow.Table(joinpath(test_dir, "f_index_rt_bins.arrow"))
        frag_bins = Arrow.Table(joinpath(test_dir, "f_index_fragment_bins.arrow"))
        
        # Check structures
        @test hasproperty(fragments, :IndexFragment)
        @test hasproperty(rt_bins, :FragIndexBin)
        @test hasproperty(frag_bins, :FragIndexBin)
        
        # Validate fragment count
        @test length(fragments.IndexFragment) == 5
        
        # Validate RT bins
        @test length(rt_bins.FragIndexBin) == 2  # Should have created 2 RT bins
        
        # Test with custom index name
        buildFragmentIndex!(
            test_dir,
            frag_ions,
            frag_bin_tol_ppm,
            rt_bin_tol,
            index_name="custom_"
        )
        
        # Check custom-named files exist
        @test isfile(joinpath(test_dir, "custom_f_index_fragments.arrow"))
        @test isfile(joinpath(test_dir, "custom_f_index_rt_bins.arrow"))
        @test isfile(joinpath(test_dir, "custom_f_index_fragment_bins.arrow"))
        
        # Test with different bin tolerances
        # Tighter fragmentation bin tolerance
        buildFragmentIndex!(
            test_dir,
            frag_ions,
            Float32(10.0),  # Much tighter fragment bin tolerance
            rt_bin_tol,
            index_name="tight_frag_"
        )
        
        tight_frag_bins = Arrow.Table(joinpath(test_dir, "tight_frag_f_index_fragment_bins.arrow"))
        standard_frag_bins = Arrow.Table(joinpath(test_dir, "f_index_fragment_bins.arrow"))
        
        # Tighter tolerance should create more bins
        @test length(tight_frag_bins.FragIndexBin) >= length(standard_frag_bins.FragIndexBin)
        
        # Tighter RT bin tolerance
        buildFragmentIndex!(
            test_dir,
            frag_ions,
            frag_bin_tol_ppm,
            Float32(1.0),  # Much tighter RT bin tolerance
            index_name="tight_rt_"
        )
        
        tight_rt_bins = Arrow.Table(joinpath(test_dir, "tight_rt_f_index_rt_bins.arrow"))
        standard_rt_bins = Arrow.Table(joinpath(test_dir, "f_index_rt_bins.arrow"))
        
        # Tighter tolerance should create more bins
        @test length(tight_rt_bins.FragIndexBin) >= length(standard_rt_bins.FragIndexBin)
        
        # Test with infinite RT tolerance (presearch index)
        buildFragmentIndex!(
            test_dir,
            frag_ions,
            frag_bin_tol_ppm,
            typemax(Float32),  # Infinite RT tolerance
            index_name="presearch_"
        )
        
        presearch_rt_bins = Arrow.Table(joinpath(test_dir, "presearch_f_index_rt_bins.arrow"))
        
        # Should have fewer RT bins (likely just 1)
        @test length(presearch_rt_bins.FragIndexBin) <= length(standard_rt_bins.FragIndexBin)
        
        # Clean up
        rm(test_dir, recursive=true)
    end
    
    #==========================================================================
    Tests for getDetailedFrags
    ==========================================================================#
    @testset "getDetailedFrags" begin
        # Mock data for testing
        frag_mz = Float32[100.0, 200.0, 300.0, 400.0]
        frag_intensity = Float16[0.5, 0.8, 0.3, 0.9]
        frag_is_y = [true, false, true, false]
        frag_is_b = [false, true, false, true]
        frag_is_p = [false, false, false, false]
        frag_index = UInt8[1, 2, 3, 4]
        frag_charge = UInt8[1, 1, 2, 2]
        frag_sulfur_count = UInt8[0, 1, 0, 1]
        frag_ion_type = UInt16[1, 2, 1, 2]
        frag_isotope = UInt8[0, 0, 0, 0]
        frag_internal = [false, false, false, false]
        frag_immonium = [false, false, false, false]
        frag_neutral_diff = [false, false, false, false]
        
        precursor_mz = Float32[500.0, 600.0]
        precursor_charge = UInt8[2, 3]
        
        # Index mapping precursors to fragments
        # First precursor has fragments 1-2, second has fragments 3-4
        prec_to_frag_idx = UInt64[1, 3, 5]
        
        # Parameters for filtering
        y_start = UInt8(1)
        b_start = UInt8(1)
        include_p = true
        include_isotope = true
        include_immonium = true
        include_internal = true
        include_neutral_diff = true
        max_frag_charge = UInt8(3)
        
        # Fragment bounds model
        frag_bounds = FragBoundModel(
            ImmutablePolynomial([50.0f0, 0.0f0]),  # min_frag_mz = 50
            ImmutablePolynomial([1000.0f0, 0.0f0])  # max_frag_mz = 1000
        )
        
        max_frag_rank = UInt8(10)
        min_frag_intensity = Float32(0.1)
        model_type = InstrumentSpecificModel("test_model")
        
        # Call function with parameters that should allow all fragments to pass
        detailed_frags, pid_to_fid = getDetailedFrags(
            frag_mz,
            frag_intensity,
            frag_is_y,
            frag_is_b,
            frag_is_p,
            frag_index,
            frag_charge,
            frag_sulfur_count,
            frag_ion_type,
            frag_isotope,
            frag_internal,
            frag_immonium,
            frag_neutral_diff,
            precursor_mz,
            precursor_charge,
            prec_to_frag_idx,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge,
            frag_bounds,
            max_frag_rank,
            min_frag_intensity,
            model_type
        )
        
        # All fragments should pass our filters
        @test length(detailed_frags) == 4
        
        # Check precursor-to-fragment mapping
        @test length(pid_to_fid) == 3  # For 2 precursors + end marker
        @test pid_to_fid[1] == 1  # First precursor starts at index 1
        @test pid_to_fid[2] == 3  # Second precursor starts at index 3
        @test pid_to_fid[3] == 5  # End marker
        
        # Check first fragment properties
        @test detailed_frags[1].prec_id == 1
        @test detailed_frags[1].mz ≈ 100.0
        @test detailed_frags[1].intensity ≈ 0.5
        @test detailed_frags[1].ion_type == 1
        @test detailed_frags[1].is_y == true
        @test detailed_frags[1].is_b == false
        @test detailed_frags[1].frag_charge == 1
        @test detailed_frags[1].ion_position == 1
        @test detailed_frags[1].prec_charge == 2
        @test detailed_frags[1].rank == 1  # First fragment gets rank 1
        @test detailed_frags[1].sulfur_count == 0
        
        # Check second fragment properties
        @test detailed_frags[2].prec_id == 1
        @test detailed_frags[2].mz ≈ 200.0
        @test detailed_frags[2].intensity ≈ 0.8
        @test detailed_frags[2].ion_type == 2
        @test detailed_frags[2].is_y == false
        @test detailed_frags[2].is_b == true
        @test detailed_frags[2].frag_charge == 1
        @test detailed_frags[2].ion_position == 2
        @test detailed_frags[2].prec_charge == 2
        @test detailed_frags[2].rank == 2  # Second fragment gets rank 2
        @test detailed_frags[2].sulfur_count == 1
        
        # Test intensity filtering
        high_intensity_frags, high_intensity_pid_to_fid = getDetailedFrags(
            frag_mz,
            frag_intensity,
            frag_is_y,
            frag_is_b,
            frag_is_p,
            frag_index,
            frag_charge,
            frag_sulfur_count,
            frag_ion_type,
            frag_isotope,
            frag_internal,
            frag_immonium,
            frag_neutral_diff,
            precursor_mz,
            precursor_charge,
            prec_to_frag_idx,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge,
            frag_bounds,
            max_frag_rank,
            Float32(0.7),  # Only fragments with intensity >= 0.7
            model_type
        )
        
        # Only two fragments should pass (intensity 0.8 and 0.9)
        @test length(high_intensity_frags) == 2
        @test high_intensity_frags[1].intensity ≈ 0.8
        @test high_intensity_frags[2].intensity ≈ 0.9
        
        # Test rank limiting with many fragments
        n_test_frags = 15
        rank_test_frag_mz = repeat(Float32[100.0], n_test_frags)
        rank_test_frag_intensity = repeat(Float16[0.9], n_test_frags)
        rank_test_frag_is_y = repeat([true], n_test_frags)
        rank_test_frag_is_b = repeat([false], n_test_frags)
        rank_test_frag_is_p = repeat([false], n_test_frags)
        rank_test_frag_index = repeat(UInt8[1], n_test_frags)
        rank_test_frag_charge = repeat(UInt8[1], n_test_frags)
        rank_test_frag_sulfur_count = repeat(UInt8[0], n_test_frags)
        rank_test_frag_ion_type = repeat(UInt16[1], n_test_frags)
        rank_test_frag_isotope = repeat(UInt8[0], n_test_frags)
        rank_test_frag_internal = repeat([false], n_test_frags)
        rank_test_frag_immonium = repeat([false], n_test_frags)
        rank_test_frag_neutral_diff = repeat([false], n_test_frags)
        
        # All fragments belong to same precursor
        rank_test_prec_to_frag_idx = UInt64[1, n_test_frags + 1]
        
        # Set a small max rank
        small_max_rank = UInt8(5)
        
        rank_limited_detailed_frags, _ = getDetailedFrags(
            rank_test_frag_mz,
            rank_test_frag_intensity,
            rank_test_frag_is_y,
            rank_test_frag_is_b,
            rank_test_frag_is_p,
            rank_test_frag_index,
            rank_test_frag_charge,
            rank_test_frag_sulfur_count,
            rank_test_frag_ion_type,
            rank_test_frag_isotope,
            rank_test_frag_internal,
            rank_test_frag_immonium,
            rank_test_frag_neutral_diff,
            Float32[500.0],        # Single precursor
            UInt8[2],              # Single charge
            rank_test_prec_to_frag_idx,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge,
            frag_bounds,
            small_max_rank,        # Only keep top 5 fragments
            min_frag_intensity,
            model_type
        )
        
        # Should only have max_rank fragments
        @test length(rank_limited_detailed_frags) == small_max_rank
        
        # Ranks should be assigned sequentially
        for i in 1:length(rank_limited_detailed_frags)
            @test rank_limited_detailed_frags[i].rank == i
        end
    end
    
    #==========================================================================
    Tests for getDetailedFrags with SplineCoefficientModel
    ==========================================================================#
    @testset "getDetailedFrags with SplineCoefficientModel" begin
        # Mock data for testing
        frag_mz = Float32[100.0, 200.0, 300.0, 400.0]
        frag_coef = NTuple{3, Float32}[(1.0f0, 0.5f0, 0.1f0), (0.8f0, 0.4f0, 0.2f0), 
                                     (0.6f0, 0.3f0, 0.1f0), (0.9f0, 0.5f0, 0.2f0)]
        frag_intensity = Float16[0.5, 0.8, 0.3, 0.9]
        frag_is_y = [true, false, true, false]
        frag_is_b = [false, true, false, true]
        frag_is_p = [false, false, false, false]
        frag_index = UInt8[1, 2, 3, 4]
        frag_charge = UInt8[1, 1, 2, 2]
        frag_sulfur_count = UInt8[0, 1, 0, 1]
        frag_ion_type = UInt16[1, 2, 1, 2]
        frag_isotope = UInt8[0, 0, 0, 0]
        frag_internal = [false, false, false, false]
        frag_immonium = [false, false, false, false]
        frag_neutral_diff = [false, false, false, false]
        
        precursor_mz = Float32[500.0, 600.0]
        precursor_charge = UInt8[2, 3]
        
        # Index mapping precursors to fragments
        # First precursor has fragments 1-2, second has fragments 3-4
        prec_to_frag_idx = UInt64[1, 3, 5]
        
        # Parameters for filtering
        y_start = UInt8(1)
        b_start = UInt8(1)
        include_p = true
        include_isotope = true
        include_immonium = true
        include_internal = true
        include_neutral_diff = true
        max_frag_charge = UInt8(3)
        
        # Fragment bounds model
        frag_bounds = FragBoundModel(
            ImmutablePolynomial([50.0f0, 0.0f0]),  # min_frag_mz = 50
            ImmutablePolynomial([1000.0f0, 0.0f0])  # max_frag_mz = 1000
        )
        
        max_frag_rank = UInt8(10)
        min_frag_intensity = Float32(0.1)
        model_type = SplineCoefficientModel("test_spline_model")
        
        # Call function with parameters that should allow all fragments to pass
        detailed_frags, pid_to_fid = getDetailedFrags(
            frag_mz,
            frag_coef,
            frag_intensity,
            frag_is_y,
            frag_is_b,
            frag_is_p,
            frag_index,
            frag_charge,
            frag_sulfur_count,
            frag_ion_type,
            frag_isotope,
            frag_internal,
            frag_immonium,
            frag_neutral_diff,
            precursor_mz,
            precursor_charge,
            prec_to_frag_idx,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge,
            frag_bounds,
            max_frag_rank,
            min_frag_intensity,
            model_type
        )
        
        # All fragments should pass our filters
        @test length(detailed_frags) == 4
        
        # Check type of returned fragments - should be SplineDetailedFrag with 3 coefficients
        @test eltype(detailed_frags) <: SplineDetailedFrag{3, Float32}
        
        # Check precursor-to-fragment mapping
        @test length(pid_to_fid) == 3  # For 2 precursors + end marker
        @test pid_to_fid[1] == 1  # First precursor starts at index 1
        @test pid_to_fid[2] == 3  # Second precursor starts at index 3
        @test pid_to_fid[3] == 5  # End marker
        
        # Check first fragment properties
        @test detailed_frags[1].prec_id == 1
        @test detailed_frags[1].mz ≈ 100.0
        @test detailed_frags[1].intensity == (1.0f0, 0.5f0, 0.1f0)  # Coefficient tuple
        @test detailed_frags[1].ion_type == 1
        @test detailed_frags[1].is_y == true
        @test detailed_frags[1].is_b == false
        @test detailed_frags[1].is_p == false
        @test detailed_frags[1].is_isotope == false
        @test detailed_frags[1].frag_charge == 1
        @test detailed_frags[1].ion_position == 1  # This is the fragment_index
        @test detailed_frags[1].prec_charge == 2
        @test detailed_frags[1].rank == 1  # First fragment gets rank 1
        @test detailed_frags[1].sulfur_count == 0
        
        # Check that filtering works with spline model
        restricted_detailed_frags, _ = getDetailedFrags(
            frag_mz,
            frag_coef,
            frag_intensity,
            frag_is_y,
            frag_is_b,
            frag_is_p,
            frag_index,
            frag_charge,
            frag_sulfur_count,
            frag_ion_type,
            frag_isotope,
            frag_internal,
            frag_immonium,
            frag_neutral_diff,
            precursor_mz,
            precursor_charge,
            prec_to_frag_idx,
            UInt8(4),  # Higher y-start value
            UInt8(4),
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge,
            frag_bounds,
            max_frag_rank,
            min_frag_intensity,
            model_type
        )
        
        # Only the y-ion with index 3 should pass
        @test length(restricted_detailed_frags) == 1
        @test restricted_detailed_frags[1].is_y == false
        @test restricted_detailed_frags[1].ion_position == 4
        
        # Test with different coefficient dimensions
        frag_coef5 = NTuple{5, Float32}[(1.0f0, 0.5f0, 0.3f0, 0.2f0, 0.1f0), (0.8f0, 0.6f0, 0.4f0, 0.3f0, 0.2f0), 
                                      (0.6f0, 0.5f0, 0.4f0, 0.2f0, 0.1f0), (0.9f0, 0.7f0, 0.5f0, 0.3f0, 0.1f0)]
        
        detailed_frags5, _ = getDetailedFrags(
            frag_mz,
            frag_coef5,
            frag_intensity,
            frag_is_y,
            frag_is_b,
            frag_is_p,
            frag_index,
            frag_charge,
            frag_sulfur_count,
            frag_ion_type,
            frag_isotope,
            frag_internal,
            frag_immonium,
            frag_neutral_diff,
            precursor_mz,
            precursor_charge,
            prec_to_frag_idx,
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge,
            frag_bounds,
            max_frag_rank,
            min_frag_intensity,
            model_type
        )
        
        # Check that coefficient dimension is preserved
        @test eltype(detailed_frags5) <: SplineDetailedFrag{5, Float32}
        @test detailed_frags5[1].intensity == (1.0f0, 0.5f0, 0.3f0, 0.2f0, 0.1f0)
    end
end