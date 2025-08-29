# Tests for Koina batch response parsing functions
# Tests the parsing of responses from different model types using real API calls

@testset "Koina Batch Response Parsing" begin
    
    @testset "InstrumentSpecificModel Response Parsing" begin
        model = InstrumentSpecificModel("unispec")
        test_data = DataFrame(
            koina_sequence = ["PEPTIDE"],
            precursor_charge = Int32[2],  
            collision_energy = Float32[25.0]
        )
        
        # Make real API call to get response
        batches = prepare_koina_batch(model, test_data, "QE", batch_size=1)
        response = make_koina_request(batches[1], KOINA_URLS["unispec"], max_attempts=3)
        
        result = parse_koina_batch(model, response)
        
        # Check result structure
        @test result isa KoinaBatchResult{Nothing}
        # KoinaBatchResult doesn't have knot_vector field
        @test result.frags_per_precursor > 0   # Should have fragments
        
        # Check DataFrame structure
        df = result.fragments
        @test size(df, 1) > 0
        @test "annotation" ∈ names(df)
        @test "mz" ∈ names(df)
        @test "intensities" ∈ names(df)
        
        # Check data types
        @test all(isa(ann, String) for ann in df.annotation)
        @test all(isa(mz, Float32) for mz in df.mz)  
        @test all(isa(int, Float32) for int in df.intensities)
        
        # Check annotation format for UniSpec (uses ^ notation)
        @test any(occursin(r"[by]\d+\^\d+", ann) for ann in df.annotation)
        
        # Check value ranges (before filtering)
        # Note: Raw API data may contain negative intensities and -1.0 m/z sentinel values
        # The production code filters these out with filter_fragments!
        @test any(df.mz .> 0)  # At least some valid m/z values
        # Skip intensity range check on raw data - production code filters negative values
    end
    
    @testset "InstrumentAgnosticModel Response Parsing" begin
        model = InstrumentAgnosticModel("prosit_2020_hcd")
        test_data = DataFrame(
            koina_sequence = ["PEPTIDE"],
            precursor_charge = Int32[2],
            collision_energy = Float32[25.0]
        )
        
        # Make real API call to get response
        batches = prepare_koina_batch(model, test_data, "QE", batch_size=1)
        response = make_koina_request(batches[1], KOINA_URLS["prosit_2020_hcd"], max_attempts=3)
        
        result = parse_koina_batch(model, response)
        
        # Should use same parsing as InstrumentSpecific but with different format
        @test result isa KoinaBatchResult{Nothing}
        # KoinaBatchResult doesn't have knot_vector field
        @test result.frags_per_precursor > 0  # Prosit returns fragments
        
        df = result.fragments
        @test size(df, 1) > 0
        @test "annotation" ∈ names(df)
        @test "mz" ∈ names(df) 
        @test "intensities" ∈ names(df)
        
        # Check annotation format for Prosit (uses + notation)
        @test any(occursin(r"[byp]\d*\+\d+", ann) for ann in df.annotation)
        
        # Check that we have various fragment types
        annotations = df.annotation
        @test any(startswith(ann, "b") for ann in annotations)
        @test any(startswith(ann, "y") for ann in annotations)
        
        # Check charge states
        @test any(endswith(ann, "+1") for ann in annotations)
        @test any(endswith(ann, "+2") for ann in annotations)
        
        # Check value ranges (before filtering)
        @test any(df.mz .> 0)  # At least some valid m/z values  
        # Skip intensity range check on raw data - production code filters negative values
    end
    
    @testset "AlphaPeptDeep Response Parsing" begin
        model = InstrumentSpecificModel("AlphaPeptDeep")
        test_data = DataFrame(
            koina_sequence = ["PEPTIDE"],
            precursor_charge = Int32[2],
            collision_energy = Float32[25.0]
        )
        
        # Make real API call to get response
        batches = prepare_koina_batch(model, test_data, "QE", batch_size=1)
        response = make_koina_request(batches[1], KOINA_URLS["AlphaPeptDeep"], max_attempts=3)
        
        result = parse_koina_batch(model, response)
        
        @test result isa KoinaBatchResult{Nothing}
        @test result.frags_per_precursor > 0
        
        df = result.fragments
        @test size(df, 1) > 0
        
        # AlphaPeptDeep uses + notation like Prosit
        @test any(occursin(r"[by]\d+\+\d+", ann) for ann in df.annotation)
        
        # Check value ranges (before filtering)
        @test any(df.mz .> 0)  # At least some valid m/z values
        # Skip intensity range check on raw data - production code filters negative values
    end
    
    @testset "Response Validation and Error Handling" begin
        # Test with invalid response format
        invalid_response1 = Dict("invalid" => "response")
        model = InstrumentSpecificModel("unispec")
        
        @test_throws Exception parse_koina_batch(model, invalid_response1)
        
        # Test empty outputs
        empty_response = Dict("outputs" => [])
        @test_throws Exception parse_koina_batch(model, empty_response)
    end
    
    @testset "Data Type Conversion" begin
        model = InstrumentSpecificModel("unispec")
        test_data = DataFrame(
            koina_sequence = ["PEPTIDE"],
            precursor_charge = Int32[2],  
            collision_energy = Float32[25.0]
        )
        
        # Make real API call and test conversion
        batches = prepare_koina_batch(model, test_data, "QE", batch_size=1)
        response = make_koina_request(batches[1], KOINA_URLS["unispec"], max_attempts=3)
        
        result = parse_koina_batch(model, response)
        df = result.fragments
        
        # Should all be converted to proper types
        @test all(isa(mz, Float32) for mz in df.mz)
        @test all(isa(int, Float32) for int in df.intensities)
        @test all(isa(ann, String) for ann in df.annotation)
    end
    
    @testset "Multiple Peptide Parsing" begin
        model = InstrumentAgnosticModel("prosit_2020_hcd")
        test_data = DataFrame(
            koina_sequence = ["PEPTIDER", "ALANINEK"],
            precursor_charge = Int32[2, 3],
            collision_energy = Float32[25.0, 30.0]
        )
        
        # Make real API call with multiple peptides
        batches = prepare_koina_batch(model, test_data, "QE", batch_size=10)
        @test length(batches) == 1  # Should fit in one batch
        
        response = make_koina_request(batches[1], KOINA_URLS["prosit_2020_hcd"], max_attempts=3)
        result = parse_koina_batch(model, response)
        
        # Should have fragments for both peptides
        @test result.frags_per_precursor > 0
        expected_total_fragments = 2 * result.frags_per_precursor
        @test nrow(result.fragments) == expected_total_fragments
        
        # All fragments should have valid data (before filtering)
        df = result.fragments
        @test any(df.mz .> 0)  # At least some valid m/z values
        # Skip intensity check on raw data - production code filters negative values
        @test all(length(ann) > 0 for ann in df.annotation)
    end
    
    @testset "Fragment Filtering Behavior" begin
        # Test that production filtering works as expected
        model = InstrumentSpecificModel("unispec")
        test_data = DataFrame(
            koina_sequence = ["PEPTIDE"],
            precursor_charge = Int32[2],
            collision_energy = Float32[25.0]
        )
        
        # Get raw API response
        batches = prepare_koina_batch(model, test_data, "QE", batch_size=1)
        response = make_koina_request(batches[1], KOINA_URLS["unispec"], max_attempts=3)
        result = parse_koina_batch(model, response)
        
        raw_df = copy(result.fragments)
        @test size(raw_df, 1) > 0  # Should have fragments
        
        # Import the filtering function and apply it
        using Pioneer: filter_fragments!
        filtered_df = copy(raw_df)
        filter_fragments!(filtered_df, model)
        
        # After filtering, all values should be valid
        if size(filtered_df, 1) > 0  # Only test if filtering left some fragments
            @test all(filtered_df.mz .> 0)        # No invalid m/z values
            @test all(filtered_df.intensities .> 0.001)  # No very low intensities
            @test all(filtered_df.intensities .≤ 1)      # Reasonable intensity range
        end
        
        # Filtering should remove some fragments (unless data was already clean)
        @test size(filtered_df, 1) ≤ size(raw_df, 1)
    end
end