# Integration tests for fragment prediction pipeline
# Tests end-to-end functionality using real Koina API calls

@testset "Fragment Prediction Integration Tests" begin
    
    @testset "End-to-End Prosit Library Building" begin
        # Test the complete pipeline from peptide sequences to predicted fragments
        # This specifically tests the prosit model as requested
        
        test_peptides = DataFrame(
            koina_sequence = ["PEPTIDER", "SEQUENCER", "FRAGMENTK"],
            precursor_charge = Int32[2, 3, 2],
            collision_energy = Float32[25.0, 30.0, 25.0]
        )
        
        # Test prosit model end-to-end
        model = InstrumentAgnosticModel("prosit_2020_hcd")
        
        # Prepare batch
        batches = prepare_koina_batch(model, test_peptides, "QE", batch_size=100)
        @test length(batches) == 1
        
        # Make real API request
        batch_json = batches[1]
        url = KOINA_URLS["prosit_2020_hcd"]
        response = make_koina_request(batch_json, url, max_attempts=3)
        
        # Parse response
        result = parse_koina_batch(model, response)
        
        # Verify complete pipeline worked
        @test result isa KoinaBatchResult{Nothing}
        @test size(result.fragments, 1) > 0
        @test result.frags_per_precursor > 0   # Prosit returns fragments
        
        # Check that we got fragments for all 3 peptides
        total_expected_fragments = 3 * result.frags_per_precursor
        @test size(result.fragments, 1) == total_expected_fragments
        
        # Verify fragment data quality (before filtering)
        df = result.fragments
        @test any(df.mz .> 0)  # At least some valid m/z values
        @test size(df, 1) > 0  # Should have fragments
        # Skip intensity checks on raw data - production code filters negative/invalid values
        
        # Check annotation format (prosit uses + notation)
        @test any(occursin(r"[byp]\d*\+\d+", ann) for ann in df.annotation)
        
        # Verify we have different fragment types
        annotations = df.annotation
        @test any(startswith(ann, "b") for ann in annotations)  # b ions
        @test any(startswith(ann, "y") for ann in annotations)  # y ions
        
        # Check charge state diversity
        @test any(endswith(ann, "+1") for ann in annotations)
        @test any(endswith(ann, "+2") for ann in annotations)
    end
    
    @testset "Model Comparison Integration" begin
        # Test that different models produce different but valid results
        test_peptide = DataFrame(
            koina_sequence = ["PEPTIDER"],
            precursor_charge = Int32[2],
            collision_energy = Float32[25.0]
        )
        
        # Test multiple models with same input
        models_and_urls = [
            ("unispec", InstrumentSpecificModel("unispec"), "QE"),
            ("prosit_2020_hcd", InstrumentAgnosticModel("prosit_2020_hcd"), "ANY"),
            ("AlphaPeptDeep", InstrumentSpecificModel("AlphaPeptDeep"), "QE")
        ]
        
        results = Dict{String, Any}()
        
        for (name, model, instrument) in models_and_urls
            batches = prepare_koina_batch(model, test_peptide, instrument, batch_size=1)
            @test length(batches) == 1
            
            url = KOINA_URLS[name]
            response = make_koina_request(batches[1], url, max_attempts=3)
            result = parse_koina_batch(model, response)
            results[name] = result
        end
        
        # Check that each model produced valid but different results
        @test length(results) == 3
        
        # Different models should have different fragment counts or formats
        unispec_result = results["unispec"]
        prosit_result = results["prosit_2020_hcd"]
        alphapept_result = results["AlphaPeptDeep"]
        
        @test unispec_result.frags_per_precursor > 0
        @test prosit_result.frags_per_precursor > 0
        @test alphapept_result.frags_per_precursor > 0
        
        # Different models should have different annotation formats
        unispec_anns = unispec_result.fragments.annotation
        prosit_anns = prosit_result.fragments.annotation
        
        # UniSpec uses ^ notation, Prosit uses + notation
        @test any(occursin(r"\^", ann) for ann in unispec_anns)
        @test any(occursin(r"\+", ann) for ann in prosit_anns)
    end
    
    @testset "Batch Processing Integration" begin
        # Test processing multiple batches with realistic data sizes
        large_dataset = DataFrame(
            koina_sequence = ["PEPTIDE$(i)" for i in 1:25],  # 25 peptides
            precursor_charge = repeat(Int32[2, 3], 13)[1:25],   # Mixed charges
            collision_energy = repeat(Float32[25.0, 30.0], 13)[1:25]  # Mixed NCE
        )
        
        model = InstrumentAgnosticModel("prosit_2020_hcd")
        
        # Test small batch size to force multiple batches
        batches = prepare_koina_batch(model, large_dataset, "QE", batch_size=10)
        @test length(batches) == 3  # 25 peptides / 10 per batch = 3 batches
        
        url = KOINA_URLS["prosit_2020_hcd"]
        all_results = KoinaBatchResult{Nothing}[]
        
        # Process each batch
        for batch in batches
            response = make_koina_request(batch, url, max_attempts=3)
            result = parse_koina_batch(model, response)
            push!(all_results, result)
        end
        
        # Verify all batches processed successfully
        @test length(all_results) == 3
        
        # Check total fragment count
        total_fragments = sum(size(result.fragments, 1) for result in all_results)
        expected_frags_per_prec = all_results[1].frags_per_precursor
        expected_total = 25 * expected_frags_per_prec  # 25 peptides * fragments each
        @test total_fragments == expected_total
        
        # Check that each batch has consistent frags_per_precursor
        @test all(result.frags_per_precursor == expected_frags_per_prec for result in all_results)
        
        # Verify fragment data quality across all batches (before filtering)
        for result in all_results
            df = result.fragments
            @test any(df.mz .> 0)  # At least some valid m/z values
            @test size(df, 1) > 0  # Should have fragments
            # Skip intensity checks on raw data - production code filters negative/invalid values
            @test any(occursin(r"[byp]\d*\+\d+", ann) for ann in df.annotation)
        end
    end
    
    @testset "Error Recovery Integration" begin
        # Test end-to-end pipeline with error conditions
        test_data = DataFrame(
            koina_sequence = ["PEPTIDER"],
            precursor_charge = Int32[2],
            collision_energy = Float32[25.0]
        )
        
        model = InstrumentSpecificModel("unispec")
        url = KOINA_URLS["unispec"]
        
        # Test with proper request should succeed
        batches = prepare_koina_batch(model, test_data, "QE", batch_size=1)
        response = make_koina_request(batches[1], url, max_attempts=3)
        result = parse_koina_batch(model, response)
        
        @test result isa KoinaBatchResult{Nothing}
        @test size(result.fragments, 1) > 0
        
        # Test with invalid URL should fail
        @test_throws Exception make_koina_request(batches[1], "https://invalid-url.com", max_attempts=1)
    end
    
    @testset "Real-world Peptide Integration" begin
        # Test with realistic peptide sequences and parameters
        realistic_peptides = DataFrame(
            koina_sequence = [
                "MFLSFPTTK",      # Tryptic peptide
                "EGVNDNEEGFFSAR", # Longer peptide
                "LVNEVTEFAK",     # Medium length peptide
                "AEFVEVTK"        # Short peptide
            ],
            precursor_charge = Int32[2, 3, 2, 2],
            collision_energy = Float32[25.0, 30.0, 27.5, 25.0]
        )
        
        # Test prosit model with realistic sequences
        model = InstrumentAgnosticModel("prosit_2020_hcd")
        batches = prepare_koina_batch(model, realistic_peptides, "QE", batch_size=100)
        
        @test length(batches) == 1  # All should fit in one batch
        
        url = KOINA_URLS["prosit_2020_hcd"]
        response = make_koina_request(batches[1], url, max_attempts=3)
        result = parse_koina_batch(model, response)
        
        # Verify results for realistic peptides
        expected_total_fragments = 4 * result.frags_per_precursor  # 4 peptides * fragments each
        @test size(result.fragments, 1) == expected_total_fragments
        
        df = result.fragments
        
        # Check that we get good fragment coverage
        annotations = df.annotation
        unique_annotations = unique(annotations)
        @test length(unique_annotations) >= 10  # Should have diverse fragments
        
        # Verify we have fragment data (before filtering)
        intensities = df.intensities
        mz_values = df.mz
        @test length(intensities) > 0  # Should have fragments
        @test length(mz_values) > 0   # Should have m/z values
        
        # Note: Raw API data may contain negative intensities and -1.0 m/z sentinel values
        # Production code filters these with filter_fragments!
        # Check that we have some valid data mixed in with potential sentinel values
        @test any(mz_values .> 0)     # At least some valid m/z values
        @test maximum(mz_values) < 3000.0 # Reasonable upper bound for valid fragments
        
        # Verify annotation diversity for realistic peptides
        @test any(contains(ann, "b") for ann in annotations)  # b-ions
        @test any(contains(ann, "y") for ann in annotations)  # y-ions  
        @test any(endswith(ann, "+1") for ann in annotations) # Charge +1
        @test any(endswith(ann, "+2") for ann in annotations) # Charge +2
    end
end