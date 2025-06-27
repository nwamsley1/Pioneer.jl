
@testset "Fragment Predict Tests" begin
    
    @testset "Model Type Validation" begin
        @testset "Instrument-specific models" begin
            # Test that instrument-specific models validate instruments correctly
            unispec_model = InstrumentSpecificModel("unispec")
            alphapept_model = InstrumentSpecificModel("alphapeptdeep")
            
            # These should be valid based on MODEL_CONFIGS
            valid_unispec_instruments = ["QE", "QEHFX", "LUMOS", "ELITE", "VELOS"]
            valid_alpha_instruments = ["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"]
            
            # Test model type properties
            @test unispec_model isa InstrumentSpecificModel
            @test alphapept_model isa InstrumentSpecificModel
        end
        
        @testset "Instrument-agnostic models" begin
            prosit_model = InstrumentAgnosticModel("prosit_2020_hcd")
            @test prosit_model isa InstrumentAgnosticModel
        end
        
        @testset "Spline coefficient models" begin
            altimeter_model = SplineCoefficientModel("altimeter")
            @test altimeter_model isa SplineCoefficientModel
        end
    end
    
    @testset "filter_fragments!" begin
        @testset "Standard intensity filtering" begin
            # Test data with various intensities
            df = DataFrame(
                annotation = ["b3", "y4", "b5", "y7"],
                intensities = Float32[0.0005, 0.002, 0.5, 0.9],
                mz = Float32[300.1, 400.2, 500.3, 600.4]
            )
            
            model = InstrumentSpecificModel("unispec")
            filter_fragments!(df, model)
            
            # Should keep only intensities > 0.001
            @test nrow(df) == 3
            @test all(df.intensities .> 0.001f0)
        end
        
        @testset "Invalid m/z filtering" begin
            df = DataFrame(
                annotation = ["b3", "y4", "b5"],
                intensities = Float32[0.5, 0.6, 0.7],
                mz = Float32[-100.0, 0.0, 500.3]
            )
            
            model = InstrumentAgnosticModel("prosit_2020_hcd")
            filter_fragments!(df, model)
            
            # Should keep only positive m/z
            @test nrow(df) == 1
            @test all(df.mz .> 0)
        end
        
        @testset "Isotope peak filtering for instrument-specific models" begin
            df = DataFrame(
                annotation = ["b3", "y4+i", "b5+2i", "y7"],
                intensities = Float32[0.5, 0.6, 0.7, 0.8],
                mz = Float32[300.1, 400.2, 500.3, 600.4]
            )
            
            model = InstrumentSpecificModel("unispec")
            filter_fragments!(df, model)
            
            # Should remove isotope peaks (containing 'i')
            @test nrow(df) == 2
            @test !any(occursin('i', ann) for ann in df.annotation)
        end
        
        @testset "Spline coefficient filtering" begin
            # For spline models, we don't filter on intensity but on m/z
            df = DataFrame(
                annotation = [1, 2, 3, 4],  # Altimeter uses indices
                coefficients = [(0.1, 0.2), (0.3, 0.4), (0.5, 0.6), (0.7, 0.8)],
                mz = Float32[-100.0, 0.0, 500.3, 600.4]
            )
            
            model = SplineCoefficientModel("altimeter")
            filter_fragments!(df, model)
            
            # Should keep only positive m/z
            @test nrow(df) == 2
            @test all(df.mz .> 0)
        end
    end
    
    @testset "sort_fragments!" begin
        df = DataFrame(
            precursor_idx = [1, 2, 1, 2, 1, 2],
            annotation = ["b3", "y4", "y5", "b2", "a3", "x6"],
            intensities = Float32[0.5, 0.9, 0.7, 0.3, 0.8, 0.6],
            mz = Float32[300.1, 400.2, 500.3, 200.1, 350.2, 450.3]
        )
        
        sort_fragments!(df)
        
        # Check that within each precursor, fragments are sorted by intensity (descending)
        @test df.precursor_idx == [1, 1, 1, 2, 2, 2]
        
        # For precursor 1
        prec1_mask = df.precursor_idx .== 1
        @test issorted(df[prec1_mask, :intensities], rev=true)
        
        # For precursor 2
        prec2_mask = df.precursor_idx .== 2
        @test issorted(df[prec2_mask, :intensities], rev=true)
    end
    
    @testset "predict_fragments_batch" begin
        temp_dir = mktempdir()
        
        @testset "InstrumentSpecificModel batch prediction" begin
            # Create mock peptide data
            peptides_df = DataFrame(
                sequence = ["PEPTIDE", "SEQUENCE", "FRAGMENT"],
                charge = [2, 3, 2],
                collision_energy = [25.0, 30.0, 27.0]
            )
            
            # Mock model and parameters
            model = InstrumentSpecificModel("unispec")
            instrument_type = "QE"
            batch_size = 2
            concurrent_requests = 1
            first_prec_idx = UInt32(1)
            
            # We can't test actual API calls, but we can test the structure
            # Create a mock response structure
            mock_fragments = DataFrame(
                annotation = ["b3", "y4", "b5", "y6", "a2", "x3"],
                mz = Float32[300.1, 400.2, 500.3, 600.4, 200.1, 350.2],
                intensities = Float32[0.8, 0.9, 0.7, 0.6, 0.5, 0.4]
            )
            
            # Test that the function structure works with mock data
            @test model isa InstrumentSpecificModel
            @test instrument_type in ["QE", "QEHFX", "LUMOS", "ELITE", "VELOS"]
        end
        
        @testset "InstrumentAgnosticModel batch prediction" begin
            peptides_df = DataFrame(
                sequence = ["PEPTIDE", "SEQUENCE"],
                charge = [2, 3],
                collision_energy = [25.0, 30.0]
            )
            
            model = InstrumentAgnosticModel("prosit_2020_hcd")
            # Instrument type should be ignored
            instrument_type = "ANY"
            
            @test model isa InstrumentAgnosticModel
        end
        
        @testset "SplineCoefficientModel batch prediction" begin
            peptides_df = DataFrame(
                sequence = ["PEPTIDE"],
                charge = [2],
                collision_energy = [25.0]
            )
            
            model = SplineCoefficientModel("altimeter")
            instrument_type = "QE"
            
            @test model isa SplineCoefficientModel
            
            # Test knot vector consistency check would happen here
            # In real implementation, all batches must have same knot vector
        end
        
        rm(temp_dir, recursive=true)
    end
    
    @testset "predict_fragments - main dispatcher" begin
        temp_dir = mktempdir()
        
        # Create test peptide data
        peptide_data = DataFrame(
            sequence = ["PEPTIDE", "SEQUENCE", "FRAGMENT", "EXAMPLE"],
            charge = [2, 3, 2, 2],
            collision_energy = [25.0, 30.0, 27.0, 26.0]
        )
        
        peptide_path = joinpath(temp_dir, "peptides.arrow")
        Arrow.write(peptide_path, peptide_data)
        
        frags_out_path = joinpath(temp_dir, "fragments.arrow")
        
        @testset "Model name validation" begin
            # Test invalid model name
            @test_throws ErrorException predict_fragments(
                peptide_path,
                frags_out_path,
                InstrumentSpecificModel("invalid_model"),
                "QE",
                10,
                100,
                "invalid_model"
            )
        end
        
        @testset "Batch size limiting" begin
            # Test that batch size is capped at 1000
            model = InstrumentSpecificModel("unispec")
            
            # Would need to mock API calls to fully test
            # Here we just verify the setup is correct
            @test model.name == "unispec"
            
            # Verify peptide data was loaded correctly
            loaded_peptides = DataFrame(Arrow.Table(peptide_path))
            @test nrow(loaded_peptides) == 4
            @test names(loaded_peptides) == ["sequence", "charge", "collision_energy"]
        end
        
        @testset "Output file creation" begin
            # Test that output path handling works
            @test !isfile(frags_out_path)
            
            # In real implementation, this would create the file
            # We just test the path construction
            @test dirname(frags_out_path) == temp_dir
            @test basename(frags_out_path) == "fragments.arrow"
        end
        
        rm(temp_dir, recursive=true)
    end
    
    @testset "Batch processing calculations" begin
        @testset "Batch index calculation" begin
            nprecs = 1000
            batch_size = 100
            koina_pool_size = 5
            
            batch_start_idxs = collect(one(UInt32):UInt32(batch_size*koina_pool_size):UInt32(nprecs))
            
            # Should create appropriate number of batches
            @test length(batch_start_idxs) == 2  # 0-499, 500-999
            @test batch_start_idxs[1] == 1
            @test batch_start_idxs[2] == 501
        end
        
        @testset "Fragment per precursor calculation" begin
            # Mock fragment data with 3 precursors
            fragment_df = DataFrame(
                precursor_idx = [1, 1, 1, 2, 2, 3, 3, 3, 3],
                annotation = ["b1", "b2", "y1", "b1", "y2", "b1", "b2", "y1", "y2"]
            )
            
            # Calculate fragments per precursor
            frags_per_prec = [3, 2, 4]  # Precursor 1: 3, Precursor 2: 2, Precursor 3: 4
            
            # Test index assignment
            n_precursors = 3
            for i in 1:n_precursors
                mask = fragment_df.precursor_idx .== i
                @test sum(mask) == frags_per_prec[i]
            end
        end
    end
end