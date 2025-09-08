# Tests for Koina batch preparation functions
# Tests the request formatting for different model types

@testset "Koina Batch Preparation" begin
    
    @testset "InstrumentSpecificModel Preparation" begin
        model = InstrumentSpecificModel("unispec")
        test_data = DataFrame(
            koina_sequence = ["PEPTIDE", "SEQUENCE"],
            precursor_charge = Int32[2, 3],  
            collision_energy = Float32[25.0, 30.0]
        )
        
        # Test with valid instrument
        batches = prepare_koina_batch(model, test_data, "QE", batch_size=1)
        @test length(batches) == 2  # Two batches for batch_size=1
        
        # Parse and validate first batch
        batch1 = JSON.parse(batches[1])
        @test haskey(batch1, "id")
        @test haskey(batch1, "inputs")
        @test length(batch1["inputs"]) == 4  # sequence, charge, nce, instrument_types
        
        # Check input structure
        inputs = batch1["inputs"]
        input_names = [inp["name"] for inp in inputs]
        @test "peptide_sequences" in input_names
        @test "precursor_charges" in input_names  
        @test "collision_energies" in input_names
        
        # Find specific inputs
        seq_input = inputs[findfirst(inp -> inp["name"] == "peptide_sequences", inputs)]
        charge_input = inputs[findfirst(inp -> inp["name"] == "precursor_charges", inputs)]
        nce_input = inputs[findfirst(inp -> inp["name"] == "collision_energies", inputs)]
        
        @test seq_input["shape"] == [1, 1]  # batch_size=1
        @test seq_input["datatype"] == "BYTES"
        @test seq_input["data"] == ["PEPTIDE"]
        
        @test charge_input["shape"] == [1, 1]
        @test charge_input["datatype"] == "INT32" 
        @test charge_input["data"] == [2]
        
        @test nce_input["shape"] == [1, 1]
        @test nce_input["datatype"] == "FP32"
        @test nce_input["data"] == [25.0f0]
        
        # Parse and validate second batch
        batch2 = JSON.parse(batches[2])
        inputs2 = batch2["inputs"]
        seq_input2 = inputs2[findfirst(inp -> inp["name"] == "peptide_sequences", inputs2)]
        @test seq_input2["data"] == ["SEQUENCE"]
        
        # Test invalid instrument type
        @test_throws ArgumentError prepare_koina_batch(model, test_data, "INVALID_INSTRUMENT")
        
        # Test with different valid instruments
        for instrument in ["QE", "QEHFX", "LUMOS", "ELITE", "VELOS"]
            batches = prepare_koina_batch(model, test_data, instrument, batch_size=100)
            @test length(batches) == 1  # All data fits in one batch
            batch = JSON.parse(batches[1])
            
            # Should include instrument_types input for unispec
            input_names = [inp["name"] for inp in batch["inputs"]]
            @test "instrument_types" in input_names
            
            instrument_input = batch["inputs"][findfirst(inp -> inp["name"] == "instrument_types", batch["inputs"])]
            @test instrument_input["data"] == [instrument, instrument]  # Repeated for each peptide
        end
    end
    
    @testset "InstrumentAgnosticModel Preparation" begin
        model = InstrumentAgnosticModel("prosit_2020_hcd")
        test_data = DataFrame(
            koina_sequence = ["PEPTIDE", "SEQUENCE"],
            precursor_charge = Int32[2, 3],
            collision_energy = Float32[25.0, 30.0]
        )
        
        # Prosit ignores instrument_type parameter
        batches = prepare_koina_batch(model, test_data, "ANY_VALUE", batch_size=100)
        @test length(batches) == 1  # All in one batch
        
        batch = JSON.parse(batches[1])
        input_names = [inp["name"] for inp in batch["inputs"]]
        
        # Prosit should NOT have instrument_types input
        @test "instrument_types" ∉ input_names
        @test "peptide_sequences" in input_names
        @test "precursor_charges" in input_names
        @test "collision_energies" in input_names
        
        # Check that all peptides are included
        seq_input = batch["inputs"][findfirst(inp -> inp["name"] == "peptide_sequences", batch["inputs"])]
        @test seq_input["data"] == ["PEPTIDE", "SEQUENCE"]
        @test seq_input["shape"] == [2, 1]  # 2 peptides
    end
    
    @testset "SplineCoefficientModel Preparation" begin
        model = SplineCoefficientModel("altimeter")
        test_data = DataFrame(
            koina_sequence = ["PEPTIDE", "SEQUENCE"],
            precursor_charge = Int32[2, 3],
            collision_energy = Float32[25.0, 30.0]  # Should be ignored for altimeter
        )
        
        batches = prepare_koina_batch(model, test_data, "QE", batch_size=100)
        @test length(batches) == 1
        
        batch = JSON.parse(batches[1])
        input_names = [inp["name"] for inp in batch["inputs"]]
        
        # Altimeter should only have sequence and charge, no NCE
        @test "peptide_sequences" in input_names
        @test "precursor_charges" in input_names
        @test "collision_energies" ∉ input_names  # Altimeter doesn't use NCE
        @test "instrument_types" ∉ input_names   # Altimeter doesn't use instruments
        
        @test length(batch["inputs"]) == 2  # Only sequence and charge
    end
    
    @testset "Batch Size Handling" begin
        model = InstrumentSpecificModel("unispec")
        
        # Test with larger dataset  
        large_data = DataFrame(
            koina_sequence = ["PEPTIDE$i" for i in 1:150],
            precursor_charge = repeat(Int32[2, 3], 75),
            collision_energy = repeat(Float32[25.0, 30.0], 75)
        )
        
        # Test normal batch size
        batches = prepare_koina_batch(model, large_data, "QE", batch_size=50)
        @test length(batches) == 3  # 150 peptides / 50 per batch = 3 batches
        
        # Test that batch size is enforced
        for (i, batch_json) in enumerate(batches)
            batch = JSON.parse(batch_json)
            seq_input = batch["inputs"][findfirst(inp -> inp["name"] == "peptide_sequences", batch["inputs"])]
            
            if i < 3  # First two batches should have 50 peptides
                @test length(seq_input["data"]) == 50
                @test seq_input["shape"] == [50, 1]
            else  # Last batch should have remaining 50 peptides
                @test length(seq_input["data"]) == 50
            end
        end
        
        # Test batch size larger than data
        batches = prepare_koina_batch(model, large_data, "QE", batch_size=200) 
        @test length(batches) == 1  # All data fits in one batch
        
        batch = JSON.parse(batches[1])
        seq_input = batch["inputs"][findfirst(inp -> inp["name"] == "peptide_sequences", batch["inputs"])]
        @test length(seq_input["data"]) == 150
        @test seq_input["shape"] == [150, 1]
    end
    
    @testset "Edge Cases" begin
        model = InstrumentSpecificModel("unispec")
        
        # Test empty dataframe
        empty_data = DataFrame(
            koina_sequence = String[],
            precursor_charge = Int32[],
            collision_energy = Float32[]
        )
        
        batches = prepare_koina_batch(model, empty_data, "QE", batch_size=10)
        @test length(batches) == 0
        
        # Test single row
        single_data = DataFrame(
            koina_sequence = ["PEPTIDE"],
            precursor_charge = Int32[2],
            collision_energy = Float32[25.0]
        )
        
        batches = prepare_koina_batch(model, single_data, "QE", batch_size=10)
        @test length(batches) == 1
        
        batch = JSON.parse(batches[1])
        seq_input = batch["inputs"][findfirst(inp -> inp["name"] == "peptide_sequences", batch["inputs"])]
        @test seq_input["data"] == ["PEPTIDE"]
        @test seq_input["shape"] == [1, 1]
    end
    
    @testset "Data Type Consistency" begin
        model = InstrumentSpecificModel("unispec")
        test_data = DataFrame(
            koina_sequence = ["PEPTIDE", "SEQUENCE"],
            precursor_charge = Int32[2, 3],
            collision_energy = Float32[25.0, 30.0]
        )
        
        batches = prepare_koina_batch(model, test_data, "QE", batch_size=100)
        batch = JSON.parse(batches[1])
        
        # Check that all inputs have correct datatypes specified
        for input in batch["inputs"]
            if input["name"] == "peptide_sequences"
                @test input["datatype"] == "BYTES"
                @test all(isa(seq, String) for seq in input["data"])
            elseif input["name"] == "precursor_charges"
                @test input["datatype"] == "INT32"
                @test all(isa(charge, Int) for charge in input["data"])
            elseif input["name"] == "collision_energies"
                @test input["datatype"] == "FP32"
                @test all(isa(nce, AbstractFloat) for nce in input["data"])
            elseif input["name"] == "instrument_types"
                @test input["datatype"] == "BYTES"
                @test all(isa(inst, String) for inst in input["data"])
            end
        end
    end
end