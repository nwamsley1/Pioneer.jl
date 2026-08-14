

@testset "Fragment Bounds Tests" begin
    
    @testset "get_fragment_bounds - from MS data" begin
        @testset "Basic fragment bounds calculation" begin
            # Create test MS2 data
            center_mass = Float32[500.0, 600.0, 700.0, 800.0, 900.0]
            isolation_width = Float32[2.0, 2.0, 4.0, 2.0, 4.0]
            ms_order = UInt8[2, 2, 2, 2, 2]
            low_frag_mass = Float32[100.0, 110.0, 120.0, 130.0, 140.0]
            high_frag_mass = Float32[1400.0, 1500.0, 1600.0, 1700.0, 1800.0]
            
            frag_bounds, prec_min, prec_max = get_fragment_bounds(
                center_mass,
                isolation_width,
                ms_order,
                low_frag_mass,
                high_frag_mass
            )
            
            # Test that we get a FragBoundModel
            @test frag_bounds isa FragBoundModel
            @test frag_bounds.low_mass isa ImmutablePolynomial
            @test frag_bounds.high_mass isa ImmutablePolynomial
            
            # Test precursor bounds
            # Min should be min(center - width/2) = min(499, 599, 698, 799, 898) = 499
            @test prec_min ≈ 499.0
            # Max should be max(center + width/2) = max(501, 601, 702, 801, 902) = 902
            @test prec_max ≈ 902.0
            
            # Test that polynomials produce reasonable values
            test_prec_mz = 650.0f0
            low_bound_val = frag_bounds.low_mass(test_prec_mz)
            high_bound_val = frag_bounds.high_mass(test_prec_mz)
            
            @test low_bound_val > 0
            @test high_bound_val > low_bound_val
        end
        
        @testset "Mixed MS orders" begin
            # Include MS1 scans that should be filtered out
            center_mass = Union{Missing, Float32}[500.0, missing, 600.0, 700.0]
            isolation_width = Union{Missing, Float32}[2.0, missing, 2.0, 4.0]
            ms_order = UInt8[2, 1, 2, 2]  # MS1 scan at position 2
            low_frag_mass = Float32[100.0, 0.0, 110.0, 120.0]
            high_frag_mass = Float32[1400.0, 0.0, 1500.0, 1600.0]
            
            frag_bounds, prec_min, prec_max = get_fragment_bounds(
                center_mass,
                isolation_width,
                ms_order,
                low_frag_mass,
                high_frag_mass
            )
            
            # Should only use MS2 scans
            @test frag_bounds isa FragBoundModel
            @test prec_min ≈ 499.0  # From first MS2 scan
            @test prec_max ≈ 702.0  # From last MS2 scan
        end
        
        @testset "Linear model fitting" begin
            # Create data that should produce a clear linear relationship
            n_points = 10
            center_mass = Float32.(range(400, 1000, n_points))
            isolation_width = fill(Float32(2.0), n_points)
            ms_order = fill(UInt8(2), n_points)
            
            # Linear relationships: low = 0.1 * precursor, high = 2.0 * precursor
            low_frag_mass = Float32.(0.1 .* center_mass)
            high_frag_mass = Float32.(2.0 .* center_mass)
            
            frag_bounds, _, _ = get_fragment_bounds(
                center_mass,
                isolation_width,
                ms_order,
                low_frag_mass,
                high_frag_mass
            )
            
            # Test predictions at various points
            for test_mz in [450.0, 700.0, 950.0]
                low_pred = frag_bounds.low_mass(test_mz)
                high_pred = frag_bounds.high_mass(test_mz)
                
                # Should approximately follow the linear relationships
                # Note: exact values depend on polynomial fitting
                @test low_pred > 0
                @test high_pred > low_pred
            end
        end
    end
    
    @testset "get_fragment_bounds - with auto-detection" begin
        temp_dir = mktempdir()
        
        @testset "Successful auto-detection from file" begin
            # Create a mock MS data file
            ms_data = DataFrame(
                centerMz = Float32[500.0, 600.0, 700.0],
                isolationWidthMz = Float32[2.0, 2.0, 4.0],
                msOrder = UInt8[2, 2, 2],
                lowMz = Float32[100.0, 110.0, 120.0],
                highMz = Float32[1400.0, 1500.0, 1600.0]
            )
            raw_file_path = joinpath(temp_dir, "test_ms_data.arrow")
            Arrow.write(raw_file_path, ms_data)
            
            result = get_fragment_bounds(
                true,  # auto_detect
                raw_file_path,
                (Float32(50.0), Float32(2000.0)),  # default bounds
                (Float32(400.0), Float32(1200.0))  # default precursor bounds
            )
            
            @test result.frag_bounds isa FragBoundModel
            # Should have adjusted bounds by ±1.0
            @test result.prec_mz_min ≈ 498.0  # 499 - 1
            @test result.prec_mz_max ≈ 703.0  # 702 + 1
        end
        
        @testset "File not found - use defaults" begin
            temp_dir = mktempdir()
            non_existent_path = joinpath(temp_dir, "non_existent.arrow")
            
            result = get_fragment_bounds(
                true,  # auto_detect
                non_existent_path,
                (Float32(50.0), Float32(2000.0)),
                (Float32(400.0), Float32(1200.0))
            )
            
            # Should fall back to defaults
            @test result.frag_bounds isa FragBoundModel
            @test result.prec_mz_min == 400.0
            @test result.prec_mz_max == 1200.0
            
            # Check that default polynomials return constant values
            @test result.frag_bounds.low_mass(500.0) == 50.0
            @test result.frag_bounds.high_mass(500.0) == 2000.0
        end
        
        @testset "Invalid file format - use defaults" begin
            # Create a file with wrong format
            bad_data = DataFrame(
                wrong_column = [1, 2, 3]
            )
            
            bad_file_path = joinpath(temp_dir, "bad_format.arrow")
            Arrow.write(bad_file_path, bad_data)
            
            result = get_fragment_bounds(
                true,  # auto_detect
                bad_file_path,
                (Float32(100.0), Float32(1800.0)),
                (Float32(350.0), Float32(1500.0))
            )
            
            # Should fall back to defaults with warning
            @test result.frag_bounds isa FragBoundModel
            @test result.prec_mz_min == 350.0
            @test result.prec_mz_max == 1500.0
        end
        
        @testset "Auto-detection disabled" begin
            # Even with valid file, should use defaults when disabled
            ms_data = DataFrame(
                centerMz = Float32[500.0],
                isolationWidthMz = Float32[2.0],
                msOrder = UInt8[2],
                lowMz = Float32[100.0],
                highMz = Float32[1400.0]
            )
            
            valid_file_path = joinpath(temp_dir, "valid_ms_data.arrow")
            Arrow.write(valid_file_path, ms_data)
            
            result = get_fragment_bounds(
                false,  # auto_detect disabled
                valid_file_path,
                (Float32(75.0), Float32(1750.0)),
                (Float32(425.0), Float32(1300.0))
            )
            
            # Should use defaults directly
            @test result.frag_bounds.low_mass(500.0) == 75.0
            @test result.frag_bounds.high_mass(500.0) == 1750.0
            @test result.prec_mz_min == 425.0
            @test result.prec_mz_max == 1300.0
        end
        
        safe_rmdir(temp_dir)
    end
    
    @testset "FragBoundModel usage" begin
        # Test the FragBoundModel interface
        low_poly = ImmutablePolynomial([50.0, 0.1])  # 50 + 0.1*x
        high_poly = ImmutablePolynomial([100.0, 2.0])  # 100 + 2.0*x
        
        frag_bounds = FragBoundModel(low_poly, high_poly)
        
        @testset "Polynomial evaluation" begin
            # Test at various precursor m/z values
            test_masses = [400.0, 600.0, 800.0, 1000.0]
            
            for mass in test_masses
                low_bound = frag_bounds.low_mass(mass)
                high_bound = frag_bounds.high_mass(mass)
                
                # Check expected values
                @test low_bound ≈ 50.0 + 0.1 * mass
                @test high_bound ≈ 100.0 + 2.0 * mass
                
                # High bound should always be greater than low bound
                @test high_bound > low_bound
            end
        end
        
        @testset "Edge cases" begin
            # Test with very low and high m/z
            @test frag_bounds.low_mass(0.0) ≈ 50.0
            @test frag_bounds.high_mass(0.0) ≈ 100.0
            
            # Very high m/z
            @test frag_bounds.low_mass(5000.0) ≈ 50.0 + 0.1 * 5000.0
            @test frag_bounds.high_mass(5000.0) ≈ 100.0 + 2.0 * 5000.0
        end
    end

    @testset "manual frag_bounds spec" begin
        manual(spec, fmin = 150.0f0, fmax = 2000.0f0) = get_fragment_bounds(
            false, "", (fmin, fmax), (340.0f0, 1250.0f0), spec,
        ).frag_bounds

        @testset "absent key is the old constant behaviour" begin
            # The test that protects every params file written before this key
            # existed: no spec must give flat bounds at the configured values.
            m = manual(nothing)
            for prec in (340.0f0, 700.0f0, 1250.0f0)
                @test m(prec) == (150.0f0, 2000.0f0)
            end
            @test m.low_mass == ImmutablePolynomial(150.0f0)
            @test m.high_mass == ImmutablePolynomial(2000.0f0)
        end

        @testset "preset names resolve to the documented coefficients" begin
            @test parse_frag_bounds_spec("thermo_auto_documented") ==
                (low_slope = 0.0f0, low_intercept = 0.0f0,
                 high_slope = 2.0f0, high_intercept = 10.0f0)
            @test parse_frag_bounds_spec("thermo_auto") ==
                (low_slope = 0.0f0, low_intercept = 0.0f0,
                 high_slope = 2.04f0, high_intercept = 24.1f0)
            # Case and surrounding space are forgiven; the rule is not.
            @test parse_frag_bounds_spec("  Thermo_Auto  ") ==
                parse_frag_bounds_spec("thermo_auto")
            # "constant" is flat, so it must agree with an absent key.
            @test manual(parse_frag_bounds_spec("constant"))(700.0f0) ==
                manual(nothing)(700.0f0)
        end

        @testset "the ceiling slopes, and the low bound does not" begin
            m = manual(parse_frag_bounds_spec("thermo_auto_documented"))
            # 2.00 * 340 + 10
            @test m(340.0f0) == (150.0f0, 690.0f0)
            @test m(900.0f0) == (150.0f0, 1810.0f0)
            # Every Thermo method measured held the low bound fixed, so the
            # preset leaves it at frag_mz_min rather than sloping it too.
            @test first(m(340.0f0)) == first(m(1250.0f0)) == 150.0f0
        end

        @testset "a sloped ceiling is clamped to frag_mz_max" begin
            # Unclamped this reaches 2.0 * 1250 + 10 = 2510 m/z, which no
            # Orbitrap records -- the failure the feature exists to prevent,
            # reintroduced at the top of the range.
            m = manual(parse_frag_bounds_spec("thermo_auto_documented"))
            @test last(m(1250.0f0)) == 2000.0f0
            # The clamp starts exactly where the line crosses the ceiling.
            @test last(m(995.0f0)) == 2000.0f0
            @test last(m(994.0f0)) < 2000.0f0
            # thermo_auto crosses at (2000 - 24.1) / 2.04.
            ma = manual(parse_frag_bounds_spec("thermo_auto"))
            @test last(ma(969.0f0)) == 2000.0f0
            @test last(ma(968.0f0)) < 2000.0f0
        end

        @testset "explicit coefficients" begin
            spec = parse_frag_bounds_spec(Dict{String, Any}(
                "low" => Dict{String, Any}("slope" => 0.0, "intercept" => 150.0),
                "high" => Dict{String, Any}("slope" => 1.5, "intercept" => 5.0),
            ))
            @test manual(spec)(400.0f0) == (150.0f0, 605.0f0)
        end

        @testset "bad specs are rejected, not guessed at" begin
            # A misspelled preset must not quietly build the wrong ceiling.
            @test_throws InvalidParametersError parse_frag_bounds_spec("thermo-auto")
            @test_throws InvalidParametersError parse_frag_bounds_spec(42)
            @test_throws InvalidParametersError parse_frag_bounds_spec(
                Dict{String, Any}("high" => Dict{String, Any}("slope" => "fast")))
            # A preset and coefficients together is ambiguous.
            @test_throws InvalidParametersError parse_frag_bounds_spec(
                Dict{String, Any}("preset" => "thermo_auto",
                                  "low" => Dict{String, Any}(),
                                  "high" => Dict{String, Any}()))
        end
    end

end