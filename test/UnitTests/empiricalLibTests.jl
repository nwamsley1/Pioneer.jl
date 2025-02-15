
@testset "BasicEmpiricalLibrary Tests" begin
    # Test constructor and column renaming
    println("dir ", @__DIR__)
    lib = BasicEmpiricalLibrary(joinpath(@__DIR__,"../../data/libParseTests/empiricalLibSortTest.csv"))
    df = getDF(lib)
    
    @test names(df)[1] == "precursor_idx"  # Test precursor_idx is first column
    @test "prec_mz" in names(df)  # Test PrecursorMz was renamed
    @test "irt" in names(df)  # Test Tr_recalibrated was renamed
    
    # Test precursor_idx assignment
    @test length(unique(df.precursor_idx)) == 3  # Should have 3 unique precursors
    @test all(df[df.ModifiedPeptide .== "PEPT[+80]IDE", :precursor_idx] .== 
              first(df[df.ModifiedPeptide .== "PEPT[+80]IDE", :precursor_idx]))
    
    # Test nested sorting with different rt_bin_tolerances
    nestedLibrarySort!(lib, rt_bin_tol = 1.5)
    df = getDF(lib)
    test_df_a = DataFrame(
        PrecursorMz = [500.5, 500.5, 500.5, 600.7, 600.7, 300.2],
        ModifiedPeptide = ["PEPT[+80]IDE", "PEPT[+80]IDE", "PEPT[+80]IDE", "SAMPL[+16]E", "SAMPL[+16]E", "TEST"],
        PrecursorCharge = [2, 2, 2, 3, 3, 2],
        Tr_recalibrated = [10.5, 10.5, 10.5, 9.2, 9.2, 20.1],
        LibraryIntensity = [1000.0, 500.0, 100.0, 800.0, 200.0, 600.0],
        ProductMz = [200.1, 300.2, 400.3, 250.1, 350.2, 450.3],
        FragmentType = ["y", "b", "y", "y", "b", "y"],
        FragmentCharge = [1, 1, 1, 1, 1, 1],
        FragmentSeriesNumber = [2, 3, 4, 2, 3, 4]
    )
    @test all(df.prec_mz .== test_df_a.PrecursorMz)
    @test all(df.ModifiedPeptide .== test_df_a.ModifiedPeptide)
    @test all(df.PrecursorCharge .== test_df_a.PrecursorCharge)
    @test all(df.irt .== test_df_a.Tr_recalibrated)
    @test all(df.LibraryIntensity .== test_df_a.LibraryIntensity)
    @test all(df.ProductMz .== test_df_a.ProductMz)
    @test all(df.FragmentType .== test_df_a.FragmentType)
    @test all(df.FragmentCharge .== test_df_a.FragmentCharge)
    @test all(df.FragmentSeriesNumber .== test_df_a.FragmentSeriesNumber)

    lib = BasicEmpiricalLibrary(joinpath(@__DIR__,"../../data/libParseTests/empiricalLibSortTest.csv"))
    nestedLibrarySort!(lib, rt_bin_tol = 1.0)
    df = getDF(lib)
    test_df_a = DataFrame(
        PrecursorMz = [600.7, 600.7, 500.5, 500.5, 500.5, 300.2],
        ModifiedPeptide = ["SAMPL[+16]E", "SAMPL[+16]E", "PEPT[+80]IDE", "PEPT[+80]IDE", "PEPT[+80]IDE", "TEST"],
        PrecursorCharge = [3, 3, 2, 2, 2, 2],
        Tr_recalibrated = [9.2, 9.2, 10.5, 10.5, 10.5, 20.1],
        LibraryIntensity = [800.0, 200.0, 1000.0, 500.0, 100.0, 600.0],
        ProductMz = [250.1, 350.2, 200.1, 300.2, 400.3, 450.3],
        FragmentType = ["y", "b", "y", "b", "y", "y"],
        FragmentCharge = [1, 1, 1, 1, 1, 1],
        FragmentSeriesNumber = [2, 3, 2, 3, 4, 4]
     )
    @test all(df.prec_mz .== test_df_a.PrecursorMz)
    @test all(df.ModifiedPeptide .== test_df_a.ModifiedPeptide)
    @test all(df.PrecursorCharge .== test_df_a.PrecursorCharge)
    @test all(df.irt .== test_df_a.Tr_recalibrated)
    @test all(df.LibraryIntensity .== test_df_a.LibraryIntensity)
    @test all(df.ProductMz .== test_df_a.ProductMz)
    @test all(df.FragmentType .== test_df_a.FragmentType)
    @test all(df.FragmentCharge .== test_df_a.FragmentCharge)
    @test all(df.FragmentSeriesNumber .== test_df_a.FragmentSeriesNumber)

    lib = BasicEmpiricalLibrary(joinpath(@__DIR__,"../../data/libParseTests/empiricalLibSortTest.csv"))
    nestedLibrarySort!(lib, rt_bin_tol = typemax(Float64))
    df = getDF(lib)
    test_df_a = DataFrame(
        PrecursorMz = [300.2, 500.5, 500.5, 500.5, 600.7, 600.7],
        ModifiedPeptide = ["TEST", "PEPT[+80]IDE", "PEPT[+80]IDE", "PEPT[+80]IDE", "SAMPL[+16]E", "SAMPL[+16]E"],
        PrecursorCharge = [2, 2, 2, 2, 3, 3],
        Tr_recalibrated = [20.1, 10.5, 10.5, 10.5, 9.2, 9.2],
        LibraryIntensity = [600.0, 100.0, 500.0, 1000.0, 200.0, 800.0],
        ProductMz = [450.3, 400.3, 300.2, 200.1, 350.2, 250.1],
        FragmentType = ["y", "y", "b", "y", "b", "y"],
        FragmentCharge = [1, 1, 1, 1, 1, 1],
        FragmentSeriesNumber = [4, 4, 3, 2, 3, 2]
     )
    @test all(df.prec_mz .== test_df_a.PrecursorMz)
    @test all(df.ModifiedPeptide .== test_df_a.ModifiedPeptide)
    @test all(df.PrecursorCharge .== test_df_a.PrecursorCharge)
    @test all(df.irt .== test_df_a.Tr_recalibrated)
    @test all(df.LibraryIntensity .== test_df_a.LibraryIntensity)
    @test all(df.ProductMz .== test_df_a.ProductMz)
    @test all(df.FragmentType .== test_df_a.FragmentType)
    @test all(df.FragmentCharge .== test_df_a.FragmentCharge)
    @test all(df.FragmentSeriesNumber .== test_df_a.FragmentSeriesNumber)
end