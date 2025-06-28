using Test
using DataFrames

# Test getModKey function
@testset "getModKey tests" begin
    @test getModKey("(5,M,x)(5,M,Unimod:4)(5,M,Unimod:1)") == "Unimod:1;Unimod:4;x"
    @test getModKey("(5,M,Unimod:4)(5,M,Unimod:35)") == "Unimod:35;Unimod:4"
    @test getModKey("(5,M,Unimod:35)(5,M,Unimod:4)") == "Unimod:35;Unimod:4"
    @test getModKey("(5,M,x)") == "x"
end

# Tests for assign_entrapment_pairs!
@testset "assign_entrapment_pairs! tests" begin
    # Test 1: Basic pairing with one target and one entrapment group
    # Expected: 1 pair with 2 members (group 0 + one from group 1)
    test_df1 = DataFrame(
        base_pep_id = [1, 1, 1],
        prec_charge = [2, 2, 2],
        is_decoy = [false, false, false],
        mod_key = ["", "", ""],
        entrapment_group_id = [0, 1, 1],
        sequence = ["PEPTIDE", "PEPTIDE", "PEPTIDE"]
    )
    assign_entrapment_pairs!(test_df1)
    
    @test !any(ismissing.(test_df1.entrap_pair_id[1:2]))  # First two should be paired
    @test test_df1.entrap_pair_id[1] == test_df1.entrap_pair_id[2]  # Same pair ID
    @test ismissing(test_df1.entrap_pair_id[3])  # Extra entrapment seq not considered (this shouldn't happen in real data)
    
    # Test 2: Multiple entrapment groups
    # Expected: 2 pairs (one for each group 0), each with 3 members (group 0 + one from group 1 + one from group 2)
    test_df2 = DataFrame(
        base_pep_id = fill(1, 6),
        prec_charge = fill(2, 6),
        is_decoy = fill(false, 6),
        mod_key = fill("", 6),
        entrapment_group_id = [0, 0, 1, 1, 2, 2],
        sequence = fill("PEPTIDE", 6)
    )
    assign_entrapment_pairs!(test_df2)
    
    # Each group 0 should be paired with one from group 1 and one from group 2
    @test length(unique(test_df2.entrap_pair_id[.!ismissing.(test_df2.entrap_pair_id)])) == 2  # Two unique pairs
    @test count(test_df2.entrap_pair_id .== test_df2.entrap_pair_id[1]) == 3  # First pair has 3 members
    @test count(test_df2.entrap_pair_id .== test_df2.entrap_pair_id[2]) == 3  # Second pair has 3 members
    
    # Test 3: No valid pairs (no group 0)
    # Expected: No pairs created
    test_df3 = DataFrame(
        base_pep_id = [1, 1, 1],
        prec_charge = [2, 2, 2],
        is_decoy = [false, false, false],
        mod_key = ["", "", ""],
        entrapment_group_id = [1, 1, 2],
        sequence = ["PEPTIDE", "PEPTIDE", "PEPTIDE"]
    )
    assign_entrapment_pairs!(test_df3)
    
    @test all(ismissing.(test_df3.entrap_pair_id))  # No pairs should be created
    
    # Test 4: Different base peptides should not be paired
    # Expected: 2 pairs (one for each group), each with 2 members
    test_df4 = DataFrame(
        base_pep_id = [1, 1, 2, 2],
        prec_charge = [2, 2, 2, 2],
        is_decoy = [false, false, false, false],
        mod_key = ["", "", "", ""],
        entrapment_group_id = [0, 1, 0, 1],
        sequence = ["PEPTIDE", "PEPTIDE", "PROTEIN", "PROTEIN"]
    )
    assign_entrapment_pairs!(test_df4)
    
    @test length(unique(skipmissing(test_df4.entrap_pair_id))) == 2  # Two separate pairs
    @test !ismissing(test_df4.entrap_pair_id[1]) && !ismissing(test_df4.entrap_pair_id[2])  # First group paired
    @test !ismissing(test_df4.entrap_pair_id[3]) && !ismissing(test_df4.entrap_pair_id[4])  # Second group paired
    @test test_df4.entrap_pair_id[1] == test_df4.entrap_pair_id[2]  # First peptide pair
    @test test_df4.entrap_pair_id[3] == test_df4.entrap_pair_id[4]  # Second peptide pair
    @test test_df4.entrap_pair_id[1] != test_df4.entrap_pair_id[3]  # Different pairs
    
    # Test 5: Round-robin distribution
    # Expected: 3 pairs (one for each group 0), each with 2 members
    # The group 1 members are distributed round-robin style
    test_df5 = DataFrame(
        base_pep_id = fill(1, 5),
        prec_charge = fill(2, 5),
        is_decoy = fill(false, 5),
        mod_key = fill("", 5),
        entrapment_group_id = [0, 0, 0, 1, 1],
        sequence = fill("PEPTIDE", 5)
    )
    assign_entrapment_pairs!(test_df5)
    
    # Three group 0 entries, two group 1 entries
    # Each group 0 gets paired with one group 1 (round-robin)
    # First group 0 pairs with first group 1
    # Second group 0 pairs with second group 1  
    # Third group 0 pairs with first group 1 again (round-robin)
    @test length(unique(skipmissing(test_df5.entrap_pair_id))) == 3  # Three unique pairs
    @test all(test_df5.entrap_pair_id .== UInt32[1, 2, 3, 3, 2])
    
    println("All assign_entrapment_pairs! tests passed!")
end