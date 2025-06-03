@testset "protein_inference_namedtuple" begin
    function sort_by_key(dict::Dictionary)
        sorted_dict = Dictionary{keytype(dict), valtype(dict)}()
        for key in sort(collect(keys(dict)), by=k -> (k.peptide, k.decoy, k.entrap_id))
            insert!(sorted_dict, key, dict[key])
        end
        return sorted_dict
    end
    
    # Test Case A: Distinct proteins
    # Protein A has peptides 1, 2
    # Protein B has peptides 3, 4
    proteins = [(protein_name = "A", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "B", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "B", decoy = false, entrap_id = UInt8(1))]
    peptides = ["pep1", "pep2", "pep3", "pep4"]
    
    dict_a = Dictionary{NamedTuple{(:peptide, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}, 
    NamedTuple{(:protein_name, :decoy, :entrap_id, :retain), Tuple{String, Bool, UInt8, Bool}}}()
    insert!(dict_a, (peptide = "pep1", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_a, (peptide = "pep2", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_a, (peptide = "pep3", decoy = false, entrap_id = UInt8(1)), (protein_name = "B", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_a, (peptide = "pep4", decoy = false, entrap_id = UInt8(1)), (protein_name = "B", decoy = false, entrap_id = UInt8(1), retain = true))
    
    result_a = infer_proteins(proteins, peptides)
    @test dict_a == sort_by_key(result_a)
    
    # Test Case B: Differentiable proteins
    # Protein A has peptides 1, 2, 3
    # Protein B has peptides 3, 4, 5
    proteins = [(protein_name = "A", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A;B", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A;B", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "B", decoy = false, entrap_id = UInt8(1))]
    peptides = ["pep1", "pep2", "pep3", "pep4"]
    
    dict_b = Dictionary{NamedTuple{(:peptide, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}, 
                        NamedTuple{(:protein_name, :decoy, :entrap_id, :retain), Tuple{String, Bool, UInt8, Bool}}}()
    insert!(dict_b, (peptide = "pep1", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_b, (peptide = "pep2", decoy = false, entrap_id = UInt8(1)), (protein_name = "A;B", decoy = false, entrap_id = UInt8(1), retain = false))
    insert!(dict_b, (peptide = "pep3", decoy = false, entrap_id = UInt8(1)), (protein_name = "A;B", decoy = false, entrap_id = UInt8(1), retain = false))
    insert!(dict_b, (peptide = "pep4", decoy = false, entrap_id = UInt8(1)), (protein_name = "B", decoy = false, entrap_id = UInt8(1), retain = true))
    
    result_b = infer_proteins(proteins, peptides)
    @test dict_b == sort_by_key(result_b)
    
    # Test Case C: Indistinguishable proteins
    # Protein A has peptides 1, 2, 3, 4
    # Protein B has peptides 1, 2, 3, 4
    proteins = [(protein_name = "A;B", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A;B", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A;B", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A;B", decoy = false, entrap_id = UInt8(1))]
    peptides = ["pep1", "pep2", "pep3", "pep4"]
    
    dict_c = Dictionary{NamedTuple{(:peptide, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}, 
                        NamedTuple{(:protein_name, :decoy, :entrap_id, :retain), Tuple{String, Bool, UInt8, Bool}}}()
    insert!(dict_c, (peptide = "pep1", decoy = false, entrap_id = UInt8(1)), (protein_name = "A;B", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_c, (peptide = "pep2", decoy = false, entrap_id = UInt8(1)), (protein_name = "A;B", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_c, (peptide = "pep3", decoy = false, entrap_id = UInt8(1)), (protein_name = "A;B", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_c, (peptide = "pep4", decoy = false, entrap_id = UInt8(1)), (protein_name = "A;B", decoy = false, entrap_id = UInt8(1), retain = true))
    
    result_c = infer_proteins(proteins, peptides)
    @test dict_c == sort_by_key(result_c)
    
    # Test Case D: Subset proteins (MODIFIED)
    # With maximum parsimony, only protein A is needed
    proteins = [(protein_name = "A", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A;B", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A;B", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A;B", decoy = false, entrap_id = UInt8(1))]
    peptides = ["pep1", "pep2", "pep3", "pep4"]
    
    dict_d = Dictionary{NamedTuple{(:peptide, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}, 
    NamedTuple{(:protein_name, :decoy, :entrap_id, :retain), Tuple{String, Bool, UInt8, Bool}}}()
    insert!(dict_d, (peptide = "pep1", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_d, (peptide = "pep2", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_d, (peptide = "pep3", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_d, (peptide = "pep4", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    
    result_d = infer_proteins(proteins, peptides)
    @test dict_d == sort_by_key(result_d)
    
    # Test Case E: Subsumable proteins (MODIFIED)
    # Maximum parsimony results in proteins A and C
    proteins = [(protein_name = "A", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A;B", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "B;C", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "C", decoy = false, entrap_id = UInt8(1))]
    peptides = ["pep1", "pep2", "pep3", "pep4"]
    
    dict_e = Dictionary{NamedTuple{(:peptide, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}, 
    NamedTuple{(:protein_name, :decoy, :entrap_id, :retain), Tuple{String, Bool, UInt8, Bool}}}()
    insert!(dict_e, (peptide = "pep1", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_e, (peptide = "pep2", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_e, (peptide = "pep3", decoy = false, entrap_id = UInt8(1)), (protein_name = "C", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_e, (peptide = "pep4", decoy = false, entrap_id = UInt8(1)), (protein_name = "C", decoy = false, entrap_id = UInt8(1), retain = true))
    
    result_e = infer_proteins(proteins, peptides)
    @test dict_e == sort_by_key(result_e)
    
    # Test Case F: Protein groups with shared peptides only
    proteins = [(protein_name = "A;B", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A;B;C", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A;B;C", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A;C", decoy = false, entrap_id = UInt8(1))]
    peptides = ["pep1", "pep2", "pep3", "pep4"]
    
    dict_f = Dictionary{NamedTuple{(:peptide, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}, 
    NamedTuple{(:protein_name, :decoy, :entrap_id, :retain), Tuple{String, Bool, UInt8, Bool}}}()
    insert!(dict_f, (peptide = "pep1", decoy = false, entrap_id = UInt8(1)), (protein_name = "A;B;C", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_f, (peptide = "pep2", decoy = false, entrap_id = UInt8(1)), (protein_name = "A;B;C", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_f, (peptide = "pep3", decoy = false, entrap_id = UInt8(1)), (protein_name = "A;B;C", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_f, (peptide = "pep4", decoy = false, entrap_id = UInt8(1)), (protein_name = "A;B;C", decoy = false, entrap_id = UInt8(1), retain = true))
    
    result_f = infer_proteins(proteins, peptides)
    @test dict_f == sort_by_key(result_f)
    
    # Test Case G: Complex case (MODIFIED)
    # Maximum parsimony results in proteins A and C
    proteins = [(protein_name = "A", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A;B", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "B;C", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "B;C", decoy = false, entrap_id = UInt8(1))]
    peptides = ["pep1", "pep2", "pep3", "pep4"]
    
    dict_g = Dictionary{NamedTuple{(:peptide, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}, 
    NamedTuple{(:protein_name, :decoy, :entrap_id, :retain), Tuple{String, Bool, UInt8, Bool}}}()
    insert!(dict_g, (peptide = "pep1", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_g, (peptide = "pep2", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_g, (peptide = "pep3", decoy = false, entrap_id = UInt8(1)), (protein_name = "C", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_g, (peptide = "pep4", decoy = false, entrap_id = UInt8(1)), (protein_name = "C", decoy = false, entrap_id = UInt8(1), retain = true))
    
    result_g = infer_proteins(proteins, peptides)
    @test dict_g == sort_by_key(result_g)
    
    # Test Case H: Combination of all previous test cases (MODIFIED to reflect changes in D, E, G)
    proteins = [(protein_name = "A", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "A", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "B", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "B", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "C", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "C;D", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "C;D", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "D", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "E;F", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "E;F", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "E;F", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "E;F", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "G", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "G;H", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "G;H", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "G;H", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "I", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "I;J", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "J;K", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "K", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "L;M", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "L;M;N", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "L;M;N", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "L;N", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "O", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "O;P", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "P;Q", decoy = false, entrap_id = UInt8(1)),
                (protein_name = "P;Q", decoy = false, entrap_id = UInt8(1))]
    
    peptides = ["pep1", "pep2", "pep3", "pep4",
                "pep5", "pep6", "pep7", "pep8",
                "pep9", "pep10", "pep11", "pep12",
                "pep13", "pep14", "pep15", "pep16",
                "pep17", "pep18", "pep19", "pep20",
                "pep21", "pep22", "pep23", "pep24",
                "pep25", "pep26", "pep27", "pep28"]
    
    dict_h = Dictionary{NamedTuple{(:peptide, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}, 
    NamedTuple{(:protein_name, :decoy, :entrap_id, :retain), Tuple{String, Bool, UInt8, Bool}}}()
    
    insert!(dict_h, (peptide = "pep1", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep2", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep3", decoy = false, entrap_id = UInt8(1)), (protein_name = "B", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep4", decoy = false, entrap_id = UInt8(1)), (protein_name = "B", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep5", decoy = false, entrap_id = UInt8(1)), (protein_name = "C", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep6", decoy = false, entrap_id = UInt8(1)), (protein_name = "C;D", decoy = false, entrap_id = UInt8(1), retain = false))
    insert!(dict_h, (peptide = "pep7", decoy = false, entrap_id = UInt8(1)), (protein_name = "C;D", decoy = false, entrap_id = UInt8(1), retain = false))
    insert!(dict_h, (peptide = "pep8", decoy = false, entrap_id = UInt8(1)), (protein_name = "D", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep9", decoy = false, entrap_id = UInt8(1)), (protein_name = "E;F", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep10", decoy = false, entrap_id = UInt8(1)), (protein_name = "E;F", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep11", decoy = false, entrap_id = UInt8(1)), (protein_name = "E;F", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep12", decoy = false, entrap_id = UInt8(1)), (protein_name = "E;F", decoy = false, entrap_id = UInt8(1), retain = true))
    
    # Case D related (peptides 13-16) - ALL attributed to G
    insert!(dict_h, (peptide = "pep13", decoy = false, entrap_id = UInt8(1)), (protein_name = "G", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep14", decoy = false, entrap_id = UInt8(1)), (protein_name = "G", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep15", decoy = false, entrap_id = UInt8(1)), (protein_name = "G", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep16", decoy = false, entrap_id = UInt8(1)), (protein_name = "G", decoy = false, entrap_id = UInt8(1), retain = true))
    
    # Case E related (peptides 17-20) - Split between I and K
    insert!(dict_h, (peptide = "pep17", decoy = false, entrap_id = UInt8(1)), (protein_name = "I", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep18", decoy = false, entrap_id = UInt8(1)), (protein_name = "I", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep19", decoy = false, entrap_id = UInt8(1)), (protein_name = "K", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep20", decoy = false, entrap_id = UInt8(1)), (protein_name = "K", decoy = false, entrap_id = UInt8(1), retain = true))
    
    insert!(dict_h, (peptide = "pep21", decoy = false, entrap_id = UInt8(1)), (protein_name = "L;M;N", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep22", decoy = false, entrap_id = UInt8(1)), (protein_name = "L;M;N", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep23", decoy = false, entrap_id = UInt8(1)), (protein_name = "L;M;N", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep24", decoy = false, entrap_id = UInt8(1)), (protein_name = "L;M;N", decoy = false, entrap_id = UInt8(1), retain = true))
    
    # Case G related (peptides 25-28) - Split between O and P;Q
    insert!(dict_h, (peptide = "pep25", decoy = false, entrap_id = UInt8(1)), (protein_name = "O", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep26", decoy = false, entrap_id = UInt8(1)), (protein_name = "O", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep27", decoy = false, entrap_id = UInt8(1)), (protein_name = "Q", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_h, (peptide = "pep28", decoy = false, entrap_id = UInt8(1)), (protein_name = "Q", decoy = false, entrap_id = UInt8(1), retain = true))
    
    result_h = infer_proteins(proteins, peptides)
    expected_keys = Set(keys(dict_h))
    result_keys = Set(keys(result_h))
    
    # Check if all expected keys are present
    @test expected_keys == result_keys
    
    # Check if values match for all keys
    for key in keys(dict_h)
        @test dict_h[key] == result_h[key] #"Value mismatch for key $key: Expected $(dict_h[key]), got $(result_h[key])"
    end
    
    # Test a case with decoy proteins
    proteins_with_decoy = [
        (protein_name = "A", decoy = false, entrap_id = UInt8(1)),
        (protein_name = "B", decoy = true, entrap_id = UInt8(1)),  # Decoy protein
        (protein_name = "C", decoy = false, entrap_id = UInt8(1))
    ]
    peptides_with_decoy = ["pep1", "pep2", "pep3"]
    
    dict_decoy = Dictionary{NamedTuple{(:peptide, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}, 
    NamedTuple{(:protein_name, :decoy, :entrap_id, :retain), Tuple{String, Bool, UInt8, Bool}}}()
    
    insert!(dict_decoy, (peptide = "pep1", decoy = false, entrap_id = UInt8(1)), (protein_name = "A", decoy = false, entrap_id = UInt8(1), retain = true))
    insert!(dict_decoy, (peptide = "pep2", decoy = true, entrap_id = UInt8(1)), (protein_name = "B", decoy = true, entrap_id = UInt8(1), retain = true))
    insert!(dict_decoy, (peptide = "pep3", decoy = false, entrap_id = UInt8(1)), (protein_name = "C", decoy = false, entrap_id = UInt8(1), retain = true))
    
    result_decoy = infer_proteins(proteins_with_decoy, peptides_with_decoy)
    @test dict_decoy == sort_by_key(result_decoy)
    
    println("All protein inference tests for NamedTuple version passed!")
end