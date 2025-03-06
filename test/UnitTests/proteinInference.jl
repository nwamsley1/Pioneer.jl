@testset "protein_inference" begin

    function sort_by_key(dict::Dictionary)
        sorted_dict = Dictionary{keytype(dict), valtype(dict)}()
        for key in sort(collect(keys(dict)))
            insert!(sorted_dict, key, dict[key])
        end
        return sorted_dict
    end
    # Test Case A: Distinct proteins
    # Protein A has peptides 1, 2
    # Protein B has peptides 3, 4
    proteins = ["A","A","B","B"]
    peptides = ["pep1", "pep2","pep3","pep4"]
    dict_a = Dictionary{String, Tuple{String, Bool}}()
    insert!(dict_a, "pep1", ("A", true))
    insert!(dict_a, "pep2", ("A", true))
    insert!(dict_a, "pep3", ("B", true))
    insert!(dict_a, "pep4", ("B", true))
    @assert dict_a == sort_by_key(infer_proteins(proteins, peptides))
    

    # Test Case B: Differentiable proteins
    # Protein A has peptides 1, 2, 3
    # Protein B has peptides 3, 4, 5
    proteins = ["A","A;B","A;B","B"]
    peptides = ["pep1", "pep2","pep3","pep4"]
    dict_b = Dictionary{String, Tuple{String, Bool}}()
    insert!(dict_b, "pep1", ("A", true))
    insert!(dict_b, "pep2", ("A;B", false))
    insert!(dict_b, "pep3", ("A;B", false))
    insert!(dict_b, "pep4", ("B", true))
    @assert dict_b == sort_by_key(infer_proteins(proteins, peptides))
    
    # Test Case C: Indistinguishable proteins
    # Protein A has peptides 1, 2, 3, 4
    # Protein B has peptides 1, 2, 3, 4
    proteins = ["A;B","A;B","A;B","A;B"]
    peptides = ["pep1", "pep2","pep3","pep4"]
    dict_c = Dictionary{String, Tuple{String, Bool}}()
    insert!(dict_c, "pep1", ("A;B", true))
    insert!(dict_c, "pep2", ("A;B", true))
    insert!(dict_c, "pep3", ("A;B", true))
    insert!(dict_c, "pep4", ("A;B", true))
    @assert dict_c == sort_by_key(infer_proteins(proteins, peptides))

    # Test Case D: Subset proteins
    # Protein A has peptides 1, 2, 3, 4
    # Protein B has peptides 1, 2
    proteins = ["A","A;B","A;B","A;B"]
    peptides = ["pep1", "pep2","pep3","pep4"]
    dict_d = Dictionary{String, Tuple{String, Bool}}()
    insert!(dict_d, "pep1", ("A", true))
    insert!(dict_d, "pep2", ("A;B", false))
    insert!(dict_d, "pep3", ("A;B", false))
    insert!(dict_d, "pep4", ("A;B", false))
    @assert dict_d == sort_by_key(infer_proteins(proteins, peptides))

    # Test Case E: Subsumable proteins
    # Protein A has peptides 1, 2
    # Protein B has peptides 2, 3
    # Protein C has peptides 3, 4
    proteins = ["A","A;B","B;C","C"]
    peptides = ["pep1", "pep2","pep3","pep4"]
    dict_e = Dictionary{String, Tuple{String, Bool}}()
    insert!(dict_e, "pep1", ("A", true))
    insert!(dict_e, "pep2", ("A;B", false))
    insert!(dict_e, "pep3", ("B;C", false))
    insert!(dict_e, "pep4", ("C", true))
    @assert dict_e == sort_by_key(infer_proteins(proteins, peptides))
    
    # Test Case F: Protein groups with shared peptides only
    proteins = ["A;B","A;B;C","A;B;C","A;C"]
    peptides = ["pep1", "pep2","pep3","pep4"]
    dict_f = Dictionary{String, Tuple{String, Bool}}()
    insert!(dict_f, "pep1", ("A;B;C", true))
    insert!(dict_f, "pep2", ("A;B;C", true))
    insert!(dict_f, "pep3", ("A;B;C", true))
    insert!(dict_f, "pep4", ("A;B;C", true))
    @assert dict_f == sort_by_key(infer_proteins(proteins, peptides))


    # Test Case G:
    # Protein A has peptides 1, 2
    # Protein B has peptides 2, 3, 4
    # Protein C has peptides 3, 4
    proteins = ["A","A;B","B;C","B;C"]
    peptides = ["pep1", "pep2","pep3","pep4"]
    dict_g = Dictionary{String, Tuple{String, Bool}}()
    insert!(dict_g, "pep1", ("A", true))
    insert!(dict_g, "pep2", ("A;B", false))
    insert!(dict_g, "pep3", ("B;C", false))
    insert!(dict_g, "pep4", ("B;C", false))
    @assert dict_g == sort_by_key(infer_proteins(proteins, peptides))

    # Test Case H:
    #Combination of all previous test cases 

    proteins = ["A","A","B","B",
    "C","C;D","C;D","D",
    "E;F","E;F","E;F","E;F",
    "G","G;H","G;H","G;H",
    "I","I;J","J;K","K",
    "L;M","L;M;N","L;M;N","L;N",
    "O","O;P","P;Q","P;Q"]
    peptides = ["pep1","pep2","pep3","pep4",
                "pep5","pep6","pep7","pep8",
                "pep9","pep10","pep11","pep12",
                "pep13","pep14","pep15","pep16",
                "pep17","pep18","pep19","pep20",
                "pep21","pep22","pep23","pep24",
                "pep25","pep26","pep27","pep28"]
    dict_h = Dictionary{String, Tuple{String, Bool}}()
    insert!(dict_h, "pep1", ("A", true))
    insert!(dict_h, "pep2", ("A", true))
    insert!(dict_h, "pep3", ("B", true))
    insert!(dict_h, "pep4", ("B", true))
    insert!(dict_h, "pep5", ("C", true))
    insert!(dict_h, "pep6", ("C;D", false))
    insert!(dict_h, "pep7", ("C;D", false))
    insert!(dict_h, "pep8", ("D", true))
    insert!(dict_h, "pep9", ("E;F", true))
    insert!(dict_h, "pep10", ("E;F", true))
    insert!(dict_h, "pep11", ("E;F", true))
    insert!(dict_h, "pep12", ("E;F", true))
    insert!(dict_h, "pep13", ("G", true))
    insert!(dict_h, "pep14", ("G;H", false))
    insert!(dict_h, "pep15", ("G;H", false))
    insert!(dict_h, "pep16", ("G;H", false))
    insert!(dict_h, "pep17", ("I", true))
    insert!(dict_h, "pep18", ("I;J", false))
    insert!(dict_h, "pep19", ("J;K", false))
    insert!(dict_h, "pep20", ("K", true))
    insert!(dict_h, "pep21", ("L;M;N", true))
    insert!(dict_h, "pep22", ("L;M;N", true))
    insert!(dict_h, "pep23", ("L;M;N", true))
    insert!(dict_h, "pep24", ("L;M;N", true))
    insert!(dict_h, "pep25", ("O", true))
    insert!(dict_h, "pep26", ("O;P", false))
    insert!(dict_h, "pep27", ("P;Q", false))
    insert!(dict_h, "pep28", ("P;Q", false))
    solution = infer_proteins(proteins, peptides)
    @test all([dict_h[key] == solution[key] for (key, val) in pairs(solution)])
    println("All protein inference tests passed!")
end

