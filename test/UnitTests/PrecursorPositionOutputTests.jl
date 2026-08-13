@testset "precursor protein-position output" begin
    mktempdir() do temp_dir
        chunk_path = joinpath(temp_dir, "precursors.arrow")
        Arrow.write(chunk_path, DataFrame(
            file_name = ["run1", "run1"],
            species = ["human", "human"],
            inferred_protein_group = ["P1", "P2"],
            accession_numbers = ["P1;P3", "P2"],
            peptide_start_positions = ["17;42", "9"],
            sequence = ["PEPTIDEK", "ANOTHERK"],
            charge = UInt8[2, 2],
            structural_mods = Union{Missing,String}[missing, missing],
            isotopic_mods = Union{Missing,String}[missing, missing],
            prec_mz = Float32[500, 600],
            missed_cleavage = UInt8[0, 1],
            num_enzymatic_termini = UInt8[1, 2],
            global_prob = Float32[0.99, 0.98],
            global_qval = Float32[0.001, 0.002],
            use_for_protein_quant = [true, true],
            precursor_idx = UInt32[1, 2],
            target = [true, true],
            entrapment_group_id = UInt8[0, 0],
            peak_area = Float32[1000, 2000],
        ))

        proteins_path = joinpath(temp_dir, "proteins.arrow")
        Arrow.write(proteins_path, (
            accession = ["P1", "P2", "P3"],
            gene_name = ["G1", "G2", "G3"],
            protein_name = ["Protein 1", "Protein 2", "Protein 3"],
        ))
        proteins = Pioneer.SetProteins(Arrow.Table(proteins_path))
        refs = [Pioneer.PSMFileReference(chunk_path)]

        Pioneer.writePrecursorCSV_chunked(
            refs,
            temp_dir,
            ["run1"],
            false,
            proteins;
            write_csv = true,
        )

        long = CSV.read(
            joinpath(temp_dir, "precursors_long.tsv"),
            DataFrame;
            delim = '\t',
            types = Dict(:peptide_start_positions => String),
        )
        wide = CSV.read(
            joinpath(temp_dir, "precursors_wide.tsv"),
            DataFrame;
            delim = '\t',
            types = Dict(:peptide_start_positions => String),
        )
        wide_arrow = Arrow.Table(joinpath(temp_dir, "precursors_wide.arrow"))

        expected_order = [
            :prec_mz,
            :peptide_start_positions,
            :missed_cleavage,
            :num_enzymatic_termini,
        ]
        long_names = Symbol.(names(long))
        wide_names = Symbol.(names(wide))
        @test long.peptide_start_positions == ["17;42", "9"]
        @test wide.peptide_start_positions == ["17;42", "9"]
        @test collect(wide_arrow.peptide_start_positions) == ["17;42", "9"]
        @test filter(in(expected_order), long_names) == expected_order
        @test filter(in(expected_order), wide_names) == expected_order
        @test Symbol.(propertynames(wide_arrow))[
            findall(in(expected_order), Symbol.(propertynames(wide_arrow)))
        ] == expected_order
        @test wide.missed_cleavage == UInt8[0, 1]
    end
end
