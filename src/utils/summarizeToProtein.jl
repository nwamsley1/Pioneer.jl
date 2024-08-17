



function getProteinGroups(
    passing_psms_folder::String,
    protein_groups_folder::String;
    min_peptides = 2
)

    function getProteinGroupsDict(
        psms_path::String;
        min_peptides::Int64 = 2)
        psms_table = Arrow.Table(psms_path)
        protein_groups = Dictionary{@NamedTuple{accession_numbers::String, target::Bool},
        @NamedTuple{
            max_pg_score::Float32, 
            peptides::Set{String}}}()
        for i in range(1, length(psms_table[1]))
            precursor_idx = psms_table[:precursor_idx][i]
            sequence = precursors[:sequence][precursor_idx]
            score = psms_table[:prob][i]
            accession_numbers = precursors[:accession_numbers][precursor_idx]
            keyname = (accession_numbers = accession_numbers, target = psms_table[:target][i])
            if haskey(protein_groups, keyname)
                max_pg_score, peptides = protein_groups[keyname]
                if score > max_pg_score
                    max_pg_score = score
                end
                push!(peptides, sequence)
                protein_groups[keyname] = (max_pg_score = max_pg_score, peptides = peptides)
            else
                sequences = Set{String}((sequence,))
                insert!(protein_groups,
                keyname,
                        (max_pg_score = score,
                        peptides = sequences)
                )
            end
        end
        filter!(x->length(x[:peptides])>=min_peptides, protein_groups)
        return protein_groups
    end

    function writeProteinGroups(
                                    protein_groups::Dictionary{
                                    @NamedTuple{accession_numbers::String, target::Bool},
                                    @NamedTuple{max_pg_score::Float32,  peptides::Set{String}}
                                    },
                                    protein_groups_path::String)
        # Extract keys and values
        keys_array = keys(protein_groups)
        values_array = values(protein_groups)

        # Create vectors for each column
        accession_numbers = [k[:accession_numbers] for k in keys_array]
        target = [k[:target] for k in keys_array]
        max_pg_score = [v[:max_pg_score] for v in values_array]
        #peptides = [join(v[:peptides], ";") for v in values_array]  # Convert Set to String

        # Create DataFrame
        df = DataFrame((
            accession_numbers = accession_numbers,
            target = target,
            max_pg_score = max_pg_score,
        )
        )

        # Convert DataFrame to Arrow.Table
        Arrow.write(protein_groups_path, df)
        return nothing
    end

    for file_path in readdir(passing_psms_folder, join=true)
        protein_groups_path = joinpath(protein_groups_folder, basename(file_path))
        protein_groups =  getProteinGroupsDict(
            file_path;
            min_peptides = min_peptides
        )
        writeProteinGroups(
            protein_groups,
            protein_groups_path
        )
    end
end

#test_table = Arrow.Table("/Users/n.t.wamsley/Desktop/testresults/RESULTS/Search/temp/passing_psms/E10H50Y40_02.arrow")




