



function getProteinGroups(
    passing_psms_folder::String,
    protein_groups_folder::String,
    accession_numbers::AbstractVector{String},
    accession_number_to_id::Dict{Sring, UInt32},
    precursor_sequence::AbstractVector{String};
    min_peptides = 2,
    protein_q_val_threshold::Float32 = 0.01f0
)

    function getProteinGroupsDict(
        psms_path::String,
        accession_numbers::AbstractVector{String},
        accession_number_to_id::Dict{Sring, UInt32},
        precursor_sequence::AbstractVector{String};
        min_peptides::Int64 = 2)


        psms_table = Arrow.Table(psms_path)

        protein_groups = Dictionary{@NamedTuple{protein_idx::UInt32, target::Bool},
        @NamedTuple{
            max_pg_score::Float32, 
            peptides::Set{String}}
        }()

        for i in range(1, length(psms_table[1]))
            precursor_idx = psms_table[:precursor_idx][i]
            sequence = precursor_sequence[precursor_idx]
            score = psms_table[:prob][i]
            protein_idx = accession_number_to_id[accession_numbers[precursor_idx]]
            keyname = (protein_idx, protein_idx, target = psms_table[:target][i])
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
                                    @NamedTuple{protein_idx::UInt32, target::Bool},
                                    @NamedTuple{max_pg_score::Float32,  peptides::Set{String}}
                                    },
                                    protein_groups_path::String)
        # Extract keys and values
        keys_array = keys(protein_groups)
        values_array = values(protein_groups)

        # Create vectors for each column
        protein_idx = [k[:protein_idx] for k in keys_array]
        target = [k[:target] for k in keys_array]
        max_pg_score = [v[:max_pg_score] for v in values_array]
        #peptides = [join(v[:peptides], ";") for v in values_array]  # Convert Set to String

        # Create DataFrame
        df = DataFrame((
            protein_idx = protein_idx,
            target = target,
            max_pg_score = max_pg_score,
        )
        )

        # Convert DataFrame to Arrow.Table
        Arrow.write(protein_groups_path, df)
        return size(df, 1)
    end


    pg_count = 0
    for file_path in readdir(passing_psms_folder, join=true)
        _, extention = splitext(file_path)
        if extention != ".arrow"
            continue
        end
        protein_groups_path = joinpath(protein_groups_folder, basename(file_path))
        protein_groups =  getProteinGroupsDict(
            file_path,
            accession_numbers,
            accession_number_to_id,
            precursor_sequence;
            min_peptides = min_peptides
        )
        pg_count += writeProteinGroups(
            protein_groups,
            protein_groups_path
        )
    end

    max_pg_scores = Vector{@NamedTuple{score::Float32, target::Bool}}(undef, pg_count)
    n = 1
    protein_groups_paths = readdir(protein_groups_folder, join=true)
    for file_path in protein_groups_paths
        pg_table = Arrow.Table(file_path)
        for i in range(1, length(pg_table[:max_pg_score]))
            max_pg_scores[n] = (score = pg_table[:max_pg_score][i], target = pg_table[:target][i])
            n += 1
        end
    end
    sort!(max_pg_scores, alg=QuickSort, by=x->x[:score],rev = true)
    #return max_pg_scores
    # Second pass: Find the probability threshold
    pg_score_threshold = find_score_threshold(max_pg_scores, protein_q_val_threshold)
    for file_path in protein_groups_paths
        pg_table = Arrow.Table(file_path)
        pg_df = DataFrame(pg_table)[(pg_table[:max_pg_score].>=pg_score_threshold),:]
        Arrow.write(joinpath(protein_groups_folder, basename(file_path)),
        pg_df)
    end

end

#test_table = Arrow.Table("/Users/n.t.wamsley/Desktop/testresults/RESULTS/Search/temp/passing_psms/E10H50Y40_02.arrow")




