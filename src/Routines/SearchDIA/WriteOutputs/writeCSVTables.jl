function parseMods(mods_string::AbstractString)::Base.RegexMatchIterator
    #Example: "1(7,M,Oxidation)(13,K,AnExampleMod)"
    mods_regex = r"(?<=\().*?(?=\))"
    return eachmatch(mods_regex, mods_string)
end

function insert_at_indices(original::S, insertions::Vector{Tuple{String, UInt8}}) where {S<:AbstractString}
    # Convert the original string into an array of characters for easier manipulation
    char_array = collect(original)

    # Sort the insertions by index in ascending order
    sorted_insertions = sort(insertions, by = x -> x[2])

    # Adjust the index for each insertion
    offset = 0
    for (substr, idx) in sorted_insertions
        # Adjust the index with the current offset
        insertion_index = idx + offset
        # Insert each character of the substring at the specified index
        for (i, char) in enumerate(substr)
            insert!(char_array, insertion_index + i, char)
        end
        # Update the offset by the length of the inserted substring
        offset += length(substr)
    end

    # Join the array of characters back into a single string
    return join(char_array)
end
function getModifiedSequence(
    sequence::S,
    isotope_mods::String,
    structural_mods::String,
    charge::UInt8) where {S<:AbstractString}

    mods = structural_mods*isotope_mods
    mods = [("("*getModName(mod.match)*")", getModIndex(mod.match)) for mod in parseMods(mods)]
    return "_"*insert_at_indices(sequence, mods)*"_."*string(charge)
end
function getModifiedSequence(
    sequence::S,
    isotope_mods::String,
    structural_mods::Missing,
    charge::UInt8) where {S<:AbstractString}

    mods = isotope_mods
    mods = [("("*getModName(mod.match)*")", getModIndex(mod.match)) for mod in parseMods(mods)]
    return "_"*insert_at_indices(sequence, mods)*"_."*string(charge)
end

function getModifiedSequence(
    sequence::S,
    isotope_mods::Missing,
    structural_mods::String,
    charge::UInt8) where {S<:AbstractString}

    mods = structural_mods
    mods = [("("*getModName(mod.match)*")", getModIndex(mod.match)) for mod in parseMods(mods)]
    return "_"*insert_at_indices(sequence, mods)*"_."*string(charge)
end
function getModifiedSequence(
    sequence::S,
    isotope_mods::Missing,
    structural_mods::Missing,
    charge::UInt8) where {S<:AbstractString}

    mods = ""
    mods = [("("*getModName(mod.match)*")", getModIndex(mod.match)) for mod in parseMods(mods)]
    return "_"*insert_at_indices(sequence, mods)*"_."*string(charge)
end
#Assume sorted by protein,peptide keys. Do this in batches and write a long and wide form .csv without
#loading the entire table into memory. 
function writePrecursorCSV(
    long_precursors_path::String,
    file_names::Vector{String},
    normalized::Bool;
    write_csv::Bool = true,
    batch_size::Int64 = 2000000)

    function makeWideFormat(
        longdf::DataFrame,
        normalized::Bool)
        if normalized
            return unstack(longdf,
            [:species,:accession_numbers,:sequence,:structural_mods,:isotopic_mods,:precursor_idx,:target],
            :file_name,:peak_area_normalized)
        else
            return unstack(longdf,
            [:species,:accession_numbers,:sequence,:structural_mods,:isotopic_mods,:precursor_idx,:target],
            :file_name,:peak_area)
        end
    end
    precursors_long = DataFrame(Arrow.Table(long_precursors_path))
    n_rows = size(precursors_long, 1)

    out_dir, arrow_path = splitdir(long_precursors_path)
    long_precursors_path = joinpath(out_dir,"precursors_long.tsv")
    wide_precursors_path = joinpath(out_dir,"precursors_wide.tsv")
    wide_precursors_arrow_path = joinpath(out_dir,"precursors_wide.arrow")
    if isfile(wide_precursors_arrow_path)
        rm(wide_precursors_arrow_path)
    end
    wide_columns = ["species"
    "accession_numbers"
    "sequence"
    "structural_mods"
    "isotopic_mods"
    "precursor_idx"
    "target"]

    sorted_columns = vcat(wide_columns, file_names)
    open(long_precursors_path,"w") do io1
        open(wide_precursors_path, "w") do io2
            #Make file headers 
            if write_csv
                println(io1,join(names(precursors_long),"\t"))
                println(io2,join(sorted_columns,"\t"))
            end
            batch_start_idx, batch_end_idx = 1, min(batch_size+1,n_rows)
            n_writes = 0
            while batch_start_idx <= n_rows
                #For the wide format, can't split a precursor between two batches.
                last_pid = precursors_long[batch_end_idx,:precursor_idx]
                while batch_end_idx < n_rows
                    if precursors_long[batch_end_idx + 1,:precursor_idx] != last_pid
                        break
                    end
                    batch_end_idx += 1
                end

                subdf =  precursors_long[range(batch_start_idx, batch_end_idx),:]
                batch_start_idx = batch_end_idx + 1
                batch_end_idx = min(batch_start_idx + batch_size, n_rows)
                if write_csv 
                    CSV.write(io1, subdf, append=true, header=false, delim="\t")
                end
                subunstack = makeWideFormat(subdf, normalized)
                col_names = names(subunstack)
                for fname in file_names
                    if fname ∉ col_names
                        subunstack[!,fname] .= missing
                    end
                end
                if write_csv
                    CSV.write(io2, subunstack[!,sorted_columns], append=true,header=false,delim="\t")
                end
                if iszero(n_writes)
                    if isfile(wide_precursors_arrow_path)
                        rm(wide_precursors_arrow_path)
                    end
                    open(wide_precursors_arrow_path, "w") do io
                        Arrow.write(io, subunstack[!,sorted_columns]; file=false)  # file=false creates stream format
                    end
                else
                    Arrow.append(wide_precursors_arrow_path,subunstack[!,sorted_columns])
                end
                n_writes += 1
            end
        end
    end
    if write_csv == false
        rm(long_precursors_path, force = true)
        rm(wide_precursors_path, force = true)
    end
    return wide_precursors_arrow_path
end
function writeProteinGroupsCSV(
    long_pg_path::String,
    sequences::AbstractVector{<:AbstractString},
    isotope_mods::AbstractVector{Union{Missing, String}},
    structural_mods::AbstractVector{Union{String,Missing}},
    precursor_charge::AbstractVector{UInt8},
    file_names::Vector{String};
    write_csv::Bool = true,
    batch_size::Int64 = 2000000)

    function makeWideFormat(
        longdf::DataFrame)
        # First create a DataFrame with the non-abundance columns we want to keep
        metadata_df = unique(longdf[:, [:species, :protein, :target, :n_peptides]])
        
        # Create the abundance wide format
        abundance_df = unstack(longdf,
            [:species, :protein, :target],
            :file_name, :abundance)
            
        # Join the metadata with the abundance data
        return leftjoin(metadata_df, abundance_df, on=[:species, :protein, :target])
    end

    protein_groups_long = DataFrame(Arrow.Table(long_pg_path))
    n_rows = size(protein_groups_long, 1)

    out_dir, arrow_path = splitdir(long_pg_path)
    long_protein_groups_path = joinpath(out_dir,"protein_groups_long.tsv")
    wide_protein_groups_path = joinpath(out_dir,"protein_groups_wide.tsv")
    wide_protein_groups_arrow = joinpath(out_dir,"protein_groups_wide.arrow")
    if isfile(wide_protein_groups_arrow)
        rm(wide_protein_groups_arrow)
    end
    # Update wide columns to include n_peptides
    wide_columns = ["species", "protein", "target", "n_peptides"]

    sorted_columns = vcat(wide_columns, file_names)
    open(long_protein_groups_path,"w") do io1
        open(wide_protein_groups_path, "w") do io2
            #Make file headers 
            if write_csv 
                println(io1,join(names(protein_groups_long),"\t"))
                println(io2,join(sorted_columns,"\t"))
            end
            batch_start_idx, batch_end_idx = 1, min(batch_size+1,n_rows)
            n_writes = 0
            while batch_start_idx <= n_rows
                #For the wide format, can't split a precursor between two batches.
                last_protein_group = protein_groups_long[batch_end_idx,:protein]
                while batch_end_idx < n_rows
                    if  protein_groups_long[batch_end_idx + 1,:protein] != last_protein_group
                        break
                    end
                    batch_end_idx += 1
                end
                subdf = protein_groups_long[range(batch_start_idx, batch_end_idx),:]
                batch_start_idx = batch_end_idx + 1
                batch_end_idx = min(batch_start_idx + batch_size, n_rows)
                subdf[!,:modified_sequence] = Vector{String}(undef, size(subdf, 1))
                for i in range(1, size(subdf, 1))
                    peptides = subdf[i,:peptides]
                    modified_sequences = Vector{String}(undef, length(peptides))
                    for (j, pid) in enumerate(peptides)
                        if ismissing(pid)
                            modified_sequences[j] = ""
                            continue
                        end
                        modified_sequences[j] = getModifiedSequence(
                            sequences[pid],
                            isotope_mods[pid],
                            structural_mods[pid],
                            precursor_charge[pid])
                    end
                    subdf[i,:modified_sequence] = join(modified_sequences,';')
                end
                subdf[!,:peptides] = subdf[!,:modified_sequence]
                select!(subdf, Not([:modified_sequence]))
                if write_csv 
                    CSV.write(io1, subdf, append=true, header=false, delim="\t")
                end
                subunstack = makeWideFormat(subdf)
                col_names = names(subunstack)
                for fname in file_names
                    if fname ∉ col_names
                        subunstack[!,fname] .= missing
                    end
                end
                if write_csv
                    CSV.write(io2, subunstack[!,sorted_columns], append=true,header=false,delim="\t")
                end

                if iszero(n_writes)
                    if isfile(wide_protein_groups_arrow)
                        rm(wide_protein_groups_arrow)
                    end
                    open(wide_protein_groups_arrow, "w") do io
                        Arrow.write(io, subunstack[!,sorted_columns]; file=false)  # file=false creates stream format
                    end
                else
                    Arrow.append(wide_protein_groups_arrow,subunstack[!,sorted_columns])
                end
                n_writes += 1
            end
        end
        if write_csv == false
            rm(long_protein_groups_path, force = true)
            rm(wide_protein_groups_path, force = true)
        end
    end
    return wide_protein_groups_arrow
end