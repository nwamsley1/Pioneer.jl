#Assume sorted by protein,peptide keys. Do this in batches and write a long and wide form .csv without
#loading the entire table into memory. 
function writePrecursorCSV(
    long_precursors_path::String,
    file_names::Vector{String};
    batch_size::Int64 = 2000000)

    function makeWideFormat(
        longdf::DataFrame)
        return unstack(longdf,
        [:species,:accession_numbers,:sequence,:structural_mods,:isotopic_mods,:precursor_idx,:target],
        :file_name,:peak_area_normalized)
    end
    precursors_long = DataFrame(Arrow.Table(long_precursors_path))
    n_rows = size(precursors_long, 1)

    out_dir, arrow_path = splitdir(long_precursors_path)
    long_precursors_path = joinpath(out_dir,"precursors_long.tsv")
    wide_precursors_path = joinpath(out_dir,"precursors_wide.tsv")

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
            println(io1,join(names(precursors_long),"\t"))
            println(io2,join(sorted_columns,"\t"))
            for start_idx in 1:batch_size:n_rows
                end_idx = min(start_idx + batch_size - 1, n_rows)

                #For the wide format, can't split a precursor between two batches.
                last_pid = precursors_long[end_idx,:precursor_idx]
                while end_idx < n_rows
                    if precursors_long[end_idx + 1,:precursor_idx] != last_pid
                        break
                    end
                    end_idx += 1
                end

                subdf =  precursors_long[start_idx:end_idx,:]
                CSV.write(io1, subdf, append=true, header=false, delim="\t")
                subunstack = makeWideFormat(subdf)
                col_names = names(subunstack)
                for fname in file_names
                    if fname ∉ col_names
                        subunstack[!,fname] .= missing
                    end
                end
                CSV.write(io2, subunstack[!,sorted_columns], append=true,header=false,delim="\t")
            end
        end
    end
end
