function writeTransitionList(best_psms::DataFrame, f_out::String)
    open(f_out, "w") do io
        # loop over data and write each line
        for row in eachrow(best_psms)

            #data = join(append!([row[:proteinNames]*","*row[:sequence]], row[:names]),",")
            data = append!([row[:proteinNames], 
                            row[:sequence]], #replace(row[:sequence], r"\[(.*?)\]" => "")
                            row[:names])
            write(io, join(data,",")*"\n")
        end
    end
end

function writeIAPIMethod(best_psms::DataFrame, f_out::String)
    open(f_out, "w") do io
        # loop over data and write each line
        write(io, join(["protein_name","sequence","precursor_mz","precursor_intensity","transition_mz"],",")*"\n")
        for row in eachrow(best_psms)

            #data = join(append!([row[:proteinNames]*","*row[:sequence]], row[:names]),",")
            data = append!([row[:proteinNames], row[:sequence], row[:precursor_mz], row[:MS1_PEAK_HEIGHT]], row[:transition_mzs])
            write(io, join(data,",")*"\n")
        end
    end
end