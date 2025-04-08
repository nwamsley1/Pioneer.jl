test_lib = DataFrame(CSV.File("/Users/nathanwamsley/Data/Apr_2025/Kevin_Tag6/MSFRAGGER_tag6_HY.tsv"))
test_lib[!,:ModifiedPeptideSequence] = replace.(test_lib[!,:ModifiedPeptideSequence], ".[308.1161]" => "")
test_lib[!,:ModifiedPeptideSequence] = replace.(test_lib[!,:ModifiedPeptideSequence], "[436.21106]" => "")
test_lib[!,:ModifiedPeptideSequence] = replace.(test_lib[!,:ModifiedPeptideSequence], ".[350.1267]" => "")
rename!(test_lib,:NormalizedRetentionTime=>:Tr_recalibrated)
rename!(test_lib,:ModifiedPeptideSequence=>:ModifiedPeptide)
CSV.write("/Users/nathanwamsley/Data/Apr_2025/Kevin_Tag6/MSFRAGGER_tag6_HY_ntw.tsv", test_lib, delim = '\t')
ParseSpecLib("/Users/nathanwamsley/Data/Apr_2025/Kevin_Tag6/parse_lib_params_channel_decoys.json")


using CSV, DataFrames
test_lib = DataFrame(CSV.File("/Users/nathanwamsley/Data/Apr_2025/Kevin_Tag6/Astral_tag6_SageLibrary_aligned.tsv"))

function has_invalid_tag6(sequence::String)
    # Find all positions where "(tag6)" occurs
    tag_positions = findall("(tag6)", sequence)
    
    for pos in tag_positions
        start_pos = pos.start
        
        # Check if this position is valid:
        # 1. It comes right after 'K'
        if start_pos > 1 && sequence[start_pos-1] == 'K'
            continue
        # 2. It comes after the first letter of the sequence
        elseif start_pos == 2
            continue
        # Otherwise, it's invalid
        else
            return true
        end
    end
    
    # If we made it here, all tag6 positions are valid
    return false
end

# Apply to your dataframe
test_lib = filter(row -> has_invalid_tag6(row.ModifiedPeptide)==false, test_lib)
test_lib[!,:ModifiedPeptide] = replace.(test_lib[!,:ModifiedPeptide], "(tag6)" => "")
#test_lib[!,:ModifiedPeptideSequence] = replace.(test_lib[!,:ModifiedPeptideSequence], ".[308.1161]" => "")
#test_lib[!,:ModifiedPeptideSequence] = replace.(test_lib[!,:ModifiedPeptideSequence], "[436.21106]" => "")
#test_lib[!,:ModifiedPeptideSequence] = replace.(test_lib[!,:ModifiedPeptideSequence], ".[350.1267]" => "")
#rename!(test_lib,:NormalizedRetentionTime=>:Tr_recalibrated)
CSV.write("/Users/nathanwamsley/Data/Apr_2025/Kevin_Tag6/Astral_tag6_SageLibrary_aligned_ntw.tsv", test_lib, delim = '\t')
ParseSpecLib("/Users/nathanwamsley/Data/Apr_2025/Kevin_Tag6/parse_lib_params_channel_decoys_sage.json")

using CSV, DataFrames
using FilePathsBase
using Parquet2: Dataset

ds = Dataset("/Users/nathanwamsley/Data/Apr_2025/Kevin_Tag6/KMD_sage_tag6_040825/transfer-HY-full.parquet")  # this only loads metadata
df = DataFrame(ds; copycols=false)  # load the entire table as a DataFrame
