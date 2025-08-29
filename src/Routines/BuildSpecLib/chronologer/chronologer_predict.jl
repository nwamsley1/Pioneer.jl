# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# src/chronologer/chronologer_predict.jl

"""
    predict_retention_times(chronologer_out_path::String)

Predict retention times for peptides using either Koina's Chronologer service
or local Chronologer installation as fallback.

Parameters:
- chronologer_out_path::String: Path to Arrow file containing peptide data.
                               Must have 'chronologer_sequence' column.
                               Will be updated in-place with predictions.

Notes:
- First attempts prediction through Koina API
- Falls back to local Chronologer if Koina fails
- Handles UniMod code conversion for local Chronologer
- Updates the input file in place with RT predictions
"""
function predict_retention_times(chronologer_in_path::String, chronologer_out_path::String)
    # Try Koina service first
    try
        chronologer_table = DataFrame(Tables.columntable(Arrow.Table(chronologer_in_path)))
        predictions = predict_rt_koina(chronologer_table)
        chronologer_table[!, :rt] = predictions
        Arrow.write(chronologer_out_path, chronologer_table)
        return
    catch e
        @user_warn "Chronologer failed through Koina. Falling back to local installation..." exception=e
        rethrow(e)
    end
    # Fall back to local Chronologer
    # no longer included. See commits before 
    #predict_rt_local(chronologer_out_path)
end

"""
Helper function to predict RTs using Koina service.
"""
function predict_rt_koina(chronologer_table::DataFrame)::Vector{Float32}
    model = RetentionTimeModel("chronologer")
    
    # Prepare batches
    batches = prepare_koina_batch(
        model,
        chronologer_table,
        batch_size=1000
    )
    
    # Make requests
    results = make_koina_batch_requests(
        batches,
        KOINA_URLS["chronologer"]
    )
    
    # Parse results
    rt_predictions = Float32[]
    for result in results
        batch_result = parse_koina_batch(model, result)
        append!(rt_predictions, batch_result.fragments.rt)
    end
    
    return rt_predictions
end

"""
Helper function to predict RTs using local Chronologer installation.
"""
function predict_rt_local(chronologer_out_path::String)
    # Convert Arrow to TSV for Chronologer
    chronologer_out_tsv = replace(chronologer_out_path, r"\.arrow$" => ".tsv")
    
    # Read data and convert UniMod codes
    chronologer_table = DataFrame(Tables.columntable(Arrow.Table(chronologer_out_path)))
    unimod_df = DataFrame(CSV.File(joinpath(@__DIR__, "../../../../chronologer/data/UniModToMass.txt")))
    unimod_dict = Dict(zip(unimod_df[!, :name], unimod_df[!, :mz]))
    
    # Replace UniMod codes with mass values
    chronologer_table[!,"chronologer_sequence"] = replace_unimod_codes(chronologer_table[!, "koina_sequence"], unimod_dict)
    # Write TSV for Chronologer
    CSV.write(chronologer_out_tsv, chronologer_table, delim='\t')

    # Run local Chronologer
    chronologer_script = joinpath(@__DIR__, "../../../../chronologer/Predict_RT.py")
    run(`python3.9 $chronologer_script $chronologer_out_tsv $chronologer_out_tsv`)
    
    # Read results back and convert to Arrow
    Arrow.write(chronologer_out_path, DataFrame(CSV.File(chronologer_out_tsv)))
    
    # Clean up
    rm(chronologer_out_tsv)
end

"""
Helper function to replace UniMod codes with mass values.
"""
function replace_unimod_codes(
    sequences::Vector{String},
    unimod_dict::Dict{String15, Float64}
)
    new_sequences = copy(sequences)
    for (code, mass) in unimod_dict
        for (i, sequence) in enumerate(new_sequences)
            pattern = "[" * code * "]"
            replacement = "[+" * string(round(mass, digits=8)) * "]"
            new_sequences[i] = replace(sequence, pattern => replacement)
        end
    end
    return new_sequences
end