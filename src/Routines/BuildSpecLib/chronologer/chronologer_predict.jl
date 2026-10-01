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
    predict_ion_mobility(in_path, out_path, im_model)

Read the precursor table at `in_path`, add `ccs` (Å², from the Koina model
`im_model`) and `inv_ion_mobility` (1/K0, Vs/cm²) columns, and write it to
`out_path`. Written to a new file rather than in place for the same reason
`predict_retention_times` is: the input Arrow file may still be mmap-locked.
"""
function predict_ion_mobility(in_path::String, out_path::String, im_model::String)
    table = DataFrame(Tables.columntable(Arrow.Table(in_path)))
    ccs = predict_ccs_koina(table, im_model)
    table[!, :ccs] = ccs
    table[!, :inv_ion_mobility] =
        ccs_to_inv_ion_mobility.(ccs, table.precursor_charge, table.mz)
    Arrow.write(out_path, table)
    return
end

"""
Helper function to predict CCS values using a Koina ion-mobility model.
"""
function predict_ccs_koina(table::DataFrame, im_model::String)::Vector{Float32}
    model = IonMobilityModel(im_model)
    batches = prepare_koina_batch(model, table, batch_size=1000)
    results = make_koina_batch_requests(batches, KOINA_URLS[im_model])
    ccs = Float32[]
    for result in results
        append!(ccs, parse_koina_batch(model, result).fragments.ccs)
    end
    length(ccs) == nrow(table) || error(
        "ion-mobility model $im_model returned $(length(ccs)) values for $(nrow(table)) precursors")
    return ccs
end

"""
    ccs_to_inv_ion_mobility(ccs, charge, mz)

Mason–Schamp conversion as used by AlphaPeptDeep for Bruker timsTOF data
(N2 drift gas, 28 Da): `1/K0 = CCS · sqrt(μ) / (z · 1059.62245)`, where
`μ = M·28/(M+28)` is the reduced mass of the ion (`M = mz·z`) and N2.
"""
function ccs_to_inv_ion_mobility(ccs::Real, charge::Integer, mz::Real)::Float32
    M = Float64(mz) * charge
    μ = M * 28.0 / (M + 28.0)
    return Float32(Float64(ccs) * sqrt(μ) / (charge * 1059.62245))
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