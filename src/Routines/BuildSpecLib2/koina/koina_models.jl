"""
Prepare batch request for instrument-specific models.
"""
function prepare_koina_batch(model::InstrumentSpecificModel,
                           data::DataFrame,
                           instrument_type::String;
                           batch_size::Int = 1000)::Vector{String}
    
    if instrument_type ∉ MODEL_CONFIGS[model.name].instruments
        throw(ArgumentError("Invalid instrument type '$instrument_type' for model '$(model.name)'"))
    end
    
    n_rows = nrow(data)
    n_batches = ceil(Int, n_rows/batch_size)
    json_batches = Vector{String}(undef, n_batches)
    
    for i in 1:n_batches
        start_idx = (i-1) * batch_size + 1
        end_idx = min(i * batch_size, n_rows)
        batch_data = data[start_idx:end_idx, :]
        batch_dict = Dict(
            "id" => "0",
            "inputs" => [
                Dict(
                    "name" => "peptide_sequences",
                    "shape" => [nrow(batch_data), 1],
                    "datatype" => "BYTES",
                    "data" => batch_data.koina_sequence
                ),
                Dict(
                    "name" => "precursor_charges",
                    "shape" => [nrow(batch_data), 1],
                    "datatype" => "INT32",
                    "data" => batch_data.precursor_charge
                ),
                Dict(
                    "name" => "collision_energies",
                    "shape" => [nrow(batch_data), 1],
                    "datatype" => "FP32",
                    "data" => batch_data.collision_energy
                ),
                Dict(
                    "name" => "instrument_types",
                    "shape" => [nrow(batch_data), 1],
                    "datatype" => "BYTES",
                    "data" => fill(instrument_type, nrow(batch_data))
                )
            ]
        )
        json_batches[i] = JSON.json(batch_dict)
    end
    
    return json_batches
end

"""
Prepare batch request for instrument-agnostic models.
"""
function prepare_koina_batch(model::InstrumentAgnosticModel,
                           data::DataFrame;
                           batch_size::Int = 1000)::Vector{String}
    n_rows = nrow(data)
    n_batches = ceil(Int, n_rows/batch_size)
    json_batches = Vector{String}(undef, n_batches)
    
    for i in 1:n_batches
        start_idx = (i-1) * batch_size + 1
        end_idx = min(i * batch_size, n_rows)
        batch_data = data[start_idx:end_idx, :]
        batch_dict = Dict(
            "id" => "0",
            "inputs" => [
                Dict(
                    "name" => "peptide_sequences",
                    "shape" => [nrow(batch_data), 1],
                    "datatype" => "BYTES",
                    "data" => batch_data.koina_sequence
                ),
                Dict(
                    "name" => "precursor_charges",
                    "shape" => [nrow(batch_data), 1],
                    "datatype" => "INT32",
                    "data" => batch_data.precursor_charge
                ),
                Dict(
                    "name" => "collision_energies",
                    "shape" => [nrow(batch_data), 1],
                    "datatype" => "FP32",
                    "data" => batch_data.collision_energy
                )
            ]
        )
        json_batches[i] = JSON.json(batch_dict)
    end
    
    return json_batches
end

"""
Prepare batch request for retention time predictions.
"""
function prepare_koina_batch(model::RetentionTimeModel,
                           data::DataFrame;
                           batch_size::Int = 1000)::Vector{String}
    
    n_rows = nrow(data)
    n_batches = ceil(Int, n_rows/batch_size)
    json_batches = Vector{String}(undef, n_batches)
    
    for i in 1:n_batches
        start_idx = (i-1) * batch_size + 1
        end_idx = min(i * batch_size, n_rows)
        batch_data = data[start_idx:end_idx, :]
        batch_dict = Dict(
            "id" => "0",
            "inputs" => [
                Dict(
                    "name" => "peptide_sequences",
                    "shape" => [nrow(batch_data), 1],
                    "datatype" => "BYTES",
                    "data" => strip.(batch_data.chronologer_sequence)
                )
            ]
        )
        json_batches[i] = JSON.json(batch_dict)
    end
    
    return json_batches
end