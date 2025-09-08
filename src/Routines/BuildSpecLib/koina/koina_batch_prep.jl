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
Prepare batch request for instrument-agnostic models (with ignored instrument parameter).
"""
function prepare_koina_batch(model::InstrumentAgnosticModel,
                           data::DataFrame,
                           instrument_type::String;  # ignored for agnostic models
                           batch_size::Int = 1000)::Vector{String}
    # Call the generic method that doesn't need instrument_type
    return prepare_koina_batch(model, data; batch_size=batch_size)
end

"""
Prepare batch request for instrument-agnostic models.
"""
function prepare_koina_batch(model::KoinaModelType,
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
Prepare batch request for instrument-agnostic models.
"""
function prepare_koina_batch(model::SplineCoefficientModel,
                           data::DataFrame,
                           ::String;
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
                )
            ]
        )
        json_batches[i] = JSON.json(batch_dict)
    end
    
    return json_batches
end

"""
Prepare batch request for retention time predictions (with ignored instrument parameter).
"""
function prepare_koina_batch(model::RetentionTimeModel,
                           data::DataFrame,
                           instrument_type::String;  # ignored for RT models
                           batch_size::Int = 1000)::Vector{String}
    # Call the generic method that doesn't need instrument_type
    return prepare_koina_batch(model, data; batch_size=batch_size)
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
                    "data" => strip.(batch_data.koina_sequence)
                )
            ]
        )
        json_batches[i] = JSON.json(batch_dict)
    end
    
    return json_batches
end