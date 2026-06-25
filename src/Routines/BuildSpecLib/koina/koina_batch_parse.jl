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

# Reshape the flat (n_precs × K × F) coefficient array into one NTuple{K} per
# fragment. Function barrier specialized on K (Val) so the hot loop is type
# stable: `flat` is a concrete Vector{Float32} and `ntuple(_, Val(K))` returns a
# concrete NTuple{K,Float32}, eliminating the per-element boxing that occurs when
# the loop runs over the `Any`-typed response data directly.
@inline function _build_coefs(flat::Vector{Float32}, n_precs::Int, n_frags::Int,
                              ::Val{K}) where {K}
    coefs = Vector{NTuple{K, Float32}}(undef, n_precs * n_frags)
    idx = 1
    @inbounds for i in 1:n_precs
        base_i = (i - 1) * K * n_frags
        for j in 1:n_frags
            coefs[idx] = ntuple(k -> flat[base_i + j + (k - 1) * n_frags], Val(K))
            idx += 1
        end
    end
    return coefs
end

"""
Parse results for spline coefficient models.
"""
function parse_koina_batch(model::SplineCoefficientModel,
                          response::Dict{String,Any})::KoinaBatchResult{Vector{Float32}}
    df = DataFrame()
    # Look up the coefficients output by name rather than relying on a fixed
    # array index — Koina is free to reorder outputs and synthetic clients
    # may emit them in a different order than the live service.
    coefficients_output = nothing
    for o in response["outputs"]
        if o["name"] == "coefficients"
            coefficients_output = o
            break
        end
    end
    coefficients_output === nothing && throw(ArgumentError(
        "SplineCoefficientModel response is missing a 'coefficients' output"))
    # shape comes from a Dict{String,Any}; make the dims concrete Ints so the
    # loops/typing below are stable.
    n_precs, n_coef_per_frag, n_frags = Int.(coefficients_output["shape"])
    knot_vector = Float32[]

    for col in response["outputs"]
        col_name = Symbol(col["name"])
        if col_name == :coefficients
            # `col["data"]` is Any; assert the converted result is Vector{Float32}
            # so the reshape loop (in _build_coefs) is allocation-free.
            flat_coeffs = Float32.(col["data"])::Vector{Float32}
            df[!, :coefficients] = _build_coefs(flat_coeffs, n_precs, n_frags,
                                                Val(n_coef_per_frag))
        elseif col_name == :knots
            knot_vector = Float32.(col["data"])::Vector{Float32}
        elseif col_name == :mz
            df[!, :mz] = Float32.(col["data"])::Vector{Float32}
        elseif col_name == :annotations
            df[!, :annotation] = Int32.(col["data"])::Vector{Int32}
        end
    end

    return KoinaBatchResult(df, Int64(n_frags), knot_vector)
end

"""
Parse results for retention time prediction models.
"""
function parse_koina_batch(model::RetentionTimeModel,
                          response::Dict{String,Any})::KoinaBatchResult{Nothing}
    df = DataFrame()

    for col in response["outputs"]
        col_name = Symbol(col["name"])
        if col_name == :rt
            df[!, col_name] = Float32.(col["data"])
        end
    end

    return KoinaBatchResult(df, 1, nothing)
end
