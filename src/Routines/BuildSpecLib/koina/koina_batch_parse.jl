"""
Parse results for standard intensity prediction models.
"""
function parse_koina_batch(model::InstrumentSpecificModel,
                          response::Dict{String,Any})::KoinaBatchResult{Nothing}
    df = DataFrame()
    n_precs, n_frags = first(response["outputs"])["shape"]
    
    for col in response["outputs"]
        col_name = Symbol(col["name"])
        if col_name ∈ [:intensities, :mz]
            df[!, col_name] = Float32.(col["data"])
        else
            df[!, :annotation] = string.(col["data"])
        end
    end
    
    return KoinaBatchResult(df, Int64(n_frags), nothing)
end

"""
Parse results for instrument-agnostic models.
"""
function parse_koina_batch(model::InstrumentAgnosticModel,
                          response::Dict{String,Any})::KoinaBatchResult{Nothing}
    # Currently same as InstrumentSpecificModel
    parse_koina_batch(InstrumentSpecificModel(model.name), response)
end

"""
Parse results for spline coefficient models.
"""
function parse_koina_batch(model::SplineCoefficientModel,
                          response::Dict{String,Any})::KoinaBatchResult{Vector{Float32}}
    df = DataFrame()
    n_precs, n_coef_per_frag, n_frags = response["outputs"][1]["shape"]
    knot_vector = Float32[]
    
    for col in response["outputs"]
        col_name = Symbol(col["name"])
        if col_name == :coefficients
            flat_coeffs = Float32.(col["data"])
            coefs = Vector{NTuple{n_coef_per_frag, Float32}}(undef, n_precs * n_frags)
            
            idx = 1
            for i in 1:n_precs
                for j in 1:n_frags
                    coefs[idx] = ntuple(k -> 
                        flat_coeffs[(i-1)*n_coef_per_frag*n_frags + j + (k-1)*n_frags],
                        n_coef_per_frag)
                    idx += 1
                end
            end
            df[!, :coefficients] = coefs
            
        elseif col_name == :knots
            knot_vector = Float32.(col["data"])
        elseif col_name == :mz
            df[!, col_name] = Float32.(col["data"])
        elseif col_name == :annotation
            df[!, col_name] = string.(col["data"])
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