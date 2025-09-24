# Utility helpers for working with LightGBM directly through its native API.

struct LightGBMModel
    booster::Union{LightGBM.LGBMClassification, Nothing}
    features::Vector{Symbol}
    constant_prediction::Union{Float32, Nothing}
end

const LightGBMModelVector = Vector{LightGBMModel}
const LightGBMModelWrapper = LightGBMModel

"""
    feature_matrix(df, features) -> Matrix{Float32}

Construct a dense matrix with the columns in `features` converted to `Float32`.
Missing values are imputed with sensible defaults per type.
"""
function feature_matrix(df::AbstractDataFrame, features::Vector{Symbol})
    n = nrow(df)
    m = length(features)
    matrix = Matrix{Float32}(undef, n, m)

    for (j, feat) in enumerate(features)
        column = df[!, feat]
        T = nonmissingtype(eltype(column))

        if T <: AbstractFloat
            if eltype(column) <: Union{Missing, T}
                matrix[:, j] = Float32.(coalesce.(column, zero(T)))
            else
                matrix[:, j] = Float32.(column)
            end
        elseif T <: Integer
            if eltype(column) <: Union{Missing, T}
                matrix[:, j] = Float32.(coalesce.(column, zero(T)))
            else
                matrix[:, j] = Float32.(column)
            end
        elseif T <: Bool
            if eltype(column) <: Union{Missing, Bool}
                matrix[:, j] = Float32.(coalesce.(column, false))
            else
                matrix[:, j] = Float32.(column)
            end
        else
            throw(ArgumentError("Unsupported feature type $(eltype(column)) for LightGBM"))
        end
    end

    return matrix
end

function build_lightgbm_classifier(; num_iterations::Integer = 100,
                                    max_depth::Integer = -1,
                                    num_leaves::Integer = 31,
                                    learning_rate::Real = 0.05,
                                    feature_fraction::Real = 0.5,
                                    bagging_fraction::Real = 1.0,
                                    bagging_freq::Integer = 0,
                                    min_child_weight::Integer = 500,
                                    min_gain_to_split::Real = 1.0,
                                    lambda_l2::Real = 0.0,
                                    max_bin::Integer = 255,
                                    num_threads::Integer = Threads.nthreads(),
                                    metric = ["binary_logloss"],
                                    objective::AbstractString = "binary",
                                    verbosity::Integer = -1)
    return LightGBM.LGBMClassification(
        objective = objective,
        metric = metric,
        learning_rate = float(learning_rate),
        num_iterations = Int(num_iterations),
        num_leaves = Int(num_leaves),
        max_depth = Int(max_depth),
        feature_fraction = float(feature_fraction),
        bagging_fraction = float(bagging_fraction),
        bagging_freq = Int(bagging_freq),
        min_data_in_leaf = Int(min_child_weight),
        min_gain_to_split = float(min_gain_to_split),
        lambda_l2 = float(lambda_l2),
        max_bin = Int(max_bin),
        num_threads = Int(num_threads),
        num_class = 1,
        verbosity = Int(verbosity),
    )
end

function _prepare_labels(labels)
    label_vec = collect(labels)
    if isempty(label_vec)
        throw(ArgumentError("LightGBM requires at least one training example"))
    end

    if eltype(label_vec) <: Bool
        label_vec = Int.(label_vec)
    elseif eltype(label_vec) <: Integer
        label_vec = Int.(label_vec)
    elseif eltype(label_vec) <: AbstractFloat
        label_vec = Int.(round.(label_vec))
    else
        throw(ArgumentError("Unsupported label type $(eltype(label_vec)) for LightGBM"))
    end

    if any(x -> x ∉ (0, 1), label_vec)
        throw(ArgumentError("LightGBM requires binary labels encoded as 0 or 1."))
    end

    return label_vec
end

function fit_lightgbm_model(model::LightGBM.LGBMClassification,
                            feature_data::AbstractDataFrame,
                            labels::AbstractVector;
                            positive_label = true)
    features = Symbol.(names(feature_data))
    X = feature_matrix(feature_data, features)
    y_int = _prepare_labels(labels)

    unique_labels = unique(y_int)
    if length(unique_labels) == 1
        constant_prob = unique_labels[1] == 0 ? 0.0f0 : 1.0f0
        return LightGBMModel(nothing, features, constant_prob)
    end

    LightGBM.fit!(model, X, y_int; verbosity = -1)
    return LightGBMModel(model, features, nothing)
end

function lightgbm_predict(model::LightGBMModel,
                          feature_data::AbstractDataFrame;
                          output_type = Float64)
    if model.booster === nothing
        n = nrow(feature_data)
        prob = model.constant_prediction === nothing ? 0.0f0 : model.constant_prediction
        return fill(convert(output_type, prob), n)
    end

    X = feature_matrix(feature_data, model.features)
    raw = LightGBM.predict(model.booster, X)
    ŷ = ndims(raw) == 2 ? dropdims(raw; dims = 2) : raw
    return convert.(output_type, ŷ)
end

function lightgbm_feature_importances(model::LightGBMModel)
    if model.booster === nothing
        return nothing
    end

    try
        return LightGBM.gain_importance(model.booster)
    catch
        return nothing
    end
end

predict(model::LightGBMModel, df::AbstractDataFrame) =
    lightgbm_predict(model, df; output_type = Float32)

function importance(model::LightGBMModel)
    gains = lightgbm_feature_importances(model)
    return gains === nothing ? nothing : collect(zip(model.features, gains))
end
