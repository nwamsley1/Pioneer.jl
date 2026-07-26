# Utility helpers for working with LightGBM directly through its native API.

struct LightGBMModel
    booster::Union{LightGBM.LGBMClassification, LightGBM.LGBMRanking, Nothing}
    features::Vector{Symbol}
    constant_prediction::Union{Float32, Nothing}
end

const LightGBMModelVector = Vector{LightGBMModel}
const LightGBMModelWrapper = LightGBMModel

# Serializes every entry into LightGBM's C API. LightGBM's bundled OpenMP
# (libomp) aborts in __kmpc_end_serialized_parallel when predict is entered
# concurrently from multiple Julia threads (pass1_oom.jl's parallel_foreach!).
const LGBM_C_LOCK = ReentrantLock()

"""
    feature_matrix(df, features) -> Matrix{Float32}

Construct a dense matrix with the columns in `features` converted to `Float32`.
Missing values are imputed with `0.0f0`. Columns are filled in parallel via a
typed inner helper to avoid Float64 intermediates and per-column broadcast
temporaries.
"""
function _fill_column!(M::Matrix{Float32}, j::Int, col::AbstractVector{T}) where {T<:Real}
    @inbounds @simd for i in eachindex(col)
        M[i, j] = Float32(col[i])
    end
end

function _fill_column!(M::Matrix{Float32}, j::Int, col::AbstractVector{<:Union{Missing, T}}) where {T<:Real}
    @inbounds for i in eachindex(col)
        v = col[i]
        M[i, j] = v === missing ? 0.0f0 : Float32(v)
    end
end

_fill_column!(::Matrix{Float32}, ::Int, col::AbstractVector) =
    throw(ArgumentError("Unsupported feature type $(eltype(col)) for LightGBM"))

function feature_matrix(df::AbstractDataFrame, features::Vector{Symbol})
    n = nrow(df)
    m = length(features)
    matrix = Matrix{Float32}(undef, n, m)
    cols = AbstractVector[df[!, f] for f in features]
    Threads.@threads for j in 1:m
        _fill_column!(matrix, j, cols[j])
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
                                    min_data_in_leaf::Integer = 500,
                                    min_gain_to_split::Real = 1.0,
                                    lambda_l1::Real = 0.0,
                                    lambda_l2::Real = 0.0,
                                    max_bin::Integer = 255,
                                    monotone_constraints::AbstractVector{<:Integer} = Int[],
                                    monotone_constraints_method::AbstractString = "basic",
                                    num_threads::Integer = Threads.nthreads(),
                                    metric = ["binary_logloss"],
                                    objective::AbstractString = "binary",
                                    is_unbalance = false,
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
        min_data_in_leaf = Int(min_data_in_leaf),
        min_gain_to_split = float(min_gain_to_split),
        lambda_l1 = float(lambda_l1),
        lambda_l2 = float(lambda_l2),
        max_bin = Int(max_bin),
        monotone_constraints = Int.(monotone_constraints),
        monotone_constraints_method = String(monotone_constraints_method),
        num_threads = Int(num_threads),
        num_class = 1,
        verbosity = Int(verbosity),
        is_unbalance = is_unbalance,
        seed = 1776, # potentialy needed for stable results
        deterministic = true, # potentialy needed for stable results
        force_row_wise = true # potentialy needed for stable results
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
