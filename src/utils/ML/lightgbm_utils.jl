# Utility helpers for working with LightGBM through the MLJ model interface.

struct LightGBMModelWrapper
    model::LightGBM.LGBMClassifier
    fitresult::Any
    cache::Any
    report::Any
    feature_names::Vector{Symbol}
    positive_label::Any
end

function build_lightgbm_classifier(; num_iterations::Integer = 100,
                                    max_depth::Integer = -1,
                                    learning_rate::Real = 0.1,
                                    feature_fraction::Real = 1.0,
                                    bagging_fraction::Real = 1.0,
                                    bagging_freq::Integer = 0,
                                    min_child_weight::Integer = 1,
                                    min_gain_to_split::Real = 0.0)
    return LightGBM.LGBMClassifier(
        objective = "binary",
        num_iterations = Int(num_iterations),
        max_depth = Int(max_depth),
        learning_rate = float(learning_rate),
        feature_fraction = float(feature_fraction),
        bagging_fraction = float(bagging_fraction),
        bagging_freq = Int(bagging_freq),
        min_child_weight = Int(min_child_weight),
        min_gain_to_split = float(min_gain_to_split),
    )
end

function fit_lightgbm_model(model::LightGBM.LGBMClassifier,
                            feature_data::AbstractDataFrame,
                            labels::AbstractVector;
                            positive_label = true)
    label_vec = collect(Bool.(labels))
    fitresult, cache, report = fit(model, 0, feature_data, label_vec)
    return LightGBMModelWrapper(model, fitresult, cache, report, Symbol.(names(feature_data)), positive_label)
end

function lightgbm_predict(wrapper::LightGBMModelWrapper,
                          feature_data::AbstractDataFrame;
                          output_type = Float64)
    data_view = feature_data[:, wrapper.feature_names]
    raw_predictions = predict(wrapper.model, wrapper.fitresult, data_view)
    probabilities = _lightgbm_probabilities(raw_predictions, wrapper.positive_label)
    return convert.(output_type, probabilities)
end

function _lightgbm_probabilities(predictions::AbstractVector{<:UnivariateFinite}, positive_label)
    return [pdf(p, positive_label) for p in predictions]
end

function _lightgbm_probabilities(predictions::AbstractVector{<:Real}, _)
    return Float64.(predictions)
end

function _lightgbm_probabilities(predictions::AbstractMatrix{<:Real}, _)
    if size(predictions, 2) == 1
        return vec(predictions)
    else
        return vec(predictions[:, end])
    end
end

function _lightgbm_probabilities(predictions, _)
    return Float64.(collect(predictions))
end

function lightgbm_feature_importances(wrapper::LightGBMModelWrapper)
    try
        if hasmethod(LightGBM.feature_importances, Tuple{typeof(wrapper.fitresult)})
            return LightGBM.feature_importances(wrapper.fitresult)
        elseif hasmethod(LightGBM.feature_importances, Tuple{typeof(wrapper.model), typeof(wrapper.fitresult)})
            return LightGBM.feature_importances(wrapper.model, wrapper.fitresult)
        end
    catch
        return nothing
    end
    return nothing
end
