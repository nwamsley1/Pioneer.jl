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

#=
Model training dispatch functions.
Extracted from percolatorSortOf.jl lines 807-830.
Uses AbstractPSMContainer for data access.
=#

"""
    train_model(model::PSMScoringModel, psms::AbstractPSMContainer, features, num_rounds) -> TrainedModel

Train a PSM scoring model on the given data.
Returns a trained model object that can be used for prediction.
"""
function train_model(model::LightGBMScorer, psms::AbstractPSMContainer,
                     features::Vector{Symbol}, num_rounds::Int)
    hp = model.hyperparams

    classifier = build_lightgbm_classifier(
        num_iterations = num_rounds,
        max_depth = get(hp, :max_depth, 10),
        learning_rate = get(hp, :learning_rate, 0.05),
        num_leaves = get(hp, :num_leaves, 63),
        feature_fraction = get(hp, :feature_fraction, 0.5),
        bagging_fraction = get(hp, :bagging_fraction, 0.25),
        bagging_freq = get(hp, :bagging_fraction, 0.25) < 1 ? 1 : 0,
        min_data_in_leaf = get(hp, :min_data_in_leaf, 500),
        min_gain_to_split = get(hp, :min_gain_to_split, 0.5),
        is_unbalance = true
    )

    # Convert to DataFrame for LightGBM
    df = to_dataframe(psms)
    feature_frame = df[:, features]
    labels = collect(get_column(psms, :target))

    return fit_lightgbm_model(classifier, feature_frame, labels; positive_label=true)
end

function train_model(model::ProbitScorer, psms::AbstractPSMContainer,
                     features::Vector{Symbol}, ::Int)
    y = get_labels(psms)
    X = get_feature_matrix(psms, features)

    n = nrows(psms)
    chunk_size = max(1, n ÷ (10 * Threads.nthreads()))
    data_chunks = Iterators.partition(1:n, chunk_size)

    beta = zeros(Float64, length(features))
    beta = ProbitRegression(beta, Matrix{Float64}(X), y, data_chunks;
                           max_iter=model.max_iter)

    return ProbitModel(beta, features)
end

"""
    predict_scores(model, psms::AbstractPSMContainer) -> Vector{Float32}

Get probability predictions from a trained model.
"""
function predict_scores(model::LightGBMModel, psms::AbstractPSMContainer)
    df = to_dataframe(psms)
    return lightgbm_predict(model, df; output_type=Float32)
end

function predict_scores(model::ProbitModel, psms::AbstractPSMContainer)
    n = nrows(psms)
    scores = zeros(Float32, n)
    chunk_size = max(1, n ÷ (10 * Threads.nthreads()))
    data_chunks = Iterators.partition(1:n, chunk_size)
    X = get_feature_matrix(psms, model.features)
    ModelPredictProbs!(scores, Matrix{Float64}(X), model.beta, data_chunks)
    return scores
end

# DataFrame convenience wrappers for backward compatibility
function train_model(model::PSMScoringModel, psms::DataFrame,
                     features::Vector{Symbol}, num_rounds::Int)
    container = DataFramePSMContainer(psms, Val(:unsafe))
    return train_model(model, container, features, num_rounds)
end

function predict_scores(model, psms::DataFrame)
    container = DataFramePSMContainer(psms, Val(:unsafe))
    return predict_scores(model, container)
end
