# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Pioneer.jl is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

const GLOBAL_PROTEIN_SCORE_FEATURES = Symbol[
    :empirical_global_score,
    :top1_pg_score,
    :top2_pg_score,
    :top3_pg_score,
    :top2_logodds_score,
    :top3_logodds_score,
    :mean_pg_score,
    :median_pg_score,
    :std_pg_score,
    :min_pg_score,
    :top1_top2_gap,
    :top2_top3_gap,
    :n_unique_peptides_observed,
    :max_n_peptides_observed,
    :global_peptide_coverage,
    :max_peptide_coverage,
    :n_passing_runs,
    :n_score_gt_0_5,
    :n_score_gt_0_9,
    :n_score_gt_0_99,
    :observed_run_centrality_mean,
    :observed_run_centrality_max,
    :missing_run_similarity_mass_approx,
]

const GLOBAL_PROTEIN_LGBM_HP = (
    num_iterations = 100,
    learning_rate = 0.05,
    max_depth = 4,
    num_leaves = 15,
    min_data_in_leaf = 10,
    min_gain_to_split = 0.0,
    feature_fraction = 1.0,
    bagging_fraction = 0.8,
    bagging_freq = 1,
    is_unbalance = true,
    max_bin = 255,
    lambda_l1 = 1.0,
    lambda_l2 = 1.0,
)

const GLOBAL_PROTEIN_MIN_TRAINING_CLASS_COUNT = 100
const GLOBAL_PROTEIN_MAX_TRAIN = 1_000_000
const GlobalProteinKey = Tuple{String, Bool, UInt8}

struct GlobalProteinRunScore
    ms_file_idx::UInt32
    score::Float32
end

struct GlobalProteinInputs
    run_scores::Dict{GlobalProteinKey, Vector{GlobalProteinRunScore}}
    observed_peptides::Dict{GlobalProteinKey, Set{String}}
    max_n_peptides::Dict{GlobalProteinKey, Int}
    n_possible_unique_peptides::Dict{GlobalProteinKey, Int}
    folds::Dict{GlobalProteinKey, UInt8}
end

function _collect_global_protein_inputs(
    pg_refs::Vector{ProteinGroupFileReference},
    protein_to_cv_fold::Dictionary{
        String,
        @NamedTuple{best_score::Float32, cv_fold::UInt8},
    },
    n_proteins::Int,
)
    run_scores = Dict{GlobalProteinKey, Vector{GlobalProteinRunScore}}()
    observed_peptides = Dict{GlobalProteinKey, Set{String}}()
    max_n_peptides = Dict{GlobalProteinKey, Int}()
    n_possible_unique_peptides = Dict{GlobalProteinKey, Int}()
    folds = Dict{GlobalProteinKey, UInt8}()
    sizehint!(run_scores, n_proteins)
    sizehint!(observed_peptides, n_proteins)
    sizehint!(max_n_peptides, n_proteins)
    sizehint!(n_possible_unique_peptides, n_proteins)
    sizehint!(folds, n_proteins)

    for ref in pg_refs
        table = Arrow.Table(file_path(ref))
        @inbounds for row in eachindex(table.protein_name)
            protein_name = String(table.protein_name[row])
            key = (
                protein_name,
                Bool(table.target[row]),
                UInt8(table.entrap_id[row]),
            )
            protein_run_scores = get!(run_scores, key) do
                GlobalProteinRunScore[]
            end
            push!(
                protein_run_scores,
                GlobalProteinRunScore(
                    UInt32(table.file_idx[row]),
                    Float32(table.pg_score[row]),
                ),
            )

            protein_peptides = get!(observed_peptides, key) do
                Set{String}()
            end
            for peptide in split(table.peptide_list[row], ';')
                push!(protein_peptides, String(peptide))
            end

            n_peptides = Int(table.n_peptides[row])
            max_n_peptides[key] = max(get(max_n_peptides, key, 0), n_peptides)
            n_possible_unique_peptides[key] =
                Int(table.n_possible_unique_peptides[row])
            folds[key] = protein_to_cv_fold[protein_name].cv_fold
        end
    end

    isempty(run_scores) &&
        throw(ArgumentError("Global protein scoring requires observed protein groups"))
    return GlobalProteinInputs(
        run_scores,
        observed_peptides,
        max_n_peptides,
        n_possible_unique_peptides,
        folds,
    )
end

function _build_global_protein_feature_table(
    inputs::GlobalProteinInputs,
    n_runs_total::Int,
    run_score_floor::Float32,
    ;
    n_experiment_runs::Int = n_runs_total,
    run_similarity::Union{Nothing, RunSimilarityAtlas} = nothing,
)
    protein_keys = sort!(collect(keys(inputs.run_scores)))
    n_proteins = length(protein_keys)
    top_run_count = isqrt(n_runs_total)
    feature_columns = Dict(
        feature => Vector{Float32}(undef, n_proteins)
        for feature in GLOBAL_PROTEIN_SCORE_FEATURES
    )

    @inbounds for (row, key) in enumerate(protein_keys)
        protein_run_scores = inputs.run_scores[key]
        scores = Float32[
            run_score.score
            for run_score in protein_run_scores
        ]
        sorted_scores = sort(scores; rev = true)
        n_observed = length(sorted_scores)
        top1 = sorted_scores[1]
        top2 = n_observed >= 2 ? sorted_scores[2] : 0.0f0
        top3 = n_observed >= 3 ? sorted_scores[3] : 0.0f0
        n_unique_peptides = length(inputs.observed_peptides[key])
        max_n_peptides = inputs.max_n_peptides[key]
        n_possible_unique_peptides =
            inputs.n_possible_unique_peptides[key]

        feature_columns[:empirical_global_score][row] =
            _logodds_from_sorted(sorted_scores, top_run_count)
        feature_columns[:top1_pg_score][row] = top1
        feature_columns[:top2_pg_score][row] = top2
        feature_columns[:top3_pg_score][row] = top3
        feature_columns[:top2_logodds_score][row] =
            _logodds_from_sorted(sorted_scores, 2)
        feature_columns[:top3_logodds_score][row] =
            _logodds_from_sorted(sorted_scores, 3)
        feature_columns[:mean_pg_score][row] = Float32(mean(scores))
        feature_columns[:median_pg_score][row] = Float32(median(scores))
        feature_columns[:std_pg_score][row] =
            n_observed > 1 ? Float32(std(scores)) : 0.0f0
        feature_columns[:min_pg_score][row] = minimum(scores)
        feature_columns[:top1_top2_gap][row] =
            n_observed >= 2 ? top1 - top2 : 0.0f0
        feature_columns[:top2_top3_gap][row] =
            n_observed >= 3 ? top2 - top3 : 0.0f0
        feature_columns[:n_unique_peptides_observed][row] = Float32(n_unique_peptides)
        feature_columns[:max_n_peptides_observed][row] = Float32(max_n_peptides)
        feature_columns[:global_peptide_coverage][row] =
            Float32(n_unique_peptides) / Float32(n_possible_unique_peptides)
        feature_columns[:max_peptide_coverage][row] =
            Float32(max_n_peptides) / Float32(n_possible_unique_peptides)
        feature_columns[:n_score_gt_0_5][row] = Float32(count(>(0.5f0), scores))
        feature_columns[:n_score_gt_0_9][row] = Float32(count(>(0.9f0), scores))
        feature_columns[:n_score_gt_0_99][row] = Float32(count(>(0.99f0), scores))

        centrality_sum = 0.0f0
        centrality_max = 0.0f0
        n_passing_runs = 0
        for run_score in protein_run_scores
            run_score.score >= run_score_floor || continue
            n_passing_runs += 1
            centrality = run_similarity === nothing ?
                0.0f0 :
                run_centrality(run_similarity, run_score.ms_file_idx)
            centrality_sum += centrality
            centrality_max = max(centrality_max, centrality)
        end
        centrality_mean = n_passing_runs > 0 ?
            centrality_sum / Float32(n_passing_runs) : 0.0f0
        n_missing_runs = max(n_experiment_runs - n_passing_runs, 0)
        feature_columns[:n_passing_runs][row] = Float32(n_passing_runs)
        feature_columns[:observed_run_centrality_mean][row] = centrality_mean
        feature_columns[:observed_run_centrality_max][row] = centrality_max
        feature_columns[:missing_run_similarity_mass_approx][row] =
            centrality_mean * Float32(n_missing_runs)
    end

    table = DataFrame(
        protein_name = String[key[1] for key in protein_keys],
        target = Bool[key[2] for key in protein_keys],
        entrap_id = UInt8[key[3] for key in protein_keys],
        cv_fold = UInt8[inputs.folds[key] for key in protein_keys],
    )
    for feature in GLOBAL_PROTEIN_SCORE_FEATURES
        table[!, feature] = feature_columns[feature]
    end
    return table
end

function _build_protein_score_dicts(
    protein_keys::Vector{GlobalProteinKey},
    scores::AbstractVector{<:Real},
)
    global_score_dict = Dict{GlobalProteinKey, Float32}()
    protein_key_score_dict = Dict{ProteinKey, Float32}()
    sizehint!(global_score_dict, length(protein_keys))
    sizehint!(protein_key_score_dict, length(protein_keys))

    @inbounds for row in eachindex(protein_keys)
        key = protein_keys[row]
        score = Float32(scores[row])
        global_score_dict[key] = score
        protein_key_score_dict[ProteinKey(key...)] = score
    end
    return global_score_dict, protein_key_score_dict
end

function _build_empirical_global_protein_score_dicts(
    inputs::GlobalProteinInputs,
    top_run_count::Int,
)
    protein_keys = sort!(collect(keys(inputs.run_scores)))
    scores = Float32[
        logodds(
            Float32[
                run_score.score
                for run_score in inputs.run_scores[key]
            ],
            top_run_count,
        )
        for key in protein_keys
    ]
    return _build_protein_score_dicts(protein_keys, scores)
end

function _global_protein_training_data_sufficient(
    inputs::GlobalProteinInputs,
    cv_folds::Vector{UInt8},
)
    for test_fold in cv_folds
        n_targets = 0
        n_decoys = 0
        for key in keys(inputs.run_scores)
            inputs.folds[key] == test_fold && continue
            if key[2]
                n_targets += 1
            else
                n_decoys += 1
            end
        end
        min(n_targets, n_decoys) >= GLOBAL_PROTEIN_MIN_TRAINING_CLASS_COUNT ||
            return false
    end
    return true
end

function _global_protein_model_eligible(
    inputs::GlobalProteinInputs,
    cv_folds::Vector{UInt8},
    protein_probit_scoring_succeeded::Bool,
)
    return protein_probit_scoring_succeeded &&
           length(cv_folds) > 1 &&
           _global_protein_training_data_sufficient(inputs, cv_folds)
end

"""
    build_global_protein_score_dicts(
        pg_refs,
        protein_to_cv_fold,
        n_proteins,
        n_runs_total,
        protein_probit_scoring_succeeded,
        run_score_floor,
    )

Build one score-distribution, peptide-support, and run-context feature row per
protein group and select out-of-fold global LightGBM scores only when they
retain more total protein-run target IDs after applying both the configured
global and experiment-wide q-value thresholds than the empirical top-run
log-odds scores. Use the empirical score without training when probit scoring
was not applied, the observed protein groups belong to one CV fold, or any
training split has fewer than 100 targets or 100 decoys.
"""
function build_global_protein_score_dicts(
    pg_refs::Vector{ProteinGroupFileReference},
    protein_to_cv_fold::Dictionary{
        String,
        @NamedTuple{best_score::Float32, cv_fold::UInt8},
    },
    n_proteins::Int,
    n_runs_total::Int,
    protein_probit_scoring_succeeded::Bool,
    run_score_floor::Float32,
    ;
    n_experiment_runs::Int = n_runs_total,
    run_similarity::Union{Nothing, RunSimilarityAtlas} = nothing,
    q_value_threshold::Float32 = SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD,
)
    inputs = _collect_global_protein_inputs(
        pg_refs,
        protein_to_cv_fold,
        n_proteins,
    )
    cv_folds = sort!(unique(collect(values(inputs.folds))))
    if !_global_protein_model_eligible(
        inputs,
        cv_folds,
        protein_probit_scoring_succeeded,
    )
        return _build_empirical_global_protein_score_dicts(
            inputs,
            isqrt(n_runs_total),
        )
    end

    @user_info "Training global protein scoring model..."
    table = _build_global_protein_feature_table(
        inputs,
        n_runs_total,
        run_score_floor,
        n_experiment_runs = n_experiment_runs,
        run_similarity = run_similarity,
    )
    scored = _score_global_features_oof(
        table,
        GLOBAL_PROTEIN_SCORE_FEATURES,
        cv_folds;
        scoring_name = "Global protein",
        lgbm_hp = GLOBAL_PROTEIN_LGBM_HP,
        min_training_class_count = GLOBAL_PROTEIN_MIN_TRAINING_CLASS_COUNT,
        max_train = GLOBAL_PROTEIN_MAX_TRAIN,
    )
    selected = _select_global_scores(
        scored.scores,
        table.empirical_global_score,
        table.target;
        scoring_name = "Global protein",
        experiment_wide_id_counts = Int.(table.n_passing_runs),
        q_threshold = q_value_threshold,
    )

    protein_keys = GlobalProteinKey[
        (table.protein_name[row], table.target[row], table.entrap_id[row])
        for row in axes(table, 1)
    ]
    score_dicts = _build_protein_score_dicts(protein_keys, selected.scores)
    @debug_l1 "Global protein scoring scored $(length(protein_keys)) protein groups " *
              "with $(length(GLOBAL_PROTEIN_SCORE_FEATURES)) features; " *
              "selected $(selected.source) after semi-supervised iteration $(scored.iter)"
    return score_dicts
end
