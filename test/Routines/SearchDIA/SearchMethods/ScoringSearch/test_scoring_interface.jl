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

using Test
using DataFrames
using Statistics
using Pioneer  # Load the Pioneer module

# Import the actual types and functions from Pioneer
using Pioneer: MBRFilterMethod, ThresholdFilter, ProbitFilter, LightGBMFilter
using Pioneer: FilterResult, train_and_evaluate, apply_filtering
using Pioneer: select_mbr_features, prepare_mbr_features, logodds
using Pioneer: ProbitRegression, ModelPredict!


@testset "Scoring Interface Tests" begin

    @testset "MBR Running Statistics" begin
        # Test the update_pair_statistics function
        using Pioneer: update_pair_statistics

        # Test first probability
        initial_stats = (
            best_prob_1 = 0.0f0,
            best_prob_2 = 0.0f0,
            worst_prob_1 = 0.0f0,
            worst_prob_2 = 0.0f0,
            mean_prob = 0.0f0,
            count_pairs = Int32(0),
            dummy_field = 42  # Test that other fields are preserved
        )

        updated = update_pair_statistics(initial_stats, 0.8f0)
        @test updated.best_prob_1 ≈ 0.8f0
        @test updated.worst_prob_1 ≈ 0.8f0
        @test updated.mean_prob ≈ 0.8f0
        @test updated.count_pairs == 1
        @test updated.dummy_field == 42  # Preserved

        # Test second probability
        updated2 = update_pair_statistics(updated, 0.9f0)
        @test updated2.best_prob_1 ≈ 0.9f0
        @test updated2.best_prob_2 ≈ 0.8f0
        @test updated2.worst_prob_1 ≈ 0.8f0
        @test updated2.worst_prob_2 ≈ 0.9f0
        @test updated2.mean_prob ≈ 0.85f0
        @test updated2.count_pairs == 2

        # Test third probability (lower than both)
        updated3 = update_pair_statistics(updated2, 0.6f0)
        @test updated3.best_prob_1 ≈ 0.9f0
        @test updated3.best_prob_2 ≈ 0.8f0
        @test updated3.worst_prob_1 ≈ 0.6f0
        @test updated3.worst_prob_2 ≈ 0.8f0
        @test updated3.mean_prob ≈ (0.9f0 + 0.8f0 + 0.6f0) / 3.0f0
        @test updated3.count_pairs == 3

        # Test fourth probability (middle value)
        updated4 = update_pair_statistics(updated3, 0.75f0)
        @test updated4.best_prob_1 ≈ 0.9f0
        @test updated4.best_prob_2 ≈ 0.8f0
        @test updated4.worst_prob_1 ≈ 0.6f0
        @test updated4.worst_prob_2 ≈ 0.75f0
        @test updated4.count_pairs == 4
    end
    
    @testset "MBR Filter Method Traits" begin
        # Test that the trait types can be instantiated
        threshold_method = ThresholdFilter()
        probit_method = ProbitFilter()
        lightgbm_method = LightGBMFilter()

        @test threshold_method isa MBRFilterMethod
        @test probit_method isa MBRFilterMethod
        @test lightgbm_method isa MBRFilterMethod
        
        # Test FilterResult construction
        result = FilterResult("TestMethod", [0.1, 0.2, 0.3], 0.15, 2)
        @test result.method_name == "TestMethod"
        @test result.scores == [0.1, 0.2, 0.3]
        @test result.threshold == 0.15
        @test result.n_passing == 2
    end
    
    @testset "Feature Selection" begin
        # Create test DataFrame with various features
        df = DataFrame(
            prob = [0.9, 0.8, 0.7, 0.6],
            irt_error = [0.1, -0.2, 0.3, -0.1],
            MBR_max_pair_prob = [0.8, 0.9, 0.7, 0.6],
            MBR_best_irt_diff = [0.5, 0.2, 0.8, 0.3],
            some_other_feature = [1, 2, 3, 4],
            missing_feature = [missing, missing, missing, missing]
        )
        
        features = select_mbr_features(df)
        
        # Should include prob and irt_error
        @test :prob in features
        @test :irt_error in features
        @test :MBR_max_pair_prob in features
        @test :MBR_best_irt_diff in features

        # Should not include features not in candidate list
        @test :some_other_feature ∉ features

        # Should not include all-missing features
        @test :missing_feature ∉ features
    end
    
    @testset "Feature Preparation" begin
        # Create test DataFrame with missing values
        df = DataFrame(
            prob = [0.9, 0.8, missing, 0.6],
            irt_error = [0.1, -0.2, 0.3, missing]
        )
        
        X, feature_names = prepare_mbr_features(df)
        
        # Check dimensions
        @test size(X, 1) == 4  # Same number of rows
        @test size(X, 2) == 3  # intercept + 2 features
        
        # Check feature names
        @test feature_names[1] == :intercept
        @test :prob in feature_names
        @test :irt_error in feature_names
        
        # Check intercept column (should be all ones)
        @test all(X[:, 1] .== 1.0)
        
        # Check no missing values remain
        @test !any(isnan.(X))
        @test !any(ismissing.(X))
    end
    
    @testset "Threshold Method Training" begin
        # Create minimal test data
        candidate_data = DataFrame(
            prob = [0.9, 0.8, 0.7, 0.6, 0.5],
            target = [true, true, true, false, false]
        )
        candidate_labels = [false, false, true, true, false]  # bad transfer flags
        
        # Mock parameters
        params = (max_MBR_false_transfer_rate = 0.2,)
        
        method = ThresholdFilter()
        result = train_and_evaluate(method, candidate_data, candidate_labels, params)
        
        @test result !== nothing
        @test result.method_name == "Threshold"
        @test length(result.scores) == 5
        @test result.scores == candidate_data.prob
        @test result.threshold isa Real
        @test result.n_passing isa Int
        @test result.n_passing >= 0
        @test result.n_passing <= 5
    end
    
    @testset "Filter Application" begin
        # Create test data
        merged_df = DataFrame(
            prob = Float32[0.9, 0.8, 0.7, 0.6, 0.5, 0.4],
            other_col = [1, 2, 3, 4, 5, 6]
        )
        candidate_mask = [true, true, false, true, true, false]
        
        # Test threshold filtering
        threshold_result = FilterResult("Threshold", merged_df.prob, 0.65, 3)
        filtered_probs = apply_filtering(threshold_result, merged_df, candidate_mask, nothing)
        
        @test length(filtered_probs) == 6
        @test filtered_probs[1] == 0.9f0  # Above threshold, not filtered
        @test filtered_probs[2] == 0.8f0  # Above threshold, not filtered
        @test filtered_probs[3] == 0.7f0  # Not a candidate, not filtered
        @test filtered_probs[4] == 0.0f0  # Below threshold, filtered
        @test filtered_probs[5] == 0.0f0  # Below threshold, filtered
        @test filtered_probs[6] == 0.4f0  # Not a candidate, not filtered
        
        # Test ML filtering (scores correspond to candidate rows only)
        ml_scores = [0.1, 0.8, 0.2, 0.9]  # Only for candidates (indices 1,2,4,5)
        ml_result = FilterResult("Probit", ml_scores, 0.5, 2)
        filtered_probs_ml = apply_filtering(ml_result, merged_df, candidate_mask, nothing)
        
        @test length(filtered_probs_ml) == 6
        @test filtered_probs_ml[1] == 0.0f0  # Score 0.1 < 0.5, filtered (bad transfer)
        @test filtered_probs_ml[2] == 0.8f0  # Score 0.8 >= 0.5, not filtered (good transfer)
        @test filtered_probs_ml[3] == 0.7f0  # Not a candidate, not filtered
        @test filtered_probs_ml[4] == 0.0f0  # Score 0.2 < 0.5, filtered (bad transfer)
        @test filtered_probs_ml[5] == 0.5f0  # Score 0.9 >= 0.5, not filtered (good transfer)
        @test filtered_probs_ml[6] == 0.4f0  # Not a candidate, not filtered
    end
    
    @testset "Logodds Function" begin
        # Test single probability
        @test logodds([0.8], 1) == 0.8
        
        # Test multiple probabilities
        probs = [0.9, 0.8, 0.7, 0.6]
        result = logodds(probs, 2)  # Take top 2
        @test result isa Float64
        @test 0.0 <= result <= 1.0
        
        # Test edge cases
        @test 0.0 < logodds([0.0, 0.0], 2) < 1.0  # Should handle zeros
        @test 0.0 < logodds([1.0, 1.0], 2) < 1.0  # Should handle ones
        
        # Test empty case
        @test logodds(Float64[], 1) == 0.0f0  # Real function returns Float32
    end
    
    @testset "Method Training Interface" begin
        # Test that training functions exist and work
        candidate_data = DataFrame(
            prob = [0.9, 0.8, 0.7, 0.6, 0.5],
            cv_fold = [1, 2, 1, 2, 1],
            irt_error = [0.1, -0.2, 0.3, -0.1, 0.2],
            MBR_max_pair_prob = [0.8, 0.9, 0.7, 0.6, 0.5]
        )
        candidate_labels = [false, true, false, true, false]
        params = (max_MBR_false_transfer_rate = 0.2,)
        
        # Test threshold method (should always work)
        threshold_method = ThresholdFilter()
        threshold_result = train_and_evaluate(threshold_method, candidate_data, candidate_labels, params)
        @test threshold_result !== nothing
        @test threshold_result isa FilterResult
        @test threshold_result.method_name == "Threshold"
        @test length(threshold_result.scores) == nrow(candidate_data)
        
        # Test ML methods (these might return nothing if features/cv_fold are insufficient)
        ml_methods = [ProbitFilter(), LightGBMFilter()]
        
        for method in ml_methods
            result = train_and_evaluate(method, candidate_data, candidate_labels, params)
            if result !== nothing  # Some methods might return nothing due to missing requirements
                @test result isa FilterResult
                @test result.method_name isa String
                @test length(result.scores) == nrow(candidate_data)
            end
        end
    end
    
    @testset "Error Handling" begin
        # Test with problematic data
        empty_df = DataFrame()
        candidate_labels = Bool[]
        params = (max_MBR_false_transfer_rate = 0.1,)
        
        # Test threshold method with empty data
        threshold_method = ThresholdFilter()
        result = train_and_evaluate(threshold_method, empty_df, candidate_labels, params)
        @test result === nothing
        
        # Test probit method without cv_fold
        df_no_cv = DataFrame(prob = [0.8, 0.7])
        labels_no_cv = [true, false]
        probit_method = ProbitFilter()
        result_probit = train_and_evaluate(probit_method, df_no_cv, labels_no_cv, params)
        @test result_probit === nothing
        
        # Test LightGBM method without cv_fold
        lightgbm_method = LightGBMFilter()
        result_lightgbm = train_and_evaluate(lightgbm_method, df_no_cv, labels_no_cv, params)
        @test result_lightgbm === nothing
    end
    
    @testset "Integration Test Components" begin
        # Test individual components work together
        n_samples = 20
        merged_df = DataFrame(
            prob = Float32[0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05,
                          0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.08],
            irt_error = randn(Float32, n_samples) .* 0.5,
            MBR_max_pair_prob = rand(Float32, n_samples) .* 0.7 .+ 0.2
        )
        
        # Create candidate mask (first 15 are candidates)
        candidate_mask = [trues(15); falses(5)]
        
        # Create bad transfer labels (first 3 candidates are bad transfers)
        is_bad_transfer = [trues(3); falses(17)]
        
        # Test feature selection
        features = select_mbr_features(merged_df)
        @test :prob in features
        @test length(features) >= 1
        
        # Test feature preparation
        X, feature_names = prepare_mbr_features(merged_df[:, features])
        @test size(X, 1) == n_samples
        @test size(X, 2) == length(features) + 1  # +1 for intercept
        @test feature_names[1] == :intercept
        
        # Test threshold method on candidates only
        candidate_data = merged_df[candidate_mask, :]
        candidate_labels = is_bad_transfer[candidate_mask]
        params = (max_MBR_false_transfer_rate = 0.2,)
        
        threshold_method = ThresholdFilter()
        result = train_and_evaluate(threshold_method, candidate_data, candidate_labels, params)
        
        @test result !== nothing
        @test result.method_name == "Threshold"
        @test length(result.scores) == 15  # Number of candidates
        
        # Test filter application
        filtered_probs = apply_filtering(result, merged_df, candidate_mask, params)
        @test length(filtered_probs) == n_samples
        @test all(filtered_probs .>= 0.0f0)
        @test all(filtered_probs .<= 1.0f0)
    end
end