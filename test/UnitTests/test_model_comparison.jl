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
using Random

# Include the model comparison framework
include("../../src/Routines/SearchDIA/SearchMethods/ScoringSearch/model_comparison.jl")

@testset "Model Comparison Framework Tests" begin
    
    @testset "Model Configuration" begin
        configs = create_model_configurations()
        
        @test length(configs) == 3
        @test any(c -> c.name == "SimpleXGBoost", configs)
        @test any(c -> c.name == "ProbitRegression", configs)
        @test any(c -> c.name == "SuperSimplified", configs)
        
        simple_xgb = configs[findfirst(c -> c.name == "SimpleXGBoost", configs)]
        @test simple_xgb.model_type == :xgboost
        @test length(simple_xgb.features) > 30  # Should have many features
        @test haskey(simple_xgb.hyperparams, :eta)
        
        probit = configs[findfirst(c -> c.name == "ProbitRegression", configs)]
        @test probit.model_type == :probit
        @test haskey(probit.hyperparams, :n_folds)
        
        simplified = configs[findfirst(c -> c.name == "SuperSimplified", configs)]
        @test simplified.model_type == :xgboost
        @test length(simplified.features) == 5  # Minimal feature set
    end
    
    @testset "Data Splitting" begin
        # Create a mock PSM DataFrame
        n_psms = 1000
        psms = DataFrame(
            target = vcat(fill(true, 600), fill(false, 400)),  # 60% targets, 40% decoys
            score = rand(n_psms),
            cv_fold = rand(1:3, n_psms)
        )
        
        train_indices, val_indices = create_train_validation_split(psms, 0.2, seed=42)
        
        @test length(train_indices) + length(val_indices) == n_psms
        @test length(val_indices) ≈ n_psms * 0.2 atol=50  # Allow some tolerance
        
        # Check stratification
        train_targets = sum(psms.target[train_indices])
        val_targets = sum(psms.target[val_indices])
        total_targets = sum(psms.target)
        
        # Both should have similar target proportions
        train_prop = train_targets / length(train_indices)
        val_prop = val_targets / length(val_indices)
        @test abs(train_prop - val_prop) < 0.1  # Within 10%
    end
    
    @testset "Performance Metrics" begin
        # Test AUC calculation
        predictions = Float32[0.9, 0.8, 0.7, 0.6, 0.4, 0.3, 0.2, 0.1]
        targets = Bool[true, true, true, false, false, false, false, true]
        
        auc = compute_auc(predictions, targets)
        @test 0.0 <= auc <= 1.0
        @test auc > 0.5  # Should be better than random
        
        # Test accuracy
        accuracy = compute_accuracy(predictions, targets, 0.5f0)
        @test 0.0 <= accuracy <= 1.0
        
        # Test sensitivity and specificity
        sensitivity, specificity = compute_sensitivity_specificity(predictions, targets, 0.5f0)
        @test 0.0 <= sensitivity <= 1.0
        @test 0.0 <= specificity <= 1.0
    end
    
    @testset "Model Selection Logic" begin
        # Create mock performance results
        performances = [
            ModelPerformance("ModelA", 100, 0.85, 0.8, 0.75, 0.85, 10.0, 30),
            ModelPerformance("ModelB", 120, 0.80, 0.82, 0.78, 0.86, 15.0, 25),  # Best by targets
            ModelPerformance("ModelC", 100, 0.90, 0.85, 0.80, 0.90, 20.0, 35)   # Same targets as A, higher AUC
        ]
        
        best_model = select_best_model(performances)
        @test best_model == "ModelB"  # Should select the one with most targets passing
        
        # Test tie-breaking by AUC
        performances_tie = [
            ModelPerformance("ModelA", 100, 0.85, 0.8, 0.75, 0.85, 10.0, 30),
            ModelPerformance("ModelB", 100, 0.90, 0.82, 0.78, 0.86, 15.0, 25)   # Same targets, higher AUC
        ]
        
        best_model_tie = select_best_model(performances_tie)
        @test best_model_tie == "ModelB"  # Should select higher AUC
    end
    
    @testset "Edge Cases" begin
        # Test empty performance list (should error)
        @test_throws BoundsError select_best_model(ModelPerformance[])
        
        # Test AUC with no targets or no decoys
        all_targets = fill(true, 10)
        all_decoys = fill(false, 10)
        predictions = rand(Float32, 10)
        
        auc_all_targets = compute_auc(predictions, all_targets)
        auc_all_decoys = compute_auc(predictions, all_decoys)
        
        @test auc_all_targets == 0.5  # No discrimination possible
        @test auc_all_decoys == 0.5   # No discrimination possible
        
        # Test data splitting with very small dataset
        small_psms = DataFrame(
            target = [true, false, true, false],
            cv_fold = [1, 2, 1, 2]
        )
        
        train_small, val_small = create_train_validation_split(small_psms, 0.5)
        @test length(train_small) + length(val_small) == 4
    end
end

println("Model comparison framework unit tests completed!")