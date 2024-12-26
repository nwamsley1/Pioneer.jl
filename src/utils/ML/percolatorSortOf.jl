function getBestScorePerPrec(psms::SubDataFrame)
    # Create dictionary to store scores per precursor
    prec_to_best_score = Dictionary{
        @NamedTuple{prec_idx::UInt32, isotopes::Tuple{Int8,Int8}}, 
        @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}
    }()
    
    # Process each row in the dataframe
    for (i, row) in enumerate(eachrow(psms))
        key = (prec_idx = row.precursor_idx, isotopes = row.isotopes_captured)
        prob = row.prob
        
        if haskey(prec_to_best_score, key)
            max_prob, mean_prob, min_prob, n = prec_to_best_score[key]
            prec_to_best_score[key] = (
                max_prob = max(max_prob, prob),
                mean_prob = mean_prob + prob,
                min_prob = min(min_prob, prob),
                n = n + one(UInt16)
            )
        else
            insert!(prec_to_best_score, key, (
                max_prob = prob,
                mean_prob = prob,
                min_prob = prob,
                n = one(UInt16)
            ))
        end
    end
    
    # Calculate final mean probabilities
    for (key, value) in pairs(prec_to_best_score)
        max_prob, mean_prob, min_prob, n = value
        prec_to_best_score[key] = (
            max_prob = max_prob,
            mean_prob = mean_prob / n,
            min_prob = min_prob,
            n = n
        )
    end
    
    return prec_to_best_score
end

function summarize_precursors!(
    psms::SubDataFrame)
    function set_column!(vals::AbstractVector{Float32}, value::Float32) 
        for i in range(1, length(vals))
            vals[i] = value
        end
    end
    function summarize_prob(probs::AbstractVector{Float32})
        minimum(probs), maximum(probs), mean(probs)
    end
    for (key, sub_psms) in pairs(groupby(psms, [:precursor_idx, :isotopes_captured]))
        min_prob, max_prob, mean_prob = summarize_prob(sub_psms[!,:prob])
        set_column!(sub_psms[!,:min_prob], min_prob)
        set_column!(sub_psms[!,:max_prob], max_prob)
        set_column!(sub_psms[!,:mean_prob], mean_prob)
    end
end

function rankPSMs!(psms::DataFrame, 
                  file_paths::Vector{String},
                  features::Vector{Symbol};
                  colsample_bytree::Float64 = 0.5,
                  colsample_bynode::Float64 = 0.5,
                  eta::Float64 = 0.15,
                  min_child_weight::Int = 1,
                  subsample::Float64 = 0.5,
                  gamma::Int = 0,
                  max_depth::Int = 10,
                  iter_scheme::Vector{Int} = [100, 100, 200],
                  print_importance::Bool = true)
    
    # Helper functions for fold selection
    selectTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool = cv_fold == test_fold
    excludeTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool = cv_fold != test_fold
    
    # Initialize probability columns
    psms[!, :prob] = zeros(Float32, size(psms, 1))
    psms[!, :max_prob] = zeros(Float32, size(psms, 1))
    psms[!, :mean_prob] = zeros(Float32, size(psms, 1))
    psms[!, :min_prob] = zeros(Float32, size(psms, 1))
    #Faster if sorted first 
    @info "sort time..."
    @time sort!(psms, [:precursor_idx, :isotopes_captured])
    #Final prob estimates
    prob_estimates = zeros(Float32, size(psms, 1))

    #Subsample psms for training. 
    psms_training_subset = @view(psms[
                                sample(MersenneTwister(1776), #seed for reproducibility
                                1:size(psms, 1), 
                                min(10000000, size(psms, 1)), #Number of psms to use for training 
                                replace=false),:])

    unique_cv_folds = unique(psms[!, :cv_fold])
    models = Dict{UInt8, Vector{Booster}}()
    
    # Train models for each fold
    Random.seed!(1776)
    for test_fold_idx in unique_cv_folds
        #Clear prob stats 
        psms[!, :prob] .= zero(Float32)
        psms[!, :max_prob] .= zero(Float32)
        psms[!, :mean_prob] .= zero(Float32)
        psms[!, :min_prob] .= zero(Float32)
        # Get training data
        psms_train = @view(psms_training_subset[findall(x -> x != test_fold_idx, psms_training_subset[!, :cv_fold]), :])
        # Train models for each iteration
        fold_models = Vector{Booster}()
        #Each round updates, max_prob, mean_prob, and min_prob, and uses these as training features in the next round 
        for num_round in iter_scheme
            bst = xgboost(
                (psms_train[!, features], psms_train[!, :target]),
                num_round = num_round,
                colsample_bytree = colsample_bytree,
                colsample_bynode = colsample_bynode,
                gamma = gamma,
                max_depth = max_depth,
                eta = eta,
                min_child_weight = min_child_weight,
                subsample = subsample,
                objective = "binary:logistic",
                seed = rand(UInt32),
                watchlist = (;)
            )
            # Store feature names and print importance if requested
            bst.feature_names = string.(features)
            if print_importance
                println(collect(zip(importance(bst))))
            end
            
            push!(fold_models, bst)
            test_fold_idxs = findall(x -> x == test_fold_idx, psms[!, :cv_fold])
            test_fold_psms = @view(psms[test_fold_idxs,:])

            test_fold_psms[!,:prob] = XGBoost.predict(bst, test_fold_psms[!,features])
            summarize_precursors!(test_fold_psms)
            psms_train[!,:prob] =  XGBoost.predict(bst, psms_train[!, features])
            summarize_precursors!(psms_train)
        end
        # Make predictions on hold out data.
        test_fold_idxs = findall(x -> x == test_fold_idx, psms[!, :cv_fold])
        prob_estimates[test_fold_idxs] = psms[test_fold_idxs,:prob]
        # Store models for this fold
        models[test_fold_idx] = fold_models
    end
    psms[!,:prob] = prob_estimates
    for (ms_file_idx, gpsms) in pairs(groupby(psms, :ms_file_idx))
        Arrow.write(file_paths[ms_file_idx[:ms_file_idx]], gpsms)
    end
    return models
end