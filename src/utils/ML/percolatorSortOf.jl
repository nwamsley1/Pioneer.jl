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

function sort_of_percolator_in_memory!(psms::DataFrame, 
                  file_paths::Vector{String},
                  features::Vector{Symbol},
                  match_between_runs::Bool = true;
                  colsample_bytree::Float64 = 0.5,
                  colsample_bynode::Float64 = 0.5,
                  eta::Float64 = 0.15,
                  min_child_weight::Int = 1,
                  subsample::Float64 = 0.5,
                  gamma::Int = 0,
                  max_depth::Int = 10,
                  iter_scheme::Vector{Int} = [100, 100, 200],
                  print_importance::Bool = true)
    
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
    # Helper functions for fold selection
    selectTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool = cv_fold == test_fold
    excludeTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool = cv_fold != test_fold
    
    # Initialize probability columns
    psms[!, :prob] = zeros(Float32, size(psms, 1))
    if match_between_runs
        psms[!, :max_prob] = zeros(Float32, size(psms, 1))
        psms[!, :mean_prob] = zeros(Float32, size(psms, 1))
        psms[!, :min_prob] = zeros(Float32, size(psms, 1))
    end
    #Faster if sorted first 
    @info "sort time..."
    @time sort!(psms, [:precursor_idx, :isotopes_captured])
    #Final prob estimates
    prob_estimates = zeros(Float32, size(psms, 1))

    unique_cv_folds = unique(psms[!, :cv_fold])
    models = Dict{UInt8, Vector{Booster}}()
    
    # Train models for each fold
    Random.seed!(1776)
    for test_fold_idx in unique_cv_folds
        #Clear prob stats 
        psms[!, :prob] .= zero(Float32)
        if match_between_runs
            psms[!, :max_prob] .= zero(Float32)
            psms[!, :mean_prob] .= zero(Float32)
            psms[!, :min_prob] .= zero(Float32)
        end
        # Get training data
        psms_train = @view(psms[findall(x -> x != test_fold_idx, psms[!, :cv_fold]), :])
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
            psms_train[!,:prob] =  XGBoost.predict(bst, psms_train[!, features])

            if match_between_runs
                summarize_precursors!(test_fold_psms)
                summarize_precursors!(psms_train)
            end
        end
        # Make predictions on hold out data.
        test_fold_idxs = findall(x -> x == test_fold_idx, psms[!, :cv_fold])
        prob_estimates[test_fold_idxs] = psms[test_fold_idxs,:prob]
        # Store models for this fold
        models[test_fold_idx] = fold_models
    end
    psms[!,:prob] = prob_estimates
    for (ms_file_idx, gpsms) in pairs(groupby(psms, :ms_file_idx))
        fpath = file_paths[ms_file_idx[:ms_file_idx]]
        writeArrow(fpath, gpsms)
    end
    return models
end

function sort_of_percolator_out_of_memory!(psms::DataFrame, 
                    file_paths::Vector{String},
                    features::Vector{Symbol},
                    match_between_runs::Bool = true; 
                    colsample_bytree::Float64 = 0.5, 
                    colsample_bynode::Float64 = 0.5,
                    eta::Float64 = 0.15, 
                    min_child_weight::Int = 1, 
                    subsample::Float64 = 0.5, 
                    gamma::Int = 0, 
                    max_depth::Int = 10,
                    iter_scheme::Vector{Int} = [100, 100, 200],
                    print_importance::Bool = true)

    function getBestScorePerPrec!(
        prec_to_best_score_new::Dictionary{@NamedTuple{prec_idx::UInt32,isotopes::Tuple{Int8,Int8}}, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}},
        test_fold_filter::Function,
        test_fold_idx::UInt8,
        file_paths::Vector{String},
        bst::Booster,
        features::Vector{Symbol})
    
        #Reset counts for new scores
        for (key, value) in pairs(prec_to_best_score_new)
            max_prob, mean_prob, min_prob, n = value
            prec_to_best_score_new[key] = (max_prob = zero(Float32), mean_prob = zero(Float32), min_prob = one(Float32), n = zero(UInt16))
        end
            
        for file_path in file_paths
            psms_subset = DataFrame(Arrow.Table(file_path))
            probs = XGBoost.predict(bst, psms_subset[!,features])
            #Update maximum probabilities for tracked precursors 
            for (i, prec_idx) in enumerate(psms_subset[!,:precursor_idx])
                if !test_fold_filter(psms_subset[i,:cv_fold], test_fold_idx)::Bool continue end
                prob = probs[i]
                key = (prec_idx = prec_idx, isotopes = psms_subset[i,:isotopes_captured])
                if haskey(prec_to_best_score_new, key)
                    max_prob, mean_prob, min_prob, n = prec_to_best_score_new[key]
                    if max_prob < prob
                        max_prob = prob
                    end
                    if min_prob > prob
                        min_prob = prob
                    end
                    mean_prob += prob
                    n += one(UInt16)
                    prec_to_best_score_new[key] = (max_prob = max_prob, mean_prob = mean_prob, min_prob = min_prob, n = n)
                else
                    insert!(prec_to_best_score_new, key, 
                    (max_prob = prob, mean_prob = prob, min_prob = prob, n = one(UInt16)))
                end
            end
        end
    
        for file_path in file_paths
            psms_subset = DataFrame(Tables.columntable(Arrow.Table(file_path)))
            probs = XGBoost.predict(bst, psms_subset[!,features])
            for (i, prec_idx) in enumerate(psms_subset[!,:precursor_idx])
                if !test_fold_filter(psms_subset[i,:cv_fold], test_fold_idx)::Bool continue end
                psms_subset[i,:prob] = probs[i]
                key = (prec_idx = prec_idx, isotopes = psms_subset[i,:isotopes_captured])
                if haskey(prec_to_best_score_new, key)
                    max_prob, mean_prob, min_prob, n = prec_to_best_score_new[key]
                    psms_subset[i,:max_prob] = max_prob
                    psms_subset[i,:min_prob] = min_prob
                    if n > 0
                        psms_subset[i,:mean_prob] = mean_prob/n
                    else
                        psms_subset[i,:mean_prob] = zero(Float32)
                    end
                end
            end
            Arrow.write(file_path, psms_subset)
        end
        return prec_to_best_score_new
    end
    
    function getBestScorePerPrec!(
        prec_to_best_score_new::Dictionary{@NamedTuple{prec_idx::UInt32,isotopes::Tuple{Int8,Int8}}, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}},
        psms::DataFrame)
        m = 0
        for (i, prec_idx) in enumerate(psms[!,:precursor_idx])
            key = (prec_idx = prec_idx, isotopes = psms[i,:isotopes_captured])
            if haskey(prec_to_best_score_new, key)
                max_prob, mean_prob, min_prob, n = prec_to_best_score_new[key]
                psms[i,:max_prob] = max_prob
                psms[i,:min_prob] = min_prob
                psms[i,:mean_prob] = mean_prob/n
                m += 1
            end
        end
    end
                    
    function selectTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool
        return cv_fold == test_fold
    end

    function excludeTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool
        return cv_fold != test_fold
    end

    if match_between_runs
        psms[!,:max_prob], psms[!,:mean_prob], psms[!,:min_prob] = zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1)), zeros(Float32, size(psms, 1))
    end
    
    psms[!,:prob] = zeros(Float32, size(psms, 1))
    folds = psms[!,:cv_fold]
    unique_cv_folds = unique(folds)
    #Train the model for 1:K-1 cross validation folds and apply to the held-out fold
    models = Dictionary{UInt8, Vector{Booster}}()
    pbar = ProgressBar(total=length(unique_cv_folds)*length(iter_scheme))
    Random.seed!(1776);
    for test_fold_idx in unique_cv_folds#(0, 1)#range(1, n_folds)
        train_fold_idxs_all = findall(x->x!=test_fold_idx, folds)
        psms_train = psms[train_fold_idxs_all,:]
        prec_to_best_score = Dictionary{@NamedTuple{prec_idx::UInt32,isotopes::Tuple{Int8,Int8}}, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}}()
        for (train_iter, num_round) in enumerate(iter_scheme)
            println("training: ", test_fold_idx, " ", train_iter, " ", size(psms,1), "\n")
            ###################
            #Train a model on the n-1 training folds.
            _seed_ = rand(UInt32)
            bst = xgboost((psms_train[!,features], psms_train[!,:target]), 
                            num_round=num_round, 
                            #monotone_constraints = monotone_constraints,
                            colsample_bytree = colsample_bytree, 
                            colsample_bynode = colsample_bynode,
                            gamma = gamma, 
                            max_depth=max_depth, 
                            eta = eta, 
                            min_child_weight = min_child_weight, 
                            subsample = subsample, 
                            objective="binary:logistic",
                            seed = _seed_,
                            #max_bin = 128,
                            watchlist=(;)
                            )
            if !haskey(models, test_fold_idx)
                insert!(
                    models,
                    test_fold_idx,
                    Vector{Booster}([bst])
                )
            else
                push!(models[test_fold_idx], bst)
            end
            bst.feature_names = [string(x) for x in features]
            print_importance ? println(collect(zip(importance(bst)))[1:30]) : nothing
            #println(collect(zip(importance(bst)))[1:5])
            prec_to_best_score = getBestScorePerPrec!(
                prec_to_best_score,
                excludeTestFold,
                test_fold_idx,
                file_paths,
                bst,
                features
            )
            getBestScorePerPrec!(
                prec_to_best_score,
                psms_train
            )
            update(pbar)
        end
    end
    pbar = ProgressBar(total=length(unique_cv_folds)*length(iter_scheme))
    for test_fold_idx in unique_cv_folds
        prec_to_best_score = Dictionary{@NamedTuple{prec_idx::UInt32,isotopes::Tuple{Int8,Int8}}, @NamedTuple{max_prob::Float32, mean_prob::Float32, min_prob::Float32, n::UInt16}}()
        for (train_iter, num_round) in enumerate(iter_scheme)
            model = models[test_fold_idx][train_iter]
            prec_to_best_score = getBestScorePerPrec!(
                prec_to_best_score,
                selectTestFold,
                test_fold_idx,
                file_paths,
                model,
                features
            )
            update(pbar)
        end
    end
    
    return models#bst, vcat(psms_test_folds...)
end
