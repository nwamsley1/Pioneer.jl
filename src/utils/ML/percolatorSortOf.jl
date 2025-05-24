function sort_of_percolator_in_memory!(psms::DataFrame, 
                  file_paths::Vector{String},
                  features::Vector{Symbol},
                  match_between_runs::Bool = true;
                  max_q_value_xgboost_rescore::Float32 = 0.01f0,
                  max_q_value_xgboost_mbr_rescore::Float32 = 0.20f0,
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

    #Faster if sorted first 
    @info "sort time..."
    sort!(psms, [:pair_id, :isotopes_captured])

    #Final prob estimates
    prob_estimates = zeros(Float32, size(psms, 1))

    #println(sum(psms.decoy), " ", sum(psms.target), " ", size(psms, 1), "\n")

    unique_cv_folds = unique(psms[!, :cv_fold])
    models = Dict{UInt8, Vector{Booster}}()
    
    # Train models for each fold
    Random.seed!(1776)
    pbar = ProgressBar(total=length(unique_cv_folds)*length(iter_scheme))
    for test_fold_idx in unique_cv_folds
        #Clear prob stats 
        clear_prob_and_group_features!(psms, match_between_runs)
        # Get training data
        psms_train = @view(psms[findall(x -> x != test_fold_idx, psms[!, :cv_fold]), :])
        # Train models for each iteration
        fold_models = Vector{Booster}()
        #Each round updates, max_prob, mean_prob, and min_prob, and uses these as training features in the next round 
        for (itr, num_round) in enumerate(iter_scheme)
            
            

            psms_train_itr = get_training_data_for_iteration!(psms_train, 
                                                                itr,
                                                                match_between_runs, 
                                                                max_q_value_xgboost_rescore,
                                                                max_q_value_xgboost_mbr_rescore)
                                         
            bst = xgboost(
                (psms_train_itr[!, features], psms_train_itr[!, :target]),
                num_round = num_round,
                colsample_bytree = colsample_bytree,
                colsample_bynode = colsample_bynode,
                scale_pos_weight = sum(psms_train_itr.decoy) / sum(psms_train_itr.target),
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
            #print_importance = true
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

            update(pbar)
        end
        # Make predictions on hold out data.
        test_fold_idxs = findall(x -> x == test_fold_idx, psms[!, :cv_fold])
        prob_estimates[test_fold_idxs] = psms[test_fold_idxs,:prob]
        # Store models for this fold
        models[test_fold_idx] = fold_models
    end
    psms[!,:prob] = prob_estimates
    dropVectorColumns!(psms) # avoids writing issues

    transform!(
        groupby(psms, :precursor_idx),
        :prob => (p -> maximum(p)) => :global_prob
    )

    transform!(
        groupby(psms, [:precursor_idx, :ms_file_idx]),
        :prob => (p -> 1-exp(sum(log1p.(-p)))) => :prec_prob
    )
    
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
                    max_q_value_xgboost_rescore::Float32 = 0.01f0,
                    max_q_value_xgboost_mbr_rescore::Float32 = 0.20f0,
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
        prec_to_best_score_new::Dictionary,
        test_fold_filter::Function,
        test_fold_idx::UInt8,
        file_paths::Vector{String},
        bst::Booster,
        features::Vector{Symbol},
        match_between_runs::Bool;
        dropVectorCols::Bool = false)

        prec_to_global_score = Dictionary{UInt32, Float32}()
    
        #Reset counts for new scores
        for (key, value) in pairs(prec_to_best_score_new)
            max_prob, mean_prob, min_prob, n, best_log2_weights, best_irts, is_best_decoy = value
            prec_to_best_score_new[key] = (max_prob = zero(Float32), 
                                            mean_prob = zero(Float32), 
                                            min_prob = one(Float32), 
                                            n = zero(UInt16),
                                            best_log2_weights = Vector{Float32}(),
                                            best_irts = Vector{Float32}(),
                                            is_best_decoy = false)
        end
            
        for file_path in file_paths
            psms_subset = DataFrame(Arrow.Table(file_path))
            probs = XGBoost.predict(bst, psms_subset[!,features])
            
            if match_between_runs
                #Update maximum probabilities for tracked precursors 
                for (i, pair_id) in enumerate(psms_subset[!,:pair_id])
                    if !test_fold_filter(psms_subset[i,:cv_fold], test_fold_idx)::Bool continue end
                    prob = probs[i]
                    key = (pair_id = pair_id, isotopes = psms_subset[i,:isotopes_captured])
                    if haskey(prec_to_best_score_new, key)
                        max_prob, mean_prob, min_prob, n, best_log2_weights, best_irts, is_best_decoy = prec_to_best_score_new[key]
                        if max_prob < prob
                            max_prob = prob
                            best_log2_weights = log2.(psms_subset.weights[i])
                            best_irts = psms_subset.irts[i]
                            is_best_decoy = psms_subset.decoy[i]
                        end
                        if min_prob > prob
                            min_prob = prob
                        end
                        mean_prob += prob
                        n += one(UInt16)
                        prec_to_best_score_new[key] = (max_prob = max_prob, 
                                                        mean_prob = mean_prob, 
                                                        min_prob = min_prob, 
                                                        n = n,
                                                        best_log2_weights = best_log2_weights,
                                                        best_irts = best_irts,
                                                        is_best_decoy = is_best_decoy)
                    else
                        insert!(prec_to_best_score_new, key, (max_prob = prob, 
                                                                mean_prob = prob, 
                                                                min_prob = prob, 
                                                                n = one(UInt16),
                                                                best_log2_weights = log2.(psms_subset.weights[i]),
                                                                best_irts = psms_subset.irts[i],
                                                                is_best_decoy = psms_subset.decoy[i]))
                    end
                end
            end

            # update global prob
            for (i, precursor_idx) in enumerate(psms_subset[!,:precursor_idx])
                if !test_fold_filter(psms_subset[i,:cv_fold], test_fold_idx)::Bool continue end
                prob = probs[i]
                if haskey(prec_to_global_score, precursor_idx)
                    prec_to_global_score[precursor_idx] = max(prob, prec_to_global_score[precursor_idx])
                else
                    insert!(prec_to_global_score, precursor_idx, prob)
                end
            end
        end
    
        for file_path in file_paths
            psms_subset = DataFrame(Tables.columntable(Arrow.Table(file_path)))
            probs = XGBoost.predict(bst, psms_subset[!,features])
            
            for (i, pair_id) in enumerate(psms_subset[!,:pair_id])
                if !test_fold_filter(psms_subset[i,:cv_fold], test_fold_idx)::Bool continue end
                psms_subset[i,:prob] = probs[i]
                if match_between_runs
                    key = (pair_id = pair_id, isotopes = psms_subset[i,:isotopes_captured])
                    if haskey(prec_to_best_score_new, key)
                        max_prob, mean_prob, min_prob, n, best_log2_weights, best_irts, is_best_decoy = prec_to_best_score_new[key]
                        psms_subset[i,:max_prob] = max_prob
                        psms_subset[i,:min_prob] = min_prob
                        if n > 0
                            psms_subset[i,:mean_prob] = mean_prob/n
                        else
                            psms_subset[i,:mean_prob] = zero(Float32)
                        end
                        
                        best_log2_weights_padded, weights_padded = pad_equal_length(best_log2_weights, log2.(psms_subset.weights[i]))
                        best_iRTs_padded, iRTs_padded = pad_rt_equal_length(best_irts, psms_subset.irts[i])

                        best_irt_at_apex = best_irts[argmax(best_log2_weights)]
                        psms_subset.best_irt_diff[i] = abs(best_irt_at_apex - psms_subset.irts[i][argmax(psms_subset.weights[i])])
                        psms_subset.rv_coefficient[i] = rv_coefficient(best_log2_weights_padded, best_iRTs_padded, weights_padded, iRTs_padded)
                        psms_subset.is_best_decoy[i] = is_best_decoy
                        psms_subset.num_runs[i] = n
                    end
                end
            end

            # update global prob
            if !("global_prob" in names(psms_subset))
                psms_subset.global_prob = zeros(Float32, size(psms_subset, 1))
            end

            for (i, precursor_idx) in enumerate(psms_subset[!,:precursor_idx])
                if !test_fold_filter(psms_subset[i,:cv_fold], test_fold_idx)::Bool continue end
                if haskey(prec_to_global_score, precursor_idx)
                    psms_subset[i,:global_prob] = prec_to_global_score[precursor_idx]
                end
            end

            if dropVectorCols # last iteration
                Arrow.write(file_path, dropVectorColumns!(psms_subset))
            else
                Arrow.write(file_path, convert_subarrays(psms_subset))
            end
        end
        return prec_to_best_score_new
    end
                    
    function selectTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool
        return cv_fold == test_fold
    end

    function excludeTestFold(cv_fold::UInt8, test_fold::UInt8)::Bool
        return cv_fold != test_fold
    end

    unique_cv_folds = unique(psms[!, :cv_fold])
    #Train the model for 1:K-1 cross validation folds and apply to the held-out fold
    models = Dictionary{UInt8, Vector{Booster}}()
    pbar = ProgressBar(total=length(unique_cv_folds)*length(iter_scheme))
    Random.seed!(1776);
    for test_fold_idx in unique_cv_folds#(0, 1)#range(1, n_folds)
        #Clear prob stats 
        clear_prob_and_group_features!(psms, match_between_runs)
        # Get training data
        psms_train = @view(psms[findall(x -> x != test_fold_idx, psms[!, :cv_fold]), :])

        for (itr, num_round) in enumerate(iter_scheme)

            psms_train_itr = get_training_data_for_iteration!(psms_train, 
                                                                itr,
                                                                match_between_runs, 
                                                                max_q_value_xgboost_rescore,
                                                                max_q_value_xgboost_mbr_rescore)
            ###################
            #Train a model on the n-1 training folds.
            _seed_ = rand(UInt32)
            bst = xgboost((psms_train_itr[!,features], psms_train_itr[!,:target]), 
                            num_round=num_round, 
                            #monotone_constraints = monotone_constraints,
                            colsample_bytree = colsample_bytree, 
                            colsample_bynode = colsample_bynode,
                            scale_pos_weight = sum(psms_train_itr.decoy) / sum(psms_train_itr.target),
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
            #print_importance = true
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
            
            # Get probabilities for training sample so we can get q-values
            psms_train[!,:prob] = XGBoost.predict(bst, psms_train[!, features])
            
            if match_between_runs
                summarize_precursors!(psms_train)
            end

            update(pbar)
        end
    end
    pbar = ProgressBar(total=length(unique_cv_folds)*length(iter_scheme))
    for test_fold_idx in unique_cv_folds
        prec_to_best_score = Dictionary{@NamedTuple{pair_id::UInt32,
                                                    isotopes::Tuple{Int8,Int8}}, 
                                        @NamedTuple{max_prob::Float32,
                                                    mean_prob::Float32, 
                                                    min_prob::Float32, 
                                                    n::UInt16,
                                                    best_log2_weights::Vector{Float32},
                                                    best_irts::Vector{Float32},
                                                    is_best_decoy::Bool}}()

        for (train_iter, num_round) in enumerate(iter_scheme)
            model = models[test_fold_idx][train_iter]
            prec_to_best_score = getBestScorePerPrec!(
                prec_to_best_score,
                selectTestFold,
                test_fold_idx,
                file_paths,
                model,
                features,
                match_between_runs;
                dropVectorCols = (train_iter == length(iter_scheme)) && (test_fold_idx == unique_cv_folds[end]) # remove on last iteration
            )
            update(pbar)
        end
    end
    
    return models#bst, vcat(psms_test_folds...)
end














function summarize_precursors!(psms::AbstractDataFrame)
    # Collect the grouped pairs so we can iterate in a threaded for loop.
    grouped_data = collect(pairs(groupby(psms, [:pair_id, :isotopes_captured])))

    Threads.@threads for idx in eachindex(grouped_data)
        key, sub_psms = grouped_data[idx]
        
        min_prob, max_prob, mean_prob = summarize_prob(sub_psms[!,:prob])
        best_idx = argmax(sub_psms.prob)
        best_log2_weights = log2.(sub_psms.weights[best_idx])
        best_iRTs = sub_psms.irts[best_idx]

        for i in 1:nrow(sub_psms)
            best_log2_weights_padded, weights_padded = pad_equal_length(best_log2_weights, log2.(sub_psms.weights[i]))
            best_iRTs_padded, iRTs_padded = pad_rt_equal_length(best_iRTs, sub_psms.irts[i])
                    
            best_irt_at_apex = sub_psms.irts[best_idx][argmax(best_log2_weights)]
            sub_psms.best_irt_diff[i] = abs(best_irt_at_apex - sub_psms.irts[i][argmax(sub_psms.weights[i])])
            sub_psms.rv_coefficient[i] = rv_coefficient(best_log2_weights_padded, best_iRTs_padded, weights_padded, iRTs_padded)
            sub_psms.min_prob[i] = min_prob
            sub_psms.max_prob[i] = max_prob
            sub_psms.mean_prob[i] = mean_prob
            sub_psms.is_best_decoy[i] = sub_psms.decoy[best_idx]
            sub_psms.num_runs[i] = nrow(sub_psms)
        end
    end

end

function clear_prob_and_group_features!(
    psms::AbstractDataFrame,
    match_between_runs::Bool
)
    # Clear main :prob column
    psms[!, :prob] .= zero(Float32)

    if match_between_runs
        psms[!, :max_prob]       .= zero(Float32)
        psms[!, :mean_prob]      .= zero(Float32)
        psms[!, :min_prob]       .= zero(Float32)
        psms[!, :rv_coefficient] .= zero(Float32)
        psms[!, :best_irt_diff]  .= zero(Float32)
        psms[!, :num_runs]       .= zero(Int32)
        psms[!, :is_best_decoy]  .= zero(Bool)
        psms[!, :q_value]        .= zero(Float64)
    end

    return psms
end

function initialize_group_features!(
    psms::AbstractDataFrame,
    match_between_runs::Bool
)
    n = nrow(psms)
    psms[!, :prob] = zeros(Float32, n)

    if match_between_runs
        psms[!, :max_prob]       = zeros(Float32, n)
        psms[!, :mean_prob]      = zeros(Float32, n)
        psms[!, :min_prob]       = zeros(Float32, n)
        psms[!, :rv_coefficient] = zeros(Float32, n)
        psms[!, :best_irt_diff]  = zeros(Float32, n)
        psms[!, :num_runs]       = zeros(Int32, n)
        psms[!, :is_best_decoy]  = zeros(Bool, n)
        psms[!, :q_value]        = zeros(Float64, n)
    end

    return psms
end

function get_training_data_for_iteration!(
    psms_train::AbstractDataFrame,
    itr::Int,
    match_between_runs::Bool,
    max_q_value_xgboost_rescore::Float32,
    max_q_value_xgboost_mbr_rescore::Float32
)
    # Train on all precursors during first iteration
    psms_train_itr = psms_train

    if itr > 1
        # For subsequent iterations, train on top scoring precursors
        get_qvalues!(
            psms_train_itr[!,:prob],
            psms_train_itr[!,:target],
            psms_train_itr[!,:q_value]
        )

        # Also train on top scoring MBR candidates if requested
        if match_between_runs
            # Determine prob threshold for precursors passing the q-value threshold
            max_prob_threshold = minimum(
                psms_train_itr.prob[
                    psms_train_itr.target .& (psms_train_itr.q_value .<= max_q_value_xgboost_rescore)
                ]
            )

            # Hacky way to ensure anything passing the initial q-value threshold
            # will pass the next q-value threshold
            psms_train.q_value[psms_train_itr.q_value .<= max_q_value_xgboost_rescore] .= 0.0
            psms_train.q_value[psms_train_itr.q_value .> max_q_value_xgboost_rescore]  .= 1.0

            # Must have at least one precursor passing the q-value threshold,
            # and the best precursor can't be a decoy
            psms_train_mbr = subset(
                psms_train_itr,
                [:is_best_decoy, :max_prob, :prob] => ByRow((d, mp, p) ->
                    (!d) && (mp >= max_prob_threshold) && (p < max_prob_threshold)
                );
                view = true
            )

            # Compute MBR q-values.
            get_qvalues!(
                psms_train_mbr[!,:prob],
                psms_train_mbr[!,:target],
                psms_train_mbr[!,:q_value]
            )
            # Take all decoys and targets passing q_thresh (all 0's now) or mbr_q_thresh
            psms_train_itr = subset(
                psms_train,
                [:target, :q_value] => ByRow((t,q) -> (!t) || (t && q <= max_q_value_xgboost_mbr_rescore))
            )
        else
            # Take all decoys and targets passing q_thresh
            psms_train_itr = subset(
                psms_train,
                [:target, :q_value] => ByRow((t,q) -> (!t) || (t && q <= max_q_value_xgboost_rescore))
            )
        end
    end

    return psms_train_itr
end



function dropVectorColumns!(df)
    to_drop = String[]
    for col in names(df)
        if eltype(df[!, col]) <: AbstractVector
            push!(to_drop, col)
        end
    end
    # 2) Drop those columns in place
    select!(df, Not(to_drop))
end


function rv_coefficient(weights_A::AbstractVector{<:Real},
    times_A::AbstractVector{<:Real},
    weights_B::AbstractVector{<:Real},
    times_B::AbstractVector{<:Real})

    # Construct two Nx2 matrices, each row is (weight, time)
    X = hcat(collect(weights_A), collect(times_A))
    Y = hcat(collect(weights_B), collect(times_B))

    # Compute cross-products (Gram matrices)
    Sx = X' * X
    Sy = Y' * Y

    # Numerator: trace(Sx * Sy)
    numerator = tr(Sx * Sy)

    # Denominator: sqrt( trace(Sx*Sx)* trace(Sy*Sy) )
    denominator = sqrt(tr(Sx * Sx) * tr(Sy * Sy))

    # Protect against zero in denominator (e.g. if X or Y is all zeros)
    if denominator == 0
        return 0.0
    end

    return numerator / denominator
end

function pad_equal_length(x::AbstractVector, y::AbstractVector)

    function pad(data, d)
        left  = div(d, 2)
        right = d - left        # put extra on the right if d is odd
        padded_data = vcat(zeros(eltype(data), left), data, zeros(eltype(data), right))
        return padded_data
    end

    d = length(y) - length(x)
    if d == 0
        # No padding needed
        return (x, y)
    elseif d > 0
        # Pad x
        return (pad(x, d), y)
    else
        # Pad y
        return (x, pad(y, abs(d)))
    end
end


function pad_rt_equal_length(x::AbstractVector, y::AbstractVector)

    function pad(data, n, d, rt_step)
        left  = div(d, 2) 
        right = d - left        # put extra on the right if d is odd
        padded_data = vcat(zeros(eltype(data), left), data, zeros(eltype(data), right))
        @inbounds @fastmath for i in range(1, left)
            padded_data[i] = padded_data[left+1] - (rt_step * ((left+1) - i))
        end 
        @inbounds @fastmath for i in range(1, right)
            padded_data[i+left+n] = padded_data[left+n] + (rt_step * i)
        end 
        return padded_data
    end

    nx = length(x)
    ny = length(y)
    d = ny - nx
    rt_step_x = (x[end] - x[1]) / nx
    rt_step_y = (y[end] - y[1]) / ny
    rt_step = nx > 1 ? rt_step_x : rt_step_y

    if d == 0
        # No padding needed
        return (x, y)
    elseif d > 0
        # Pad x
        return (pad(x, nx, d, rt_step), y)
    else
        # Pad y
        return (x, pad(y, ny, abs(d), rt_step))
    end
end

function summarize_prob(probs::AbstractVector{Float32})
    minimum(probs), maximum(probs), mean(probs)
end

function convert_subarrays(df::DataFrame)
    for col_name in names(df)
        col_data = df[!, col_name]
        # If it's a SubArray, let's convert it to a plain vector:
        if eltype(col_data) <: AbstractVector
           df[!, col_name] = [Vector(col_data[i]) for i in eachindex(col_data)]
        end
    end
    return df
end