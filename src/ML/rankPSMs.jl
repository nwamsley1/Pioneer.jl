function rankPSMs!(PSMs::DataFrame, features::Vector{Symbol}; n_folds::Int = 3, n_trees::Int = 500, n_features::Int = 10, max_depth::Int = 10, min_sample_size::Int = 5000)
    
    X = Matrix(PSMs[:,features])
    X_labels = PSMs[:, :decoy]

    function getCVFolds(PSMs::DataFrame, n_folds::Int)
        permutation = randperm(size(PSMs)[1])
        fold_size = length(permutation)÷n_folds
        return [((n-1)*fold_size + 1):(n*fold_size) for n in range(1, n_folds)]
    end

    function getTrainFold(folds::Vector{UnitRange{Int64}}, test_fold_idx::Int, X::Matrix{Float64}, X_labels::Vector{Bool})
        train_fold_idxs = vcat([folds[fold] for fold in range(1, length(folds)) if fold != test_fold_idx]...)
        return X[train_fold_idxs,:], X_labels[train_fold_idxs,1]
    end
    for test_fold_idx in range(1, n_folds)
        train_features, train_classes = getTrainFold(folds, test_fold_idx, X, X_labels)
        fraction = min(min_sample_size/length(train_classes), 0.5)

        model = build_forest(train_classes, train_features, features, n_trees, fraction, max_depth)
        probs = apply_forest_proba(model, X[folds[test_fold_idx],:],[true, false])

        PSMs[folds[test_fold_idx],:target_prob] = probs[:,2]
        println(features[sortperm(split_importance(model))])
    end
    return nothing
end

function getQvalues!(PSMs::DataFrame)
    #Could bootstratp to get more reliable values. 
    q_values = zeros(Float64, (size(PSMs)[1],))
    order = reverse(sortperm(PSMs[:,:target_prob]))
    targets = 0
    decoys = 0
    for i in order
        if PSMs[i,:decoy] == true
            decoys += 1
        else
            targets += 1
        end
        q_values[i] = decoys/(targets + decoys)
    end
    PSMs[:,:q_value] = q_values;
end

combine(sdf -> sdf[argmax(sdf.hyperscore),:], groupby(PSMs[PSMs[:,:q_values].<0.1,:], :precursor_idx))