using DecisionTree
function rankPSMs!(PSMs::DataFrame, features::Vector{Symbol}; n_folds::Int = 3, n_trees::Int = 500, n_features::Int = 10, max_depth::Int = 10, min_sample_size::Int = 5000)
    
    X = Matrix(PSMs[:,features])
    X_labels = PSMs[:, :decoy]

    function getCVFolds(PSMs::DataFrame, n_folds::Int)
        permutation = randperm(size(PSMs)[1])
        fold_size = length(permutation)÷n_folds
        return [((n-1)*fold_size + 1):(n*fold_size) for n in range(1, n_folds)]
    end

    PSMs[:,:prob] = zeros(Float64, size(PSMs)[1])
    model = ""
    for test_fold_idx in range(1, n_folds)
        println(test_fold_idx)
        train_fold_idxs = vcat([folds[fold] for fold in range(1, length(folds)) if fold != test_fold_idx]...)
        #println("train_fold_idxs ", train_fold_idxs)
        train_features = X[train_fold_idxs,:]
        train_classes = X_labels[train_fold_idxs,1]
        fraction = min(min_sample_sie/length(train_classes), 0.5)
        model = build_forest(train_classes, train_features, features, n_trees, fraction, max_depth)
        probs = apply_forest_proba(model, X[folds[test_fold_idx],:],[true, false])
        PSMs[folds[test_fold_idx],:prob] = probs[:,2]
        println(features[sortperm(split_importance(model))])
    end
    #model = build_forest(X_labels, X', 4, 2000, 0.5, 3)
    #probs = apply_forest_proba(model, X',[true, false])
    #PSMs[:,:prob] = probs[:,2]
    return model
end