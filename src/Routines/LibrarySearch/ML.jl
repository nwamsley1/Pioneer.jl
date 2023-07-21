function rankPSMs!(PSMs::DataFrame, features::Vector{Symbol}; n_folds::Int = 3, n_trees::Int = 500, n_features::Int = 10, max_depth::Int = 10, fraction::AbstractFloat = 0.1, print_importance::Bool = false)
   
    #[:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:diff_scribe,:charge,:city_block,:matched_ratio,:weight,:intensity,:count,:SN]
    X = Matrix(PSMs[:,features])
    X_labels = PSMs[:, :decoy]

    permutation = randperm(size(PSMs)[1])
    fold_size = length(permutation)÷n_folds

    folds = [((n-1)*fold_size + 1):(n*fold_size) for n in range(1, n_folds)]

    PSMs[:,:prob] = zeros(Float64, size(PSMs)[1])
    model = ""
    for test_fold_idx in range(1, n_folds)
        train_fold_idxs = vcat([folds[fold] for fold in range(1, length(folds)) if fold != test_fold_idx]...)
        train_features = X[train_fold_idxs,:]
        train_classes = X_labels[train_fold_idxs,1]
        model = build_forest(train_classes, train_features, n_features, n_trees, fraction, max_depth)
        probs = apply_forest_proba(model, X[folds[test_fold_idx],:],[true, false])
        PSMs[folds[test_fold_idx],:prob] = probs[:,2]
        if print_importance
            println(features[sortperm(split_importance(model))])
        end
    end
    return model
end

function getQvalues!(PSMs::DataFrame, probs::Vector{Union{Missing, Float64}}, labels::Vector{Union{Missing, Bool}})
    #Could bootstratp to get more reliable values. 
    q_values = zeros(Float64, (length(probs),))
    order = reverse(sortperm(probs))
    targets = 0
    decoys = 0
    for i in order
        if labels[i] == true
            decoys += 1
        else
            targets += 1
        end
        q_values[i] = decoys/(targets + decoys)
    end
    PSMs[:,:q_value] = q_values;
end

CSV.write("/Users/n.t.wamsley/Projects/TEST_DATA/PSMs_072023_05.csv", PSMs)


PSMs = DataFrame(CSV.File("/Users/n.t.wamsley/Projects/TEST_DATA/PSMs_072023_04.csv"))

PSMs[isnan.(PSMs[:,:matched_ratio]),:matched_ratio] .= Inf
PSMs[(PSMs[:,:matched_ratio]).==Inf,:matched_ratio] .= 416119.4

#Try XGBoost
function rankPSMs!(PSMs::DataFrame, features::Vector{Symbol}; n_folds::Int = 3, colsample_bytree::Float64 = 0.5, num_round::Int = 25, eta::Float64 = 0.15, min_child_weight::Int = 1, subsample::Float64 = 0.5, gamma::Int = 0, max_depth::Int = 10, print_importance::Bool = false)
   
    #[:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:diff_scribe,:charge,:city_block,:matched_ratio,:weight,:intensity,:count,:SN]
    X = Matrix(PSMs[:,features])
    X_labels = PSMs[:, :decoy]

    permutation = randperm(size(PSMs)[1])
    fold_size = length(permutation)÷n_folds

    folds = [((n-1)*fold_size + 1):(n*fold_size) for n in range(1, n_folds)]

    PSMs[:,:prob] = zeros(Float64, size(PSMs)[1])
    model = ""
    for test_fold_idx in range(1, n_folds)
        train_fold_idxs = vcat([folds[fold] for fold in range(1, length(folds)) if fold != test_fold_idx]...)
        train_features = X[train_fold_idxs,:]
        train_classes = X_labels[train_fold_idxs,1]
        #train_features = DMatrix(train_features, feature_names=[string(x) for x in features])
        #XGBoost.setfeaturenames!(train_features, [string(x) for x in features])
        bst = xgboost((train_features, train_classes), num_round=num_round, colsample_bytree = colsample_bytree, gamma = gamma, max_depth=max_depth, eta = eta, min_child_weight = min_child_weight, subsample = subsample, objective="binary:logistic")
        #model = build_forest(train_classes, train_features, n_features, n_trees, fraction, max_depth)
        ŷ = XGBoost.predict(bst, X[folds[test_fold_idx],:])
        #probs = apply_forest_proba(model, X[folds[test_fold_idx],:],[true, false])
        PSMs[folds[test_fold_idx],:prob] = (1 .- ŷ) #probs[:,2]
        #if print_importance
        #    println(features[sortperm(split_importance(model))])
        #end
    end
    bst.feature_names = [string(x) for x in features]
    return bst
end

features = [:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:charge,:city_block,:matched_ratio,:weight,:kendall,:rank_hyper,:rank_poisson, :rank_scribe,:rank_total,:len]
@time bst = rankPSMs!(PSMs, features, colsample_bytree = 1.0, min_child_weight = 10, gamma = 10, n_folds = 3, num_round = 50, eta = 0.15)
@time getQvalues!(PSMs, PSMs[:,:prob], PSMs[:,:decoy]);
length(unique(PSMs[(PSMs[:,:q_value].<=0.01).&(PSMs[:,:decoy].==false),:precursor_idx]))
length(unique(PSMs[(PSMs[:,:q_value].<=0.1).&(PSMs[:,:decoy].==false),:precursor_idx]))


features = [:hyperscore,:total_ions,:intensity_explained,:error,:poisson,:spectral_contrast_all, :spectral_contrast_matched,:RT_error,:scribe_score,:y_ladder,:b_ladder,:RT,:diff_hyper,:median_ions,:n_obs,:charge,:city_block,:matched_ratio,:weight,:kendall,:rank_hyper,:rank_poisson, :rank_scribe,:rank_total,:len,:intensity,:count,:SN]
best_psms = best_psms[(best_psms[:,:intensity].>0).&(best_psms[:,:count].>=5),:]
@time bst = rankPSMs!(best_psms, features, colsample_bytree = 1.0, min_child_weight = 10, gamma = 10, subsample = 1.0, n_folds = 5, num_round = 200, eta = 0.075)
@time getQvalues!(best_psms, best_psms[:,:prob], best_psms[:,:decoy]);
length(unique(best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:precursor_idx]))
#length(unique(best_psms[(best_psms[:,:q_value].<=0.1).&(PSMs[:,:best_psms].==false),:precursor_idx]))


best_psms = DataFrame(CSV.File("/Users/n.t.wamsley/Projects/TEST_DATA/best_psms_072123_01.csv"))
#["Targets ≤ 0.1% FDR",""
function plotStepHist(PSMs::DataFrame, group_a::BitVector, group_b::BitVector, column::Symbol, b_range::Any = nothing; normalize::Bool = true, transform::Any = x->x, title::String = "TITLE", label_a::String="Y1", label_b::String="Y2", f_out::String = "test.pdf")
    theme(:wong)
    p = plot(title=title, legend =:topleft);
    #b_range = range(0, 1, length=1000)
    stephist(p, transform.(PSMs[group_a,column]), bins = b_range, alpha = 1, normalize=normalize, labels = label_a)
    stephist!(transform.(PSMs[group_b,column]), bins = b_range, alpha = 1, normalize=normalize, labels =label_b)
    savefig(f_out)
end

targets = PSMs[:,:decoy].==false
decoys = PSMs[:,:decoy].==true
targets_01fdr = (PSMs[:,:decoy].==false) .& (PSMs[:,:q_value].<0.01)

plotStepHist(PSMs, targets, decoys, :prob, range(0, 1, length=1000), normalize = false, title = "Discriminant Score", label_a = "Targets", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/discriminant_score.pdf")
plotStepHist(PSMs, targets, decoys, :q_value, range(0, 0.5, length=100), normalize = false, title = "Q-Value", label_a = "Targets", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/q_value.pdf")
plotStepHist(PSMs, targets_01fdr, decoys, :spectral_contrast_all, range(0.5, 1, length=250), title = "Cosine Similarity Score", label_a = "Targets ≤ 0.1% FDR", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/CSS.pdf")
plotStepHist(PSMs, targets_01fdr, decoys, :matched_ratio, range(-5, 10, length=1000), transform = x->log2(x), title = "Log2 Predicted Spectra Explained Ratio", label_a = "Targets ≤ 0.1% FDR", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/matched_ratio.pdf")
plotStepHist(PSMs, targets_01fdr, decoys, :scribe_score, range(0, 20, length=1000), title = "Scribe Score", label_a = "Targets ≤ 0.1% FDR", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/scribe_score.pdf")
plotStepHist(PSMs, targets_01fdr, decoys, :kendall, range(-20, 0, length=100), title = "Kendall Correlation log2(p-val)", label_a = "Targets ≤ 0.1% FDR", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/kendall.pdf")
plotStepHist(PSMs, targets_01fdr, decoys, :hyperscore, range(1,100, length=1000), title = "XTandem HyperScore", label_a = "Targets ≤ 0.1% FDR", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/hyperscore.pdf")
plotStepHist(PSMs, targets_01fdr, decoys, :rank_total, range(1, 50, length=50), title = "Matched Ions Rank", label_a = "Targets ≤ 0.1% FDR", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/matched_ions_rank.pdf")
plotStepHist(PSMs, targets_01fdr, decoys, :rank_hyper, range(1, 50, length=50), title = "Hyperscore Rank", label_a = "Targets ≤ 0.1% FDR", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/hyperscore_rank.pdf")
plotStepHist(PSMs, targets_01fdr, decoys, :rank_scribe, range(1, 50, length=50), title = "Scribe Score Rank", label_a = "Targets ≤ 0.1% FDR", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/scribe_score_rank.pdf")
plotStepHist(PSMs, targets_01fdr, decoys, :poisson, range(-50, 0, length=100), title = "Poisson", label_a = "Targets ≤ 0.1% FDR", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/poisson.pdf")
plotStepHist(PSMs, targets_01fdr, decoys, :RT_error, range(0, 40, length=100), title = "RT Error", label_a = "Targets ≤ 0.1% FDR", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/rt_error.pdf")
plotStepHist(PSMs, targets_01fdr, decoys, :city_block,  range(-5, 0, length=1000), title = "City Block Distance", label_a = "Targets ≤ 0.1% FDR", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/city_block.pdf")
plotStepHist(PSMs, targets_01fdr, decoys, :total_ions,   range(1, 40, length=40), title = "Total Ions", label_a = "Targets ≤ 0.1% FDR", label_b = "Decoys", f_out = "/Users/n.t.wamsley/Projects/TEST_DATA/figs/total_ions.pdf")
merge_pdfs(readdir("/Users/n.t.wamsley/Projects/TEST_DATA/figs/"; join=true), "/Users/n.t.wamsley/Projects/TEST_DATA/figs/discriminant_scores.pdf")

b_range = range(1, 50, length=50)
stephist((PSMs[(PSMs[:,:decoy].==false) .& (PSMs[:,:q_value].<0.01),:rank_hyper]), bins = b_range,alpha = 0.5, normalize=true)
stephist!((PSMs[(PSMs[:,:decoy].==true),:rank_hyper]), bins = b_range, alpha = 0.5, normalize=true)
#stephist!((PSMs[(PSMs[:,:decoy].==false),:rank_hyper]), bins = b_range,alpha = 0.5, normalize=true)

b_range = range(1, 50, length=50)
stephist((PSMs[(PSMs[:,:decoy].==false) .& (PSMs[:,:q_value].<0.01),:rank_scribe]), bins = b_range,alpha = 0.5, normalize=true)
stephist!((PSMs[(PSMs[:,:decoy].==true),:rank_scribe]), bins = b_range, alpha = 0.5, normalize=true)
#stephist!((PSMs[(PSMs[:,:decoy].==false),:rank_scribe]), bins = b_range,alpha = 0.5, normalize=true)

b_range = range(-50, 0, length=100)
stephist((PSMs[(PSMs[:,:decoy].==false) .& (PSMs[:,:q_value].<0.01),:poisson]), bins = b_range,alpha = 0.5, normalize=true)
stephist!((PSMs[(PSMs[:,:decoy].==true),:poisson]), bins = b_range, alpha = 0.5, normalize=true)
#stephist!((PSMs[(PSMs[:,:decoy].==false),:poisson]), bins = b_range,alpha = 0.5, normalize=true)

b_range = range(0, 40, length=100)
stephist((PSMs[(PSMs[:,:decoy].==false) .& (PSMs[:,:q_value].<0.01),:RT_error]), bins = b_range,alpha = 0.5, normalize=true)
stephist!((PSMs[(PSMs[:,:decoy].==true),:RT_error]), bins = b_range, alpha = 0.5, normalize=true)
#stephist!((PSMs[(PSMs[:,:decoy].==false),:RT_error]), bins = b_range,alpha = 0.5, normalize=true)

b_range = range(-5, 0, length=1000)
stephist((PSMs[(PSMs[:,:decoy].==false) .& (PSMs[:,:q_value].<0.01),:city_block]), bins = b_range,alpha = 0.5, normalize=true)
stephist!((PSMs[(PSMs[:,:decoy].==true),:city_block]), bins = b_range, alpha = 0.5, normalize=true)
#stephist!((PSMs[(PSMs[:,:decoy].==false),:city_block]), bins = b_range,alpha = 0.5, normalize=true)


b_range = range(1, 40, length=40)
stephist((PSMs[(PSMs[:,:decoy].==false) .& (PSMs[:,:q_value].<0.01),:total_ions]), bins = b_range,alpha = 0.5, normalize=true)
stephist!((PSMs[(PSMs[:,:decoy].==true),:total_ions]), bins = b_range, alpha = 0.5, normalize=true)