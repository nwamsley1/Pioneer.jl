# allocate separate vectors
prob_test   = zeros(Float32, nrow(psms))  # final CV predictions
prob_train  = zeros(Float32, nrow(psms))  # temporary, used during training
MBR_estimates = zeros(Float32, nrow(psms)) # optional MBR layer

for test_fold_idx in unique_cv_folds
    # split data
    train_idx = train_indices[test_fold_idx]
    test_idx  = fold_indices[test_fold_idx]

    psms_train = @view psms[train_idx, :]
    psms_test  = @view psms[test_idx, :]

    for (itr, num_round) in enumerate(iter_scheme)
        # get filtered training data for this iteration
        psms_train_itr = get_training_data_for_iteration!(psms_train, ...)

        # features for this iteration
        train_feats = itr < mbr_start_iter ? non_mbr_features : features

        # train booster
        bst = train_booster(psms_train_itr, train_feats, num_round; ...)

        # **temporary predictions for training only**
        prob_train[train_idx] = predict(bst, psms_train)
        psms_train[!,:prob] = prob_train[train_idx]

        # **predict held-out fold**
        prob_test[test_idx] = predict(bst, psms_test)
        psms_test[!,:prob] = prob_test[test_idx]

        # optionally update MBR features for next iteration
        if match_between_runs
            update_mbr_features!(psms_train, psms_test, prob_train, test_idx, itr, mbr_start_iter, max_q_value_xgboost_rescore)
        end
    end

    # after all iterations of this fold
    if match_between_runs
        MBR_estimates[test_idx] = psms_test.prob
    else
        prob_test[test_idx] = psms_test.prob
    end
end

# assign final probabilities
if match_between_runs
    psms[!,:prob] = MBR_estimates
else
    psms[!,:prob] = prob_test
end
