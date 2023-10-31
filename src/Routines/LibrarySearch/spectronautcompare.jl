combined_set = load("/Users/n.t.wamsley/Desktop/good_precursors_a.jld2")
combined_set  = combined_set["combined_set"]
combined_set = Set(replace.(combined_set, "[Carbamidomethyl (C)]" => ""))

MS_TABLE = Arrow.Table(MS_TABLE_PATH) 

main_search_params[:min_spectral_contrast] = Float32(-Inf)
main_search_params[:matched_ratio] = Float32(-Inf)
main_search_params[matched_ratio] = Float32(-Inf)
combined_set = load("C://Users//n.t.wamsley//Desktop//good_precursors_a.jld2")
combined_set = collect(combined_set["combined_set"])
combined_set = Set(replace.(combined_set, "[Carbamidomethyl (C)]" => ""))

sub_search_time = @timed PSMs = mainLibrarySearch(
    MS_TABLE,
    prosit_lib["f_index"],
    prosit_lib["f_det"],
    RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
    UInt32(ms_file_idx), #MS_FILE_IDX
    frag_err_dist_dict[ms_file_idx],
    16.1,
    main_search_params,
    #scan_range = (201389, 204389),
    #scan_range = (55710, 55710),
    #scan_range = (50426, 51000),
    scan_range = (1, length(MS_TABLE[:masses]))
);

refinePSMs!(PSMs, MS_TABLE, prosit_lib["precursors"]);
psms_file_name = "nOF5_test_PSMs.jld2"

jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", psms_file_name); PSMs);
PSMs = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", psms_file_name));
PSMs = PSMs["PSMs"]
###########
#All passing seed 
refinePSMs!(PSMs, MS_TABLE, prosit_lib["precursors"]);
PSMs[:,:stripped_sequence] = replace.(PSMs[:,:sequence], "M(ox)" => "M");
pioneer_passing = Set("_".*PSMs[!,:stripped_sequence].*"_.".*string.(PSMs[!,:charge]))
setdiff(combined_set, pioneer_passing) #6524


###########
#All passing 25% FDR 
psms_passing = PSMs[PSMs[:,:q_value].<=0.25,:]
pioneer_passing = Set("_".*psms_passing[!,:stripped_sequence].*"_.".*string.(psms_passing[!,:charge]))
setdiff(combined_set, pioneer_passing) #8433

PSMs[!,:AIE] = log2.((2 .^ Float32.(PSMs[!,:TIC])) .* (2 .^ Float32.(PSMs[!,:log2_intensity_explained])))
FORM = FormulaTerm(
        (Term(:target),),
        (
        Term(:scribe_corrected),
         Term(:spectral_contrast_corrected),
         Term(:spectral_contrast),
         Term(:scribe),
         Term(:RT_error),
         Term(:missed_cleavage),
         Term(:Mox),
         Term(:TIC),
         Term(:city_block),
         Term(:entropy_score),
         #Term(:adjusted_intensity_explained),
         Term(:total_ions),
         #Term(:error),
         Term(:err_norm),
         #Term(:b_y_overlap),
         Term(:charge),
         #Term(:total_ions) & Term(:entropy_score)
    ))
    model_fit = glm(FORM, PSMs, 
                            Binomial(), 
                            ProbitLink())
    Y′ = Float16.(GLM.predict(model_fit, PSMs));
    getQvalues!(PSMs, allowmissing(Y′),  allowmissing(PSMs[:,:decoy]));
    
psms_passing = PSMs[PSMs[:,:q_value].<=0.25,:]
pioneer_passing = Set("_".*psms_passing[!,:stripped_sequence].*"_.".*string.(psms_passing[!,:charge]))
setdiff(combined_set, pioneer_passing) #7745

filter!(x->x.q_value<=0.25, PSMs);
setdiff(combined_set, pioneer_passing) #7745
filter!(x->x.matched_ratio>0, PSMs)
pioneer_passing = Set("_".*PSMs[!,:stripped_sequence].*"_.".*string.(PSMs[!,:charge]))
setdiff(combined_set, pioneer_passing) #7833
filter!(x->x.spectral_contrast_corrected>0.5, PSMs)
pioneer_passing = Set("_".*PSMs[!,:stripped_sequence].*"_.".*string.(PSMs[!,:charge]))
setdiff(combined_set, pioneer_passing) #7879
filter!(x->x.RT_error<10.0, PSMs)
pioneer_passing = Set("_".*PSMs[!,:stripped_sequence].*"_.".*string.(PSMs[!,:charge]))
setdiff(combined_set, pioneer_passing) #8363

################
#Enforce iRT bound. 
sub_search_time = @timed PSMs = mainLibrarySearch(
    MS_TABLE,
    prosit_lib["f_index"],
    prosit_lib["precursors"],
    prosit_lib["f_det"],
    RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
    UInt32(ms_file_idx), #MS_FILE_IDX
    frag_err_dist_dict[ms_file_idx],
    16.1,
    main_search_params,
    #scan_range = (201389, 204389),
    #scan_range = (55710, 55710),
    #scan_range = (50426, 51000),
    scan_range = (1, length(MS_TABLE[:masses]))
);

refinePSMs!(PSMs, MS_TABLE, prosit_lib["precursors"]);
psms_file_name = "nOF5_test_PSMs_r.jld2"

jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", psms_file_name); PSMs);
PSMs = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", psms_file_name));
PSMs = PSMs["PSMs"]
###########
#All passing seed 
refinePSMs!(PSMs, MS_TABLE, prosit_lib["precursors"]);
PSMs[:,:stripped_sequence] = replace.(PSMs[:,:sequence], "M(ox)" => "M");
pioneer_passing = Set("_".*PSMs[!,:stripped_sequence].*"_.".*string.(PSMs[!,:charge]))
setdiff(combined_set, pioneer_passing) #6633

filter!(x->x.RT_error<10.0, PSMs)
pioneer_passing = Set("_".*PSMs[!,:stripped_sequence].*"_.".*string.(PSMs[!,:charge]))
setdiff(combined_set, pioneer_passing) #6955

refinePSMs!(PSMs, MS_TABLE, prosit_lib["precursors"]);
psms_passing = PSMs[PSMs[:,:q_value].<=0.25,:]
pioneer_passing = Set("_".*psms_passing[!,:stripped_sequence].*"_.".*string.(psms_passing[!,:charge]))
setdiff(combined_set, pioneer_passing) #8601

PSMs[!,:AIE] = log2.((2 .^ Float32.(PSMs[!,:TIC])) .* (2 .^ Float32.(PSMs[!,:log2_intensity_explained])))
FORM = FormulaTerm(
        (Term(:target),),
        (
        Term(:scribe_corrected),
         Term(:spectral_contrast_corrected),
         Term(:spectral_contrast),
         Term(:scribe),
         Term(:RT_error),
         Term(:missed_cleavage),
         Term(:Mox),
         Term(:TIC),
         Term(:city_block),
         Term(:entropy_score),
         #Term(:adjusted_intensity_explained),
         Term(:total_ions),
         #Term(:error),
         Term(:err_norm),
         #Term(:b_y_overlap),
         Term(:charge),
         #Term(:total_ions) & Term(:entropy_score)
    ))
    model_fit = glm(FORM, PSMs, 
                            Binomial(), 
                            ProbitLink())
    Y′ = Float16.(GLM.predict(model_fit, PSMs));
    getQvalues!(PSMs, allowmissing(Y′),  allowmissing(PSMs[:,:decoy]));
    
psms_passing = PSMs[PSMs[:,:q_value].<=0.25,:]
pioneer_passing = Set("_".*psms_passing[!,:stripped_sequence].*"_.".*string.(psms_passing[!,:charge]))
setdiff(combined_set, pioneer_passing) #8102


filter!(x->x.q_value<=0.25, PSMs);
setdiff(combined_set, pioneer_passing) #8102
filter!(x->x.matched_ratio>0, PSMs)
pioneer_passing = Set("_".*PSMs[!,:stripped_sequence].*"_.".*string.(PSMs[!,:charge])) #289K
setdiff(combined_set, pioneer_passing) #8195
filter!(x->x.spectral_contrast_corrected>0.5, PSMs)
pioneer_passing = Set("_".*PSMs[!,:stripped_sequence].*"_.".*string.(PSMs[!,:charge])) #279K
setdiff(combined_set, pioneer_passing) #8244

sort!(PSMs,:RT); #Sorting before grouping is critical. 
PSMs[!,:max_weight] .= zero(Float32)
test_chroms = groupby(PSMs, :precursor_idx);
Threads.@threads for i in ProgressBar(range(1, length(test_chroms)))
    #for i in range(1, length(grouped_precursor_df))
    best_scan = argmax(test_chroms[i].matched_ratio)
    test_chroms[i].best_scan[best_scan] = true
end
filter!(x->x.best_scan, PSMs)

bins = LinRange(0, 1, 100)
histogram(PSMs[PSMs[!,:decoy].==false,:prob], bins = bins, alpha = 0.5, norm =:probability)
histogram!(PSMs[PSMs[!,:decoy].==true,:prob], bins = bins, alpha = 0.5, norm =:probability)
transform!(PSMs, AsTable(:) => ByRow(psm -> 
prosit_lib["precursors"][psm[:precursor_idx]].mz
) => :prec_mz
);
sort!(PSMs,:RT, rev = false);
#Build RT index of precursors to integrate
rt_index = buildRTIndex(PSMs);

main_search_params[:max_peak_width] = 1.0

ms2_params =  copy(main_search_params)
ms2_params[:min_frag_count] = 1
ms2_params[:min_spectral_contrast] = Float32(-Inf)
ms2_params[:min_matched_ratio] = Float32(-Inf)
ms2_params[:min_topn] = 0
ms2_params[:max_peak_width] = 0.75
test = integrateMS2(MS_TABLE, 
                    prosit_lib["f_det"],
                    rt_index,
                    UInt32(ms_file_idx), 
                    frag_err_dist_dict[ms_file_idx],
                    ms2_params, 
                    scan_range = (1, length(MS_TABLE[:scanNumber]))
                    #scan_range = (101357, 110357)
                    );

    refinePSMs!(test, MS_TABLE, prosit_lib["precursors"]);
    sort!(test,:RT); #Sorting before grouping is critical. 
    test_chroms = groupby(test, :precursor_idx);

    time_test = @timed integratePrecursors(test_chroms, 
                                            n_quadrature_nodes = params_[:n_quadrature_nodes],
                                            intensity_filter_fraction = params_[:intensity_filter_fraction],
                                            LsqFit_tol = params_[:LsqFit_tol],
                                            Lsq_max_iter = params_[:Lsq_max_iter],
                                            tail_distance = params_[:tail_distance])
    filter!(x -> x.best_scan, test);
    test_chroms[104][:,[:sequence,:decoy,:q_value,:RT,:weight,:entropy_score,:matched_ratio]]

features = [ :FWHM,
    :FWHM_01,
    :GOF,
    :H,
    :Mox,
    :RT,
    :RT_error,
    :longest_y,
    :y_count,
    :b_count,
    #:b_ladder,
    :base_width_min,
    :best_rank,
    #:best_over_precs,
    :charge,
    :city_block,
    :data_points,
    :entropy_score,
    :err_norm,
    :error,
    :hyperscore,
    #:intensity_explained,
    :ions_sum,
    :log_sum_of_weights,
    :matched_ratio,
    :mean_log_entropy,
    :mean_log_probability,
    :mean_log_spectral_contrast,
    :mean_matched_ratio,
    :mean_scribe_score,
    :missed_cleavage,
    #:ms1_ms2_diff,
    :n_obs,
    :peak_area,
    #:peak_area_ms1,
    :points_above_FWHM,
    :points_above_FWHM_01,
    :poisson,
    :prec_mz,
    :scribe,
    :scribe_corrected,
    :sequence_length,
    :spectral_contrast,
    :spectral_contrast_corrected,
    #:spectrum_peaks,
    :topn,
    :total_ions,
    :weight,
    #:y_ladder,
    #:ρ,
    :log2_intensity_explained,
    :TIC,
    :adjusted_intensity_explained
    ]
test[:,"sequence_length"] = length.(replace.(test[:,"sequence"], "M(ox)" => "M"));
transform!(test, AsTable(:) => ByRow(psm -> 
prosit_lib["precursors"][psm[:precursor_idx]].mz
) => :prec_mz
);
filter!(x->x.weight>0, test);
for i in range(1, size(test)[1])
    if isinf(test[i,:mean_matched_ratio])
        test[i,:mean_matched_ratio] = Float16(6000)
    end
end
xgboost_time = @timed bst = rankPSMs!(test, 
                        features,
                        colsample_bytree = 1.0, 
                        min_child_weight = 5, 
                        gamma = 1, 
                        subsample = 0.5, 
                        n_folds = 2, 
                        num_round = 200, 
                        max_depth = 10, 
                        eta = 0.05, 
                        #eta = 0.0175,
                        train_fraction = 9.0/9.0,
                        n_iters = 2);

for name in names(test)
    if eltype(test[!,name]) == String
        continue
    else
        if any(isinf.(test[!,name]))
            println(name)
        end
    end
end

best_psms = test
getQvalues!(best_psms, allowmissing(best_psms[:,:prob]), allowmissing(best_psms[:,:decoy]));
best_psms[:,:stripped_sequence] = replace.(best_psms[:,:sequence], "M(ox)" => "M");
best_psms_passing = best_psms[(best_psms[!,:q_value].<=0.01).&(best_psms[!,:decoy].==false),:]
pioneer_passing_fdr = Set("_".*best_psms_passing[!,:stripped_sequence].*"_.".*string.(best_psms_passing[!,:charge]))
setdiff(combined_set, pioneer_passing_fdr)

selectRTIndexedTransitions!