jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_54M_110623.jld2"); PSMs)
combined_set = load("/Users/n.t.wamsley/Desktop/good_precursors_a.jld2")
combined_set  = combined_set["combined_set"]
combined_set = Set(replace.(combined_set, "[Carbamidomethyl (C)]" => ""))

PSMs[!,:sequence] .= ""
PSMs[!,:charge] .= zero(UInt8)
for i in ProgressBar(range(1, size(PSMs)[1]))
    PSMs[i,:sequence] = prosit_lib["precursors"][PSMs[i,:precursor_idx]].sequence;
    PSMs[i,:charge] = UInt8(getCharge(prosit_lib["precursors"][PSMs[i,:precursor_idx]]));
end
PSMs[:,:stripped_sequence] = replace.(PSMs[:,:sequence], "M(ox)" => "M");
#jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_54M_110623.jld2"); PSMs)
pioneer_passing = Set("_".*PSMs[!,:stripped_sequence].*"_.".*string.(PSMs[!,:charge]))
setdiff(combined_set, pioneer_passing)

PSMs_SUB = PSMs[PSMs[!,:best_rank].==1,:]
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)


PSMs_SUB = PSMs[PSMs[!,:topn].>1,:]
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)


PSMs_SUB = PSMs[(PSMs[!,:topn].>1).&(PSMs[!,:best_rank].==1),:]
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)

PSMs_SUB = PSMs[(PSMs[!,:matched_ratio].>=(-1.0)),:]
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)

PSMs_SUB = PSMs[(PSMs[!,:matched_ratio].>=(-1.0)).&((PSMs[!,:topn].>1).&(PSMs[!,:best_rank].==1)),:]
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)

PSMs_SUB = PSMs[((PSMs[!,:y_count].+PSMs[!,:b_count]).>=(3)),:]
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)

PSMs_SUB = PSMs[((PSMs[!,:y_count].+PSMs[!,:b_count]).>=(4)),:]
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)

PSMs_SUB = PSMs[((PSMs[!,:y_count].+PSMs[!,:b_count]).>=(5)),:]
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)

PSMs_SUB = PSMs[((PSMs[!,:y_count].+PSMs[!,:b_count]).>=(4)).&(PSMs[!,:matched_ratio].>=(-1.0)).&((PSMs[!,:topn].>1).&(PSMs[!,:best_rank].==1)),:]
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)

PSMs_SUB = PSMs[((PSMs[!,:y_count].+PSMs[!,:b_count]).>=(5)).&(PSMs[!,:matched_ratio].>=(-1.0)).&((PSMs[!,:topn].>1).&(PSMs[!,:best_rank].==1)),:]
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)

PSMs_SUB = PSMs[((PSMs[!,:y_count].+PSMs[!,:b_count]).>=(4)).&(PSMs[!,:matched_ratio].>=(-1.0)).&((PSMs[!,:topn].>1).&(PSMs[!,:best_rank].==1)).&(PSMs[!,:spectral_contrast].>=(0.5)),:]
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)

PSMs_SUB = PSMs[((PSMs[!,:y_count].+PSMs[!,:b_count]).>=(5)).&(PSMs[!,:matched_ratio].>=(-1.0)).&((PSMs[!,:topn].>1).&(PSMs[!,:best_rank].==1)).&(PSMs[!,:spectral_contrast].>=(0.5)),:]
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)

PSMs_SUB = copy(PSMs[((PSMs[!,:y_count].+PSMs[!,:b_count]).>=(5)).&(PSMs[!,:matched_ratio].>=(-1.0)).&((PSMs[!,:topn].>1).&(PSMs[!,:best_rank].==1)).&(PSMs[!,:spectral_contrast].>=(0.5)),:])
refinePSMs!(PSMs_SUB, MS_TABLE, prosit_lib["precursors"], max_rt_error = 10.0);
filter!(x->x.q_value<=0.25, PSMs_SUB);
sort!(PSMs_SUB,:RT); #Sorting before grouping is critical. 
PSMs_SUB[!,:max_weight] .= zero(Float32)
grouped_PSMs = groupby(PSMs_SUB, :precursor_idx);
for i in range(1, length(grouped_PSMs))
    best_scan = argmax(grouped_PSMs[i].matched_ratio)
    grouped_PSMs[i].best_scan[best_scan] = true
end
filter!(x->x.best_scan, PSMs_SUB);
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)

PSMs_SUB = copy(PSMs[((PSMs[!,:y_count].+PSMs[!,:b_count]).>=(5)).&(PSMs[!,:matched_ratio].>=(-1.0)).&((PSMs[!,:topn].>1).&(PSMs[!,:best_rank].==1)).&(PSMs[!,:spectral_contrast].>=(0.5)),:])
refinePSMs!(PSMs_SUB, MS_TABLE, prosit_lib["precursors"], max_rt_error = 10.0);
filter!(x->x.q_value<=0.25, PSMs_SUB);
sort!(PSMs_SUB,:RT); #Sorting before grouping is critical. 
PSMs_SUB[!,:max_weight] .= zero(Float32)
grouped_PSMs = groupby(PSMs_SUB, :precursor_idx);
for i in range(1, length(grouped_PSMs))
    best_scan = argmax(grouped_PSMs[i].matched_ratio)
    grouped_PSMs[i].best_scan[best_scan] = true
end
filter!(x->x.best_scan, PSMs_SUB);
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)

PSMs_SUB = copy(PSMs[((PSMs[!,:y_count].+PSMs[!,:b_count]).>=(5)).&(PSMs[!,:matched_ratio].>=(-1.0)).&((PSMs[!,:topn].>1).&(PSMs[!,:best_rank].==1)).&(PSMs[!,:spectral_contrast].>=(0.5)),:])
refinePSMs!(PSMs_SUB, MS_TABLE, prosit_lib["precursors"], max_rt_error = 10.0);
PSMs_SUB[:,:q_value] .= zero(Float16);
FORM = FormulaTerm(
    (Term(:target),),
    (
     Term(:scribe_corrected),
     Term(:spectral_contrast_corrected),
     Term(:spectral_contrast),
     Term(:hyperscore),
     Term(:poisson),
     Term(:scribe),
     Term(:RT_error),
     Term(:missed_cleavage),
     Term(:Mox),
     Term(:TIC),
     Term(:city_block),
     Term(:entropy_score),
     #Term(:entropy_score) & Term(:adjusted_intensity_explained),
     Term(:total_ions),
     Term(:err_norm),
     Term(:charge),
     Term(:spectrum_peak_count),
     Term(:adjusted_intensity_explained)
))

#=
for name in names(PSMs_SUB)
    if eltype(PSMs_SUB[!,name]) == String
        continue
    else
        if any(isinf.(PSMs_SUB[!,name]))
            println(name)
            println(argmax(isinf.(PSMs_SUB[!,name])))
        end
    end
end
=#
model_fit = glm(FORM,PSMs_SUB, 
                        Binomial(), 
                        ProbitLink())
Y′ = Float16.(GLM.predict(model_fit, PSMs_SUB));
getQvalues!(PSMs_SUB, allowmissing(Y′),  allowmissing(PSMs_SUB[:,:decoy]));
model_fit = glm(FORM,PSMs_SUB[(PSMs_SUB.decoy) .| ((PSMs_SUB.q_value.<=0.01)),:], 
                        Binomial(), 
                        ProbitLink())
Y′ = Float16.(GLM.predict(model_fit, PSMs_SUB));
getQvalues!(PSMs_SUB, allowmissing(Y′),  allowmissing(PSMs_SUB[:,:decoy]));
model_fit = glm(FORM,PSMs_SUB[(PSMs_SUB.decoy) .| ((PSMs_SUB.q_value.<=0.01)),:], 
                        Binomial(), 
                        ProbitLink())
Y′ = Float16.(GLM.predict(model_fit, PSMs_SUB));
getQvalues!(PSMs_SUB, allowmissing(Y′),  allowmissing(PSMs_SUB[:,:decoy]));
println("Target PSMs at 25% FDR: ", sum((PSMs_SUB.q_value.<=0.25).&(PSMs_SUB.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs_SUB.q_value.<=0.01).&(PSMs_SUB.decoy.==false)))
pioneer_passing = Set("_".*PSMs_SUB[(PSMs_SUB.q_value.<=0.25).&(PSMs_SUB.decoy.==false),:stripped_sequence].*"_.".*string.(PSMs_SUB[(PSMs_SUB.q_value.<=0.25).&(PSMs_SUB.decoy.==false),:charge]))
setdiff(combined_set, pioneer_passing)

features = [:RT_error,:entropy_score]
lda = fit(MulticlassLDA, Float32.(Matrix(PSMs_SUB[!,features]))', PSMs_SUB[!,:decoy]; outdim = 2)
Ylda = MultivariateStats.predict(lda,  Matrix(PSMs_SUB[!,features])')[2,:]
getQvalues!(PSMs_SUB, allowmissing(Ylda),  allowmissing(PSMs_SUB[:,:decoy]));

X = Float32.(Matrix(PSMs_SUB[1:1:end,features]))'
X_labels = Vector(PSMs_SUB[1:1:end,:decoy]).==false
#X_labels = [Dict(false=>0.0, true=>1.0)[x] for x in X_labels]
lda = fit(MulticlassLDA, X, X_labels; outdim=2)
Ylda =  MultivariateStats.predict(lda, X)



filter!(x->x.q_value<=0.25, PSMs_SUB);
sort!(PSMs_SUB,:RT); #Sorting before grouping is critical. 
PSMs_SUB[!,:max_weight] .= zero(Float32)
grouped_PSMs = groupby(PSMs_SUB, :precursor_idx);
for i in range(1, length(grouped_PSMs))
    best_scan = argmax(grouped_PSMs[i].matched_ratio)
    grouped_PSMs[i].best_scan[best_scan] = true
end
filter!(x->x.best_scan, PSMs_SUB);
pioneer_passing = Set("_".*PSMs_SUB[!,:stripped_sequence].*"_.".*string.(PSMs_SUB[!,:charge]))
setdiff(combined_set, pioneer_passing)

Term(:spectral_contrast),
Term(:scribe),
Term(:RT_error),
Term(:missed_cleavage),
Term(:Mox),
Term(:TIC),
Term(:city_block),
Term(:entropy_score),
Term(:total_ions),
Term(:err_norm),
Term(:charge),
Term(:spectrum_peak_count)
features = [
    :spectral_contrast,
    :spectral_contrast_corrected,
    :RT_error,
    :scribe,
    :scribe_corrected,
    :Mox,
    :TIC,
    :city_block,
    :entropy_score,
    :total_ions,
    :err_norm,
    :charge,
    :spectrum_peak_count
]
X = Matrix(PSMs_SUB[!,features])
y = ones(Float32, length(PSMs_SUB[!,:decoy]))
for i in eachindex(y)
    if PSMs_SUB[i,:decoy]
        y[i] = Float32(-1)
    else
        y[i] = Float32(1)
    end
end
λ = 0.5
logistic = LogisticRegression(λ)
theta = MLJLinearModels.fit(logistic,X, PSMs_SUB[!,:decoy])
ypred =X*reshape(theta, length(theta), 1)

model = LogisticClassifier(lambda = 1.0, penalty = :l2, scale_penalty_with_samples = false)
mach = MLJ.machine(model, PSMs_SUB[!,features], categorical(PSMs_SUB[!,:decoy].==false))
fit!(mach)
ypred = MLJ.predict(mach, PSMs_SUB[!,features])
Y′ = zeros(Float64, length(ypred))
for i in eachindex(ypred) 
    Y′[i] = ypred[i].prob_given_ref[2]
end
getQvalues!(PSMs_SUB, allowmissing(Y′),  allowmissing(PSMs_SUB[:,:decoy].==false));
println("Target PSMs at 25% FDR: ", sum((PSMs_SUB.q_value.<=0.25).&(PSMs_SUB.decoy.==false)))
println("Target PSMs at 1% FDR: ", sum((PSMs_SUB.q_value.<=0.01).&(PSMs_SUB.decoy.==false)))









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
###########pioneer_passing_fdr = Set("_".*best_psms_passing[!,:stripped_sequence].*"_.".*string.(best_psms_passing[!,:charge]))
setdiff(combined_set, pioneer_passing_fdr)
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

sort!(PSMs,:RT); #Sorting before grouping is critical. 
PSMs[!,:max_weight] .= zero(Float32)
test_chroms = groupby(PSMs, :precursor_idx);
Threads.@threads for i in ProgressBar(range(1, length(test_chroms)))
    #for i in range(1, length(grouped_precursor_df))
    best_scan = argmax(test_chroms[i].matched_ratio)
    test_chroms[i].best_scan[best_scan] = true
end
filter!(x->x.best_scan, PSMs)
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

for name in names(best_psms)
    if eltype(best_psms[!,name]) == String
        continue
    else
        if any(isnan.(best_psms[!,name]))
            println(name)
        end
    end
end


best_psms = test
getQvalues!(best_psms, allowmissing(best_psms[:,:prob]), allowmissing(best_psms[:,:decoy]));
best_psms_passing = best_psms[(best_psms[!,:q_value].<=0.01).&(best_psms[!,:decoy].==false),:]
best_psms[:,:stripped_sequence] = replace.(best_psms[:,:sequence], "M(ox)" => "M");

pioneer_passing_fdr = Set("_".*best_psms_passing[!,:stripped_sequence].*"_.".*string.(best_psms_passing[!,:charge]))
setdiff(combined_set, pioneer_passing_fdr)

pioneer_passing_fdr = Set("_".*best_psms[!,:stripped_sequence].*"_.".*string.(best_psms[!,:charge]))
setdiff(combined_set, pioneer_passing_fdr)


selectRTIndexedTransitions!

MS_TABLE = Arrow.Table(MS_TABLE_PATH) 

main_search_params[:min_spectral_contrast] = Float32(-Inf)
main_search_params[:min_matched_ratio] = Float32(0.0)

combined_set = load("/Users/n.t.wamsley/Desktop/good_precursors_a.jld2")
combined_set = collect(combined_set["combined_set"])
combined_set = Set(replace.(combined_set, "[Carbamidomethyl (C)]" => ""))


################
###############3
#filter!(x->x.matched_ratio>0, PSMs)
sub_search_time = @timed PSMs = mainLibrarySearch(
    MS_TABLE,
    prosit_lib["f_index"],
    prosit_lib["f_det"],
    RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
    UInt32(ms_file_idx), #MS_FILE_IDX
    frag_err_dist_dict[ms_file_idx],
    16.1,
    main_search_params,
    scan_range = (1, length(MS_TABLE[:masses]))
);

refinePSMs!(PSMs, MS_TABLE, prosit_lib["precursors"], max_rt_error = 10.0);
filter!(x->x.q_value<=0.25, PSMs);
pioneer_passing = Set("_".*PSMs[!,:stripped_sequence].*"_.".*string.(PSMs[!,:charge]))
setdiff(combined_set, pioneer_passing) #8102

sort!(PSMs,:RT); #Sorting before grouping is critical. 
PSMs[!,:max_weight] .= zero(Float32)
test_chroms = groupby(PSMs, :precursor_idx);
Threads.@threads for i in ProgressBar(range(1, length(test_chroms)))
    #for i in range(1, length(grouped_precursor_df))
    best_scan = argmax(test_chroms[i].matched_ratio)
    test_chroms[i].best_scan[best_scan] = true
end
filter!(x->x.best_scan, PSMs)
transform!(PSMs, AsTable(:) => ByRow(psm -> 
prosit_lib["precursors"][psm[:precursor_idx]].mz
) => :prec_mz
);
sort!(PSMs,:RT, rev = false);
#Build RT index of precursors to integrate
rt_index = buildRTIndex(PSMs);



ms2_params = copy(main_search_params)
ms2_params[:min_frag_count] = 1
ms2_params[:min_spectral_contrast] = Float32(-Inf)
ms2_params[:min_matched_ratio] = Float32(-Inf)
ms2_params[:min_topn] = 0
ms2_params[:max_peak_width] = 0.75
MS2_CHROMS = integrateMS2(MS_TABLE, 
                    prosit_lib["f_det"],
                    rt_index,
                    UInt32(ms_file_idx), 
                    frag_err_dist_dict[ms_file_idx],
                    16.1,
                    ms2_params, 
                    scan_range = (1, length(MS_TABLE[:scanNumber]))
                    #scan_range = (101357, 110357)
                    );

_refinePSMs!(MS2_CHROMS, MS_TABLE, prosit_lib["precursors"]);
sort!(MS2_CHROMS,:RT); #Sorting before grouping is critical. 
MS2_CHROMS_GROUPED = groupby(MS2_CHROMS, :precursor_idx);

time_test = @timed integratePrecursors(MS2_CHROMS_GROUPED, 
                                        n_quadrature_nodes = params_[:n_quadrature_nodes],
                                        intensity_filter_fraction = params_[:intensity_filter_fraction],
                                        LsqFit_tol = params_[:LsqFit_tol],
                                        Lsq_max_iter = params_[:Lsq_max_iter],
                                        tail_distance = params_[:tail_distance])
filter!(x -> x.best_scan, MS2_CHROMS);

filter!(x->x.weight>0, MS2_CHROMS);
for i in range(1, size(MS2_CHROMS)[1])
    if isinf(MS2_CHROMS[i,:mean_matched_ratio])
        MS2_CHROMS[i,:mean_matched_ratio] = Float16(6000)
    end
end
transform!( MS2_CHROMS, AsTable(:) => ByRow(psm -> 
prosit_lib["precursors"][psm[:precursor_idx]].mz
) => :prec_mz
);
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

xgboost_time = @timed bst = rankPSMs!(MS2_CHROMS, 
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

getQvalues!(MS2_CHROMS, allowmissing(MS2_CHROMS[:,:prob]), allowmissing(MS2_CHROMS[:,:decoy]));
MS2_CHROMS[:,:stripped_sequence] = replace.(MS2_CHROMS[:,:sequence], "M(ox)" => "M");
best_psms_passing = MS2_CHROMS[(MS2_CHROMS[!,:q_value].<=0.01).&(MS2_CHROMS[!,:decoy].==false),:]
pioneer_passing_fdr = Set("_".*best_psms_passing[!,:stripped_sequence].*"_.".*string.(best_psms_passing[!,:charge]))
setdiff(combined_set, pioneer_passing_fdr)

best_psms_passing = MS2_CHROMS[(MS2_CHROMS[!,:q_value].<=0.01).&(MS2_CHROMS[!,:decoy].==false).&(MS2_CHROMS[!,:data_points].>4),:]
pioneer_passing_fdr = Set("_".*best_psms_passing[!,:stripped_sequence].*"_.".*string.(best_psms_passing[!,:charge]))
setdiff(combined_set, pioneer_passing_fdr)


MS2_CHROMS = MS2_CHROMS["best_psms"]
MS2_CHROMS_GROUPED = groupby(MS2_CHROMS, [:file_path,:sequence,:charge]);

MS2_CHROMS_GROUPED[( file_path = "/Users/n.t.wamsley/TEST_DATA/PXD028735/LFQ_Orbitrap_AIF_Condition_B_Sample_Beta_03.arrow", 
                        sequence = "VTAEVVLAHLGGGSTSR", 
                        charge = 0x03
                    )][!,[:sequence,:weight,:entropy_score,:spectral_contrast,:matched_ratio,:RT]]

gx, gw = gausslegendre( params_[:n_quadrature_nodes])
integratePrecursorMS2(MS2_CHROMS_GROUPED[( file_path = "/Users/n.t.wamsley/TEST_DATA/PXD028735/LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_03.arrow", 
                    sequence = "VTAEVVLAHLGGGSTSR", 
                    charge = 0x03
                    )],
                    gx::Vector{Float64},
                    gw::Vector{Float64},
                    intensity_filter_fraction = params_[:intensity_filter_fraction],
                    α = 0.01f0,
                    half_width_at_α = 0.15f0,
                    LsqFit_tol = 1e-5,
                    Lsq_max_iter = 1000,
                    tail_distance = tail_distance = params_[:tail_distance],
                    isplot = true
                    )

integratePrecursorMS2(MS2_CHROMS_GROUPED[( file_path = "/Users/n.t.wamsley/TEST_DATA/PXD028735/LFQ_Orbitrap_AIF_Condition_B_Sample_Beta_03.arrow", 
                    sequence = "VTAEVVLAHLGGGSTSR", 
                    charge = 0x03
                    )],
                    gx::Vector{Float64},
                    gw::Vector{Float64},
                    intensity_filter_fraction = params_[:intensity_filter_fraction],
                    α = 0.01f0,
                    half_width_at_α = 0.15f0,
                    LsqFit_tol = 1e-5,
                    Lsq_max_iter = 1000,
                    tail_distance = tail_distance = params_[:tail_distance],
                    isplot = true
                    )
histogram(best_psms[(best_psms[!,:decoy].==false).&(best_psms[!,:q_value].<=0.01),:scribe_fitted], bins = LinRange(0, 10, 100), normalize = :probability, alpha = 0.5)
histogram!(best_psms[(best_psms[!,:decoy].==true),:scribe_fitted], bins = LinRange(0, 10, 100), normalize = :probability, alpha = 0.5)

histogram(best_psms[(best_psms[!,:decoy].==false).&(best_psms[!,:q_value].<=0.01),:city_block_fitted], bins = LinRange(0, 4, 100), normalize = :probability, alpha = 0.5)
histogram!(best_psms[(best_psms[!,:decoy].==true),:city_block_fitted], bins = LinRange(0, 4, 100), normalize = :probability, alpha = 0.5)

histogram(best_psms[(best_psms[!,:decoy].==false).&(best_psms[!,:q_value].<=0.01),:entropy_score], bins = LinRange(0, 1, 100), normalize = :probability, alpha = 0.5)
histogram!(best_psms[(best_psms[!,:decoy].==true),:entropy_score], bins = LinRange(0, 1, 100), normalize = :probability, alpha = 0.5)

histogram(best_psms[(best_psms[!,:decoy].==false).&(best_psms[!,:q_value].<=0.01),:mean_log_entropy], bins = LinRange(-2.5, 0, 100), normalize = :probability, alpha = 0.5)
histogram!(best_psms[(best_psms[!,:decoy].==true),:mean_log_entropy], bins = LinRange(-2.5, 0, 100), normalize = :probability, alpha = 0.5)