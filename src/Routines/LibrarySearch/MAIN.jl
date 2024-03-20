println("pwd ", pwd())
const methods_path = joinpath(pwd(), "src","Routines","LibrarySearch")
include(joinpath(methods_path,"loadParamsAndData.jl"))
###########
#Pre-Search
#Need to Estimate the following from a random sample of high-confidence targets
#1) Fragment Mass error/correction
#2) Fragment Mass tolerance
#3) iRT to RT conversion spline
###########
println("Begining Presearch")
presearch_time = @timed begin
    include(joinpath(methods_path, "parameterTuningSearch.jl"))
end
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "irt_errs.jld2"); irt_errs)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "RT_to_iRT_map_dict.jld2"); RT_to_iRT_map_dict)
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "frag_err_dist_dict.jld2"); frag_err_dist_dict)

irt_errs = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "irt_errs.jld2"))["irt_errs"]
RT_to_iRT_map_dict = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "RT_to_iRT_map_dict.jld2"))["RT_to_iRT_map_dict"]
frag_err_dist_dict = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "frag_err_dist_dict.jld2"))["frag_err_dist_dict"]

println("Finished presearch in ", presearch_time.time, " seconds")

###########
#Main PSM Search
###########
println("Begining Main Search...")
main_search_time = @timed begin
    include(joinpath(methods_path,"firstSearch.jl"))
end
jldsave(joinpath(results_folder, "PSMs_Dict.jld2"); PSMs_Dict)
println("Finished main search in ", main_search_time.time, " seconds")
#PSMs_Dict = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_020824_M0.jld2"))["PSMs_Dict"]
############
#Build Retention Time Index
println("Combining Main Search Results...")
combine_results_time = @timed begin
    include(joinpath(methods_path,"combineFirstSearchResults.jl"))
end
jldsave(joinpath(results_folder, "iRT_RT_spline.jld2"); iRT_RT)
jldsave(joinpath(results_folder,  "RT_iRT_spline.jld2"); RT_iRT)
jldsave(joinpath(results_folder, "precID_to_iRT.jld2"); precID_to_iRT)
jldsave(joinpath(results_folder, "RT_INDICES.jld2"); RT_INDICES)
jldsave(joinpath(results_folder, "precID_to_cv_fold.jld2"); precID_to_cv_fold)
println("Combined main search results in ", combine_results_time.time, " seconds")
############
#New Inplace Arrays for Integration
println("Begining Quantitative Search...")
quant_search_time = @timed begin
    include(joinpath(methods_path,"quantitativeSearch.jl"))
end
println("Combined main search results in ", quant_search_time.time, " seconds")

###########
#XGBoost
###########
println("Begining XGBoost...")
score_traces_time = @timed begin
    include(joinpath(methods_path,"scoreTraces.jl"))
end
jldsave(joinpath(results_folder, "best_psms_scored.jld2"); best_psms)
best_psms_passing = best_psms[best_psms[!,:q_value].<=0.01,:]
jldsave(joinpath(results_folder, "best_psms_passing.jld2"); best_psms_passing)
println("Scored Traces In ", score_traces_time.time, " seconds")

###########
#QC Plots
###########
println("Generating QC Plots...")
score_traces_time = @timed begin
    include(joinpath(methods_path,"qcPlots.jl"))
end