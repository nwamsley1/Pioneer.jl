include("loadParamsAndData.jl")


###########
#Pre-Search
#Need to Estimate the following from a random sample of high-confidence targets
#1) Fragment Mass error/correction
#2) Fragment Mass tolerance
#3) iRT to RT conversion spline
###########
println("Begining Presearch")
presearch_time = @timed begin
    include("parameterTuningSearch.jl")
end
println("Finished presearch in ", presearch_time.time, " seconds")

###########kkkk
#Main PSM Search
###########
println("Begining Main Search...")
main_search_time = @timed begin
    include("mainSearch.jl")
end
println("Finished main search in ", main_search_time.time, " seconds")
#PSMs_Dict = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_020824_M0.jld2"))["PSMs_Dict"]
############
#Build Retention Time Index
println("Combining Main Search Results...")
combine_results_time = @timed begin
    include("combineMainSearchResults.jl")
end
println("Combined main search results in ", combine_results_time.time, " seconds")
############
#New Inplace Arrays for Integration
println("Begining Quantitative Search...")
quant_search_time = @timed begin
    include("quantitativeSearch.jl")
end
println("Combined main search results in ", quant_search_time.time, " seconds")

println("Begining XGBoost...")
score_traces_time = @timed begin
    include("scoreTraces.jl")
end
println("Scored Traces In ", score_traces_time.time, " seconds")



