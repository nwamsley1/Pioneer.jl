function pioneer_main()
println("pwd ", pwd())
#total_time = @timed begin
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
println("Finished presearch in ", presearch_time.time, " seconds")

###########
#Main PSM Search
###########
println("Begining Main Search...")
main_search_time = @timed begin
    include(joinpath(methods_path,"firstSearch.jl"))
end

############
#Build Retention Time Index
println("Combining Main Search Results...")
combine_results_time = @timed begin
    include(joinpath(methods_path,"combineFirstSearchResults.jl"))
end
#jldsave(joinpath(results_folder, "iRT_RT_spline.jld2"); iRT_RT)

println("Combined main search results in ", combine_results_time.time, " seconds")
############
#New Inplace Arrays for Integration
println("Begining Quantitative Search...")
BPSMS = Dict{Int64, DataFrame}()
quant_search_time = @timed begin
    include(joinpath(methods_path,"quantitativeSearch.jl"))
end
best_psms = vcat(values(BPSMS)...)
println("Combined main search results in ", quant_search_time.time, " seconds")

###########
#XGBoost
##########
println("Begining XGBoost...")
score_traces_time = @timed begin
    include(joinpath(methods_path,"scoreTraces.jl"))
end
traces_passing = Set(best_psms[(best_psms[!,:q_value].<=0.01).&(best_psms[!,:target]),:precursor_idx])
println("Scored Traces In ", score_traces_time.time, " seconds")

###########
#Score Protein Groups
scored_proteins = scoreProteinGroups!(best_psms)
protein_to_q_value = Dict{Tuple{UInt32, String}, Float32}()
for i in ProgressBar(range(1, size(scored_proteins, 1)))
    protein_to_q_value[
        (
            UInt32(scored_proteins[i,:ms_file_idx]),
            scored_proteins[i,:accession_numbers]
        )
    ] = scored_proteins[i,:q_value]
end

###########
#Re-quantify with 1% fdr precursors 
best_psms[!,:peak_area] = zeros(Float32, size(best_psms, 1))
best_psms[!,:new_best_scan] = zeros(UInt32, size(best_psms, 1))
grouped_best_psms = groupby(best_psms, :file_name)
println("Begining Quantitative Search...")
quant_search_time = @timed begin
    include(joinpath(methods_path,"secondQuant.jl"))
end

###########
#Ungroup precursors and get the best trace for each 
best_psms = DataFrame(grouped_best_psms)
#Must be based on weight and not peak_area because peak_area was not calculated for decoys
getBestTrace!(best_psms, 0.01, :weight)
filter!(x->x.best_trace, best_psms)
#precursor level q_value 
getQvalues!(best_psms[!,:prob], best_psms[:,:target], best_psms[!,:q_value]);
#filter out decoys and low-scoring targets
filter!(x->(x.q_value<=0.01)&(x.target)&(!isnan(x.peak_area)), best_psms)
#IDs_PER_FILE = value_counts(best_psms, [:file_name])
best_psms[!,:species] = [precursors[:proteome_identifiers][pid] for pid in best_psms[!,:precursor_idx]]
###########
#Normalize Quant 
###########
println("Cross-Run Normalization of Precursor Level Quant...")
score_traces_time = @timed begin
    include(joinpath(methods_path,"normalizeQuant.jl"))
end

############
#Prep for Protein Inference 
best_psms[!,:species] = [precursors[:proteome_identifiers][pid] for pid in best_psms[!,:precursor_idx]]
best_psms[!,:ms_file_idx] =  UInt32.(best_psms[!,:ms_file_idx])
best_psms[!,:peak_area] =  allowmissing(best_psms[!,:peak_area])
best_psms[!,:peak_area_normalized] =  allowmissing(best_psms[!,:peak_area_normalized])
gbpsms = groupby(best_psms,:ms_file_idx)
file_idx_dicts = [Dict{UInt32, @NamedTuple{prob::Float32, qvalue::Float32, peak_area::Float32}}() for _ in range(1, length(gbpsms))]
for (file_idx, bpsms) in ProgressBar(pairs(gbpsms))
    for i in range(1, size(bpsms, 1))
            file_idx_dicts[file_idx[:ms_file_idx]][bpsms[i,:precursor_idx]] = (
              prob = bpsms[i,:prob],
             qvalue = bpsms[i,:q_value],
             peak_area = bpsms[i,:peak_area])
    end
end
best_psms[!,:structural_mods] = [precursors[:structural_mods][pid] for pid in best_psms[!,:precursor_idx]]
best_psms[!,:isotopic_mods] = [precursors[:isotopic_mods][pid] for pid in best_psms[!,:precursor_idx]]

features = [:species,:accession_numbers,:sequence,:structural_mods,
             :isotopic_mods,:charge,:precursor_idx,:target,:weight,:peak_area,:peak_area_normalized,
             :scan_idx,:prob,:q_value,:prec_mz,:RT,:irt_obs,:irt_pred,:best_rank,:best_rank_iso,:topn,:topn_iso,:longest_y,:longest_b,:b_count,
             :y_count,:p_count,:non_cannonical_count,:isotope_count
             ,:matched_ratio]
sort!(best_psms,[:species,:accession_numbers,:sequence,:target])
Arrow.write(joinpath(results_folder,"best_psms.arrow"),best_psms[!,features])
wide_psms_quant = unstack(best_psms,[:species,:accession_numbers,:sequence,:structural_mods,:isotopic_mods,:precursor_idx,:target],:file_name,:peak_area_normalized)
sort!(wide_psms_quant,[:species,:accession_numbers,:sequence,:target])
CSV.write(joinpath(results_folder,"best_psms_wide.arrow"),wide_psms_quant)

#Summarize Precursor ID's
value_counts(df, col) = combine(groupby(df, col), nrow)

precursor_id_table = value_counts(best_psms,:file_name)
#CSV.write("/Users/n.t.wamsley/Desktop/precursor_ids_table.csv")



best_psms_old = copy(best_psms)
score_traces_time = @timed begin
    include(joinpath(methods_path,"proteinQuant.jl"))
end
protein_quant[!,:experiments] = UInt32.(protein_quant[!,:experiments])
protein_quant[!,:file_name] = [file_id_to_parsed_name[ms_file_idx] for ms_file_idx in protein_quant[!,:experiments]]

#Wide DataFormat 
wide_protein_quant = unstack(protein_quant,[:species,:protein,:target],:experiments,:log2_abundance)
sort!(wide_protein_quant,[:species,:protein,:target])
CSV.write(joinpath(results_folder,"proteins_wide.csv"),wide_protein_quant)
###########
#QC Plots
###########
println("Generating QC Plots...")
score_traces_time = @timed begin
    include(joinpath(methods_path,"qcPlots.jl"))
end

end