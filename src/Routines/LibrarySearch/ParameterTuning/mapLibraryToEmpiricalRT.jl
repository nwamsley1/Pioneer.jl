function mapLibraryToEmpiricalRT(psms_paths::Dict{String, String},
                                 rt_to_irt_dict_presearch::Dict{Int64, R}; 
                        min_prob::AbstractFloat = 0.9
                    ) where {R<:RtConversionModel}

    #Dictionaries mapping fild_id names to data
    irt_rt = Dictionary{String, Any}() #File name => KDEmapping from irt to rt 
    rt_irt = Dictionary{String, Any}() #File name => KDEmapping from irt to rt 
    for (key, psms_path) in pairs(psms_paths) #For each data frame 
        psms = Arrow.Table(psms_path)
        best_hits = psms[:prob].>min_prob#Map rts using only the best psms
        if sum(best_hits) > 100
            best_rts = psms[:rt][best_hits]
            best_irts = psms[:irt_predicted][best_hits]
            irt_to_rt_spline = UniformSpline(
                                        best_rts,
                                        best_irts,
                                        3, 
                                        5
            )
            rt_to_irt_spline = UniformSpline(
                best_irts,
                best_rts,
                3, 
                5
            )
            #Build rt=>irt and irt=> rt mappings for the file and add to the dictionaries 
            insert!(irt_rt, key, irt_to_rt_spline)
            insert!(rt_irt, key, rt_to_irt_spline)
        else
            ms_file_idx = first(psms[:ms_file_idx])
            rt_to_irt_spline = rt_to_irt_dict_presearch[ms_file_idx]
            low_rt = min(minimum(psms[:rt]), rt_to_irt_spline.first)
            max_rt = min(maximum(psms[:rt]), rt_to_irt_spline.last)
            bins_ = LinRange(low_rt, max_rt, 100)
            irt_to_rt_spline = UniformSpline(
                rt_to_irt_spline.(bins_),
                collect(bins_),
                3, 
                5
                )
            insert!(irt_rt, key, irt_to_rt_spline)
            insert!(rt_irt, key, rt_to_irt_spline)
        end
    end
    return irt_rt, rt_irt
end

#=
ecoli_test = DataFrame(Tables.columntable(
Arrow.Table("/Users/n.t.wamsley/TEST_DATA/ECOLI_TEST/arrow_out/old_garbage/LFQ_Orbitrap_AIF_Ecoli_01.arrow")))
filter!(x->coalesce(x.centerMz, -1.0)<600, ecoli_test)
filter!(x->coalesce(x.centerMz, -1.0)>500, ecoli_test)
filter!(x->coalesce(x.retentionTime, -1.0)<100, ecoli_test)
filter!(x->coalesce(x.retentionTime, -1.0)<75, ecoli_test)
 Arrow.write("/Users/n.t.wamsley/TEST_DATA/ECOLI_TEST/arrow_out/ecoli_filtered_01.arrow", ecoli_test)

ecoli_test = DataFrame(Tables.columntable(
Arrow.Table("/Users/n.t.wamsley/TEST_DATA/ECOLI_TEST/arrow_out/old_garbage/LFQ_Orbitrap_AIF_Ecoli_02.arrow")))
filter!(x->coalesce(x.centerMz, -1.0)<600, ecoli_test)
filter!(x->coalesce(x.centerMz, -1.0)>500, ecoli_test)
filter!(x->coalesce(x.retentionTime, -1.0)<100, ecoli_test)
filter!(x->coalesce(x.retentionTime, -1.0)<75, ecoli_test)
 Arrow.write("/Users/n.t.wamsley/TEST_DATA/ECOLI_TEST/arrow_out/ecoli_filtered_02.arrow", ecoli_test)

 ecoli_test = DataFrame(Tables.columntable(
Arrow.Table("/Users/n.t.wamsley/TEST_DATA/ECOLI_TEST/arrow_out/old_garbage/LFQ_Orbitrap_AIF_Ecoli_03.arrow")))
filter!(x->coalesce(x.centerMz, -1.0)<600, ecoli_test)
filter!(x->coalesce(x.centerMz, -1.0)>500, ecoli_test)
filter!(x->coalesce(x.retentionTime, -1.0)<100, ecoli_test)
filter!(x->coalesce(x.retentionTime, -1.0)<75, ecoli_test)
 Arrow.write("/Users/n.t.wamsley/TEST_DATA/ECOLI_TEST/arrow_out/ecoli_filtered_03.arrow", ecoli_test)



=#