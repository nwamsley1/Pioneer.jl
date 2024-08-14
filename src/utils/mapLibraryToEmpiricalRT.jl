function mapLibraryToEmpiricalRT(psms_paths::Dictionary{String, String}; 
                        min_prob::AbstractFloat = 0.9
                    )

    #Dictionaries mapping fild_id names to data
    iRT_RT = Dictionary{String, Any}() #File name => KDEmapping from iRT to RT 
    RT_iRT = Dictionary{String, Any}() #File name => KDEmapping from iRT to RT 
    for (key, psms_path) in pairs(psms_paths) #For each data frame 
        psms = Arrow.Table(psms_path)

        best_hits = psms[:prob].>min_prob#Map RTs using only the best psms
        best_rts = psms[:RT][best_hits]
        best_irts = psms[:iRT_predicted][best_hits]
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
        #Build RT=>iRT and iRT=> RT mappings for the file and add to the dictionaries 
        insert!(iRT_RT, key, irt_to_rt_spline)
        insert!(RT_iRT, key, rt_to_irt_spline)
    end
    return iRT_RT, RT_iRT
end