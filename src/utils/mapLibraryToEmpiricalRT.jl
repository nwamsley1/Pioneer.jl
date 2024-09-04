function mapLibraryToEmpiricalRT(psms_paths::Dictionary{String, String}; 
                        min_prob::AbstractFloat = 0.9
                    )

    #Dictionaries mapping fild_id names to data
    irt_rt = Dictionary{String, Any}() #File name => KDEmapping from irt to rt 
    rt_irt = Dictionary{String, Any}() #File name => KDEmapping from irt to rt 
    for (key, psms_path) in pairs(psms_paths) #For each data frame 
        psms = Arrow.Table(psms_path)

        best_hits = psms[:prob].>min_prob#Map rts using only the best psms
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
    end
    return irt_rt, rt_irt
end