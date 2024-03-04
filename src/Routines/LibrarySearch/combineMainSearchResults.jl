combinetime = @time begin
iRT_RT, RT_iRT = mapRTandiRT(PSMs_Dict)
precID_to_iRT = getPrecIDtoiRT(PSMs_Dict, RT_iRT, max_precursors = 250000)
RT_INDICES = makeRTIndices(PSMs_Dict,precID_to_iRT,iRT_RT)
precID_to_cv_fold = getCVFolds(precID_to_iRT)
end
println("time ", combinetime.time)


rt_align_mat = zeros((length(RT_INDICES), 100))

for (i, rt) in enumerate(LinRange(0, 100, 100))
    j = 1
    for (key, rtmap) in pairs(RT_iRT)
        rt_align_mat[j, i] = rtmap(rt)
        j += 1
    end
end