
no_missing = any.(eachrow(ismissing.(wide_psms_quant))).==false;
single_species = length.(unique.(split.(wide_psms_quant[!,:species], ';'))).==1;
no_mox= count.("Unimod:35", wide_psms_quant[!,:structural_mods]).==0;
wide_psms_quant_p = wide_psms_quant[single_species.&no_mox.&no_missing,:]

A = 100 .*map(x->(1 + 1/(4*3))*std(x)/mean(x), eachrow(wide_psms_quant_p[!,[8, 9, 10]]))
B = 100 .*map(x->(1 + 1/(4*3))*std(x)/mean(x), eachrow(wide_psms_quant_p[!,[11, 12, 13]]))

sum((A.<20).&(B.<20))