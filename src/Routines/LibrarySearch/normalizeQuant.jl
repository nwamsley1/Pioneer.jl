
@time begin
    [normalizeQuant(
        best_psms,
        col_name,
        N = params_[:normalization_params]["n_rt_bins"],
        spline_n_knots = params_[:normalization_params]["spline_n_knots"],
        max_q_value = params_[:normalization_params]["max_q_value"],
        min_points_above_FWHM = params_[:normalization_params]["min_points_above_FWHM"]
    ) for col_name in [:trapezoid_area,:peak_area,:weight]]
end

#=
gpsms = groupby(best_psms[(best_psms[:,:q_value].<=0.01) .& (best_psms[:,:decoy].==false),:],[:precursor_idx,:isotopes_captured])

CV = Vector{Union{Missing,  Float32}}(undef, length(gpsms))
for i in range(1, length(gpsms))
    precgroup = gpsms[i]
    if size(precgroup, 1) < 3
        CV[i] = missing
        continue
    end
    CV[i] = std(precgroup[!,:trapezoid_area_normalized])/mean(precgroup[!,:trapezoid_area_normalized])
end
CV = skipmissing(CV)
sum((CV.<=0.2).&(CV.>1e-6))

CV = Vector{Union{Missing,  Float32}}(undef, length(gpsms))
for i in range(1, length(gpsms))
    precgroup = gpsms[i]
    if size(precgroup, 1) < 3
        CV[i] = missing
        continue
    end
    CV[i] = std(precgroup[!,:trapezoid_area])/mean(precgroup[!,:trapezoid_area])
end
CV = skipmissing(CV)
sum((CV.<=0.2).&(CV.>1e-6))

CV = Vector{Union{Missing,  Float32}}(undef, length(gpsms))
for i in range(1, length(gpsms))
    precgroup = gpsms[i]
    if size(precgroup, 1) < 3
        CV[i] = missing
        continue
    end
    CV[i] = std(precgroup[!,:weight])/mean(precgroup[!,:weight])
end
CV = skipmissing(CV)
sum((CV.<=0.2).&(CV.>1e-6))

CV = Vector{Union{Missing,  Float32}}(undef, length(gpsms))
for i in range(1, length(gpsms))
    precgroup = gpsms[i]
    if size(precgroup, 1) < 3
        CV[i] = missing
        continue
    end
    CV[i] = std(precgroup[!,:peak_area])/mean(precgroup[!,:peak_area])
end
CV = skipmissing(CV)
sum((CV.<=0.2).&(CV.>1e-6))
=#

