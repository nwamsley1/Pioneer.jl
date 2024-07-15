@time begin
    [normalizeQuant(
        best_psms,
        col_name,
        N = params_[:normalization_params]["n_rt_bins"],
        spline_n_knots = params_[:normalization_params]["spline_n_knots"],
        max_q_value = params_[:normalization_params]["max_q_value"],
        min_points_above_FWHM = params_[:normalization_params]["min_points_above_FWHM"]
    ) for col_name in [:peak_area]]
end

