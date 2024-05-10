best_psms[(best_psms[!,:q_value].<=0.01).&(best_psms[!,:file_name].=="01"),:]

best_precursors = unique(best_psms[(best_psms[!,:q_value].<=0.01).&(best_psms[!,:file_name].=="01"),[:precursor_idx,:iso_rank]])
jldsave("/Users/n.t.wamsley/Desktop/best_precursors.jld2"; best_precursors)
N = 10000

PSMS[!,:iRT] = RT_iRT["01"].(PSMS[!,:RT])

PSMS[!,:center_mass] = MS_TABLE[:centerMass][PSMS[!,:scan_idx]]

dtype = Float32;
gx, gw = gausslegendre(100);
state = GD_state(
    HuberParams(zero(dtype), zero(dtype),zero(dtype),zero(dtype)), #Initial params
    zeros(dtype, N), #t
    zeros(dtype, N), #y
    zeros(dtype, N), #data
    falses(N), #mask
    0, #number of iterations
    N #max index
    );

integratePrecursorMS2(MS2_CHROMS[(precursor_idx = best_precursors[N,1],iso_rank = best_precursors[N,2])],
state,
gx::Vector{Float64},
gw::Vector{Float64},
intensity_filter_fraction =  Float32(params_[:integration_params]["intensity_filter_threshold"]),
α = 0.001f0,
half_width_at_α = 0.15f0,
isplot = true
);

MS2_CHROMS[(precursor_idx = best_precursors[N,1],iso_rank = best_precursors[N,2])][!,
[:best_rank,:topn,:b_count,:y_count,:isotope_count,:scribe,:scribe_fitted,
:city_block_fitted,:spectral_contrast,:matched_ratio,:weight,:RT,:center_mass,:scan_idx]]


N += 1

diff(MS2_CHROMS[(precursor_idx = best_precursors[N,1],iso_rank = best_precursors[N,2])][!,:scan_idx])