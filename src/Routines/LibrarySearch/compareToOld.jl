tpold = load("/Users/n.t.wamsley/Desktop/asms_comparison_070924/three_proteome_filtered_asms.jld2")["three_proteome_results"]
tpold_set = Set(tpold[!,:precursor_idx])

length(setdiff(tpold_set, precursors_passing))
#=
julia> length(setdiff(tpold_set, precursors_passing))
2484
julia> length(tpold_set)
110556
=#
#Most old precursors passed scoring. So some difficulty with precision/quant. 
length(setdiff(tpold_set, Set(best_psms_passing[!,:precursor_idx])))


missing_precs = setdiff(Set(tpold[!,:modified_sequence]), Set(three_proteome_results[!,:modified_sequence]))

sort!(three_proteome_results[[x∈missing_precs for x in three_proteome_results[!,:modified_sequence]],:],:modified_sequence)
sort!(tpold[[x∈missing_precs for x in tpold[!,:modified_sequence]],:],:modified_sequence)
precs_to_test = Set(tpold[[x∈missing_precs for x in tpold[!,:modified_sequence]],:precursor_idx])
precs_to_test_list = collect(precs_to_test)

N = 1000


bpsms_old = load("/Users/n.t.wamsley/Desktop/asms_comparison_070924/best_psms_passing_final_asms_070924.jld2")["best_psms_passing"]


sort!(best_psms_passing[best_psms_passing[!,:precursor_idx].==precs_to_test_list[N],
[:file_name,:precursor_idx,:best_rank,:topn,:b_count,:y_count,:q_value,:weight,:peak_area,:scan_idx,:RT,:isotopes_captured]],
:file_name)
sort!(bpsms_old[bpsms_old[!,:precursor_idx].==precs_to_test_list[N],
[:file_name,:precursor_idx,:best_rank,:topn,:b_count,:y_count,:q_value,:weight,:peak_area,:scan_idx,:RT]],
:file_name)
N += 1
#sum(isnan.(gchroms[N][!,:intensity]))
#N += 1

TEST = vcat(getChromatograms(
            MS_TABLE, 
            params_;
            precursors = prosit_lib["precursors"],
            fragment_lookup_table = library_fragment_lookup_table,
            rt_index = RT_INDICES[file_id_to_parsed_name[ms_file_idx]],
            ms_file_idx = UInt32(ms_file_idx), 
            rt_to_irt_spline = RT_iRT[file_id_to_parsed_name[ms_file_idx]],
            mass_err_model = frag_err_dist_dict[ms_file_idx],
            irt_err = irt_err,#irt_errs[ms_file_idx]/3,
            ion_matches = ionMatches,
            ion_misses = ionMisses,
            id_to_col = IDtoCOL,
            ion_templates = ionTemplates,
            iso_splines = iso_splines,
            chromatograms = chromatograms,
            scored_psms = complex_scored_PSMs,
            unscored_psms = complex_unscored_PSMs,
            spectral_scores = complex_spectral_scores,
            precursor_weights = precursor_weights,
            precursors_passing = precursors_passing,
            quad_transmission_func = QuadTransmission(1.0f0, 1000.0f0)
            )...);

Hs, _residuals_, _weights_ = TEST[1]
solveHuber!(Hs, _residuals_, _weights_, 
Float32(params_[:deconvolution_params]["huber_delta"]), 
0.0f0, 
params_[:deconvolution_params]["max_iter_newton"], 
params_[:deconvolution_params]["max_iter_bisection"],
params_[:deconvolution_params]["max_iter_outer"],
Float32(params_[:deconvolution_params]["accuracy_newton"]),
Float32(params_[:deconvolution_params]["accuracy_bisection"]),
10.0,#Hs.n/10.0,
Float32(params_[:deconvolution_params]["max_diff"])
);


frag_err_dist_dict = Dict(
5 => MassErrorModel{Float32}(5.24532, (17.1673, 5.53177)),
  4 => MassErrorModel{Float32}(5.23162, (17.3318, 5.88154)),
  6 => MassErrorModel{Float32}(5.19919, (15.9143, 4.94332)),
  2 => MassErrorModel{Float32}(3.48961, (14.2485, 6.18278)),
  3 => MassErrorModel{Float32}(3.8312, (15.0107, 6.20125)),
  1 => MassErrorModel{Float32}(3.33906, (15.6146, 6.41751)))

N = 100000

nkey = keys(gchroms)[N]
plot(gchroms[nkey][!,:rt], gchroms[nkey][!,:intensity])
N += 1



for (key, psms) in pairs(gchroms)
    sort!(psms, :rt)
end

nkey = (precursor_idx = 724316, isotopes_captured = (0, 1))
plot(gchroms[nkey][!,:rt], gchroms[nkey][!,:intensity])


scan_idx = 152308 

@time chroms = vcat(getChromatograms(
    MS_TABLE, 
    params_;
    precursors = prosit_lib["precursors"],
    fragment_lookup_table = library_fragment_lookup_table,
    rt_index = RT_INDICES[file_id_to_parsed_name[ms_file_idx]],
    ms_file_idx = UInt32(ms_file_idx), 
    rt_to_irt_spline = RT_iRT[file_id_to_parsed_name[ms_file_idx]],
    mass_err_model = frag_err_dist_dict[ms_file_idx],
    irt_err = irt_err,#irt_errs[ms_file_idx]/3,
    ion_matches = ionMatches,
    ion_misses = ionMisses,
    id_to_col = IDtoCOL,
    ion_templates = ionTemplates,
    iso_splines = iso_splines,
    chromatograms = chromatograms,
    scored_psms = complex_scored_PSMs,
    unscored_psms = complex_unscored_PSMs,
    spectral_scores = complex_spectral_scores,
    precursor_weights = precursor_weights,
    precursors_passing = precursors_passing,
    quad_transmission_func = QuadTransmission(1.0f0, 1000.0f0)
    )...);

MATCHES, MISSES, TEMPLATES = chroms[13]
[x for x in MATCHES if x.prec_id == 724316]
[x for x in MISSES if x.prec_id == 724316]
[x for x in TEMPLATES if x.prec_id == 724316]
scan_idx = 152308 
spectra = MS_TABLE 
matchPeaks!(MATCHES, 
                                        MISSES, 
                                        TEMPLATES, 
                                        length(TEMPLATES), 
                                        spectra[:masses][scan_idx], 
                                        spectra[:intensities][scan_idx], 
                                        MassErrorModel{Float32}(3.33906, (15.6146, 6.41751)),#frag_err_dist_dict[1],
                                        spectra[:highMass][scan_idx],
                                        UInt32(scan_idx), 
                                        UInt32(ms_file_idx))
