
function quantitationSearch(
    #Mandatory Args
    spectra::Arrow.Table,
    precursors::Vector{LibraryPrecursorIon{Float32}},
    ion_list::LibraryFragmentLookup{Float32},
    rt_index::retentionTimeIndex{Float32, Float32},
    ms_file_idx::UInt32,
    err_dist::MassErrorModel{Float32},
    irt_tol::Float64,
    params::Dict,
    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
    ionTemplates::Vector{Vector{DetailedFrag{Float32}}},
    iso_splines::IsotopeSplineModel{Float64},
    scored_PSMs::Vector{Vector{S}},
    unscored_PSMs::Vector{Vector{Q}},
    spectral_scores::Vector{Vector{R}},
    precursor_weights::Vector{Vector{Float32}}
    ) where {S<:ScoredPSM{Float32, Float16},
                                            Q<:UnscoredPSM{Float32},
                                            R<:SpectralScores{Float16}}
    frag_ppm_err = Float32(getLocation(err_dist))
    
    #fragment_tolerance = quantile(err_dist, params[:frag_tol_quantile])
    return searchRAW(
        spectra, 
        missing,
        precursors,
        ion_list, 
        x->x,
        ms_file_idx,
        err_dist,
        missing,
        
        ionMatches,
        ionMisses,
        all_fmatches,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        scored_PSMs,
        unscored_PSMs,
        spectral_scores,
        precursor_weights,
        missing,
        
        frag_ppm_err = Float32(frag_ppm_err),
        unmatched_penalty_factor = params[:unmatched_penalty_factor],
        isotope_err_bounds = params[:isotope_err_bounds],
        max_peak_width = params[:max_peak_width],
        min_topn_of_m = params[:min_topn_of_m],
        filter_by_rank = true,
        huber_δ = params[:huber_δ],
        min_frag_count = params[:min_frag_count],
        min_log2_matched_ratio = params[:min_log2_matched_ratio],
        min_index_search_score = zero(UInt8),#params[:min_index_search_score],
        min_weight = params[:min_weight],
        min_max_ppm = (15.0f0, 40.0f0),
        n_frag_isotopes = params[:n_frag_isotopes],
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        rt_index = rt_index,
        rt_bounds = params[:rt_bounds],
        irt_tol = irt_tol
    )
end


BPSMS = Dict{Int64, DataFrame}()
PSMS_DIR = joinpath(MS_DATA_DIR,"Search","RESULTS")
PSM_PATHS = [joinpath(PSMS_DIR, file) for file in filter(file -> isfile(joinpath(PSMS_DIR, file)) && match(r".jld2$", file) != nothing, readdir(PSMS_DIR))];

features = [:intercept, :charge, :total_ions, :err_norm, 
:scribe, :city_block, :city_block_fitted, 
:spectral_contrast, :entropy_score, :weight]

quantitation_time = @timed for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    println("starting file $ms_file_idx")
    PSMS = vcat(quantitationSearch(MS_TABLE, 
                    prosit_lib["precursors"],
                    prosit_lib["f_det"],
                    RT_INDICES[MS_TABLE_PATH],
                    UInt32(ms_file_idx), 
                    frag_err_dist_dict[ms_file_idx],
                    irt_errs[ms_file_idx],
                    ms2_integration_params,  
                    ionMatches,
                    ionMisses,
                    all_fmatches,
                    IDtoCOL,
                    ionTemplates,
                    iso_splines,
                    complex_scored_PSMs,
                    complex_unscored_PSMs,
                    complex_spectral_scores,
                    precursor_weights,
                    )...);
        addSecondSearchColumns!(PSMS, MS_TABLE, prosit_lib["precursors"], precID_to_cv_fold);
        addIntegrationFeatures!(PSMS);
        getIsoRanks!(PSMS, MS_TABLE, ms2_integration_params[:quadrupole_isolation_width]);
        PSMS[!,:prob] = zeros(Float32, size(PSMS, 1));
        scoreSecondSearchPSMs!(PSMS,features);
        MS2_CHROMS = groupby(PSMS, [:precursor_idx,:iso_rank]);
        integratePrecursors(MS2_CHROMS, 
                            n_quadrature_nodes = params_[:n_quadrature_nodes],
                            intensity_filter_fraction = params_[:intensity_filter_fraction],
                            α = 0.001f0);
        addPostIntegrationFeatures!(PSMS, 
                                    MS_TABLE, prosit_lib["precursors"],
                                    ms_file_idx,
                                    MS_TABLE_ID_TO_PATH,
                                    RT_iRT,
                                    precID_to_iRT
                                    );
        PSMS[!,:file_path].=MS_TABLE_PATH
        BPSMS[ms_file_idx] = PSMS;
        GC.gc()
end
best_psms = vcat(values(BPSMS)...)

#best_psms_old = copy(best_psms)
#filter!(x->x.y_count>2, best_psms)
features = [ 
    #:max_prob,
    :median_prob,
    :q90_prob,
    :iso_rank,
    :assymetry,
    :fraction_censored,
    :FWHM,
    :FWHM_01,
    :GOF,
    :H,
    :base_width_min,
    :peak_area,
    :points_above_FWHM,
    :points_above_FWHM_01,

    :missed_cleavage,
    :Mox,
    :prec_mz,
    :sequence_length,
    :charge,

    :irt_pred,
    :irt_error,
    :irt_obs,
    :RT,
    :irt_diff,

    :max_y_ions,
    :y_ions_sum,
    :longest_y,
    :y_count,
    :b_count,
    :isotope_count,
    :total_ions,
    :best_rank,
    :topn,

    :max_score,
    :mean_score,
    :max_city_fitted,
    :mean_city_fitted,
    :city_block,
    :city_block_fitted,
    :entropy_score,
    :max_entropy,
    :scribe,
    :scribe_corrected,
    :scribe_fitted,
    :spectral_contrast,
    :spectral_contrast_corrected,
    :max_matched_ratio,
    :max_scribe_score,
    
    :data_points,
    :err_norm,
    :error,
    :matched_ratio,
    :poisson,

    :weight,
    :log2_intensity_explained,
    :TIC,
    :adjusted_intensity_explained
];

#best_psms[!,features]
best_psms[!,:q_value] = zeros(Float32, size(best_psms, 1));
best_psms[!,:decoy] = best_psms[!,:target].==false;
xgboost_time = @timed bst = rankPSMs!(best_psms, 
                        features,
                        colsample_bytree = 0.5, 
                        colsample_bynode = 0.5,
                        min_child_weight = 5, 
                        gamma = 1, 
                        subsample = 0.25, 
                        max_depth = 10, 
                        eta = 0.05, 
                        #eta = 0.0175,
                        train_fraction = 9.0/9.0,
                        iter_scheme = [100, 100, 200],
                        print_importance = false);
best_psms = bst[2];
best_psms[!,:prob] =Float32.(best_psms[!,:prob]);
getQvalues!(best_psms[!,:prob], best_psms[:,:target], best_psms[!,:q_value]);
value_counts(df, col) = combine(groupby(df, col), nrow);
IDs_PER_FILE = value_counts(best_psms[(best_psms[:,:q_value].<=0.01) .& (best_psms[:,:decoy].==false),:], [:file_path])
#jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "best_psms_scored_M0M1_022124_K_alltrace.jld2"); best_psms)

length(unique(best_psms[(best_psms[:,:q_value].<=0.01) .& (best_psms[:,:decoy].==false),:precursor_idx]))
println("file specific ids ", sum(IDs_PER_FILE[!,:nrow]))
transform!(best_psms, AsTable(:) => ByRow(psm -> 
prosit_lib["precursors"][psm[:precursor_idx]].accession_numbers
) => :accession_numbers
);
getBestTrace!(best_psms)
IDs_PER_FILE = value_counts(best_psms[(best_psms[:,:q_value].<=0.01) .& (best_psms[:,:decoy].==false),:], [:file_path])

jldsave(joinpath(MS_DATA_DIR,"Search", "RESULTS", "best_psms_scored_HUPO_huber1e3_lasso1e1_030624.jld2"); best_psms)
println("TEST")



#best_psms = load("/Users/n.t.wamsley/TEST_DATA/HEIL_2023/TEST_y4b3_nOf5/Search/RESULTS/best_psms_scored_T14_022924_besttrace.jld2");
best_psms = best_psms["best_psms"]


passing_gdf = groupby(best_psms, [:precursor_idx])


ssing = best_psms[best_psms[!,:q_value].<=0.01,:]

best_psms_passing[!,:condition] = [split(x,"_")[end - 1] for x in best_psms_passing[!,:file_path]]
passing_gdf = groupby(best_psms_passing, [:precursor_idx,:condition])

cvs = zeros(length(passing_gdf))
i = 1
for (key, value) in ProgressBar(pairs(passing_gdf))
    cvs[i] = std(value[!,:weight])/mean(value[!,:weight])
    i += 1
end


cvs = zeros(length(passing_gdf))
i = 1
for (key, value) in ProgressBar(pairs(passing_gdf))
    cvs[i] = std(value[!,:peak_area])/mean(value[!,:peak_area])
    i += 1
end
cvs[cvs.<0.2]

best_psms[!,:SILAC] = [split(split(x,';')[1],'_')[end] == "SILAC" for x in best_psms[!,:accession_numbers]];

BADID = best_psms[(best_psms[!,:SILAC].==false).&(best_psms[!,:condition].=="A").&(best_psms[!,:q_value].<=0.01),:]

GOODID = best_psms[(best_psms[!,:SILAC].==false).&(best_psms[!,:condition].=="K").&(best_psms[!,:q_value].<=0.01),:]

BADID[!,:matched_ratio]