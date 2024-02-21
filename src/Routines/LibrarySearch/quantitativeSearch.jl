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
    println("frag_ppm_err $frag_ppm_err")
    
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

quantitation_time = @timed for (ms_file_idx, MS_TABLE_PATH) in collect(enumerate(MS_TABLE_PATHS))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    println("starting file $ms_file_idx")
    @time PSMS = vcat(quantitationSearch(MS_TABLE, 
                    prosit_lib["precursors"],
                    prosit_lib["f_det"],
                    RT_INDICES_many[MS_TABLE_PATH],
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
    @time begin
        addSecondSearchColumns!(PSMS, MS_TABLE, prosit_lib["precursors"], precID_to_cv_fold_many);
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
        println("size(PSMS) ", size(PSMS))
        GC.gc()
    end
end
best_psms = vcat(values(BPSMS)...)
passing_gdf = groupby(best_psms, [:precursor_idx])
sort!(passing_gdf[(precursor_idx = 0x005102bb,)][!,[:prob,:sequence,:b_count,:y_count,:isotope_count,:weight,:H,:matched_ratio,:entropy_score,:scribe,:city_block_fitted,:file_path]],:file_path)

BPSMS = nothing
GC.gc()
println("Finished quant search in ", quantitation_time, "seconds")
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "BPSMS_K_022124.jld2"); BPSMS)

best_psms_passing = best_psms[best_psms[!,:q_value].<=0.01,:]

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