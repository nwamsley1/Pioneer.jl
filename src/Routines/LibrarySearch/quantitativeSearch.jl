
function quantitationSearch(
    #Mandatory Args
    spectra::Arrow.Table,
    precursors::Vector{LibraryPrecursorIon{Float32}},
    ion_list::LibraryFragmentLookup{Float32},
    rt_index::retentionTimeIndex{Float32, Float32},
    ms_file_idx::UInt32,
    err_dist::MassErrorModel{Float32},
    irt_tol::Float64,
    params::NamedTuple,
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
        
        isotope_err_bounds = params[:isotope_err_bounds],
        filter_by_rank = Bool(params[:quant_search_params]["filter_by_rank"]),
        filter_by_count = Bool(params[:quant_search_params]["filter_on_frag_count"]),
        min_index_search_score = zero(UInt8),

        δ = Float32(params[:deconvolution_params]["huber_delta"]),
        λ = Float32(params[:deconvolution_params]["lambda"]),
        max_iter_newton = Int64(params[:deconvolution_params]["max_iter_newton"]),
        max_iter_bisection = Int64(params[:deconvolution_params]["max_iter_bisection"]),
        max_iter_outer = Int64(params[:deconvolution_params]["max_iter_outer"]),
        accuracy_newton = Float32(params[:deconvolution_params]["accuracy_newton"]),
        accuracy_bisection = Float32(params[:deconvolution_params]["accuracy_bisection"]),
        max_diff = Float32(params[:deconvolution_params]["max_diff"]),

        min_topn_of_m = Tuple([Int64(x) for x in params[:quant_search_params]["min_topn_of_m"]]),
        min_frag_count = Int64(params[:quant_search_params]["min_frag_count"]),
        min_log2_matched_ratio = Float32(params[:quant_search_params]["min_log2_matched_ratio"]),
        min_max_ppm = Tuple([Float32(x) for x in params[:frag_tol_params]["frag_tol_bounds"]]),
        n_frag_isotopes = Int64(params[:quant_search_params]["n_frag_isotopes"]),
        max_best_rank = Int64(params[:quant_search_params]["max_best_rank"]),

        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        rt_index = rt_index,
        irt_tol = irt_tol,
    )
end


BPSMS = Dict{Int64, DataFrame}()
#PSMS_DIR = joinpath(MS_DATA_DIR,"Search","RESULTS")
#PSM_PATHS = [joinpath(PSMS_DIR, file) for file in filter(file -> isfile(joinpath(PSMS_DIR, file)) && match(r".jld2$", file) != nothing, readdir(PSMS_DIR))];

features = [:intercept, :charge, :total_ions, :err_norm, 
:scribe, :city_block, :city_block_fitted, 
:spectral_contrast, :entropy_score, :weight]

quantitation_time = @timed for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    @time PSMS = vcat(quantitationSearch(MS_TABLE, 
                    prosit_lib["precursors"],
                    prosit_lib["f_det"],
                    RT_INDICES[MS_TABLE_PATH],
                    UInt32(ms_file_idx), 
                    frag_err_dist_dict[ms_file_idx],
                    irt_errs[ms_file_idx],
                    params_,  
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
        getIsoRanks!(PSMS, MS_TABLE, params_[:quadrupole_isolation_width]);
        PSMS[!,:prob] = zeros(Float32, size(PSMS, 1));
        scoreSecondSearchPSMs!(PSMS,features);
        MS2_CHROMS = groupby(PSMS, [:precursor_idx,:iso_rank]);
        integratePrecursors(MS2_CHROMS, 
                            n_quadrature_nodes = Int64(params_[:integration_params]["n_quadrature_nodes"]),
                            intensity_filter_fraction = Float32(params_[:integration_params]["intensity_filter_threshold"]),
                            α = 0.001f0);
        addPostIntegrationFeatures!(PSMS, 
                                    MS_TABLE, prosit_lib["precursors"],
                                    ms_file_idx,
                                    MS_TABLE_ID_TO_PATH,
                                    RT_iRT,
                                    precID_to_iRT
                                    );
        PSMS[!,:file_path].=file_id_to_parsed_name[ms_file_idx]
        BPSMS[ms_file_idx] = PSMS;
        GC.gc()
end

best_psms = vcat(values(BPSMS)...)