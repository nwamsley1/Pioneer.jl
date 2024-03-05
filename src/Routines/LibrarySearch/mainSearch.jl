
function mainLibrarySearch(
    #Mandatory Args
    spectra::Arrow.Table, 
    frag_index::FragmentIndex{Float32},
    precursors::Vector{LibraryPrecursorIon{Float32}},
    ion_list::LibraryFragmentLookup{Float32},
    iRT_to_RT_spline::Any,
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
    precursor_weights::Vector{Vector{Float32}},
    precs::Vector{Counter{UInt32, UInt8}}) where {S<:ScoredPSM{Float32, Float16},
                                                            Q<:UnscoredPSM{Float32},
                                                            R<:SpectralScores{Float16}}#where {S<:ScoredPSM{Float32, Float16}, LibraryIon{Float32}}

    frag_ppm_err = getLocation(err_dist)
    
    return searchRAW(
        spectra, 
        frag_index, 
        precursors,
        ion_list,
        iRT_to_RT_spline,
        ms_file_idx,
        err_dist,
        searchScan!,

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
        precs,

        expected_matches = params[:expected_matches],
        frag_ppm_err = frag_ppm_err,
        isotope_err_bounds = params[:isotope_err_bounds],
        min_frag_count = params[:min_frag_count],
        min_log2_matched_ratio = params[:min_log2_matched_ratio_main_search],
        min_index_search_score = params[:min_index_search_score],
        min_topn_of_m = params[:min_topn_of_m],
        min_max_ppm = (10.0f0, 30.0f0),
        n_frag_isotopes = params[:n_frag_isotopes],
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        rt_bounds = params[:rt_bounds],
        irt_tol = irt_tol,
    )
end

PSMs_Dict = Dictionary{String, DataFrame}()
main_search_time = @timed for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
        println("starting file $ms_file_idx")
        @time begin
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)  
        @time PSMs = vcat(mainLibrarySearch(
            MS_TABLE,
            prosit_lib["f_index"],
            prosit_lib["precursors"],
            #prosit_lib["f_det"],
            library_fragment_lookup_table,
            RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
            UInt32(ms_file_idx), #MS_FILE_IDX
            frag_err_dist_dict[ms_file_idx],
            irt_errs[ms_file_idx],
            main_search_params,
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
            precs
        #scan_range = (100000, 100010)
        )...);

        addMainSearchColumns!(PSMs, MS_TABLE, prosit_lib["precursors"]);
        getRTErrs!(PSMs);

        column_names = [:spectral_contrast,:scribe,:city_block,:entropy_score,
                        :iRT_error,:missed_cleavage,:Mox,:charge,:TIC,
                        :y_count,:err_norm,:spectrum_peak_count,:intercept]

        scoreMainSearchPSMs!(PSMs,
                                    column_names,
                                    n_train_rounds = 2,
                                    max_iter_per_round = 20,
                                    max_q_value = 0.01);

        getProbs!(PSMs);
                                    
        getBestPSMs!(PSMs,
                        prosit_lib["precursors"],
                        max_q_value = 0.25,
                        max_psms = 250000);

        println("retained ", size(PSMs, 1), " psms")

        insert!(PSMs_Dict, 
            MS_TABLE_PATH, 
            PSMs#[!,
                #[:precursor_idx,:RT,:iRT_predicted,:prec_mz,:q_value,:score]
                #]
        );
        end
end


println("Finished main search in ", main_search_time.time, "seconds")
println("Finished main search in ", main_search_time, "seconds")


jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_max250000_q25_minfrag3ynoiso_030424_M0.jld2"); PSMs_Dict)

PSMs_Dict = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_max250000_q25_030424_M0.jld2"))["PSMs_Dict"]
#jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_030424_M0.jld2"); PSMs_Dict)

jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_max250000_q25_030424_M0.jld2"); PSMs_Dict)
PSMs_Dict = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_max250000_q25_030424_M0.jld2"))["PSMs_Dict"]
#jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_030424_M0.jld2"); PSMs_Dict)
#PSMs_Dict = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_030224_M0.jld2"))["PSMs_Dict"]


jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_030224_M0.jld2"); PSMs_Dict)
PSMs_Dict = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_030224_M0.jld2"))["PSMs_Dict"]
#=
jldsave(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_021924_M0.jld2"); PSMs_Dict)
PSMs_Dict = load(joinpath(MS_DATA_DIR, "Search", "RESULTS", "PSMs_Dict_021924_M0.jld2"))["PSMs_Dict"]
=#