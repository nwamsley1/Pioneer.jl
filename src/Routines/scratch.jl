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

N = length(MS_TABLE_PATHS)
MS_TABLE_PATH = MS_TABLE_PATHS[N]
MS_TABLE = Arrow.Table(MS_TABLE_PATH)
ms_file_idx = N

@time PSMS = vcat(quantitationSearch(MS_TABLE, 
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
test = PSMS[PSMS[!,:precursor_idx].==5309115,[:precursor_idx,:scan_idx,:y_count,:b_count,:isotope_count,
:weight,:matched_ratio,:entropy_score,:scribe,:city_block_fitted,:RT]]

plot(test[!,:RT], test[!,:weight], seriestype=:scatter)

test_df = best_psms[best_psms[!,:sequence].=="EAGVFVPR",[:precursor_idx,:scan_idx,:sequence,:charge,:y_count,:b_count,:isotope_count,
:weight,:H,:matched_ratio,:entropy_score,:scribe,:city_block_fitted]]
precursors[test_df[1,:precursor_idx]]