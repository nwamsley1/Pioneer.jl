
ms2_integration_params[:min_frag_count] = 0
ms2_integration_params[:min_topn_of_m] = (0,3)
MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    println("starting file $ms_file_idx")
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
        addIntegrationFeatures!(PSMS);
        getIsoRanks!(PSMS, MS_TABLE, ms2_integration_params[:quadrupole_isolation_width]);
        PSMS[!,:prob] = zeros(Float32, size(PSMS, 1));
        scoreSecondSearchPSMs!(PSMS,features);
        MS2_CHROMS = groupby(PSMS, [:precursor_idx,:iso_rank]);

        MS_TABLE = Arrow.Table(MS_TABLE_PATHS[end-1])
        @time PSMS2 = vcat(quantitationSearch(MS_TABLE, 
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
addSecondSearchColumns!(PSMS2, MS_TABLE, prosit_lib["precursors"], precID_to_cv_fold);
addIntegrationFeatures!(PSMS2);
getIsoRanks!(PSMS2, MS_TABLE, ms2_integration_params[:quadrupole_isolation_width]);
PSMS2[!,:prob] = zeros(Float32, size(PSMS2, 1));
scoreSecondSearchPSMs!(PSMS2,features);
MS2_CHROMS2 = groupby(PSMS2, [:precursor_idx,:iso_rank]);

MS_TABLE = Arrow.Table(MS_TABLE_PATHS[end-2])
@time PSMS3 = vcat(quantitationSearch(MS_TABLE, 
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
addSecondSearchColumns!(PSMS3, MS_TABLE, prosit_lib["precursors"], precID_to_cv_fold);
addIntegrationFeatures!(PSMS3);
getIsoRanks!(PSMS3, MS_TABLE, ms2_integration_params[:quadrupole_isolation_width]);
PSMS3[!,:prob] = zeros(Float32, size(PSMS3, 1));
scoreSecondSearchPSMs!(PSMS3,features);
MS2_CHROMS3 = groupby(PSMS3, [:precursor_idx,:iso_rank]);


#filter!(x->x.precursor_idx == prec_id,PSMS);
#PSMS[!,[:topn,:y_count,:entropy_score,:weight,:scan_idx]]
#test_chrom[!,[:topn,:y_count,:entropy_score,:RT,:weight,:scan_idx]]

    display_features = [:precursor_idx,:prob,:charge,:topn,:b_count,:y_count,:isotope_count,:matched_ratio,:entropy_score,:city_block_fitted,:RT,:weight,:scan_idx]

prec_id = 11948245 


prec_id =   11069235

prec_id =  11083883
iso_rank = 1
test_chrom = MS2_CHROMS[(precursor_idx = prec_id, iso_rank = iso_rank,)];#[!,[:precursor_idx,:topn,:b_count,:y_count,:isotope_count,:entropy_score,:city_block_fitted,:RT,:weight]]
test_chrom2 = MS2_CHROMS2[(precursor_idx = prec_id, iso_rank = iso_rank,)];#[!,[:precursor_idx,:topn,:b_count,:y_count,:isotope_count,:entropy_score,:city_block_fitted,:RT,:weight]]
test_chrom3 = MS2_CHROMS3[(precursor_idx = prec_id, iso_rank = iso_rank,)];#[!,[:precursor_idx,:topn,:b_count,:y_count,:isotope_count,:entropy_score,:city_block_fitted,:RT,:weight]]

#test_chrom[!,display_features]
plot(test_chrom[!,:RT].-0.1, test_chrom[!,:weight], alpha = 0.5, seriestype=:scatter)
plot!(test_chrom2[!,:RT], test_chrom2[!,:weight], alpha = 0.5, seriestype=:scatter)
plot!(test_chrom3[!,:RT], test_chrom3[!,:weight], alpha = 0.5, seriestype=:scatter)

test_chrom[!,display_features]
test_chrom2[!,display_features]
test_chrom3[!,display_features]

include("src/Routines/LibrarySearch/methods/integrateChroms.jl")
gx, gw = gausslegendre(100)
N = 500
dtype = eltype(test_chrom.weight)
state = GD_state(
    HuberParams(zero(dtype), zero(dtype),zero(dtype),zero(dtype)), #Initial params
    zeros(dtype, N), #t
    zeros(dtype, N), #y
    zeros(dtype, N), #data
    falses(N), #mask
    0, #number of iterations
    N #max index
    )
integratePrecursorMS2(test_chrom,
state,
gx::Vector{Float64},
gw::Vector{Float64},
intensity_filter_fraction = 0.01f0,
α = 0.01f0,
half_width_at_α = 0.15f0,
isplot = true
)

include("src/Routines/LibrarySearch/methods/integrateChroms.jl")
gx, gw = gausslegendre(100)
N = 500
dtype = eltype(test_chrom.weight)
state = GD_state(
    HuberParams(zero(dtype), zero(dtype),zero(dtype),zero(dtype)), #Initial params
    zeros(dtype, N), #t
    zeros(dtype, N), #y
    zeros(dtype, N), #data
    falses(N), #mask
    0, #number of iterations
    N #max index
    )
integratePrecursorMS2(test_chrom2,
state,
gx::Vector{Float64},
gw::Vector{Float64},
intensity_filter_fraction = 0.01f0,
α = 0.01f0,
half_width_at_α = 0.15f0,
isplot = true
)

include("src/Routines/LibrarySearch/methods/integrateChroms.jl")
gx, gw = gausslegendre(100)
N = 500
dtype = eltype(test_chrom.weight)
state = GD_state(
    HuberParams(zero(dtype), zero(dtype),zero(dtype),zero(dtype)), #Initial params
    zeros(dtype, N), #t
    zeros(dtype, N), #y
    zeros(dtype, N), #data
    falses(N), #mask
    0, #number of iterations
    N #max index
    )
integratePrecursorMS2(test_chrom3,
state,
gx::Vector{Float64},
gw::Vector{Float64},
intensity_filter_fraction = 0.01f0,
α = 0.01f0,
half_width_at_α = 0.15f0,
isplot = true
)