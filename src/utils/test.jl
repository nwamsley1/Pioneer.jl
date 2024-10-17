

function trim_top(data::AbstractVector{T}, proportion::Real) where {T<:Real}
    if !(0 <= proportion < 1)
        throw(ArgumentError("proportion must be between 0 and 1"))
    end
    
    n = length(data)
    k = floor(Int, n * (1 - proportion))  # number of elements to keep
    
    # Sort the data and take the first k elements
    return sort(data)[1:k]
end



testdf = getHuberLossParam(
    huber_δs,
    precursor_dict,
    MS_TABLE_PATHS,
    precursors[:is_decoy],
    frag_err_dist_dict,
    rt_index_paths,
    bin_rt_size,
    rt_irt,
    irt_errs,
    chromatograms,
    file_path_to_parsed_name,
    params_,
    spec_lib,
    ionMatches,
    ionMisses,
    IDtoCOL,
    ionTemplates,
    iso_splines,
    complex_scored_PSMs,
    complex_unscored_PSMs,
    complex_spectral_scores,
    precursor_weights);


sort!(testdf,:huber_δ)
gdf = groupby(testdf,[:precursor_idx,:scan_idx])
N = 100000
gdf = groupby(testdf, [:precursor_idx,:scan_idx])
x = log2.(gdf[N][!, :huber_δ])

# Create the first plot with left y-axis
p = plot(x, gdf[N][!, :tmean], 
         label="tmean", 
         xlabel="log2(huber_δ)", 
         ylabel="tmean",
         legend=:topleft)

# Create the second plot with right y-axis
plot!(twinx(), x, gdf[N][!, :mean_r2], 
      label="mean_r2", 
      ylabel="mean_r2",
      color=:red,
      legend=:topright)

# Display the plot
display(p)
N += 1

function testCombFunc(df::SubDataFrame)
    return (
        n = df.n[1],
        p = df.p[1],
        huber_δ = df.huber_δ[argmin(df.tmean)],
        iqr = df.iqr[end],
        tmean = df.tmean[end],
        mw = maximum(df.weight)
    )
end


gdf = groupby(testdf,[:scan_idx])
cgdf = combine(gdf, x->testCombFunc(x))


plot(log2.(cgdf[!,:huber_δ]), cgdf[!,:iqr], seriestype=:scatter)

plot(log2.(cgdf[!,:huber_δ]), cgdf[!,:mw], seriestype=:scatter, alpha = 0.1)


plot(log2.(cgdf[!,:huber_δ]), cgdf[!,:tmean], seriestype=:scatter, alpha = 0.1)



huber_δ=x.huber_δ[argmin(x.tmean)]




subt_psms = test_psms[!,[:b_count,:y_count,:total_ions,:log2_summed_intensity,:topn,:rt,:scan_idx,:precursor_idx,:score,:q_value,:spectrum_peak_count]]
ms_table = Arrow.Table(first(MS_TABLE_PATHS))

subt_psms[!,:spectrum_intensity] = [sum(ms_table[:intensity_array][scan_idx]) for scan_idx in subt_psms[!,:scan_idx]]
subt_psms[!,:matched_intensity] = exp2.(Float32.(subt_psms[!,:log2_summed_intensity]))

subt_psms[!,:matched_ions_fraction] = subt_psms[!,:total_ions]./subt_psms[!,:spectrum_peak_count]
subt_psms[!,:matched_intensity_fraction] = subt_psms[!,:matched_intensity]./subt_psms[!,:spectrum_intensity]

sort!(subt_psms, :matched_intensity_fraction, rev = true)


subt_psms[!,:sequence] = [precursors[:sequence][prec_id] for prec_id in subt_psms[!,:precursor_idx]]
subt_psms[!,:structural_mods] = [precursors[:structural_mods][prec_id] for prec_id in subt_psms[!,:precursor_idx]]
subt_psms[!,:prec_charge] = [precursors[:prec_charge][prec_id] for prec_id in subt_psms[!,:precursor_idx]]
subt_psms[!,:prec_mz] = [precursors[:mz][prec_id] for prec_id in subt_psms[!,:precursor_idx]]
subt_psms[!,:center_mz] = [ms_table[:centerMz][scan_idx] for scan_idx in subt_psms[!,:scan_idx]]


subt_psms[!,:modified_sequence] = getModifiedSequence.(subt_psms[!,:sequence], "", subt_psms[!,:structural_mods])
CSV.write(
subt_psms[!,[:scan_idx,:rt,:center_mz,:sequence,:structural_mods,:prec_charge,:prec_mz,:b_count,:y_count,:total_ions,:matched_ions_fraction,:matched_intensity_fraction,:score,:q_value]]


#basedir = "/Volumes/Active/Backpack/libraries/astral/Altimeter101324_MixedSpecies_OlsenAstral_NoEntrapment_101324_zCorrected/"
#altimeter_json_fnames = readdir("/Volumes/Active/Backpack/libraries/astral/Altimeter101324_MixedSpecies_OlsenAstral_NoEntrapment_101324_zCorrected/")
basedir = "/Volumes/Active/Backpack/libraries/astral/Altimeter101324_MixedSpecies_OlsenAstral_NoEntrapment_101324_zCorrected/"
ordered_altimiter_json_paths = joinpath.(
    basedir, 
    sort(readdir(basedir), by = x->parse(Int64, String(split(split(x,'_')[end], '.')[1])))
)

test_out = JSON.parse(read(altimeter_jsons[5000], String))
