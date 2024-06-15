#=
uncorrected_chroms = copy(chroms)
uncorrected_gchroms = groupby(uncorrected_chroms, [:precursor_idx])
=#
subdf = uncorrected_gchroms[(precursor_idx = precs_passing[N],)]
gsubdf = groupby(subdf,:isotopes_captured)
p = plot(layout = (1, 2))
best_rt = subdf[argmax(subdf[!,:intensity]),:rt]
for (key, chrom) in pairs(gsubdf)
    plot_range = (chrom[!,:rt].>(best_rt - 0.1)).&(chrom[!,:rt].<(best_rt + 0.1))
    plot!(chrom[plot_range,:rt], 
            chrom[plot_range,:intensity], 
            subplot = 1,
            title = "Uncorrected \n"*precursors[:sequence][precs_passing[N]],
            seriestype=:scatter, 
            label = key[:isotopes_captured], 
            show = true)
end

gchroms = groupby(chroms, [:precursor_idx])
subdf = gchroms[(precursor_idx = precs_passing[N],)]
gsubdf = groupby(subdf,:isotopes_captured)
#plot()
for (key, chrom) in pairs(gsubdf)
    plot_range = (chrom[!,:rt].>(best_rt - 0.1)).&(chrom[!,:rt].<(best_rt + 0.1))
    plot!(chrom[plot_range,:rt], 
            chrom[plot_range,:intensity],             
            subplot = 2, title = "Corrected \n"*precursors[:sequence][precs_passing[N]], seriestype=:scatter, 
            label = key[:isotopes_captured], show = true)
end



N += 1


filter!(x->first(x.isotopes_captured)<2, chroms)
filter!(x->first(x.isotopes_captured)>-1, chroms)
correctPrecursorAbundances!(chroms[!,:intensity],
                            iso_splines,
                            chroms[!,:isotopes_captured],
                            chroms[!,:precursor_idx],
                            precursors[:mz],
                            precursors[:prec_charge],
                            precursors[:sulfur_count])

correctPrecursorAbundance(100.0f0,
iso_splines,
(0, 0),
precursors[:mz][precs_passing[N]]*precursors[:prec_charge][precs_passing[N]],
precursors[:sulfur_count][precs_passing[N]])


correctPrecursorAbundance(100.0f0,
iso_splines,
(1, 4),
precursors[:mz][precs_passing[N]]*precursors[:prec_charge][precs_passing[N]],
precursors[:sulfur_count][precs_passing[N]])



to_keep = (chroms_coldstart[!,:intensity].>1e-6).&(chroms[!,:intensity].>1e-6)

plot(log2.(chroms_coldstart[to_keep,:intensity])[1:1000:end],
log2.(chroms[to_keep,:intensity])[1:1000:end],
seriestype=:scatter,
alpha = 0.1)



@time RESULT = quantitationSearch(
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
    );

Hs, id_to_col, _weights_ = RESULT[end - 1]
H = sparse(Hs.rowval[1:Hs.n_vals], Hs.colval[1:Hs.n_vals], Hs.nzval[1:Hs.n_vals])
Hx = sparse(Hs.rowval[1:Hs.n_vals], Hs.colval[1:Hs.n_vals], Hs.x[1:Hs.n_vals])

H*_weights_[1:Hs.n]

(466.2551f0 - 466.2399)/(466.2551f0/1e6)

(466.75677143554685- 466.7409)/(466.75677143554685/1e6)

(   466.2369- 466.2399)/( 466.2369/1e6)
( 466.2402- 466.2399)/( 466.2402/1e6)
(  466.2429- 466.2399)/( 466.242/1e6)

H = sparse(Hs.rowval[1:Hs.n_vals], Hs.colval[1:Hs.n_vals], Hs.nzval[1:Hs.n_vals])

sparse(Hs.rowval[1:Hs.n_vals], Hs.colval[1:Hs.n_vals], Hs.x[1:Hs.n_vals])

w = []

H*_weights_[1:Hs.n]

H[:, id_to_col[2930911]]
Int64(id_to_col[2930911])

[cor(H[1:2600, id_to_col[2930911]], H[1:2600,i]) for i in range(1, Hs.n)]

H[1:2600,50]
scan_precs ∩ precursors_passing
precursors[:sequence][argmax(id_to_col.vals.==50)]

(precursors[:sequence][argmax(id_to_col.vals.==50)], precursors[:is_decoy][argmax(id_to_col.vals.==50)],  _weights_[50])
(precursors[:sequence][argmax(id_to_col.vals.==50)], precursors[:is_decoy][argmax(id_to_col.vals.==40)],  _weights_[49])
(precursors[:sequence][argmax(id_to_col.vals.==50)], precursors[:is_decoy][argmax(id_to_col.vals.==50)],  _weights_[50])
(precursors[:sequence][argmax(id_to_col.vals.==50)], precursors[:is_decoy][argmax(id_to_col.vals.==50)],  _weights_[50])
(precursors[:sequence][argmax(id_to_col.vals.==50)], precursors[:is_decoy][argmax(id_to_col.vals.==50)],  _weights_[50])
(precursors[:sequence][argmax(id_to_col.vals.==50)], precursors[:is_decoy][argmax(id_to_col.vals.==50)],  _weights_[50])

(precursors[:sequence][argmax(id_to_col.vals.==49)], precursors[:is_decoy][argmax(id_to_col.vals.==49)])
(precursors[:sequence][argmax(id_to_col.vals.==52)], precursors[:is_decoy][argmax(id_to_col.vals.==52)])
(precursors[:sequence][argmax(id_to_col.vals.==47)], precursors[:is_decoy][argmax(id_to_col.vals.==47)])
(precursors[:sequence][argmax(id_to_col.vals.==46)], precursors[:is_decoy][argmax(id_to_col.vals.==46)])
(precursors[:sequence][argmax(id_to_col.vals.==45)], precursors[:is_decoy][argmax(id_to_col.vals.==45)])
(precursors[:sequence][argmax(id_to_col.vals.==44)], precursors[:is_decoy][argmax(id_to_col.vals.==44)])

library_fragment_lookup_table.frags[library_fragment_lookup_table.prec_frag_ranges[argmax(id_to_col.vals.==49)]]
library_fragment_lookup_table.frags[library_fragment_lookup_table.prec_frag_ranges[argmax(id_to_col.vals.==50)]]

precursors[:mz]


rts = []
masses = []
intensities = []

sub_ = (abs.(MS_TABLE[:centerMass].-467.46237f0).<1e-6).&(abs.(MS_TABLE[:retentionTime].-19.794031f0).<2.5)

mass_lists = [MS_TABLE[:masses][x] for x in range(1, length(MS_TABLE[:masses])) if coalesce(sub_[x], false)==true]
intensity_lists = [MS_TABLE[:intensities][x] for x in range(1, length(MS_TABLE[:masses])) if coalesce(sub_[x], false)==true]

list_rts = [MS_TABLE[:retentionTime][x] for x in range(1, length(MS_TABLE[:masses])) if coalesce(sub_[x], false)==true]
for i in range(1, length(list_rts))
    for j in range(1, length(mass_lists[i]))
        if (mass_lists[i][j] > 740) & (mass_lists[i][j] < 760)
        push!(masses, mass_lists[i][j])
        push!(intensities, intensity_lists[i][j])
        push!(rts, list_rts[i])
        end
    end
end

plot(masses, rts, marker_z=log2.(intensities), seriestype=:scatter, alpha = 1.0, xlim = (746.25, 746.5))



PlotlyJS.plot(PlotlyJS.surface(z = Float32.(log2.(intensities)), x = Float32.(masses), y = Float32.(rts)))
