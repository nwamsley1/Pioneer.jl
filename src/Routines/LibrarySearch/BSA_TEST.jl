BSA_TABLE = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/PIONEER_BSA_TEST_062024/BSA_PRM_PIONEER_ISOTOPE_CORRECTION_TEST_CORRECTED_062024_uPAC_50cm.arrow")
plot(BSA_TABLE[:retentionTime], BSA_TABLE[:TIC], xlim = (9.1,9.2))



hcat(BSA_TABLE[:centerMass][3640:3660], BSA_TABLE[:TIC][3640:3660])

"LVNELTEFAK"
best_scan_idx = 3650

#Get Transition List
tlist = DataFrame(CSV.File("/Users/n.t.wamsley/TEST_DATA/PIONEER_BSA_TEST_062024/TRANSITION_LIST.csv"))
seq_to_id = Dict("TCVADESHAGCEK" => 1,
"YICDNQDTISSK" => 2,
"LVNELTEFAK" => 3)
tlist = tlist[tlist[!,"Fragment Ion Ordinal"].>2,:]
bsa_precursors = DataFrame(Dict(:mz => Float32[488.5345, 722.3247, 582.3190],
:prec_charge => UInt8[3, 2, 2],
:sulfur_count => UInt8[2,1,0]))
io = IOBuffer()
Arrow.write(io, bsa_precursors)
seekstart(io)
bsa_precursors = Arrow.Table(io)

bsa_fragments = Vector{DetailedFrag{Float32}}(undef, size(tlist, 1))
for i in range(1, size(bsa_fragments, 1))
    sequence = tlist[i,"Peptide Sequence"]
    ion_position = tlist[i,"Fragment Ion Ordinal"]
    sulfur_count = 0
    if tlist[i,"Fragment Ion Type"]=="b"
        sulfur_cont = UInt8(count("C", sequence[1:ion_position]) +  count("M", sequence[1:ion_position]))
    else
        sequence = reverse(sequence)
        sulfur_cont = UInt8(count("C", sequence[1:ion_position]) +  count("M", sequence[1:ion_position]))
    end

    bsa_fragments[i] = DetailedFrag(
        UInt32(seq_to_id[tlist[i,"Peptide Sequence"]]),
        Float32(tlist[i,:"Product Mz"]),
        one(Float16),

        tlist[i,"Fragment Ion Type"]=="b" ? UInt8(1) : UInt8(2),
        false,

        UInt8(tlist[i,"Product Charge"]),
        UInt8(tlist[i,"Fragment Ion Ordinal"]),
        UInt8(count('+', tlist[i,"Precursor"])),
        zero(UInt8),
        UInt8(sulfur_count)
    )
end

LVNELTEFAK_FRAGS = sort([frag for frag in bsa_fragments if frag.prec_id == 3], by = x->x.mz)
best_scan_idx = 3651
nmatches, nmisses = matchPeaks!(ionMatches[1], 
                                ionMisses[1], 
                                LVNELTEFAK_FRAGS, 
                                length(LVNELTEFAK_FRAGS), 
                                BSA_TABLE[:masses][best_scan_idx], 
                                BSA_TABLE[:intensities][best_scan_idx], 
                                MassErrorModel{Float32}(0.0f0, (20.0f0, 20.0f0)),
                                BSA_TABLE[:highMass][best_scan_idx],
                                UInt32(best_scan_idx), 
                                UInt32(1))

matched_frags = ionMatches[1][1:nmatches]
sort!(matched_frags, by=x->x.intensity, rev = true)
norm_factor = maximum([x.intensity for x in matched_frags])
bsa_fragments = Vector{DetailedFrag{Float32}}()
prec_charge = UInt8(2)
sequence = "LVNELTEFAK"
for (i, frag) in enumerate(matched_frags)
    ion_position = frag.frag_index
    sulfur_count = 0
    if frag.ion_type==one(UInt8)#b_ions
        sulfur_cont = UInt8(count("C", sequence[1:ion_position]) +  count("M", sequence[1:ion_position]))
    else
        sequence = reverse(sequence)
        sulfur_cont = UInt8(count("C", sequence[1:ion_position]) +  count("M", sequence[1:ion_position]))
    end
    push!(bsa_fragments,
    DetailedFrag(
        frag.prec_id,
        frag.theoretical_mz,
        Float16(frag.intensity/norm_factor),
        frag.ion_type,
        false,
        frag.frag_charge,
        frag.frag_index,
        prec_charge,
        UInt8(i),    
        UInt8(sulfur_count)
    ))   
end              
bsa_fragment_lookup_table = LibraryFragmentLookup(
    bsa_fragments,
    [range(UInt32(0), UInt32(0)),
   range(UInt32(0), UInt32(0)),
    range(UInt32(1), UInt32(length(bsa_fragments)))]
)

N = 200
t = collect(LinRange(0.0, 4*π, N))
u = sin.(t) 
u .+= randn(N)./50
plot(t, u, seriestype=:scatter)
test_spline = UniformSpline(u, t, 3, 20)
BSA_TABLE = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/PIONEER_BSA_TEST_062024/BSA_PRM_PIONEER_ISOTOPE_CORRECTION_TEST_CORRECTED_062124_uPAC_50cm.arrow")

params_[:deconvolution_params]["huber_delta"] = Float32(1e9)
psms = vcat(secondSearch(
    BSA_TABLE, 
    params_;
    precursors=bsa_precursors,
    fragment_lookup_table = bsa_fragment_lookup_table,
    rt_index = retentionTimeIndex([rtIndexBin(0.0f0, 0.0f0, [(UInt32(3),  582.3190f0)])]),#RT_INDICES[file_id_to_parsed_name[ms_file_idx]],
    ms_file_idx = UInt32(ms_file_idx), 
    rt_to_irt_spline = test_spline,#RT_iRT[file_id_to_parsed_name[ms_file_idx]],
    mass_err_model = MassErrorModel{Float32}(0.0f0, (25.0f0, 25.0f0)),
    irt_err = Inf,
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
    )...);

psms[!,:rt] = [BSA_TABLE[:retentionTime][x] for x in psms[!,:scan_idx]]
sort!(psms, :rt)
getIsotopesCaptured!(psms, bsa_precursors[:prec_charge],bsa_precursors[:mz], BSA_TABLE)
psms[psms[!,:scan_idx].==3651,[:weight,:isotopes_captured]]
gpsms = groupby(psms, :isotopes_captured)
p = plot(layout = (1, 2))
for (i, iso_set) in enumerate([(0, 5), (0, 0), (1, 5), (2, 5)])
    psms_ = gpsms[(isotopes_captured = iso_set,)]
    plot!(psms_[!,:rt], psms_[!,:weight], seriestype=:scatter, show = true, xlim = (9.1, 9.3), color = i, label = string(iso_set), subplot = 1)
end
correctPrecursorAbundances!(psms[!,:weight],
iso_splines,
psms[!,:isotopes_captured],
psms[!,:precursor_idx],
bsa_precursors[:mz],
bsa_precursors[:prec_charge],
bsa_precursors[:sulfur_count])
for (i, iso_set) in enumerate([(0, 5), (0, 0), (1, 5), (2, 5)])
    psms_ = gpsms[(isotopes_captured = iso_set,)]
    plot!(psms_[!,:rt], psms_[!,:weight], seriestype=:scatter, show = true, xlim = (9.1, 9.3),color = i, label = string(iso_set), subplot = 2)
end



p = plot()
for (iso, psms) in pairs(gpsms)
    plot!(psms[!,:rt], psms[!,:entropy_score], seriestype=:scatter, show = true, xlim = (9.1, 9.3), subplot_idx = 2)
end

correctPrecursorAbundances!(psms[!,:weight],
iso_splines,
psms[!,:isotopes_captured],
psms[!,:precursor_idx],
bsa_precursors[:mz],
bsa_precursors[:prec_charge],
bsa_precursors[:sulfur_count])

@time RESULT = plotPSM(
    BSA_TABLE, 
    3650,
    UInt32(3),
    params_;
    precursors = bsa_precursors,
    fragment_lookup_table = bsa_fragment_lookup_table,
    rt_index =  retentionTimeIndex([rtIndexBin(0.0f0, 0.0f0, [(UInt32(3),  582.3190f0)])]),
    ms_file_idx = UInt32(ms_file_idx), 
    rt_to_irt_spline = test_spline,
    mass_err_model = MassErrorModel{Float32}(0.0f0, (25.0f0, 25.0f0)),
    irt_err = Inf,#irt_errs[ms_file_idx]/3,
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