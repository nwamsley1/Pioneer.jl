#BSA_TABLE = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/PIONEER_BSA_TEST_062024/BSA_PRM_PIONEER_ISOTOPE_CORRECTION_TEST_CORRECTED_062024_uPAC_50cm.arrow")
BSA_TABLE = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/PIONEER_BSA_TEST_062024/BSA_PRM_PIONEER_M0ONLY_062624_uPAC_50cm_20240626111010.arrow")
#best_scan_idx = 2749
BSA_TABLE = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/PIONEER_BSA_TEST_062024/BSA_PRM_PIONEER_ISOTOPE_CORRECTION_TEST_CORRECTED_062124_uPAC_50cm.arrow")
#best_scan_idx = 3658 best m0-m4
#BSA_TABLE = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/PIONEER_BSA_TEST_062024/BSA_PRM_PIONEER_ISOTOPE_CORRECTION_TEST_CORRECTED_20NCE_062624_uPAC_50cm.arrow")

hcat(BSA_TABLE[:centerMass][3655:3665], BSA_TABLE[:TIC][3655:3665])

"LVNELTEFAK"
best_scan_idx = 3658
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

bsa_fragments = Vector{DetailedFrag{Float32}}()
for i in range(1, size(tlist, 1))
    sequence = tlist[i,"Peptide Sequence"]
    ion_position = tlist[i,"Fragment Ion Ordinal"]
    sulfur_count = 0
    if tlist[i,"Fragment Ion Type"]=="b"
        sulfur_cont = UInt8(count("C", sequence[1:ion_position]) +  count("M", sequence[1:ion_position]))
    else
        sequence = reverse(sequence)
        sulfur_cont = UInt8(count("C", sequence[1:ion_position]) +  count("M", sequence[1:ion_position]))
    end

    frag_mz = Float32(tlist[i,:"Product Mz"])
    frag_charge = UInt8(tlist[i,"Product Charge"])
    for j in range(0, 3)
        push!(bsa_fragments, DetailedFrag(
            UInt32(seq_to_id[tlist[i,"Peptide Sequence"]]),
            Float32(frag_mz + j*NEUTRON/frag_charge),
            one(Float16),

            tlist[i,"Fragment Ion Type"]=="b" ? UInt8(1) : UInt8(2),
            #UInt8(i),
            false,

            frag_charge,
            UInt8(tlist[i,"Fragment Ion Ordinal"]),
            UInt8(count('+', tlist[i,"Precursor"])),
            UInt8(i),
            UInt8(sulfur_count)
        ))
    end
end

LVNELTEFAK_FRAGS = sort([frag for frag in bsa_fragments if frag.prec_id == 3], by = x->x.mz)

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
sort!(matched_frags, by=x->(x.predicted_rank, x.theoretical_mz), rev = false)
norm_factor = maximum([x.intensity for x in matched_frags])
bsa_fragments = Vector{DetailedFrag{Float32}}()
prec_charge = UInt8(2)
sequence = "LVNELTEFAK"
last_rank = 0
mono_mz = 0.0
j = 0
current_intensity = 0.0
for (i, frag) in enumerate(matched_frags)
    ion_position = frag.frag_index
    sulfur_count = 0
    if frag.ion_type==one(UInt8)#b_ions
        sulfur_cont = UInt8(count("C", sequence[1:ion_position]) +  count("M", sequence[1:ion_position]))
    else
        sequence = reverse(sequence)
        sulfur_cont = UInt8(count("C", sequence[1:ion_position]) +  count("M", sequence[1:ion_position]))
    end
    current_rank = frag.predicted_rank
    if current_rank != last_rank
        last_rank = current_rank
        j += 1
        current_intensity = frag.intensity/norm_factor
        mono_mz = frag.theoretical_mz
        push!(bsa_fragments,
        DetailedFrag(
            frag.prec_id,
            mono_mz,
            Float16(current_intensity),
            frag.ion_type,
            false,
            frag.frag_charge,
            frag.frag_index,
            prec_charge,
            UInt8(i),    
            UInt8(sulfur_count)
        ))   
    else
        current_intensity += frag.intensity/norm_factor
        bsa_fragments[j] = DetailedFrag(
            frag.prec_id,
            mono_mz,
            Float16(current_intensity),
            frag.ion_type,
            false,
            frag.frag_charge,
            frag.frag_index,
            prec_charge,
            UInt8(i),    
            UInt8(sulfur_count)
        )
    end
end      
sort!(bsa_fragments, by=x->x.intensity, rev = true)        
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
#BSA_TABLE = Arrow.Table("/Users/n.t.wamsley/TEST_DATA/PIONEER_BSA_TEST_062024/BSA_PRM_PIONEER_ISOTOPE_CORRECTION_TEST_CORRECTED_062024_uPAC_50cm.arrow")
params_[:deconvolution_params]["huber_delta"] = Float32(1000)
test_qtf = QuadTransmission(2.0, 10.0)
psms = vcat(secondSearch(
    BSA_TABLE, 
    params_;
    precursors=bsa_precursors,
    fragment_lookup_table = bsa_fragment_lookup_table,
    rt_index = retentionTimeIndex([rtIndexBin(0.0f0, 0.0f0, [(UInt32(3),  582.3190f0)])]),#RT_INDICES[file_id_to_parsed_name[ms_file_idx]],
    ms_file_idx = UInt32(ms_file_idx), 
    rt_to_irt_spline = test_spline,#RT_iRT[file_id_to_parsed_name[ms_file_idx]],
    mass_err_model = MassErrorModel{Float32}(0.0f0, (20.0f0, 20.0f0)),
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
    quad_transmission_func = test_qtf
    )...);

psms[!,:rt] = [BSA_TABLE[:retentionTime][x] for x in psms[!,:scan_idx]]
psms[!,:window_mz] = [BSA_TABLE[:centerMass][x] for x in psms[!,:scan_idx]]
sort!(psms, :rt)
getIsotopesCaptured!(psms, bsa_precursors[:prec_charge],bsa_precursors[:mz], BSA_TABLE)
gpsms = groupby(psms, :isotopes_captured)
p = plot(layout = (1, 2), title = "BSA_LVNELTEFAK")
for (i, iso_set) in enumerate([(0, 5), (0, 0), (1, 5), (2, 5)])
    psms_ = gpsms[(isotopes_captured = iso_set,)]
    plot!(psms_[!,:rt], psms_[!,:weight], show = true, xlim = (9.1, 9.3), color = i, label = string(iso_set), subplot = 1)
end
precursor_transmission = zeros(Float32, 5)
correctPrecursorAbundances!(
psms[!,:weight],
precursor_transmission,
test_qtf,
iso_splines,
psms[!,:window_mz],
psms[!,:isotopes_captured],
psms[!,:precursor_idx],
bsa_precursors[:mz],
bsa_precursors[:prec_charge],
bsa_precursors[:sulfur_count])
for (i, iso_set) in enumerate([(0, 5), (0, 0), (1, 5), (2, 5)])
    psms_ = gpsms[(isotopes_captured = iso_set,)]
    plot!(psms_[!,:rt], psms_[!,:weight], show = true, xlim = (9.1, 9.3),color = i, label = string(iso_set), subplot = 2)
end
#savefig("/Users/n.t.wamsley/Documents/PRESENTATIONS/lab_meetings/bsa_correction_v2_box.pdf")


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
psms[170:200,[:isotopes_captured,:weight,:scan_idx]]

test_qtf = QuadTransmission(2.0, 10.0)
@time RESULT = plotPSM(
    BSA_TABLE, 
    3652,
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
    quad_transmission_func = test_qtf
    );


test_quad = QuadTransmission(4.0f0, 5.0f0)
plot(LinRange(-6.0, 6.0, 100), [test_quad(0.0, x) for x in LinRange(-6.0, 6.0, 100)])

isotopes = zeros(Float32, 5)
test_frag = isotope(500.0f0, 0, 0)
test_prec = isotope(1000.0f0, 0, 0)
test_pset = (0, 5)
getFragAbundance!(
    isotopes,
    iso_splines,
    test_frag,
    test_prec,
    test_pset
)
println(isotopes)
julia> println(isotopes)
Float32[0.7593282, 0.19801514, 0.036832526, 0.005068335, 0.000585617]

isotopes = zeros(Float32, 5)
test_frag = isotope(500.0f0, 0, 0)
test_prec = isotope(1000.0f0, 0, 0)
test_pset = Float32[1.0, 1.0, 1.0, 1.0, 1.0]
getFragAbundance!(
    isotopes,
    test_pset,
    iso_splines,
    test_frag,
    test_prec,
)
println(isotopes)
julia> println(isotopes)
Float32[0.7592788, 0.197894, 0.036644634, 0.004880442, 0.00046448276]


############
#New subset
############

isotopes = zeros(Float32, 5)
test_frag = isotope(500.0f0, 0, 0)
test_prec = isotope(1000.0f0, 0, 0)
test_pset = (1, 3)
getFragAbundance!(
    isotopes,
    iso_splines,
    test_frag,
    test_prec,
    test_pset
)
println(isotopes)
julia> println(isotopes)
Float32[0.18222821, 0.19688448, 0.035286143, 0.0038709282, 0.0]

isotopes = zeros(Float32, 5)
test_frag = isotope(500.0f0, 0, 0)
test_prec = isotope(1000.0f0, 0, 0)
test_pset = Float32[0.0, 1.0, 1.0, 1.0, 0.0]
getFragAbundance!(
    isotopes,
    test_pset,
    iso_splines,
    test_frag,
    test_prec,
)
println(isotopes)
julia> println(isotopes)
Float32[0.18222821, 0.19688448, 0.035286143, 0.0038709282, 0.0]


test_qtf = QuadTransmission(2.0, 5.0)
getPrecursorIsotopeTransmission!(
    test_pset, 
    582.319f0,
    UInt8(2),
    584.5698f0,
    test_qtf
)
println(test_qtf)




function (qtf::QuadTransmission)(window_center::T, x::T) where {T<:AbstractFloat}
    return one(T)/(1 + abs((x - window_center)/qtf.width)^(T(2)*qtf.b))
end

test_qtf = QuadTransmission(2.0, 5.0)
getPrecursorIsotopeTransmission!(
    test_pset, 
    582.319f0,
    UInt8(2),
    580.0f0,
    test_qtf
)
println(test_pset[2]/test_pset[1])


plot(LinRange(-4.0, 4.0, 100), [test_qtf(0.0, x) for x in LinRange(-4.0, 4.0, 100)])

struct QuadTransmission{T<:AbstractFloat}
    width::T
    b::T
end

function (qtf::QuadTransmission)(window_center::T, x::T) where {T<:AbstractFloat}
    return one(T)/(1 + abs((x - (window_center))/(qtf.width))^(T(2)*qtf.b))
end

test_qtf = QuadTransmission(2.0, 5.0)
getPrecursorIsotopeTransmission!(
    test_pset, 
    582.319f0,
    UInt8(2),
    580.8f0,
    test_qtf
)
println(test_pset)
println(test_pset[2]/test_pset[1])


test_qtf = QuadTransmission(2.0, 5.0)
getPrecursorIsotopeTransmission!(
    test_pset, 
    582.319f0,
    UInt8(2),
    584.3198f0,
    test_qtf
)
println(test_pset)
println(test_pset[2]/test_pset[1])






test_qtf = QuadTransmission(2.0, 5.0)
getPrecursorIsotopeTransmission!(
    test_pset, 
    488.5345f0,
    UInt8(2),
    488.5345f0 - 2.319f0,
    test_qtf
)
println(test_pset)
println(test_pset[1]/test_pset[2])

test_qtf = QuadTransmission(2.0, 5.0)
getPrecursorIsotopeTransmission!(
    test_pset, 
    722.3247f0,
    UInt8(2),
    722.3247f0 - 2.319f0,
    test_qtf
)
println(test_pset)
println(test_pset[1]/test_pset[2])



plot(LinRange(-6.0, 6.0, 100), [test_qtf(0.0, x) for x in LinRange(-4.0, 4.0, 100)])

struct QuadTransmission{T<:AbstractFloat}
    width::T
    b::T
end

function (qtf::QuadTransmission)(window_center::T, x::T) where {T<:AbstractFloat}
    return one(T)/(1 + abs((x - window_center)/qtf.width)^(T(2)*qtf.b))
end







