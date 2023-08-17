

diannreport = DataFrame(CSV.File("/Users/n.t.wamsley/Desktop/report.tsv"))
function buildRTIndex(RTs::Vector{T}, prec_mzs::Vector{U}, prec_ids::Vector{I}, bin_rt_size::AbstractFloat) where {T,U<:AbstractFloat,I<:Integer}
    
    start_idx = 1
    start_RT =  RTs[start_idx]
    rt_index = retentionTimeIndex(T, U) #Initialize retention time index
    i = 1
    while i < length(RTs) + 1
        if ((RTs[min(i + 1, length(RTs))] - start_RT) > bin_rt_size) | (i == length(RTs))
            push!(rt_index.rt_bins, 
                    rtIndexBin(RTs[start_idx], #Retention time for first precursor in the bin
                          RTs[i],     #Retention time for last precursor in the bin
                        [(zero(UInt32), zero(Float32)) for _ in 1:(i - start_idx + 1)] #Pre-allocate precursors 
                        )
                )

            n = 1 #n'th precursor 
            for idx in start_idx:(min(i, length(RTs))) 
                rt_index.rt_bins[end].prec[n] = (prec_ids[idx], prec_mzs[idx]) #Add n'th precursor
                n += 1
            end

            sort!(rt_index.rt_bins[end].prec, by = x->last(x)) #Sort precursors by m/z
            i += 1
            start_idx = i
            start_RT = RTs[min(start_idx, length(RTs))]
            continue
        else
            i += 1
        end
    end


    function sortRTBins!(rt_index::retentionTimeIndex{T, U})
        for i in 1:length(rt_index.rt_bins)
            sort!(rt_index.rt_bins[i].prec, by = x->last(x));
        end
        return nothing
    end
    sortRTBins!(rt_index)
    return rt_index
end

buildRTIndex(PSMs::DataFrame, bin_rt_size::AbstractFloat = 0.1) = buildRTIndex(PSMs[:,:RT], PSMs[:,:prec_mz], PSMs[:,:precursor_idx], bin_rt_size)

missingPSMs = Dict{Symbol, Any}(
    :RT => Float32[],
    :prec_mz => Float32[],
    :precursor_idx => UInt32[]
)
seqs = Set{String}()
prec_to_rt = UnorderedDictionary{String, Float64}()

missingPSMs = Dict{Symbol, Any}(
    :RT => Float32[],
    :prec_mz => Float32[],
    :precursor_idx => UInt32[]
)

missing_peptides = setdiff(diann ∩ all_sequences, diann ∩ targets_all_big)
missing_seqs = diannreport[[(x ∈ missing_peptides) for x in diannreport[:,"Stripped.Sequence"]],:][:,"Stripped.Sequence"]
missing_charges = diannreport[[(x ∈ missing_peptides) for x in diannreport[:,"Stripped.Sequence"]],:][:,"Precursor.Charge"]
missing_peptides = Set([missing_seqs[i]*"_"*string(missing_charges[i]) for i in eachindex(missing_seqs)])
for i in ProgressBar(eachindex(precursors_mouse_detailed_33NCEcorrected_start1))
    if isassigned(precursors_mouse_detailed_33NCEcorrected_start1, i)
        seq = precursors_mouse_detailed_33NCEcorrected_start1[i].sequence
        irt = precursors_mouse_detailed_33NCEcorrected_start1[i].iRT
        charge = string(precursors_mouse_detailed_33NCEcorrected_start1[i].charge)
        stripped_seq = replace(seq, "(ox)" => "")
        if stripped_seq*"_"*charge ∉ missing_peptides
            continue
        end
        if stripped_seq ∉ seqs
            push!(missingPSMs[:RT], iRT_to_RT[irt])
            push!(missingPSMs[:prec_mz],  precursors_mouse_detailed_33NCEcorrected_start1[i].mz)
            push!(missingPSMs[:precursor_idx],  UInt32(i))
        end
    end
end
missingPSMs = DataFrame(missingPSMs)
sort!(missingPSMs,:RT)
rt_index = buildRTIndex(missingPSMs, 5.0)

DIANN_MISSING_PSMS = SearchRAW(MS_TABLE, rt_index, frags_mouse_detailed_33NCEcorrected_start1, one(UInt32))
DIANN_MISSING_PSMS_BEST = combine(sdf -> sdf[argmax(sdf.matched_ratio),:], groupby(DIANN_MISSING_PSMS, :precursor_idx))
DIANN_MISSING_PSMS_BEST[DIANN_MISSING_PSMS_BEST[:,:scan_idx].==  20189,:]

DIANN_MISSING_PSMS_BEST[DIANN_MISSING_PSMS_BEST[:,:precursor_idx].==  UInt32(865045),:]

PSMs[PSMs[:,:precursor_idx].==  UInt32(865045),:]

PSMs[PSMs[:,:scan_idx].==  UInt32(20189),:]

histogram(log2.(DIANN_MISSING_PSMS_BEST[:,:matched_ratio]), normalize = :pdf)
histogram!(log2.(PSMs[PSMs[:,:decoy].==true,:matched_ratio]), normalize = :pdf)
function SearchRAW(
                    spectra::Arrow.Table, 
                    #ptable::PrecursorDatabase,
                    rt_index::retentionTimeIndex{Float32, Float32},
                    fragment_list::Vector{Vector{LibraryFragment{Float32}}},
                    ms_file_idx::UInt32;
                    quadrupole_isolation_width::Float64 = 8.5,
                    #interpolation::Interpolations.Extrapolation{Float64, 1, Interpolations.GriddedInterpolation{Float64, 1, Float64, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Line{Nothing}};
                    #isolation_width::Float64 = 4.25,
                    #precursor_tolerance::Float64 = 5.0,
                    max_peak_width = 20.0,
                    fragment_tolerance::Float64 = 20.0,
                    #topN::Int64 = 20,
                    min_frag_count::Int64 = 4,
                    #min_matched_ratio::Float32 = Float32(0.8),
                    min_spectral_contrast::Float32 = Float32(0.0),
                    λ::Float32 = Float32(1e3),
                    γ::Float32 = zero(Float32),
                    scan_range::Tuple{Int64, Int64} = (0, 0), 
                    #max_peaks::Int = 200, 
                    max_iter::Int = 1000,
                    nmf_tol::Float32 = Float32(100.0),
                    #rt_tol::Float32 = Float32(30.0),
                    #fragment_match_ppm::U,
                    data_type::Type{T} = Float64
                    ) where {T,U<:Real}
    
    scored_PSMs = makePSMsDict(XTandem(data_type))
    ms2, MS1, MS1_i = 0, 0, 0
    precs = Counter(UInt32, UInt8, Float32, length(fragment_list)) #Prec counter
    times = Dict(:counter => 0.0, :reset => 0.0, :nmf => 0.0, :metrics => 0.0, :match_peaks => 0.0, :build => 0.0, :score => 0.0)
    n = 0
    kt = KendallTau(Float32)
    for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(spectra)))
    #for (i, spectrum) in enumerate(Tables.namedtupleiterator(spectra))

        if spectrum[:msOrder] == 1
            MS1 = spectrum[:masses]
            MS1_i = spectrum[:intensities]
            continue
        else
            ms2 += 1
        end

        if scan_range != (0, 0)
            i < first(scan_range) ? continue : nothing
            i > last(scan_range) ? continue : nothing
        end

        transitions = selectTransitions(fragment_list, rt_index, Float32(spectrum[:retentionTime]), Float32(max_peak_width/2.0), spectrum[:precursorMZ], Float32(quadrupole_isolation_width/2.0))
      
        #println("length(transitions) ", length(transitions))
        #return transitions
        #times[:reset] += @elapsed reset!(precs)
        reset!(precs)
        
        if length(transitions) == 0
            continue
        end
        
       

        #times[:match_peaks] += @elapsed fragmentMatches, fragmentMisses = matchPeaks(transitions, 
        fragmentMatches, fragmentMisses = matchPeaks(transitions, 
                                    spectrum[:masses], 
                                    spectrum[:intensities], 
                                    FragmentMatch{Float32},
                                    count_unmatched=true,
                                    δs = zeros(T, (1,)),
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = zero(Float32),
                                    ppm = fragment_tolerance
                                    )
        if iszero(length(fragmentMatches))
            continue
        end

        #times[:build] += @elapsed X, H, UNMATCHED, IDtoROW = buildDesignMatrix(fragmentMatches, fragmentMisses, topN)
        X, Hs, Hst, IDtoROW, matched_cols = buildDesignMatrix(fragmentMatches, fragmentMisses)
        
        times[:nmf] += @elapsed weights = sparseNMF(Hst, Hs, X; λ=λ,γ=γ, max_iter=max_iter, tol=nmf_tol)[:]

        scribe_score, city_block, matched_ratio, spectral_contrast_matched, spectral_contrast_all, kt_pval = getDistanceMetrics(Hst, X, matched_cols, kt)
        
        #For progress and debugging. 

        unscored_PSMs = UnorderedDictionary{UInt32, XTandem{T}}()

        ScoreFragmentMatches!(unscored_PSMs, fragmentMatches)

        #times[:score] += @elapsed Score!(scored_PSMs, unscored_PSMs, 
        Score!(scored_PSMs, unscored_PSMs, 
                length(spectrum[:intensities]), 
                Float64(sum(spectrum[:intensities])), 
                0.0, scribe_score, city_block, matched_ratio, spectral_contrast_matched, spectral_contrast_all, kt_pval, weights, IDtoROW,
                scan_idx = Int64(i),
                min_spectral_contrast = min_spectral_contrast,
                min_frag_count = min_frag_count
                )
        n += 1
    end

    println("processed $ms2 scans!")
    println("counter: ", times[:counter]/n)
    println("reset: ", times[:reset]/n)
    println("nmf : ", times[:nmf]/n)
    println("metrics : ", times[:metrics]/n)
    println("build : ", times[:build]/n)
    println("score : ", times[:score]/n)
    println("match_peaks : ", times[:match_peaks]/n)
    return DataFrame(scored_PSMs)
end
testPSMs = SearchRAW(MS_TABLE, prosit_mouse_33NCEcorrected_start1_5ppm_15irt,  frags_mouse_detailed_33NCEcorrected_start1, UInt32(1), linear_spline,
                        min_frag_count = 4, 
                        topN = 1000, 
                        fragment_tolerance = 15.6, 
                        λ = Float32(1e3), 
                        γ =Float32(1),
                        max_peaks = 10000, 
                        scan_range = (0, 20190), #101357 #22894
                        precursor_tolerance = 20.0,
                        min_spectral_contrast =  Float32(0.5),
                        min_matched_ratio = Float32(0.45),
                        rt_tol = Float32(20.0)
                        )
                        SearchRAW(MS_TABLE, prosit_mouse_33NCEcorrected_start1_5ppm_15irt,  frags_mouse_detailed_33NCEcorrected_start1, UInt32(1), linear_spline,
                        min_frag_count = 4, 
                        topN = 1000, 
                        fragment_tolerance = 15.6, 
                        λ = Float32(1e3), 
                        γ =Float32(1),
                        max_peaks = 10000, 
                        scan_range = (20189, 20189), #101357 #22894
                        precursor_tolerance = 20.0,
                        min_spectral_contrast =  Float32(0.5),
                        min_matched_ratio = Float32(0.45),
                        rt_tol = Float32(20.0)
                        )

                        testPSMs[testPSMs[:,:precursor_idx].==  UInt32(865045),:]

                        testPSMs[testPSMs[:,:scan_idx].==  20189,:]

A = Int64[]
for i in eachindex(test_precs.counts)
    if iszero(first(test_precs.counts[i])) == false
        push!(A, i)
    end
end


#########
#Compare diann seqs at 1% fdr
#########
diannreport = DataFrame(CSV.File("/Users/n.t.wamsley/Desktop/report.tsv"))
diannreport[:,"Unique.Sequence"] .= replace.(diannreport[:,"Modified.Sequence"], "C(UniMod:4)" => "C")
diannreport[:,"Unique.Sequence"] .= replace.(diannreport[:,"Unique.Sequence"], "M(UniMod:35)" => "M(ox)")


#=diannreport[:,"Unique.Sequence"] .= diannreport[:,"Unique.Sequence"].*"_".*string.(diannreport[:,"Precursor.Charge"])
#best_psms[:,:stripped_sequence] .= replace.(best_psms[:,:sequence], "(ox)" => "")
#PSMs[:,:stripped_sequence] .= replace.(PSMs[:,:sequence], "(ox)" => "")
#Unique (unmodified) sequences 

diann_and_titus = diann_passed ∩ titus_passed

diann_int = diannreport[[x ∈ diann_and_titus for x in diannreport[:,"Unique.Sequence"]],"Precursor.Quantity"]

best_psms[:,:unique_seq] = best_psms[:,:sequence]
best_psms[titus_passed,:unique_seq] = best_psms[titus_passed,:sequence].*"_".*string.(best_psms[titus_passed,:charge])

titus_int = best_psms[[x ∈ diann_and_titus for x in best_psms[:,:unique_seq]], :intensity_ms1]

Plots.plot(log10.(diann_int), log10.(titus_int), seriestype=:scatter)

diannreport[[x ∈ diann_and_titus for x in diannreport[:,"Unique.Sequence"]],"Precursor.Quantity"]=#

diann_passed = Set(diannreport[:,"Unique.Sequence"].*"_".*string.(diannreport[:,"Precursor.Charge"]))
titus_passed = (best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false)
titus_passed = Set(best_psms[titus_passed,:sequence].*"_".*string.(best_psms[titus_passed,:charge]))
#titus_seed= (PSMs[:,:q_value].<=0.1).&(PSMs[:,:decoy].==false)
titus_seed= (PSMs[:,:decoy].==false)
titus_seed = Set(PSMs[titus_seed,:sequence].*"_".*string.(PSMs[titus_seed,:charge]))

CSV.write("/Users/n.t.wamsley/Desktop/diann_passed_precursor.tsv",  Tables.table(collect(diann_passed)), writeheader=false)
CSV.write("/Users/n.t.wamsley/Desktop/titus_passed_precursor.tsv",  Tables.table(collect(titus_passed)), writeheader=false)
CSV.write("/Users/n.t.wamsley/Desktop/titus_seed_precursor.tsv",  Tables.table(collect(titus_seed)), writeheader=false)


Plots.histogram(best_psms[(best_psms[:,:q_value].<0.01).&(best_psms[:,:decoy].==false),:cross_cor], normalize = :probability, alpha = 0.5,bins = 100)
Plots.histogram!(best_psms[(best_psms[:,:decoy].==true),:cross_cor], normalize = :probability, alpha = 0.5,bins = 100)

Plots.histogram(best_psms[(best_psms[:,:q_value].<0.01).&(best_psms[:,:decoy].==false),:offset], normalize = :probability, alpha = 0.5,bins = 100)
Plots.histogram!(best_psms[(best_psms[:,:decoy].==true),:offset], normalize = :probability, alpha = 0.5,bins = 100)

####Pep level
diann_passed = Set(diannreport[:,"Unique.Sequence"])
titus_passed = (best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false)
titus_passed = Set(best_psms[titus_passed,:sequence])
#titus_seed= (PSMs[:,:q_value].<=0.1).&(PSMs[:,:decoy].==false)
titus_seed= (PSMs[:,:decoy].==false)
titus_seed = Set(PSMs[titus_seed,:sequence])


setdiff(diann_passed, titus_seed)
diann_passed ∩ titus_passed
pep = "GGGALVENTTTGLSR"
diannreport[diannreport[:,:"Unique.Sequence"] .== pep,["RT.Start","RT","RT.Stop"]]
AVMNSQQGIEYILSNQGYVR
"YNPGYVLAGR" ∈ setdiff(diann_passed, titus_passed)
######
length(setdiff(diann_passed, titus_passed)) #21455

length(setdiff(diann_passed, titus_seed)) #16752

length(setdiff(setdiff(diann_passed, titus_passed), setdiff(diann_passed, titus_seed))) #4703

length(setdiff(titus_passed, diann_passed)) #4236

length(titus_passed ∩ diann_passed)

######
diann_passed = diann_passed ∩ titus_seed
length(setdiff(diann_passed, titus_passed))
length(diann_passed ∩ titus_passed)
length(setdiff(titus_passed, diann_passed))
"""
    There are ~7000 sequences that didn't pass the initial 10% fdr threshold at the psm level
    There are another ~7000 sequences that pass the initial 10% but that don't achieve 1% fdr
    after peak integration.

    There are also some sequences that diann does not detect at 1% fdr but that titus does. Are 
    these good id's or false positives?
"""

missing_seqs = setdiff(setdiff(diann_passed, titus_passed), setdiff(diann_passed, titus_seed))

missing_peptides = best_psms[[x∈missing_seqs for x in best_psms[:,:sequence]],:]
titus_targets = best_psms[(best_psms[:,:q_value].<0.01).&(best_psms[:,:decoy].==false),:]
titus_decoys = best_psms[(best_psms[:,:q_value].>0.01).&(best_psms[:,:decoy].==true),:]

###########
histogram(log2.(missing_peptides[:,:matched_ratio]), normalise =:pdf, alpha = 0.5, bins = 100)
histogram!(log2.(titus_decoys[:,:matched_ratio]), normalise = :pdf, alpha = 0.5,bins = 100)
histogram!(log2.(titus_targets[:,:matched_ratio]), normalise = :pdf, alpha = 0.5, bins = 100)
"""
Missing Peptides have lower matched ratios, similar to the decoys

julia> median(missing_peptides[:,:matched_ratio])
2.134892f0
julia> median(titus_decoys[:,:matched_ratio])
1.7648417f0
median(titus_targets[:,:matched_ratio])
8.028733f0
"""

###########
histogram((missing_peptides[:,:scribe_score]), normalise = :probability, alpha = 0.5,bins = 100)
histogram!((titus_decoys[:,:scribe_score]), normalise =:probability, alpha = 0.5,bins = 100)
histogram!((titus_targets[:,:scribe_score]), normalise =:probability, alpha = 0.5,bins = 100)
"""
Missing Peptides have lower scribe scores

julia> median(missing_peptides[:,:scribe_score])
10.916118f0
julia> median(titus_decoys[:,:scribe_score])
10.210345f0
julia> median(titus_targets[:,:scribe_score])
14.558072f0
"""

############
histogram((missing_peptides[:,:spectral_contrast_all]), normalize = :probability, alpha = 0.5,bins = 100)
histogram!((titus_decoys[:,:spectral_contrast_all]), normalize=:probability, alpha = 0.5,bins = 100)
histogram!((titus_targets[:,:spectral_contrast_all]), normalize =:probability, alpha = 0.5,bins = 100)

"""
Missing Peptides have lower scribe scores

julia> median(missing_peptides[:,:spectral_contrast_all])
0.852202f0
julia> median(titus_decoys[:,:spectral_contrast_all])
0.81355864f0
julia> median(titus_targets[:,:spectral_contrast_all])
0.9612154f0
"""

############
histogram((missing_peptides[:,:RT_error]), normalize = :probability, alpha = 0.5,bins = 100)
histogram!((titus_decoys[:,:RT_error]), normalize =:probability, alpha = 0.5,bins = 100)
histogram!((titus_targets[:,:RT_error]),normalize =:probability, alpha = 0.5,bins = 100)
"""
Interesting because here the decoys have the greatest RT errors but the missing peptides
and passing targets have similair rt errors. 

julia> median(missing_peptides[:,:RT_error])
4.41528702089964
julia> median(titus_decoys[:,:RT_error])
5.4556267604548765
julia> median(titus_targets[:,:RT_error])
4.376064849485516
"""

############
histogram((missing_peptides[:,:total_ions]), normalize = :probability, alpha = 0.5,bins = 100)
histogram!((titus_decoys[:,:total_ions]), normalize =:probability, alpha = 0.5,bins = 100)
histogram!((titus_targets[:,:total_ions]), normalize =:probability, alpha = 0.5,bins = 100)
"""
julia> median(missing_peptides[:,:total_ions])
5.0

julia> median(titus_decoys[:,:total_ions])
5.0

julia> median(titus_targets[:,:total_ions])
9.0
"""

bar(countmap(missing_peptides[:,:topn])./sum(countmap(missing_peptides[:,:topn])), normalize =:probability, alpha = 0.5)
bar!(countmap(titus_decoys[:,:topn])./sum(countmap(titus_decoys[:,:topn])), normalize = :probability, alpha = 0.5)
bar!(countmap(titus_targets[:,:topn])./sum(countmap(titus_targets[:,:topn])), normalize = :probability, alpha = 0.5)
"""
Missing Peptides have lower matched ratios, similar to the decoys

julia> median(missing_peptides[:,:matched_ratio])
2.134892f0
julia> median(titus_decoys[:,:matched_ratio])
1.7648417f0
median(titus_targets[:,:matched_ratio])
8.028733f0
"""

############
b_range = range(-1, 12, length=100)
Plots.histogram(log2.(missing_peptides[:,:matched_ratio]), bins = b_range, normalize = :probability, alpha = 0.5, label = "DIA-NN Only")
Plots.histogram!(log2.(titus_decoys[:,:matched_ratio]), bins = b_range, normalize = :probability, alpha = 0.25,  label = "TITUS Decoys")
Plots.histogram!(log2.(titus_targets[:,:matched_ratio]), bins = b_range, normalize = :probability, alpha = 0.5,label = "TITUS <=1% FDR")

Plots.histogram(log2.(missing_peptides[:,:matched_ratio]), normalize = :probability, alpha = 0.5,bins = 100, label = "DIA-NN Only")
Plots.histogram!(log10.(titus_decoys[:,:intensity]), normalize =:probability, alpha = 0.5, bins = 100, label = "TITUS Decoys")
Plots.histogram!(log10.(titus_targets[:,:intensity]), normalize =:probability, alpha = 0.5, bins = 100, label = "TITUS <=1% FDR")

Plots.histogram(log10.(missing_peptides[:,:intensity]), alpha = 0.5,bins = 100, label = "DIA-NN Only")
#Plots.histogram!(log10.(titus_decoys[:,:intensity]), alpha = 0.25, bins = 100, label = "TITUS Decoys")
Plots.histogram!(log10.(titus_targets[:,:intensity]), alpha = 0.5, bins = 100, label = "TITUS <=1% FDR")

b_range = range(0, 25, length=100)
Plots.histogram((missing_peptides[:,:scribe_score]), bins = b_range, normalize = :probability, alpha = 0.5, label = "DIA-NN Only")
Plots.histogram!((titus_decoys[:,:scribe_score]), bins = b_range, normalize = :probability, alpha = 0.25,  label = "TITUS Decoys")
Plots.histogram!((titus_targets[:,:scribe_score]), bins = b_range, normalize = :probability, alpha = 0.5,label = "TITUS <=1% FDR")
"""
Missing peptides actually have lower intensities than the decoys

julia> median(missing_peptides[:,:intensity])
5.0

julia> median(titus_decoys[:,:intensity])
5.0

julia> median(titus_targets[:,:intensity])
9.0
"""
#=
###############################################################################################
#Is it possible that our RT predictions for the diann-specific peptides were particularly weak 
#and that DIA-NN had better predictions?
#1) get DIA-NN RT errors for titus_targets
#2) get DIA-NN RT errors for missing_peptides
#3) Compare these errors to the prosit errors
#Hypothesis is that the DIA-NN RT errors are greater for the missing peptides than for teh titus_targets
#even though the prosit errors are the same for both groups
###############################################################################################
=#
titus_targets_set = Set(titus_targets[:,:stripped_sequence])
missing_peptides_set = Set(missing_peptides[:,:stripped_sequence])
diann_titus_targets = diannreport[[x∈titus_targets_set for x in diannreport[:,"Stripped.Sequence"]],:]
diann_missing_peptides  = diannreport[[x∈missing_peptides_set for x in diannreport[:,"Stripped.Sequence"]],:]

############
histogram(abs.(diann_titus_targets[:,"RT"] .- diann_titus_targets[:,"Predicted.RT"]), normalize = :probability, alpha = 0.5,bins = 100)
histogram!(abs.(diann_missing_peptides[:,"RT"] .- diann_missing_peptides[:,"Predicted.RT"]), normalize = :probability, alpha = 0.5,bins = 100)
"""
RT error distributions are the same, but about half the prosit errors

julia> median(abs.(diann_titus_targets[:,"RT"] .- diann_titus_targets[:,"Predicted.RT"]))
2.5766999999999953

julia> median(abs.(diann_missing_peptides[:,"RT"] .- diann_missing_peptides[:,"Predicted.RT"]))
2.4657999999999944
"""

##########
#Are there other features in DIA-NN we could use to discriminate?
histogram(log2.(diann_titus_targets[:,"Precursor.Normalised"]), normalize = :probability, alpha = 0.5,bins = 100)
histogram!(log2.(diann_missing_peptides[:,"Precursor.Normalised"]), normalize = :probability, alpha = 0.5,bins = 100)



histogram([mean([parse(Float32, f) for f in n]) for n in [x[1:(end - 1)] for x in [split(x,';') for x in diann_titus_targets[:,"Fragment.Correlations"]]]], normalize = :probability, alpha = 0.5,bins = 100)
histogram!([mean([parse(Float32, f) for f in n]) for n in [x[1:(end - 1)] for x in [split(x,';') for x in diann_missing_peptides[:,"Fragment.Correlations"]]]], normalize = :probability, alpha = 0.5,bins = 100)



histogram(PSMs[PSMs[:,:decoy],:scribe_score], normalize = :probability, alpha = 0.5,bins = 100)
histogram!(PSMs[(PSMs[:,:decoy].==false) .& (PSMs[:,:q_value].<=0.01),:scribe_score], normalize = :probability, alpha = 0.5,bins = 100)

histogram!(log2.(diann_missing_peptides[:,"Precursor.Normalised"]), normalize = :probability, alpha = 0.5,bins = 100)

plot(best_psms[best_psms[:,:q_value].<0.01,:RT_pred], best_psms[best_psms[:,:q_value].<0.01,:RT], seriestype = :scatter, alpha = 0.1)
plot(best_psms[best_psms[:,:q_value].<0.01,:iRT], best_psms[best_psms[:,:q_value].<0.01,:RT], seriestype = :scatter, alpha = 0.1)
plot!(best_psms[best_psms[:,:decoy],:iRT], best_psms[best_psms[:,:decoy],:RT], seriestype = :scatter, alpha = 0.1)