

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