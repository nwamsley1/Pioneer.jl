isotopes = sort(collect(zip(getIsotopes(best_psms[:,:sequence], best_psms[:,:precursor_idx], QRoots(6), 6, 2), best_psms[:,:RT])), by = x->last(x))
export FragmentMatch, getNearest, matchPeaks, matchPeaks!
best_psms = DataFrame(CSV.File("/Users/n.t.wamsley/Desktop/best_psms_080423.csv"))
@time isotopes =  getIsotopes(best_psms[:,:sequence], best_psms[:,:precursor_idx], best_psms[:,:charge], QRoots(4), 4)
matches, misses = matchPeaks(isotopes[0x0006bbe9],
                            MS_TABLE[:masses][50009],
                            MS_TABLE[:intensities][50009],
                            PrecursorMatch{Float32},
                            count_unmatched=true,
                            δs = zeros(Float64, (1,)),
                            ppm = 10.0
                            )

prec_rt_table = sort(collect(zip(best_psms[:,:RT], UInt32.(best_psms[:,:precursor_idx]))), by = x->first(x))

function selectIsotopes(prec_list::Vector{Tuple{Float64, UInt32}}, isotope_dict::UnorderedDictionary{UInt32, Vector{Isotope{U}}}, rt::T, rt_tol::T) where {T,U<:AbstractFloat}
    isotopes = Vector{Isotope{U}}()
    i = 1
    rt_start = searchsortedfirst(prec_list, rt - rt_tol, lt=(r,x)->first(r)<x) #First RT bin to search
    rt_stop = searchsortedlast(prec_list, rt + rt_tol, lt=(x, r)->first(r)>x) #Last RT bin to search 
    #return rt_start, rt_stop
    for i in range(rt_start, rt_stop)
        append!(isotopes, isotope_dict[last(prec_list[i])])
    end
    return sort(isotopes, by = x->getMZ(x))
end

test_df = integrateRAW(MS_TABLE, prec_rt_table, isotopes, one(UInt32), precursor_tolerance = 6.5, scan_range = (0, 300000), λ = Float32(1e6), γ = Float32(1/2))

function integrateRAW(
                    spectra::Arrow.Table, 
                    #rt_index::retentionTimeIndex{T, U},
                    prec_list::Vector{Tuple{Float64, UInt32}},
                    isotopes::UnorderedDictionary{UInt32, Vector{Isotope{Float32}}},
                    ms_file_idx::UInt32;
                    precursor_tolerance::Float64 = 20.0,
                    quadrupole_isolation_width::Float64 = 8.5,
                    max_peak_width::Float64 = 2.0,
                    λ::Float32 = Float32(2e12),
                    γ::Float32 = Float32(1/2),
                    max_iter::Int = 1000,
                    nmf_tol::Float32 = Float32(100.0),
                    scan_range::Tuple{Int64, Int64} = (0, 0), 
                    ) where {T,U<:Real}
    
    ms1 = 0
    nmf = Dict(:precursor_idx => UInt32[], :weight => Float32[], :rt => Float32[])
    matches = ""
    misses = ""
    for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(spectra)))

        if spectrum[:msOrder] == 2
            continue
        else
            ms1 += 1
        end
        if scan_range != (0, 0)
            i < first(scan_range) ? continue : nothing
            i > last(scan_range) ? continue : nothing
        end
        #Get peptides that could be in the spectra
        #transitions = selectTransitions(fragment_list, rt_index, Float64(spectrum[:retentionTime]), max_peak_width/2.0, spectrum[:precursorMZ], Float32(quadrupole_isolation_width/2.0))
        #isotopes[0x0006bbe9]
        #Match fragments to peaks
        iso = selectIsotopes(prec_rt_table, isotopes, Float64(spectrum[:retentionTime]), 1.0)
        matches, misses = matchPeaks(iso,
                            spectrum[:masses],
                            spectrum[:intensities],
                            PrecursorMatch{Float32},
                            count_unmatched=true,
                            δs = zeros(Float64, (1,)),
                            ppm = precursor_tolerance
                            )

        #=fragmentMatches, fragmentMisses = matchPeaks(transitions, 
                                    spectrum[:masses], 
                                    spectrum[:intensities], 
                                    count_unmatched =true,
                                    δs = zeros(T, (1,)),
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = zero(Float32),
                                    ppm = fragment_tolerance
                                    )=#


        if iszero(length(matches))
            continue
        end

        #Build templates for regrssion. 
        #Do we need to remove precursors with less than N matched fragments?
        #if spectrum[:retentionTime] < 49.8483 
        #    continue
        #elseif spectrum[:retentionTime] > 49.8992 
        #end
        X, Hs, Hst, IDtoROW = buildDesignMatrix(matches, misses)
        #return X, Hs, Hst, IDtoROW
        weights = sparseNMF(Hst, Hs, X; λ=λ,γ=γ, max_iter=max_iter, tol=nmf_tol)

        for key in keys(IDtoROW)
            push!(nmf[:precursor_idx], key)
            push!(nmf[:weight], weights[IDtoROW[key]])
            push!(nmf[:rt], spectrum[:retentionTime])
        end
    end
    nmf = DataFrame(nmf)
    sort!(nmf, [:precursor_idx,:rt]);
    return groupby(nmf, :precursor_idx)
    #matches, misses
end


1928598
integratePrecursor(test_df, UInt32(1928598), isplot = true)
transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(test_df, UInt32(psm[:precursor_idx]), isplot = false)) => [:intensity_ms1, :count_ms1, :SN_ms1, :slope_ms1, :peak_error_ms1,:apex_ms1,:fwhm_ms1]);

test = best_psms[(best_psms[:,:intensity].>0.0).&(best_psms[:,:intensity_ms1].>0.0),[:precursor_idx,:q_value,:decoy,:sequence,:intensity,:intensity_ms1]]
histogram2d(log2.(test[:,:intensity]), log2.(test[:,:intensity_ms1]), seriestype = :scatter)
histogram2d(log2.(test[test[:,:decoy].==false,:intensity]), log2.(test[test[:,:decoy].==false,:intensity_ms1]), seriestype = :scatter)


histogram2d(log2.(test[test[:,:decoy].==true,:intensity]), log2.(test[test[:,:decoy].==true,:intensity_ms1]), seriestype = :scatter)
best_psms[10000:20000,[:precursor_idx,:q_value,:decoy,:intensity,:peak_error,:fwhm,:SN,:RT,:count,:slope,:apex]]

test[:,:ms1_ms2_ratio] = (log2.(test[:,:intensity_ms1])).-(log2.(test[:,:intensity]))
histogram(test[test[:,:decoy].==false,:ms1_ms2_ratio], alpha = 0.5, normalize = :pdf, bins = 40)

histogram!(test[test[:,:decoy].==true,:ms1_ms2_ratio], alpha = 0.5, normalize = :pdf, bins = 40)

best_psms[best_psms[:,:intensity].>0.0,[:precursor_idx,:q_value,:decoy,:intensity,:peak_error,:fwhm,:SN,:RT,:count,:slope,:apex]]

best_psms[(best_psms[:,:intensity].>0.0) .& ((best_psms[:,:decoy].==true)),[:precursor_idx,:sequence,:q_value,:decoy,:intensity,:peak_error,:fwhm,:SN,:RT,:count,:slope,:apex]]

histogram(best_psms[(best_psms[:,:intensity].>0.0) .& ((best_psms[:,:decoy].==true)),:offset], normalize = :pdf, alpha = 0.5)
histogram!(best_psms[(best_psms[:,:intensity].>0.0) .& ((best_psms[:,:decoy].==false)),:offset], normalize = :pdf, alpha = 0.5)


histogram2d(best_psms[(best_psms[:,:intensity].>0.0) .& ((best_psms[:,:decoy].==false)),:cross_cor], log2.(best_psms[(best_psms[:,:intensity].>0.0) .& ((best_psms[:,:decoy].==false)),:intensity]), normalize = :pdf, alpha = 0.5)
testlow = (best_psms[:,:intensity].>0.0) .& (best_psms[:,:decoy].==true) .& (best_psms[:,:offset] .== -10)
testhigh = (best_psms[:,:intensity].>0.0) .& (best_psms[:,:decoy].==true) .& (best_psms[:,:offset] .== 0)

histogram((best_psms[testlow,:cross_cor]), normalize = :pdf, alpha = 0.5)
histogram!((best_psms[testhigh,:cross_cor]), normalize = :pdf, alpha = 0.5)



decoy_seqs = best_psms[(best_psms[:,:intensity].>0.0) .& ((best_psms[:,:decoy].==true)),:sequence]
length(collect(eachmatch(r"ox", join(decoy_seqs,""))))/length(decoy_seqs)

findnearest(decoy_seqs[1],target_seqs,Jaro())

target_seqs = best_psms[(best_psms[:,:intensity].>0.0) .& ((best_psms[:,:decoy].==false)),:sequence]
length(collect(eachmatch(r"ox", join(target_seqs,""))))/length(target_seqs)

best_psms[startswith.(best_psms[:,:sequence],r"[KR]").&((best_psms[:,:decoy].==false)),[:sequence,:decoy]]

best_psms[startswith.(best_psms[:,:sequence],r"[KR]").&((best_psms[:,:decoy].==true)),[:sequence,:decoy]]
#length(collect(eachmatch(r"ox", join(best_psms[(best_psms[:,:intensity].>0.0) .& ((best_psms[:,:decoy].==true)),:sequence],""))))

integratePrecursor(test_df, UInt32(3791737), isplot = true)
integratePrecursor(test_df, UInt32(1203901), isplot = true)

integratePrecursor(test_df, UInt32(7928840), isplot = true)
integratePrecursor(test_df, UInt32( 8688408), isplot = true)
integratePrecursor(test_df, UInt32(6937987), isplot = true) #Very Abundant 
integratePrecursor(test_df, UInt32(8465737 ), isplot = true)
integratePrecursor(test_df, UInt32(7394541), isplot = true) 
integratePrecursor(test_df, UInt32(7946661), isplot = true) 
integratePrecursor(test_df, UInt32(47044), isplot = true) 
integratePrecursor(test_df, UInt32(440364), isplot = true) 
integratePrecursor(test_df, UInt32(4994390 ), isplot = true) 

integratePrecursor(test_df, UInt32(652040), isplot = true)
pid = UInt32(6937987)

pid = UInt32(8465737)

pid = UInt32(6937987)

pid = UInt32(2574823) #important 
argmax(crosscor(test_df[(precursor_idx=pid,)][:,:weight], chroms[(precursor_idx=pid,)][:,:weight][2:end - 1]))
sum(crosscor(test_df[(precursor_idx=pid,)][:,:weight], chroms[(precursor_idx=pid,)][:,:weight][2:end - 1]))
pid = UInt32(6937987)
argmax(crosscor(test_df[(precursor_idx=pid,)][:,:weight], chroms[(precursor_idx=pid,)][:,:weight][2:end - 1]))
sum(crosscor(test_df[(precursor_idx=pid,)][:,:weight], chroms[(precursor_idx=pid,)][:,:weight][2:end - 1]))

pid = UInt32(  1438385)
plot(test_df[(precursor_idx=pid,)][:,:rt], test_df[(precursor_idx=pid,)][:,:weight], seriestype=:scatter)
plot!(chroms[(precursor_idx=pid,)][:,:rt], chroms[(precursor_idx=pid,)][:,:weight]*10, seriestype=:scatter)


plot(range(1, 33), crosscor(test_df[(precursor_idx=pid,)][:,:weight], chroms[(precursor_idx=pid,)][:,:weight][2:end - 1]), seriestype = :scatter)


pairs = getScanPairs(test_df[(precursor_idx=pid,)][:,:rt], chroms[(precursor_idx=pid,)][:,:rt])
plot(range(1, 33), crosscor(test_df[(precursor_idx=pid,)][pairs[1],:weight], chroms[(precursor_idx=pid,)][pairs[2],:weight]), seriestype=:scatter)

collect(range(-5, 5))[argmax(crosscor(test_df[(precursor_idx=pid,)][pairs[1],:weight], chroms[(precursor_idx=pid,)][pairs[2],:weight], collect(range(-5, 5))))]

getCrossCorr(test_df[(precursor_idx=pid,)], chroms[(precursor_idx=pid,)])
transform!(best_psms, AsTable(:) => ByRow(psm -> getCrossCorr(test_df, chroms, UInt32(psm[:precursor_idx]))) => [:offset,:cross_cor]);

function getCrossCorr(MS1::GroupedDataFrame{DataFrame}, MS2::GroupedDataFrame{DataFrame}, precursor_idx::I; lag::Int = 10) where {I<:Integer}
    function missingID(gdf::GroupedDataFrame{DataFrame}, precursor_idx::I)
        if !((precursor_idx=precursor_idx,) in keys(gdf)) #If the precursor is not found
            return true
        end
        return false
    end
    if missingID(MS1, precursor_idx) | missingID(MS2, precursor_idx)
        return (missing, missing)
    end
    return getCrossCorr(MS1[(precursor_idx =precursor_idx,)], MS2[(precursor_idx = precursor_idx, )], lag)
end
function getCrossCorr(MS1::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, MS2::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, lag::Int64 = 10)
    scan_pairs = getScanPairs(MS1[:,:rt], MS2[:,:rt])
    lag = min(lag, length(scan_pairs[1]) - 1)
    if lag < 3
        return (missing, missing)
    end
    lag = collect(range(-lag, lag))
    cross_cor = crosscor(MS1[scan_pairs[1],:weight], MS2[scan_pairs[2],:weight], lag)
    return lag[argmax(cross_cor)], maximum(cross_cor)
end

Y1 = fft(test_df[(precursor_idx=pid,)][:,:weight])
Y2 = fft(chroms[(precursor_idx=pid,)][:,:weight][2:end - 1])

test[test[:,:decoy].==true,:][1000:1010,:]
test[test[:,:decoy].==false,:][1000:1010,:]

best_psms[(best_psms[:,:decoy].==true).&(best_psms[:,:q_value].<=0.01),[:precursor_idx,:sequence]]
chroms[(precursor_idx=UInt32(6937987),)][:,:weight]
integratePrecursor(test_df, UInt32(7566318), isplot = true)

histogram(best_psms[(best_psms[:,:decoy].==false),:peak_error], alpha = 0.5, normalize = :pdf, bins = 40)
#histogram!(best_psms[(best_psms[:,:decoy].==false) .& (best_psms[:,:q_value].>0.01),:SN], alpha = 0.5, normalize = :pdf)
histogram!(best_psms[best_psms[:,:decoy].==true,:peak_error], alpha = 0.5, normalize = :pdf, bins = 40)

histogram(log2.(best_psms[(best_psms[:,:decoy].==false),:intensity]), alpha = 0.5, normalize = :pdf, bins = 40)
#histogram!(best_psms[(best_psms[:,:decoy].==false) .& (best_psms[:,:q_value].>0.01),:SN], alpha = 0.5, normalize = :pdf)
histogram!(log2.(best_psms[best_psms[:,:decoy].==true,:intensity]), alpha = 0.5, normalize = :pdf, bins = 40)

diannreport = DataFrame(CSV.File("/Users/n.t.wamsley/Desktop/report.pr_matrix.tsv"))

diannreport = DataFrame(CSV.File("/Users/n.t.wamsley/Desktop/reportlib.tsv"))
histogram(diannreport[:,"Tr_recalibrated"])

transform!(best_psms, AsTable(:) => ByRow(psm -> replace(psm[:,:sequence], "(ox)" => "")) => :stripped_sequence);

best_psms[:,:stripped_sequence] .= replace.(best_psms[:,:sequence], "(ox)" => "")
PSMs[:,:stripped_sequence] .= replace.(PSMs[:,:sequence], "(ox)" => "")

diann = Set(diannreport[:,"Stripped.Sequence"]) #∩ Set(best_psms[:,:stripped_sequence])
Set(best_psms[:,:stripped_sequence])
all_sequences = String[]
for i in eachindex(precursors_mouse_detailed_33NCEfixed)
    if isassigned(precursors_mouse_detailed_33NCEfixed, i)
        push!(all_sequences, replace(precursors_mouse_detailed_33NCEfixed[i].sequence, "(ox)"=>""))
    end
end
all_sequenes = Set(all_sequences)
targets_best = Set(best_psms[(best_psms[:,:q_value].<=0.01) .& (best_psms[:,:decoy].==false),:stripped_sequence])
targets_all = Set(PSMs[:,:stripped_sequence])
length(diann ∩ all_sequences) - length(diann ∩ targets_all)
length(diann ∩ targets_best)
length(diann) - length(diann ∩ all_sequences)
setdiff(diann, all_sequences)
setdiff(diann, targets)
setdiff(targets, diann)
