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

function integrateMS1(
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

@time isotopes =  getIsotopes(best_psms[:,:sequence], best_psms[:,:precursor_idx], best_psms[:,:charge], QRoots(4), 4)
prec_rt_table = sort(collect(zip(best_psms[:,:RT], UInt32.(best_psms[:,:precursor_idx]))), by = x->first(x))
test_df = integrateRAW(MS_TABLE, prec_rt_table, isotopes, one(UInt32), precursor_tolerance = 6.5, scan_range = (0, 300000), λ = Float32(0), γ = Float32(0))


#getCrossCorr(test_df, chroms, pid)
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
    if all(isnan.(cross_cor))
        return (missing, missing)
    end
    return lag[argmax(cross_cor)], maximum(cross_cor)
end

transform!(best_psms, AsTable(:) => ByRow(psm -> getCrossCorr(test_df, chroms, UInt32(psm[:precursor_idx]))) => [:offset,:cross_cor]);
transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(test_df, UInt32(psm[:precursor_idx]), isplot = false)) => [:intensity_ms1, :count_ms1, :SN_ms1, :slope_ms1, :peak_error_ms1,:apex_ms1,:fwhm_ms1]);