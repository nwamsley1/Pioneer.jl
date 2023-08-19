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
                    prec_rt_table::Vector{Tuple{Float64, UInt32}},
                    isotopes::UnorderedDictionary{UInt32, Vector{Isotope{Float32}}},
                    ms_file_idx::UInt32;
                    precursor_tolerance::Float64 = 20.0,
                    λ::Float32 = Float32(2e12),
                    γ::Float32 = Float32(1/2),
                    max_iter::Int = 1000,
                    nmf_tol::Float32 = Float32(100.0),
                    scan_range::Tuple{Int64, Int64} = (0, 0), 
                    )
    
    ms1 = 0
    nmf = Dict(:precursor_idx => UInt32[], :weight => Float32[], :rt => Float32[])
    matches = ""
    misses = ""
    fragmentMatches = [Isotope{Float32}() for x in range(1, 10000)]
    fragmentMisses = [Isotope{Float32}() for x in range(1, 10000)]
    #for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(spectra)))
    for (i, spectrum) in enumerate(Tables.namedtupleiterator(spectra))


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
        nmatches, nmisses = matchPeaks(iso,
                            fragmentMatches,
                            fragmentMisses,
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
        X, Hs, IDtoROW = buildDesignMatrix(fragmentMatches, fragmentMisses, nmatches, nmisses)

        for i in range(1, nmatches)
            fragmentMatches[i] = FragmentMatch{Float32}()
        end

        for i in range(1, nmisses)
            fragmentMisses[i] = FragmentMatch{Float32}()
        end
        
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

