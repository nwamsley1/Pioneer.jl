isotopes = sort(collect(zip(getIsotopes(best_psms[:,:sequence], best_psms[:,:precursor_idx], QRoots(6), 6, 2), best_psms[:,:RT])), by = x->last(x))
export FragmentMatch, getNearest, matchPeaks, matchPeaks!
best_psms = DataFrame(CSV.File("/Users/n.t.wamsley/Desktop/best_psms_080423.csv"))
@time isotopes =  getIsotopes(best_psms[:,:sequence], best_psms[:,:precursor_idx], best_psms[:,:charge], QRoots(6), 6)
matches, misses = matchPeaks(isotopes[0x0006bbe9],
                            MS_TABLE[:masses][50009],
                            MS_TABLE[:intensities][50009],
                            PrecursorMatch{Float32},
                            count_unmatched=true,
                            δs = zeros(Float64, (1,)),
                            ppm = 10.0
                            )


function integrateRAW(
                    spectra::Arrow.Table, 
                    rt_index::retentionTimeIndex{T, U},
                    isotopes::UnorderedDictionary{UInt32, Vector{Isotope{Float32}}},
                    ms_file_idx::UInt32;
                    fragment_tolerance::Float64 = 40.0,
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

        #Match fragments to peaks
        fragmentMatches, fragmentMisses = matchPeaks(transitions, 
                                    spectrum[:masses], 
                                    spectrum[:intensities], 
                                    count_unmatched =true,
                                    δs = zeros(T, (1,)),
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = zero(Float32),
                                    ppm = fragment_tolerance
                                    )


        if iszero(length(fragmentMatches))
            continue
        end

        #Build templates for regrssion. 
        #Do we need to remove precursors with less than N matched fragments?
        #X, Hs, Hst, IDtoROW = buildDesignMatrix(fragmentMatches, fragmentMisses)
        #weights = sparseNMF(Hst, Hs, X; λ=λ,γ=γ, max_iter=max_iter, tol=nmf_tol)

        #for key in keys(IDtoROW)
        #    push!(nmf[:precursor_idx], key)
        #    push!(nmf[:weight], weights[IDtoROW[key]])
        #    push!(nmf[:rt], spectrum[:retentionTime])
        #end
    end
    #nmf = DataFrame(nmf)
    #sort!(nmf, [:precursor_idx,:rt]);
    #return groupby(nmf, :precursor_idx)
    fragmentMatches, fragmentMisses
end