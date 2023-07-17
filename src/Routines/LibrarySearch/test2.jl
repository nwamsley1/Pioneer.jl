function integrateRAW(
                    spectra::Arrow.Table, 
                    rt_index::retentionTimeIndex{T, U},
                    fragment_list::Vector{Vector{LibraryFragment{Float64}}},
                    ms_file_idx::UInt32;
                    fragment_tolerance::Float64 = 40.0,
                    quadrupole_isolation_width::Float64 = 8.5,
                    max_peak_width::Float64 = 2.0,
                    λ::Float32 = Float32(2e12),
                    γ::Float32 = Float32(1/2),
                    max_iter::Int = 1000,
                    nmf_tol::Float32 = Float32(100.0)
                    ) where {T,U<:Real}
    
    ms2 = 0
    nmf = Dict(:precursor_idx => UInt32[], :weight => Float32[], :rt => Float32[])
    for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(spectra)))

        if spectrum[:msOrder] == 1
            continue
        else
            ms2 += 1
        end

        #Get peptides that could be in the spectra
        transitions = selectTransitions(fragment_list, rt_index, Float64(spectrum[:retentionTime]), max_peak_width/2.0, spectrum[:precursorMZ], Float32(quadrupole_isolation_width/2.0))

        #Match fragments to peaks
        fragmentMatches = matchPeaks(transitions, 
                                    spectrum[:masses], 
                                    spectrum[:intensities], 
                                    δs = zeros(T, (1,)),
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = zero(Float32),
                                    ppm = fragment_tolerance
                                    )

        if iszero(length(fragmentMatches))
            continue
        end

        #Build templates for regrssion
        X, Hst, Hs, IDtoROW = buildDesignMatrix(fragmentMatches)
        
        #println("Hs.m ", Hs.m)
        #println("Hs.m ", Hs.n)
        #println(Hst)
        #println(X')
        #Fit non-negative adaptive LASSO
        weights = sparseNMF(Hs, Hst, X; λ=λ, γ=γ, max_iter=max_iter, tol=nmf_tol)

        for key in keys(IDtoROW)
            push!(nmf[:precursor_idx], key)
            push!(nmf[:weight], weights[IDtoROW[key]])
            push!(nmf[:rt], spectrum[:retentionTime])
        end
    end
    nmf = DataFrame(nmf)
    sort!(nmf, [:precursor_idx,:rt]);
    return groupby(nmf, :precursor_idx)
end