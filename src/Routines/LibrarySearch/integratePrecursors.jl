best_psms = combine(sdf -> sdf[argmax(sdf.hyperscore),:], groupby(PSMs[PSMs[:,:q_values].<=0.1,:], :precursor_idx))
transform!(best_psms, AsTable(:) => ByRow(psm -> prec_mzs[psm[:precursor_idx]]) => :prec_mz)
sort!(best_psms, :RT)



function integrateRAW(
                    spectra::Arrow.Table, 
                    rt_index::retentionTimeIndex{T, U},
                    fragment_list::Vector{Vector{LibraryFragment{Float64}}},
                    ms_file_idx::UInt32;
                    fragment_tolerance::Float64 = 40.0,
                    quadrupole_isolation_width::Float64 = 8.5,
                    lambda::Float64 = 1e5,
                    max_peak_width::Float64 = 2.0
                    ) where {T,U<:Real}
    
    ms2 = 0
    nmf = Dict(:precursor_idx => UInt32[], :weight => Float32[], :rt => Float32[])
    for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(spectra)))

        if spectrum[:msOrder] == 1
            continue
        else
            ms2 += 1
        end

        transitions = selectTransitions(fragment_list, rt_index, Float64(spectrum[:retentionTime]), max_peak_width/2.0, spectrum[:precursorMZ], Float32(quadrupole_isolation_width/2.0))

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

        X, H, IDtoROW = buildDesignMatrix(fragmentMatches)
        
        #Initialize weights for each precursor template. 
        #Should find a more sophisticated way of doing this. 
        W = reshape([Float32(1000) for x in range(1,size(H)[1])], (1, size(H)[1]))
        weights = (NMF.solve!(NMF.MultUpdate{Float32}(maxiter=50, verbose = false, 
                                                    lambda_w = lambda, 
                                                    tol = 100, #Need a reasonable way to choose lambda?
                                                    update_H = false #Important to keep H constant. 
                                                    ), X, W, H).W[1,:])

        #For progress and debugging. 
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

function plotChromatogram(chroms::GroupedDataFrame{DataFrame}, precursor_idx::UInt32)
    chrom = chroms[(precursor_idx=precursor_idx,)]
    rt = chrom[:,:rt]
    f_width = min(11, (length(rt)÷2)*2 - 1)
    #if f_width < 4
    #    return (0.0, 0)
    #end
    intensity = savitzky_golay(chrom[:,:weight], 5, 3).y
    println("Peak area ", integrate(rt, intensity, TrapezoidalFast()))
    plot(rt, intensity)
    plot!(rt, chrom[:,:weight], seriestype=:scatter)
end


function integratePrecursor(chroms::GroupedDataFrame{DataFrame}, precursor_idx::UInt32; isplot::Bool = false)
    if !((precursor_idx=precursor_idx,) in keys(chroms))
        return (0.0, 0)
    end
    chrom = chroms[(precursor_idx=precursor_idx,)]
    rt = chrom[:,:rt]
   #println(chrom[:,:weight])
    #for weight in chrom[:,:weight]
    function find_longest_streak_indices(arr::BitVector)
        longest_start = 0
        longest_stop = 0
        current_start = 0
        current_streak = 0
        longest_streak = 0
        
        for i in 1:length(arr)
            if arr[i] == 1
                if current_streak == 0
                    current_start = i
                end
                current_streak += 1
                if current_streak > longest_streak
                    longest_streak = current_streak
                    longest_start = current_start
                    longest_stop = i
                end
            else
                current_streak = 0
            end
        end
        return longest_start, longest_stop
    end
    start, stop = find_longest_streak_indices(chrom[:,:weight].>(maximum(chrom[:,:weight])*0.01))
    #println("A $start $stop")
    if (stop - start) < 2
        return (0.0, 0)
    end
    #println("B")
    #f_width = min(11, (length(rt)÷2)*2 - 1)
    #if f_width < 4
    #    return (0.0, 0)
    #end
    intensity = Float64[]
    if (stop - start) > 5
        intensity = savitzky_golay(chrom[:,:weight][max(1, start - 1):min(stop + 1, end)], 5, 3).y
    else
        intensity = chrom[:,:weight][max(1, start - 1):min(stop + 1, end)]
    end

    if isplot == false
        return (integrate(rt, intensity, TrapezoidalFast()), (stop - start + 1))
    end
    plot!(rt, chrom[:,:weight], seriestype=:scatter, show = true)
    plot!(rt[max(1, start - 1):min(stop + 1, end)], intensity, show = true)
    return (integrate(rt, intensity, TrapezoidalFast()), (stop - start + 1))
end

integratePrecursor(test_integrate_gb, UInt32(4171604), isplot = true)
integratePrecursor(test_integrate_gb, UInt32(1648748), isplot = true)
plotChromatogram(test_integrate_gb, UInt32(801239))
plotChromatogram(test_integrate_gb, UInt32(1648748)) #potential interference at tail. 
plotChromatogram(test_integrate_gb, UInt32(2537365))

integratePrecursor(test_integrate_gb, UInt32( 527551))

integratePrecursor(test_integrate_gb, UInt32(3012890))
test_integrate_gb[(precursor_idx=0x008f3a26,)]
#=prec_mzs = zeros(Float32, length(prosit_precs))
for prec_bin in ProgressBar(prosit_index_intensities.precursor_bins)
    for prec in prec_bin.precs
        prec_mzs[getPrecID(prec)] = getPrecMZ(prec)
    end
end=#
transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(test_integrate_gb, psm[:precursor_idx])) => [:intensity, :count])
    
sort(best_psms[best_psms[:,:intensity].>0.0,[:decoy,:precursor_idx,:RT,:intensity,:count,:matched_ratio,:spectral_contrast_all]], :matched_ratio)
best_psms

histogram(log2.(best_psms[best_psms[:,:decoy],:intensity] .+ 1), normalize=:pdf, alpha = 0.5)
histogram!(log2.(max.(best_psms[best_psms[:,:decoy].==false,:intensity], 0.0) .+ 1), normalize=:pdf, alpha = 0.5)

transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(test_integrate_gb, psm[:precursor_idx])) => [:intensity, :count])