
function integrateRAW(
                    spectra::Arrow.Table, 
                    rt_index::retentionTimeIndex{T, U},
                    fragment_list::Vector{Vector{LibraryFragment{Float32}}},
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
    
    ms2 = 0
    nmf = Dict(:precursor_idx => UInt32[], :weight => Float32[], :rt => Float32[])
    for (i, spectrum) in ProgressBar(enumerate(Tables.namedtupleiterator(spectra)))

        if spectrum[:msOrder] == 1
            continue
        else
            ms2 += 1
        end
        if scan_range != (0, 0)
            i < first(scan_range) ? continue : nothing
            i > last(scan_range) ? continue : nothing
        end
        #Get peptides that could be in the spectra
        transitions = selectTransitions(fragment_list, rt_index, Float64(spectrum[:retentionTime]), max_peak_width/2.0, spectrum[:precursorMZ], Float32(quadrupole_isolation_width/2.0))

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

        #Build templates for regrssion
        X, Hs, Hst, IDtoROW = buildDesignMatrix(fragmentMatches, fragmentMisses)
        #println(size(X))
        #println(Hs.n)
        #println(Hs.m)
        #Fit non-negative adaptive LASSO
        #λ = mean(X)*(Hs.n^2)
        #println(λ)
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
end
function integratePrecursor(chroms::GroupedDataFrame{DataFrame}, precursor_idx::Int64; isplot::Bool = false)
    if !((precursor_idx=precursor_idx,) in keys(chroms))
        return (0.0, 0, 0.0)
    end
    chrom = chroms[(precursor_idx=precursor_idx,)]
    rt = chrom[:,:rt]
    #println(precursor_idx)
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
        return (0.0, 0, Float64(0.0))
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
        return (Float64(integrate(rt, intensity, TrapezoidalFast())), (stop - start + 1), Float64(sum(chrom[start:stop,:weight])/sum(chrom[:,:weight])))
    end
    plot(rt, chrom[:,:weight], seriestype=:scatter, show = true)
    plot!(rt[max(1, start - 1):min(stop + 1, end)], intensity, show = true)
    return (Float64(integrate(rt, intensity, TrapezoidalFast())), (stop - start + 1), Float64(sum(chrom[start:stop,:weight])/sum(chrom[:,:weight])))
end

#=
chroms = integrateRAW(MS_TABLE, rt_index, prosit_detailed, one(UInt32), fragment_tolerance=15.6, λ=Float32(5e3) , γ=Float32(0), max_peak_width = 2.0, scan_range = (0, 300000));
transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(chroms, psm[:precursor_idx], isplot = false)) => [:intensity, :count, :SN]);
non_zero = best_psms[(best_psms[:,:intensity].>0).&(best_psms[:,:count].>5),:];
println("Signal to Noise ", mean(non_zero[:,:SN]))
@time getQvalues!(non_zero, non_zero[:,:prob], non_zero[:,:decoy]);
println("TARGETS at 10%FDR Intensity", sum(non_zero[non_zero[:,:decoy].==false,:intensity])*mean(non_zero[:,:SN]))
println("TARGETS at 10%FDR ", sum(non_zero[:,:q_values].<=0.1))
println("TARGETS at 1%FDR ", sum(non_zero[:,:q_values].<=0.01))
println("TARGETS at 1%FDR intensity ", sum(non_zero[:,:q_values].<=0.01)*mean(non_zero[:,:SN]))
println("TARGET/DECOY intensity ", sum(non_zero[non_zero[:,:decoy].==false,:][:,:intensity])/sum(non_zero[non_zero[:,:decoy].==true,:][:,:intensity]))
println("TARGET/DECOY count ", sum(non_zero[non_zero[:,:decoy].==false,:][:,:count])/sum(non_zero[non_zero[:,:decoy].==true,:][:,:count]))


integratePrecursor(chroms, UInt32(4171604), isplot = true)
integratePrecursor(chroms, UInt32(1648748), isplot = true)
integratePrecursor(chroms, UInt32(801239), isplot = true)
integratePrecursor(chroms, UInt32(679408), isplot = true)
integratePrecursor(chroms, UInt32(2537365), isplot = true)
integratePrecursor(chroms, UInt32(2257951), isplot = true)
integratePrecursor(chroms, UInt32( 1367150 ), isplot = true)

integratePrecursor(chroms, UInt32(  4259798 ), isplot = true)

integratePrecursor(chroms, UInt32(   3574807 ), isplot = true)
integratePrecursor(chroms, UInt32(  2642152), isplot = true)
integratePrecursor(chroms, UInt32(   508178 ), isplot = true)
=#



