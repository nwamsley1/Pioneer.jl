
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

function integratePrecursor(chroms::GroupedDataFrame{DataFrame}, precursor_idx::UInt32; isplot::Bool = false)
    if !((precursor_idx=precursor_idx,) in keys(chroms))
        return (0.0, 0, 0.0)
    end

    chrom = chroms[(precursor_idx=precursor_idx,)]

    non_zero = BitVector(undef, size(chrom)[1])
    for i in 1:(size(chrom)[1] - 1)
        if iszero(chrom[:,:weight][i])&(iszero(chrom[:,:weight][i+1]).==false)
            non_zero[i] = false
        else
            non_zero[i] = true
        end
    end

    chrom = chrom[non_zero,:]
    rt = chrom[:,:rt]
    intensity = chrom[:,:weight]
    #start, stop = find_longest_streak_indices(chrom[:,:weight].>(maximum(chrom[:,:weight])*0.01))
    best, start, stop =  getIntegrationBounds(chrom, 9, 3, 1.0/6.0, 5)
    #plot(rt, chrom[:,:weight], seriestype=:scatter, show = true)
    return 
    if (stop - start) < 2
        return (0.0, 0, Float64(0.0))
    end

    if isplot == false
        return (Float64(integrate(rt, intensity, TrapezoidalFast())), (stop - start + 1), Float64(sum(chrom[start:stop,:weight])/sum(chrom[:,:weight])))
    end
    println("A")
    println(integrate(chrom[:,:rt][start:stop], chrom[:,:weight][start:stop], TrapezoidalFast()))
    #plot(rt, chrom[:,:weight], seriestype=:scatter, show = true);
    #plot!(rt[max(1, start - 1):min(stop + 1, end)], intensity, show = true)
    return (Float64(integrate(rt, intensity, TrapezoidalFast())), (stop - start + 1), Float64(sum(chrom[start:stop,:weight])/sum(chrom[:,:weight])))
end


savitzky_golay(chrom[:,:weight][max(1, start - 1):min(stop + 1, end)], 5, 3).y

first_d = diff(test[:,:weight])./diff(test[:,:rt])
first_d_smooth = savitzky_golay(first_d, 9, 3).y
#plot(test[:,:rt][2:end], diff(test[:,:weight])./diff(test[:,:rt]), seriestype=:scatter)
plot(test[:,:rt][1:end-1], first_d_smooth, seriestype=:scatter)
plot!(test[:,:rt], test[:,:weight]*5, seriestype=:scatter)

function getSmooth()
end

function getIntegrationBounds(chrom::DataFrame, window::Int, order::Int, min_width_t::Float64, min_width::Int)
    zero_idx, zero_sign = getZeroCrossings(getSmoothDerivative(chrom[:,:rt], chrom[:,:weight], window, order), chrom[:,:rt][1:end - 1])
    return getPeakBounds(chrom[:,:weight][1:end - 1], chrom[:,:rt][1:end - 1], zero_idx, zero_sign, min_width_t, min_width)
end

function getSmoothDerivative(x::Vector{T}, y::Vector{U}, window::Int, order::Int) where {T,U<:AbstractFloat}
    return savitzky_golay(diff(y)./diff(x), window, order).y
end

function getZeroCrossings(series_y::Vector{T}, series_x::Vector{U}) where {T,U<:AbstractFloat}
    zero_crossings_idx = Vector{Int64}()
    zero_crossings_sign = Vector{T}()
    for i in 1:(length(series_y) - 1)
        if (sign(series_y[i + 1])*sign(series_y[i])) < 0
            push!(zero_crossings_idx, i)
            slope = (series_y[i + 1] - series_y[i])/(series_x[i + 1] - series_x[i])
            push!(zero_crossings_sign, sign(slope))
        end
    end
    return zero_crossings_idx, zero_crossings_sign
end

function getPeakBounds(series_y::Vector{T}, series_x::Vector{T}, zero_crossings_1d::Vector{Int64}, zero_crossings_sign::Vector{U}, min_width_t::Float64 = 10.0, min_width::Int = 5) where {T,U<:AbstractFloat}

    best_peak_intensity = 0
    best_peak = 0
    best_left = 0
    best_right = 0
    N = length(series_y) 
    for i in 1:length(zero_crossings_1d)
        if zero_crossings_sign[i] < 0
            if series_y[zero_crossings_1d[i]] > best_peak_intensity

                left_bound = i > 1 ? max(1, zero_crossings_1d[i - 1]) : 1
                right_bound = i < N ? min(N, zero_crossings_1d[i + 1]) : N
                if (series_x[right_bound] - series_x[left_bound]) >= min_width_t
                    if (right_bound - left_bound) >=min_width
                        best_peak_intensity = series_y[zero_crossings_1d[i]]
                        best_peak = zero_crossings_1d[i]
                        best_left = left_bound
                        best_right = right_bound
                    end
                end
            end
        end
    end
    return best_peak, best_left, best_right
end


function integratePrecursor(chroms::GroupedDataFrame{DataFrame}, precursor_idx::UInt32; isplot::Bool = false)
    if !((precursor_idx=precursor_idx,) in keys(chroms))
        return (0.0, 0, 0.0)
    end
    chrom = chroms[(precursor_idx=precursor_idx,)]
    chrom = chrom[chrom[:,:weight].!=0.0,:]
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
    #start, stop = find_longest_streak_indices(chrom[:,:weight].>(maximum(chrom[:,:weight])*0.01))
    start, stop = 1, size(chrom)[1]
    #plot(rt, chrom[:,:weight], seriestype=:scatter, show = true)
    #return 
    #if (stop - start) < 2
    #    return (0.0, 0, Float64(0.0))
    #end
    #println("B")
    #f_width = min(11, (length(rt)÷2)*2 - 1)
    #if f_width < 4
    #    return (0.0, 0)
    #end
    
    intensity = chrom[:,:weight]
    #if (stop - start) > 5
    #    intensity = savitzky_golay(chrom[:,:weight][max(1, start - 1):min(stop + 1, end)], 5, 3).y
    #else
    #    intensity = chrom[:,:weight][max(1, start - 1):min(stop + 1, end)]
    #end

    if isplot == false
        return (Float64(integrate(rt, intensity, TrapezoidalFast())), (stop - start + 1), Float64(sum(chrom[start:stop,:weight])/sum(chrom[:,:weight])))
    end
    plot(rt, chrom[:,:weight], seriestype=:scatter, show = true)
    plot!(rt[max(1, start - 1):min(stop + 1, end)], intensity, show = true)
    return (Float64(integrate(rt, intensity, TrapezoidalFast())), (stop - start + 1), Float64(sum(chrom[start:stop,:weight])/sum(chrom[:,:weight])))
end

integratePrecursor(chroms, UInt32(   3310087 ), isplot = false)
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


integratePrecursor(chroms, UInt32(3191213), isplot = false)


