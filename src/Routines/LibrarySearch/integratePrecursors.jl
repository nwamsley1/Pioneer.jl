
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
        #println(typeof(transitions))
        #Match fragments to peaks
        fragmentMatches, fragmentMisses = matchPeaks(transitions, 
                                    spectrum[:masses], 
                                    spectrum[:intensities], 
                                    FragmentMatch{Float32},
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
        X, Hs, Hst, IDtoROW = buildDesignMatrix(fragmentMatches, fragmentMisses)
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

"""
https://terpconnect.umd.edu/~toh/spectrum/Differentiation.html#PeakDetection
Another common use of differentiation is in the detection of peaks in a signal. It's clear from the basic properties described in the previous section 
that the first derivative of a peak has a downward-going zero-crossing at the peak maximum, which can be used to locate the x-value of the peak, 
as shown on the right (script). If there is no noise in the signal, then any data point that has lower values on both sides of it will be a peak maximum. 
But there is always at least a little noise in real experimental signals, and that will cause many false zero-crossings simply due to the noise. 
To avoid this problem, one popular technique smooths the first derivative of the signal first, before looking for downward-going zero-crossings, 
and then takes only those zero crossings whose slope exceeds a certain predetermined minimum (called the "slope threshold") at a point where the original 
signal amplitude exceeds a certain minimum (called the "amplitude threshold"). By carefully adjusting the smooth width, slope threshold, and amplitude 
threshold, it is possible to detect only the desired peaks over a wide range of peak widths and ignore peaks that are too small, too wide, or too narrow. 
Moreover, because smoothing can distort peak signals, reducing peak heights, and increasing peak widths, this technique can be extended to measure the
 position, height, and width of each peak by least-squares curve-fitting of a segment of original unsmoothed signal near the top of the peak 
 (where the signal-to-noise ratio is usually the best). Thus, even if heavy smoothing is necessary to provide reliable discrimination against noise peaks, 
 the peak parameters extracted by curve fitting are not distorted and the effect of random noise in the signal is reduced by curve fitting over multiple
  data points in the peak. This technique has been implemented in Matlab/Octave and in spreadsheets.
"""
function integratePrecursor(chroms::GroupedDataFrame{DataFrame}, precursor_idx::UInt32; max_smoothing_window::Int = 15, min_smoothing_order::Int = 3, isplot::Bool = false)
    if !((precursor_idx=precursor_idx,) in keys(chroms)) #If the precursor is not found
        return (0.0, 0, 0.0, 0.0, missing, missing, missing)
    end

    #Chromatogram for the precursor. 
    #Has columns "weight" and "rt". 
    chrom = chroms[(precursor_idx=precursor_idx,)]
    return integratePrecursor(chrom, max_smoothing_window = max_smoothing_window, min_smoothing_order = min_smoothing_order, isplot = isplot)
end

function integratePrecursor(chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}; max_smoothing_window::Int = 15, min_smoothing_order::Int = 3, isplot::Bool = false)


    #Remove empty rows 
    non_zero = BitVector(undef, size(chrom)[1])
    fillNonZero!(non_zero, chrom[:,:weight])
    chrom = chrom[non_zero,:]
    #Too few points to attempt integration
    if size(chrom)[1] < 3
        return (0.0, 0, 0.0, 0.0, missing, missing, missing)
    end

    #Smoothing parameters for first derivative
    window_size, order = getSmoothingParams(size(chrom)[1], max_smoothing_window, min_smoothing_order)

    #Add zero intensity points to both ends of the chromatogram
    rt = pad(chrom[:,:rt], chrom[1,:rt] - (chrom[2,:rt] - chrom[1,:rt]), chrom[end,:rt] + (chrom[end,:rt] - chrom[end - 1,:rt]))
    intensity = pad(chrom[:,:weight], zero(eltype(typeof(chrom[:,:weight]))), zero(eltype(typeof(chrom[:,:weight]))))

    #Use smoothed first derivative crossings to identify peak apex and left/right boundaries
    best_peak_slope, start, stop, mid =  getIntegrationBounds(rt, intensity, window_size, order, 1.0/6.0, 5)
    if mid == 0
        return (0.0, 0, Float64(0.0), 0.0, missing, missing, missing)
    end 
    fwhm = getFWHM(rt, intensity, start, mid, stop)

    #No sufficiently wide peak detected. 
    if (stop - start) < 2
        return (0.0, 0, Float64(0.0), 0.0, missing, missing, missing)
    end

    #Smoothing parameters for chromatogram
    window_size, order = getSmoothingParams(stop - start, max_smoothing_window, min_smoothing_order)
    intensity_smooth = savitzky_golay(intensity[start:stop], window_size, order).y
    intensity_smooth[intensity_smooth .< 0.0] .= zero(eltype(typeof(intensity_smooth)))

    if isplot
        plot(rt, intensity, show = true, seriestype=:scatter)
        plot!(rt[start:stop], intensity_smooth, fillrange = [0.0 for x in 1:length(intensity_smooth)], alpha = 0.25, color = :grey, show = true);
        vline!([rt[start]], color = :red);
        vline!([rt[stop]], color = :red);
    end

    peak_area = NumericalIntegration.integrate(rt[start:stop], intensity[start:stop], TrapezoidalFast())
    count = (stop - start + 1)
    SN = Float64(sum(intensity[start:stop])/sum(intensity))
    error = sum(abs.(intensity_smooth .- intensity[start:stop]))/(length(intensity_smooth)*maximum(intensity))

    return peak_area, count, SN, best_peak_slope, error, rt[start:stop][argmax(intensity_smooth)], fwhm
end

function pad(x::Vector{T}, head::T, tail::T) where {T<:AbstractFloat}
    push!(x, tail)
    pushfirst!(x, head)
    return x
end

function getSmoothingParams(n_data_points::Int, max_window::Int, min_order::Int)
    window_size = min(max_window, max((n_data_points)÷3, 3))
    if isodd(window_size) == false
        window_size += 1
    end
    order = min(min_order, window_size - 2)
    return window_size, order
end

function fillNonZero!(non_zero::BitVector, intensity::Vector{T}) where {T<:AbstractFloat}
    for i in 1:(length(intensity))
        if iszero(intensity[i])&(iszero(intensity[min(i+1, end)]) == false)
            non_zero[i] = false
        else
            non_zero[i] = true
        end
    end
end

function getIntegrationBounds(RTs::Vector{T}, intensities::Vector{U}, window::Int, order::Int, min_width_t::Float64, min_width::Int) where {T,U<:AbstractFloat}
    #Get indices of first derivative zero crossings and the slope of the first derivative at each crossing
    zero_idx, zero_sign = getZeroCrossings(getSmoothDerivative(RTs, intensities, window, order), RTs[1:end - 1])

    return getPeakBounds(intensities, RTs, zero_idx, zero_sign, min_width_t, min_width)
end

function getSmoothDerivative(x::Vector{T}, y::Vector{U}, window::Int, order::Int) where {T,U<:AbstractFloat}
    #Smooth the first derivative 
    return savitzky_golay(diff(y)./diff(x), window, order).y
end

function getZeroCrossings(Y::Vector{T}, X::Vector{U}) where {T,U<:AbstractFloat}
    zero_crossings_idx = Vector{Int64}() #Indices of zero-crossings of first derivative
    zero_crossings_slope = Vector{T}() #Slopes of first derivative at the zero crossings
    for i in 1:(length(Y) - 1)
        if (sign(Y[i + 1])*sign(Y[i])) < 0 #First derivative crossed zero
            push!(zero_crossings_idx, i) #Add zero crossing
            slope = (Y[i + 1] - Y[i])/( X[i + 1] - X[i]) #slope at zero crossing
            push!(zero_crossings_slope, slope)
        end
    end
    return zero_crossings_idx, zero_crossings_slope
end

function getPeakBounds(intensity::Vector{T}, rt::Vector{T}, zero_crossings_1d::Vector{Int64}, zero_crossings_slope::Vector{U}, min_width_t::Float64 = 10.0, min_width::Int = 5) where {T,U<:AbstractFloat}

    best_peak_intensity = 0
    best_peak = 0
    best_peak_slope = 0
    best_left = 0
    best_right = 0
    N = length(intensity) 
    if iszero(length(zero_crossings_1d)) #If the first derivative does not cross zero, there is no peak
        return 0, 1, length(intensity), 0
    end
    for i in 1:length(zero_crossings_1d)
        if zero_crossings_slope[i] < 0 #Peaks will always occur at first derivative crossings where the slope is negative 
            if intensity[zero_crossings_1d[i]] > best_peak_intensity #Select the peak where the raw intensity is greatest. 

                #Left peak boundary is the previous first derivative crossing or the leftmost datapoint
                left_bound = i > 1 ? max(1, zero_crossings_1d[i-1] + 1) : 1

                #Right peak boundary is the next first derivative crossing or the rightmost datapoint
                if i+1<length(zero_crossings_1d)
                    right_bound = i < N ? min(N, zero_crossings_1d[i + 1]) : N
                else
                    right_bound = N
                end

                if (right_bound - left_bound) >=min_width #If peak width (in data points) exceeds minimum
                    if (rt[right_bound] - rt[left_bound]) >= min_width_t #If peak width (in minutes) exceeds minimum
                        best_peak_intensity = intensity[zero_crossings_1d[i]]
                        best_peak = zero_crossings_1d[i]
                        best_peak_slope = abs(zero_crossings_slope[i])
                        best_left = left_bound
                        best_right = right_bound
                    end
                end
            end
        end
    end

    return best_peak_slope, best_left, best_right, best_peak
end

function getFWHM(X::Vector{T}, Y::Vector{U}, start::Int, mid::Int, stop::Int) where {T,U<:AbstractFloat}
    function getMid(y, x_start, x_stop, y_start, y_stop)
        slope = (x_start - x_stop)/(y_start - y_stop)
        return (y - y_start)*slope + x_start
    end

    left_x, right_x = zero(T), zero(T)
    if (mid == start) | (mid == stop)
        return zero(T)
    end
    halfmax = Y[mid]/2
    for i in start:(mid - 1)
        if (Y[i + 1]>halfmax) & (Y[i]<halfmax)
            left_x = getMid(halfmax, X[i], X[i + 1], Y[i], Y[i + 1])
        end
    end

    for i in stop:-1:(mid + 1)
        if (Y[i - 1]>halfmax) & (Y[i]<halfmax)
            right_x = getMid(halfmax, X[i-1], X[i], Y[i-1], Y[i])
        end
    end

    return right_x - left_x

end
#=-
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
=#
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

integratePrecursor(chroms, UInt32(1284418), isplot = true)
integratePrecursor(chroms, UInt32(7479714), isplot = true)
integratePrecursor(chroms, UInt32(4343373), isplot = true)
integratePrecursor(chroms, UInt32(232427), isplot = true)
integratePrecursor(chroms, UInt32(4034269), isplot = true)


integratePrecursor(chroms, UInt32( 8316025), isplot = true)
integratePrecursor(chroms, UInt32( 8316019), isplot = true)

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


#4847815 
integratePrecursor(chroms, UInt32(4847815), isplot = true) #decoy
integratePrecursor(chroms, UInt32(7632773), isplot = true) #decoy
integratePrecursor(chroms, UInt32(5913935), isplot = true) #decoy
integratePrecursor(chroms, UInt32(5953390), isplot = true) #decoy
integratePrecursor(chroms, UInt32(6266705 ), isplot = true) #decoy

integratePrecursor(chroms, UInt32(3355405), isplot = true) #bad target
integratePrecursor(chroms, UInt32( 4016479  ), isplot = true)
integratePrecursor(chroms, UInt32(  1807561 ), isplot = true)
integratePrecursor(chroms, UInt32(4097143  ), isplot = true)

integratePrecursor(chroms, UInt32(3191213), isplot = false)

=#



