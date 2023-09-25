
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
function integratePrecursor(chroms::GroupedDataFrame{DataFrame}, precursor_idx::UInt32, p0::NTuple{5, T}; max_smoothing_window::Int = 15, min_smoothing_order::Int = 3, min_scans::Int = 5, min_width::AbstractFloat = 1.0/6.0, integration_width::AbstractFloat = 4.0, integration_points::Int = 1000, isplot::Bool = false) where {T<:AbstractFloat}

    if !((precursor_idx=precursor_idx,) in keys(chroms)) #If the precursor is not found
        return (missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing)
    end

    #Chromatogram for the precursor. 
    #Has columns "weight" and "rt". 
    chrom = chroms[(precursor_idx=precursor_idx,)]
    return integratePrecursor(chrom, p0, eltype(chrom[:,:weight]), max_smoothing_window = max_smoothing_window, min_smoothing_order = min_smoothing_order, min_scans = min_scans, min_width = min_width, integration_width = integration_width, integration_points = integration_points, isplot = isplot)
end

function integratePrecursor(chrom::SubDataFrame{DataFrame, DataFrames.Index, Vector{Int64}}, p0::NTuple{5, U}, T::Type; max_smoothing_window::Int = 15, min_smoothing_order::Int = 3, min_scans::Int = 5, min_width::AbstractFloat = 1.0/6.0, integration_width::AbstractFloat = 4.0, integration_points::Int = 1000, isplot::Bool = false) where {U<:AbstractFloat}

    non_zero = BitVector(undef, size(chrom)[1])
    fillNonZero!(non_zero, chrom[:,:weight])

    #Same precursor may be isolated multiple times within a single cycle_idx
    #Retain only most abundant within each cycle .
    for i in range(1, size(chrom)[1]-1)
        if chrom[i,:cycle_idx] == chrom[i+1,:cycle_idx]
            #May not be appropriate criterion to choose between competing scans
            if chrom[i,:weight] > chrom[i + 1,:weight]
                non_zero[i + 1] = false
            else
                non_zero[i] = false
            end
        end
    end
    chrom = chrom[non_zero,:]
    #display(chrom)
    #Scan with the highest score 
    #Too few points to attempt integration
    if size(chrom)[1] < 3
        return (missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing)
    end
    #########
    #
    rt = chrom[:,:rt]
    intensity = collect(chrom[:,:weight])
    intensity[chrom[:,:rank].<2].=0.0
    frag_counts = collect(chrom[:,:frag_count])

    #########
    #Get Boundaries for fitting EGH curve 
    #Smoothing parameters for first derivative
    window_size, order = getSmoothingParams(length(rt), max_smoothing_window, min_smoothing_order)

    #Use smoothed first derivative crossings to identify peak apex and left/right boundaries
    start, stop = getIntegrationBounds(frag_counts, rt, intensity, window_size, order, min_width, min_scans)

    ########
    #Fit EGH Curve 
    EGH_FIT = nothing
    try
        EGH_FIT = curve_fit(EGH, 
                            JEGH, 
                            rt[start:stop],
                            intensity[start:stop],
                            #rt[start:stop],
                            #intensity[start:stop],
                            getP0(p0))
    catch
        return (missing, missing, missing, missing, missing, missing, missing, missing, missing, missing, missing)
    end

    peak_area = Integrate(EGH, EGH_FIT.param, (T(EGH_FIT.param[2] - integration_width/2), 
                                               T(EGH_FIT.param[2] + integration_width/2)),
                                               integration_points)
    ##########
    #Plots                                           
    if isplot
        #Plot Data
        Plots.plot(rt, intensity, show = true, seriestype=:scatter)
        #Plots Fitted EGH Curve
        X = LinRange(T(EGH_FIT.param[2] - 1.0), T(EGH_FIT.param[2] + 1.0), integration_points)
        Plots.plot!(X,  
                    EGH(T.(collect(X)), EGH_FIT.param), 
                    fillrange = [0.0 for x in 1:integration_points], 
                    alpha = 0.25, color = :grey, show = true
                    ); 
        Plots.vline!([rt[start]], color = :red);
        Plots.vline!([rt[stop]], color = :red);
    end


    ############
    #Calculate Features
    MODEL_INTENSITY = EGH(rt, EGH_FIT.param)
    RESID = abs.(intensity .- MODEL_INTENSITY)
    GOF_IDX = (MODEL_INTENSITY.>(EGH_FIT.param[end]*0.1)) .& (intensity.>(EGH_FIT.param[end]*0.1))
    GOF_IDX = (MODEL_INTENSITY.>(EGH_FIT.param[end]*0.0)) .& (intensity.>(EGH_FIT.param[end]*0.0))

    var = std(RESID)/EGH_FIT.param[4]
    GOF = 1 - sum(RESID[GOF_IDX])/sum(intensity[GOF_IDX])
    FWHM = getFWHM(0.5, EGH_FIT.param[3], EGH_FIT.param[1])
    FWHM_01 = getFWHM(0.01, EGH_FIT.param[3], EGH_FIT.param[1])
    asymmetry = atan(abs(EGH_FIT.param[3])/sqrt(abs(EGH_FIT.param[1]/2)))
    points_above_FWHM = sum(intensity[start:stop].>(EGH_FIT.param[end]*0.5))
    points_above_FWHM_01 = sum(intensity[start:stop].>(EGH_FIT.param[end]*0.01))

    σ = EGH_FIT.param[1]
    tᵣ = EGH_FIT.param[2]
    τ = EGH_FIT.param[3]
    H = EGH_FIT.param[4]
    return peak_area, GOF, FWHM, FWHM_01, asymmetry, points_above_FWHM, points_above_FWHM_01, σ, tᵣ, τ, H
end

function pad(x::Vector{T}, head::T, tail::T, n::Int = 1) where {T<:Real}
    for i in range(1, n)
        push!(x, tail)
        pushfirst!(x, head)
    end
    return x
end

function getSmoothingParams(n_data_points::Int, max_window::Int, min_order::Int)
    #window_size = n_data_points - min_order#min(max_window, max((n_data_points)÷3, 3))
    window_size = min(max_window, max((n_data_points)÷3, 3))
    #window_size = 7
    if isodd(window_size) == false
        window_size += 1
    end
    order = min(min_order, window_size - 2)
    return window_size, order
end

function fillNonZero!(non_zero::BitVector, intensity::Vector{T}) where {T<:AbstractFloat}
    for i in 1:(length(intensity))
        if iszero(intensity[i])&(iszero(intensity[min(i+1, end)])==false)
            non_zero[i] = false
        else
            non_zero[i] = true
        end
    end
end

function getIntegrationBounds(frag_counts::Vector{Int64}, RTs::Vector{T}, intensities::Vector{U}, window::Int, order::Int, min_width_t::Float64, min_width::Int) where {T,U<:AbstractFloat}
    #Get indices of first derivative zero crossings and the slope of the first derivative at each crossing
    zero_idx, zero_sign = getZeroCrossings(getSmoothDerivative(RTs, intensities, window, order), RTs[1:end - 1])
    return getPeakBounds(frag_counts, intensities, RTs, zero_idx, zero_sign, min_width_t, min_width)
end

function getSmoothDerivative(x::Vector{T}, y::Vector{U}, window::Int, order::Int) where {T,U<:AbstractFloat}
    smooth_y = max.(savitzky_golay(y, window, order).y, zero(T))
    return savitzky_golay(diff(smooth_y)./diff(x), window, order).y
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

function getPeakBounds(frag_counts::Vector{Int64}, intensity::Vector{T}, rt::Vector{T}, zero_crossings_1d::Vector{Int64}, zero_crossings_slope::Vector{U}, min_width_t::Float64 = 10.0, min_width::Int = 5) where {T,U<:AbstractFloat}

    best_peak_intensity = -Inf
    best_left = 1
    best_right = length(rt)
    N = length(intensity)
    
    for i in 1:length(zero_crossings_1d)
        if zero_crossings_slope[i] < 0 #Peaks will always occur at first derivative crossings where the slope is negative 

            if intensity[zero_crossings_1d[i]]*frag_counts[zero_crossings_1d[i]] > best_peak_intensity
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
                        #best_peak_dist = abs((rt[zero_crossings_1d[i]])-(rt[best_scan_idx]))#intensity[zero_crossings_1d[i]]
                        best_peak_intensity = intensity[zero_crossings_1d[i]]*frag_counts[zero_crossings_1d[i]]
                        best_left = max(1, left_bound)#max(left_bound - 2, 1)
                        best_right = min(right_bound, length(intensity))# min(right_bound + 2, length(intensity))
                    end
                end
            end
        end
    end

    if best_left == best_right
        return 1, length(rt)
    end
    return best_left, best_right
end

#=function getPeakBounds(intensity::Vector{T}, rt::Vector{T}, zero_crossings_1d::Vector{Int64}, zero_crossings_slope::Vector{U}, min_width_t::Float64 = 10.0, min_width::Int = 5) where {T,U<:AbstractFloat}

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

    #Roll back boundaries
    while intensity[best_right - 1] < intensity[best_right]
        best_right = best_right - 1
        if best_right == best_left
            break
        end
    end

    while intensity[best_left + 1] < intensity[best_left]
        best_left += 1
        if best_right == best_left
            break
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
=#
#=

function integrateRAW(
                    spectra::Arrow.Table, 
                    rt_index::retentionTimeIndex{T, U},
                    fragment_list::Vector{Vector{LibraryFragment{Float32}}},
                    ms_file_idx::UInt32;
                    frag_ppm_err::Float64 = 0.0,
                    fragment_tolerance::Float64 = 40.0,
                    quadrupole_isolation_width::Float64 = 8.5,
                    max_peak_width::Float64 = 2.0,
                    λ::Float32 = Float32(2e12),
                    γ::Float32 = Float32(1/2),
                    max_iter::Int = 1000,
                    nmf_tol::Float32 = Float32(100.0),
                    scan_range::Tuple{Int64, Int64} = (0, 0),
                    mz_range::Tuple{AbstractFloat, AbstractFloat} = (0.0, 2000.0)
                    ) where {T,U<:Real}
    
    ms2 = 0
    nmf = Dict(:precursor_idx => UInt32[], :weight => Float32[], :rt => Float32[], :frag_count => Int64[])
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
        if (spectrum[:precursorMZ] < first(mz_range)) | (spectrum[:precursorMZ] > last(mz_range))
           continue
        end
    
        transitions, prec_ids = selectTransitions(fragment_list, rt_index, Float64(spectrum[:retentionTime]), max_peak_width/2.0, spectrum[:precursorMZ], Float32(quadrupole_isolation_width/2.0))
  
        fragmentMatches, fragmentMisses = matchPeaks(transitions, 
                                    spectrum[:masses], 
                                    spectrum[:intensities], 
                                    FragmentMatch{Float32},
                                    count_unmatched =true,
                                    δs = [frag_ppm_err],
                                    scan_idx = UInt32(i),
                                    ms_file_idx = ms_file_idx,
                                    min_intensity = zero(Float32),
                                    ppm = fragment_tolerance
                                    )

        frag_counts = counter(UInt32)

        for match in fragmentMatches
            DataStructures.inc!(frag_counts, match.prec_id)
        end

        #filter!(x->(frag_counts[x.prec_id].>2), fragmentMatches)
        #filter!(x->(frag_counts[x.prec_id].>2), fragmentMisses)
        #filter!(x->x.frag_index.>2, fragmentMatches)
        #filter!(x->x.frag_index.>2, fragmentMisses)
        if !iszero(length(fragmentMatches))
            X, Hs, Hst, IDtoROW = buildDesignMatrix(fragmentMatches, fragmentMisses)
            weights = sparseNMF(Hst, Hs, X; λ=λ,γ=γ, max_iter=max_iter, tol=nmf_tol)
        else
            IDtoROW = UnorderedDictionary{UInt32, UInt32}()
        end
        for key in prec_ids
            if haskey(frag_counts, key)

                if haskey(IDtoROW, key)
            #for key in keys(IDtoROW)
                    push!(nmf[:precursor_idx], key)
                    push!(nmf[:weight], weights[IDtoROW[key]])
                    push!(nmf[:rt], spectrum[:retentionTime])
                    push!(nmf[:frag_count], frag_counts[key])
                else
                    push!(nmf[:precursor_idx], key)
                    push!(nmf[:weight], Float32(0.0))
                    push!(nmf[:rt], spectrum[:retentionTime])
                    push!(nmf[:frag_count], frag_counts[key])
                end
            else
                push!(nmf[:precursor_idx], key)
                push!(nmf[:weight], Float32(0.0))
                push!(nmf[:rt], spectrum[:retentionTime])
                push!(nmf[:frag_count], 0)
            end
        end
    end
    nmf = DataFrame(nmf)
    sort!(nmf, [:precursor_idx,:rt]);
    return groupby(nmf, :precursor_idx)
end

=#
