abstract type SpectralScores{T<:AbstractFloat} end

struct SpectralScoresComplex{T<:AbstractFloat} <: SpectralScores{T}
    spectral_contrast::T
    fitted_spectral_contrast::T
    gof::T
    max_matched_residual::T
    max_unmatched_residual::T
    fitted_manhattan_distance::T
    matched_ratio::T
    #entropy_score::T
end

struct SpectralScoresSimple{T<:AbstractFloat} <: SpectralScores{T} 
    scribe::T
    city_block::T
    spectral_contrast::T
    matched_ratio::T
    entropy_score::T
    percent_theoretical_ignored::T
end

function getDistanceMetrics(H::SparseArray{Ti,T}, 
    spectral_scores::Vector{SpectralScoresSimple{U}};
    relative_improvement_threshold::Float32 = 1.25f0,
    min_frags::Int64 = 5) where {Ti<:Integer,T,U<:AbstractFloat}

    @inbounds @fastmath for col in range(1, H.n)

        # Gather the indices of relevant peaks in this spectrum. 
        # We'll iteratively drop one at a time if they have too much interference
        included = Int[]  # We'll store the indices from H.colptr[col] : H.colptr[col+1]-1
        tot_pred_signal = 0.0
        for i in H.colptr[col]:(H.colptr[col+1]-1)
            # Skip if it's an isotope or something we don't consider
            if !iszero(H.isotope[i])
                continue
            end
            push!(included, i)
            tot_pred_signal += H.nzval[i]
        end

        scribe_best, city_best, cosine_similarity_best, matched_ratio_best, ent_best, next_worst_pos, next_worst_pred_signal, num_matching_peaks = computeMetricsFor(H, col, included)
        
        percent_theoretical_ignored = 0.0f0
        while (num_matching_peaks > min_frags) && (next_worst_pos > 0)
            deleteat!(included, next_worst_pos)
            scribe, city, cosine_similarity, matched_ratio, ent, worst_pos, worst_pred_signal, num_matching_peaks = computeMetricsFor(H, col, included)

            # If ignoring the worst peak doesn't increase the scribe score enough, then we're done
            if (scribe < (scribe_best * relative_improvement_threshold))
                break
            end

            # update scores
            percent_theoretical_ignored += next_worst_pred_signal

            scribe_best = scribe
            city_best = city
            cosine_similarity_best = cosine_similarity
            matched_ratio_best = matched_ratio
            ent_best = ent
            next_worst_pos = worst_pos
            next_worst_pred_signal = worst_pred_signal
        end

        spectral_scores[col] = SpectralScoresSimple(
                    Float16(scribe_best),
                    Float16(city_best),
                    Float16(cosine_similarity_best),
                    Float16(matched_ratio_best),
                    Float16(ent_best),
                    Float16(percent_theoretical_ignored / tot_pred_signal)
                )

    end
end


function computeMetricsFor(H::SparseArray{Ti,T}, col, included_indices) where {Ti<:Integer,T<:AbstractFloat}
    # We'll accumulate partial sums and compute the same metrics as your snippet.
    # (For clarity, we skip inlining optimizations like @inbounds, @fastmath here.)

    # We'll also find the worst match, i.e. the position with max difference between (nzval[i], x[i])
    worst_val = -Inf
    worst_pos = 0
    worst_idx = 0
    num_matching_peaks = 0

    # Sums for numerator/denominator
    # We also want the "normalized" predicted and observed, so we can correctly find the worst inteferring peak
    total_h = zero(T)
    total_x = zero(T)

    @inbounds @fastmath for i in included_indices
        total_h += H.nzval[i]
        total_x += H.x[i]
        if H.x[i] > 0
            num_matching_peaks += 1
        end
    end

    # Protect against zero total intensities
    total_h = max(total_h, eps()) 
    total_x = max(total_x, eps())

    h_sqrt_sum = zero(T)
    x_sqrt_sum = zero(T)

    h2_norm = zero(T)
    x2_norm = zero(T)

    h2_norm_v2 = zero(T)
    x2_norm_v2 = zero(T)

    dot_product = zero(T)
    scribe_score = zero(T)
    city_block_dist = zero(T)

    matched_sum = zero(T)
    unmatched_sum = zero(T)

    N = 0
    @inbounds @fastmath for (local_pos, i) in enumerate(included_indices)
        if iszero(H.isotope[i])==false
            continue
        end

        h_sqrt_sum += sqrt(H.nzval[i])
        x_sqrt_sum += sqrt(H.x[i])#/Xsum

        h2_norm += H.nzval[i]^2
        x2_norm += H.x[i]^2

        if H.matched[i]
            matched_sum += H.nzval[i]
        else
            unmatched_sum += H.nzval[i]
        end

        # "normalized" predicted and observed, so we can know which peak is the worst for spectral angle
        h_val_v2 = H.nzval[i] / total_h
        x_val_v2 = H.x[i] / total_x

        h2_norm_v2 += h_val_v2^2
        x2_norm_v2 += x_val_v2^2
        
        dot_product += h_val_v2 * x_val_v2

        diff = x_val_v2 - h_val_v2 # only look for positive difference because it implies there's interference
        if (diff > worst_val) && (x_val_v2 > 0)
            worst_val = diff
            worst_pos = local_pos
            worst_idx = i
        end

        N += 1
    end
        
    #Sqrt of sum of squares
    h2_norm = sqrt(h2_norm)
    x2_norm = sqrt(x2_norm)

    h2_norm_v2 = sqrt(h2_norm_v2)
    x2_norm_v2 = sqrt(x2_norm_v2)

    @inbounds @fastmath for i in included_indices
        scribe_score += ((sqrt(H.nzval[i])/h_sqrt_sum) - (sqrt(H.x[i])/x_sqrt_sum))^2

        city_block_dist += abs(
            (H.nzval[i]/h2_norm) - (H.x[i]/h2_norm)
        )   
    end

    
    scribe_score = -log2(scribe_score / N)
    city_block_dist   = -log2(city_block_dist / N)
    cosine_similarity = dot_product / (h2_norm_v2 * x2_norm_v2)
    unmatched_sum = max(unmatched_sum, eps())
    matched_ratio = log2(matched_sum / unmatched_sum)
    ent_val = -1.0 * getEntropy(H, col, included_indices)
    worst_intensity_ignored = worst_idx > 0 ? H.nzval[worst_idx] : 0.0

    return (scribe_score, city_block_dist, cosine_similarity, matched_ratio, ent_val, worst_pos, worst_intensity_ignored, num_matching_peaks)
end

function getDistanceMetrics(w::Vector{T}, r::Vector{T}, H::SparseArray{Ti,T}, spectral_scores::Vector{SpectralScoresComplex{U}}) where {Ti<:Integer,T,U<:AbstractFloat}


    ########
    #Get residuals
    #########
    @turbo for i in range(1, H.m)
        r[i] = zero(T)
    end

    for n in range(1, H.n_vals)
        if iszero(r[H.rowval[n]])
            r[H.rowval[n]] = -H.x[n]
        end
    end

    for col in range(1, H.n)
        start = H.colptr[col]
        stop = H.colptr[col+1] - 1
        for n in start:stop
            r[H.rowval[n]] += w[col]*H.nzval[n]
        end
    end

    for col in range(1, H.n)
        
        h2_sum = zero(T)
        x2_sum = zero(T)
        x_sum = zero(T)
        dot_product = zero(T)
        fitted_dotp = zero(T)
        matched_sum = zero(T)
        unmatched_sum = zero(T)    

        manhattan_distance = zero(T)
        max_matched_residual = zero(T)
        max_unmatched_residual = zero(T)
        sum_of_residuals = zero(T)
        fitted_dotp = zero(T)        
        fitted_dotp_norm1 = zero(T)
        sum_of_fitted_peaks_matched = zero(T)
        sum_of_fitted_peaks_unmatched = zero(T)

        @inbounds @fastmath for i in range(H.colptr[col], H.colptr[col + 1]-1)

            #Fitted Manhattan Distance
            x_sum += H.x[i]
            manhattan_distance += abs(w[col]*H.nzval[i] - H.x[i])

            #Normalized Dot Product 
            dot_product += H.nzval[i]*H.x[i]
            x2_sum += (H.x[i])^2
            h2_sum += (H.nzval[i])^2 
    

            sum_of_residuals += abs(r[H.rowval[i]])

            fitted_peak = w[col]*H.nzval[i]
            r_abs = abs(r[H.rowval[i]])
            sum_of_residuals += r_abs

            if H.matched[i]
                matched_sum += H.nzval[i]
                fitted_dotp += sqrt(fitted_peak + r_abs)*sqrt(fitted_peak)
                fitted_dotp_norm1 += fitted_peak + r_abs
                sum_of_fitted_peaks_matched += fitted_peak
                if r_abs > max_matched_residual
                    max_matched_residual = r_abs
                end
            else
                unmatched_sum += H.nzval[i]
                sum_of_fitted_peaks_unmatched += fitted_peak
                if r_abs > max_unmatched_residual
                    max_unmatched_residual = r_abs
                end
            end
        end
        sum_of_fitted_peaks =  sum_of_fitted_peaks_matched +  sum_of_fitted_peaks_unmatched

        fitted_dotp_norm = fitted_dotp/(sqrt(fitted_dotp_norm1)*sqrt(sum_of_fitted_peaks))
        fitted_spectral_contrast = fitted_dotp_norm#1 - 2*acos(fitted_dotp_norm)/π
        dot_product_norm = dot_product/(sqrt(h2_sum)*sqrt(x2_sum))
        spectral_contrast = dot_product_norm#1 - 2*acos(dot_product_norm)/π

        spectral_scores[col] = SpectralScoresComplex(
            Float16(spectral_contrast), #spectral_contrast
            Float16(fitted_spectral_contrast), #fitted_spectral_contrast
            Float16(-log2(sum_of_residuals/sum_of_fitted_peaks)), #gof
            Float16(-log2(max_matched_residual/sum_of_fitted_peaks_matched)),#max_matched_residual
            Float16(-log2(max_unmatched_residual/sum_of_fitted_peaks + 1e-10)), #max_unmatched_residual
            Float16(-log2(manhattan_distance/x_sum)), #fitted_manhattan_distance
            Float16(log2(matched_sum/unmatched_sum)), #matched_ratio
            #Float16(-1.0*getEntropy(H, r, col)) #entropy
        )
    end
end
    
function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}
    #println("col $col")
    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    #println(" range(H.colptr[col], H.colptr[col + 1]-1) ", range(H.colptr[col], H.colptr[col + 1]-1))
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        #iso = (one(UInt8)-H.isotope[i])
        hp = H.nzval[i]
        xp = H.x[i]
        Xsum += xp
        Hsum += hp
        HXentropy += (hp + xp)*log(hp + xp + Float32(1e-10))
        Xentropy += xp*log(xp + Float32(1e-10))
        Hentropy += hp*log(hp+Float32(1e-10))
    end
    HXsum = Hsum + Xsum
    Xentropy = log(Xsum) - Xentropy/Xsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    if (Xentropy < 3) & ((Hentropy >= 3) & (HXentropy >= 3))
        Xw = Xentropy = T(0.25)*(1 + Xentropy)
        Xentropy = zero(T)
        Xsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #iso = (one(UInt8)-H.isotope[i])
            xp = (H.x[i]^Xw)#/Xsum
            Xentropy += xp*log(xp + Float32(1e-10))
            Xsum += xp
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
    elseif (Xentropy < 3) | (Hentropy < 3) | (HXentropy < 3)

        Xw = Xentropy < 3 ? T(0.25)*(1 + Xentropy) : one(T)
        Hw = Hentropy < 3 ? T(0.25)*(1 + Hentropy) :  one(T)
        HXw = HXentropy < 3 ? T(0.25)*(1 + HXentropy) : one(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            #println("a H.nzval[i], ", H.nzval[i], "H.x[i], ",H.x[i])]        
            #iso = (one(UInt8)-H.isotope[i])
            hp = (H.nzval[i]^Hw) #/Hsum
            xp = (H.x[i]^Xw)#/Xsum
            hxp = ((H.nzval[i] + H.x[i])^HXw)
            Hentropy += hp*log(hp + Float32(1e-10))
            Xentropy += xp*log(xp + Float32(1e-10))
            HXentropy += hxp*log(hxp + Float32(1e-10))
            Xsum += xp
            Hsum += hp
            HXsum += hxp
        end
        #println("TRUE Hw $Hw Xw $Xw")
        #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
        Xentropy = log(Xsum) - Xentropy/Xsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXentropy = log(HXsum) - HXentropy/HXsum
    end
    #println("XSum $Xsum Hsum $Hsum HXsum $HXsum")
    #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
function getEntropy(H::SparseArray{Ti, T}, r::Vector{T}, w::Vector{T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}
    #println("col $col")
    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    #println(" range(H.colptr[col], H.colptr[col + 1]-1) ", range(H.colptr[col], H.colptr[col + 1]-1))
    @inbounds @fastmath for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        if H.matched[i]
            hp = H.nzval[i]*w[col]
            xp = H.nzval[i]*w[col] + abs(r[H.rowval[i]])
        else
            hp = zero(T)#H.nzval[i]*w[col]
            xp =  abs(r[H.rowval[i]])
        end
        Xsum += xp
        Hsum += hp
        HXentropy += (hp + xp)*log(hp + xp + Float32(1e-10))
        Xentropy += (xp)*log(xp + Float32(1e-10))
        Hentropy += (hp)*log(hp+Float32(1e-10))
    end
    HXsum = Hsum + Xsum
    Xentropy = log(Xsum) - Xentropy/Xsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    #=
    if (Xentropy < 3) & ((Hentropy >= 3) & (HXentropy >= 3))
        Xw = T(0.25)*(1 + Xentropy)
        Xentropy = zero(T)
        Xsum = zero(T)
        @inbounds @fastmath  for i in range(H.colptr[col], H.colptr[col + 1]-1)
            if H.matched[i]
                xp = ((H.nzval[i]*w[col] + abs(r[H.rowval[i]]))^Xw)#/Xsum
                Xentropy += xp*log(xp + Float32(1e-10))
                Xsum += xp
            end
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
    elseif (Xentropy < 3) | (Hentropy < 3) | (HXentropy < 3)

        Xw = Xentropy < 3 ? T(0.25)*(1 + Xentropy) : one(T)
        Hw = Hentropy < 3 ? T(0.25)*(1 + Hentropy) :  one(T)
        HXw = HXentropy < 3 ? T(0.25)*(1 + HXentropy) : one(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        @inbounds @fastmath  for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise    
            if H.matched[i]
                hp = H.nzval[i]*w[col]
                xp = H.nzval[i]*w[col] + abs(r[H.rowval[i]])
            else
                hp = H.nzval[i]*w[col]
                xp = zero(T)
            end
            hxp = (hp + xp)^HXw
            hp = (hp^Hw) #/Hsum
            xp = (xp^Xw)#/Xsum
            Hentropy += hp*log(hp + Float32(1e-10))
            Xentropy += xp*log(xp + Float32(1e-10))
            HXentropy += hxp*log(hxp + Float32(1e-10))
            Xsum += xp
            Hsum += hp
            HXsum += hxp
        end
        #println("TRUE Hw $Hw Xw $Xw")
        #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
        Xentropy = log(Xsum) - Xentropy/Xsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXentropy = log(HXsum) - HXentropy/HXsum
    end
    =#
    #println("XSum $Xsum Hsum $Hsum HXsum $HXsum")
    #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
function getEntropy(H::SparseArray{Ti, T}, r::Vector{T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}
    #println("col $col")
    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    #println(" range(H.colptr[col], H.colptr[col + 1]-1) ", range(H.colptr[col], H.colptr[col + 1]-1))
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]
        xp = H.x[i]
        Xsum += xp
        Hsum += hp
        HXentropy += (hp + xp)*log(hp + xp + Float32(1e-10))
        Xentropy += (xp)*log(xp + Float32(1e-10))
        Hentropy += (hp)*log(hp+Float32(1e-10))
    end
    HXsum = Hsum + Xsum
    Xentropy = log(Xsum) - Xentropy/Xsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    if (Xentropy < 3) & ((Hentropy >= 3) & (HXentropy >= 3))
        Xw = Xentropy = T(0.25)*(1 + Xentropy)
        Xentropy = zero(T)
        Xsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            xp = (H.x[i]^Xw)#/Xsum
            Xentropy += xp*log(xp + Float32(1e-10))
            Xsum += xp
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
    elseif (Xentropy < 3) | (Hentropy < 3) | (HXentropy < 3)

        Xw = Xentropy < 3 ? T(0.25)*(1 + Xentropy) : one(T)
        Hw = Hentropy < 3 ? T(0.25)*(1 + Hentropy) :  one(T)
        HXw = HXentropy < 3 ? T(0.25)*(1 + HXentropy) : one(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise    
            hp = (H.nzval[i]^Hw) #/Hsum
            xp = (H.x[i]^Xw)#/Xsum
            hxp = ((H.nzval[i] + H.x[i])^HXw)
            Hentropy += hp*log(hp + Float32(1e-10))
            Xentropy += xp*log(xp + Float32(1e-10))
            HXentropy += hxp*log(hxp + Float32(1e-10))
            Xsum += xp
            Hsum += hp
            HXsum += hxp
        end
        #println("TRUE Hw $Hw Xw $Xw")
        #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
        Xentropy = log(Xsum) - Xentropy/Xsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXentropy = log(HXsum) - HXentropy/HXsum
    end
    #println("XSum $Xsum Hsum $Hsum HXsum $HXsum")
    #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
#=

function getEntropy(H::SparseArray{Ti, T}, r::Vector{T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}
    #println("col $col")
    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    #println(" range(H.colptr[col], H.colptr[col + 1]-1) ", range(H.colptr[col], H.colptr[col + 1]-1))
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        iso = one(UInt8)-H.isotope[i]
        hp = H.nzval[i]
        xp = H.x[i]
        Xsum += iso*xp
        Hsum += iso*hp
        HXentropy += iso*(hp + xp)*log(hp + xp + Float32(1e-10))
        Xentropy += iso*(xp)*log(xp + Float32(1e-10))
        Hentropy += iso*(hp)*log(hp+Float32(1e-10))
    end
    HXsum = Hsum + Xsum
    Xentropy = log(Xsum) - Xentropy/Xsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    if (Xentropy < 3) & ((Hentropy >= 3) & (HXentropy >= 3))
        Xw = Xentropy = T(0.25)*(1 + Xentropy)
        Xentropy = zero(T)
        Xsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            iso = (one(UInt8)-H.isotope[i])
            xp = (H.x[i]^Xw)#/Xsum
            Xentropy += iso*xp*log(xp + Float32(1e-10))
            Xsum += iso*xp
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
    elseif (Xentropy < 3) | (Hentropy < 3) | (HXentropy < 3)

        Xw = Xentropy < 3 ? T(0.25)*(1 + Xentropy) : one(T)
        Hw = Hentropy < 3 ? T(0.25)*(1 + Hentropy) :  one(T)
        HXw = HXentropy < 3 ? T(0.25)*(1 + HXentropy) : one(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            #println("a H.nzval[i], ", H.nzval[i], "H.x[i], ",H.x[i])]        
            iso = (one(UInt8)-H.isotope[i])
            hp = (H.nzval[i]^Hw) #/Hsum
            xp = (H.x[i]^Xw)#/Xsum
            hxp = ((H.nzval[i] + H.x[i])^HXw)
            Hentropy += iso*hp*log(hp + Float32(1e-10))
            Xentropy += iso*xp*log(xp + Float32(1e-10))
            HXentropy += iso*hxp*log(hxp + Float32(1e-10))
            Xsum += iso*xp
            Hsum += iso*hp
            HXsum += iso*hxp
        end
        #println("TRUE Hw $Hw Xw $Xw")
        #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
        Xentropy = log(Xsum) - Xentropy/Xsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXentropy = log(HXsum) - HXentropy/HXsum
    end
    #println("XSum $Xsum Hsum $Hsum HXsum $HXsum")
    #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end

function getEntropy2(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}
    #println("col $col")
    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    #println(" range(H.colptr[col], H.colptr[col + 1]-1) ", range(H.colptr[col], H.colptr[col + 1]-1))
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        iso = one(UInt8)#H.isotope[i]
        hp = H.nzval[i]
        xp = H.x[i]
        Xsum += iso*xp
        Hsum += iso*hp
        HXentropy += iso*(hp + xp)*log(hp + xp + Float32(1e-10))
        Xentropy += iso*(xp)*log(xp + Float32(1e-10))
        Hentropy += iso*(hp)*log(hp+Float32(1e-10))
    end
    HXsum = Hsum + Xsum
    Xentropy = log(Xsum) - Xentropy/Xsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    if (Xentropy < 3) & ((Hentropy >= 3) & (HXentropy >= 3))
        Xw = Xentropy = T(0.25)*(1 + Xentropy)
        Xentropy = zero(T)
        Xsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            iso = one(UInt8)#H.isotope[i]
            xp = (H.x[i]^Xw)#/Xsum
            Xentropy += iso*xp*log(xp + Float32(1e-10))
            Xsum += iso*xp
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
    elseif (Xentropy < 3) | (Hentropy < 3) | (HXentropy < 3)

        Xw = Xentropy < 3 ? T(0.25)*(1 + Xentropy) : one(T)
        Hw = Hentropy < 3 ? T(0.25)*(1 + Hentropy) :  one(T)
        HXw = HXentropy < 3 ? T(0.25)*(1 + HXentropy) : one(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            #println("a H.nzval[i], ", H.nzval[i], "H.x[i], ",H.x[i])]        
            iso = one(UInt8)#H.isotope[i]
            hp = (H.nzval[i]^Hw) #/Hsum
            xp = (H.x[i]^Xw)#/Xsum
            hxp = ((H.nzval[i] + H.x[i])^HXw)
            Hentropy += iso*hp*log(hp + Float32(1e-10))
            Xentropy += iso*xp*log(xp + Float32(1e-10))
            HXentropy += iso*hxp*log(hxp + Float32(1e-10))
            Xsum += iso*xp
            Hsum += iso*hp
            HXsum += iso*hxp
        end
        #println("TRUE Hw $Hw Xw $Xw")
        #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
        Xentropy = log(Xsum) - Xentropy/Xsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXentropy = log(HXsum) - HXentropy/HXsum
    end
    #println("XSum $Xsum Hsum $Hsum HXsum $HXsum")
    #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
=#
#=
function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}

    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]#/Hsum
        xp = H.x[i]#/Xsum
        Hentropy += hp*log(hp + 1e-10)
        #HXentropy += (hp + xp)*log(hp + xp)
        Xentropy += xp*log(xp + 1e-10)
        Xsum += xp
        Hsum += hp
        #HXsum += xp + hp
    end

    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]/Hsum
        xp = H.x[i]/Xsum
        HXentropy += (hp + xp)*log(hp + xp + 1e-10)
        HXsum += xp + hp
    end
    Xentropy = log(Xsum) - Xentropy/Xsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    
    if (Xentropy < 3) | (Hentropy < 3)

        Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
        Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]^Hw#/Hsum
            xp = H.x[i]^Xw#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum 
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
        HXentropy = zero(Float32)
        HXsum = zero(Float32)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        HXentropy = log(HXsum) - HXentropy/HXsum
    end

    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}

    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]#/Hsum
        xp = H.x[i]#/Xsum
        Hentropy += hp*log(hp + 1e-10)
        #HXentropy += (hp + xp)*log(hp + xp)
        Xentropy += xp*log(xp + 1e-10)
        Xsum += xp
        Hsum += hp
        #HXsum += xp + hp
    end

    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]/Hsum
        xp = H.x[i]/Xsum
        HXentropy += (hp + xp)*log(hp + xp + 1e-10)
        HXsum += xp + hp
    end
    Xentropy = log(Xsum) - Xentropy/Xsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    
    if (Xentropy < 3) | (Hentropy < 3)

        Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
        Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]^Hw#/Hsum
            xp = H.x[i]^Xw#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum 
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
        HXentropy = zero(Float32)
        HXsum = zero(Float32)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        HXentropy = log(HXsum) - HXentropy/HXsum
    end

    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}

    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]#/Hsum
        xp = H.x[i]#/Xsum
        Hentropy += hp*log(hp + 1e-10)
        #HXentropy += (hp + xp)*log(hp + xp)
        Xentropy += xp*log(xp + 1e-10)
        Xsum += xp
        Hsum += hp
        #HXsum += xp + hp
    end

    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]/Hsum
        xp = H.x[i]/Xsum
        HXentropy += (hp + xp)*log(hp + xp + 1e-10)
        HXsum += xp + hp
    end
    Xentropy = log(Xsum) - Xentropy/Xsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    
    if (Xentropy < 3) | (Hentropy < 3)

        Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
        Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]^Hw#/Hsum
            xp = H.x[i]^Xw#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum 
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
        HXentropy = zero(Float32)
        HXsum = zero(Float32)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        HXentropy = log(HXsum) - HXentropy/HXsum
    end

    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}

    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]#/Hsum
        xp = H.x[i]#/Xsum
        Hentropy += hp*log(hp + 1e-10)
        #HXentropy += (hp + xp)*log(hp + xp)
        Xentropy += xp*log(xp + 1e-10)
        Xsum += xp
        Hsum += hp
        #HXsum += xp + hp
    end

    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]/Hsum
        xp = H.x[i]/Xsum
        HXentropy += (hp + xp)*log(hp + xp + 1e-10)
        HXsum += xp + hp
    end
    Xentropy = log(Xsum) - Xentropy/Xsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    
    if (Xentropy < 3) | (Hentropy < 3)

        Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
        Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]^Hw#/Hsum
            xp = H.x[i]^Xw#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum 
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
        HXentropy = zero(Float32)
        HXsum = zero(Float32)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        HXentropy = log(HXsum) - HXentropy/HXsum
    end

    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}

    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]#/Hsum
        xp = H.x[i]#/Xsum
        Hentropy += hp*log(hp + 1e-10)
        #HXentropy += (hp + xp)*log(hp + xp)
        Xentropy += xp*log(xp + 1e-10)
        Xsum += xp
        Hsum += hp
        #HXsum += xp + hp
    end

    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]/Hsum
        xp = H.x[i]/Xsum
        HXentropy += (hp + xp)*log(hp + xp + 1e-10)
        HXsum += xp + hp
    end
    Xentropy = log(Xsum) - Xentropy/Xsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    
    if (Xentropy < 3) | (Hentropy < 3)

        Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
        Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]^Hw#/Hsum
            xp = H.x[i]^Xw#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum 
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
        HXentropy = zero(Float32)
        HXsum = zero(Float32)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        HXentropy = log(HXsum) - HXentropy/HXsum
    end

    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end

function getEntropy(H::SparseArray{Ti, T}, col::Int64) where {Ti<:Integer,T<:AbstractFloat}

    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]#/Hsum
        xp = H.x[i]#/Xsum
        Hentropy += hp*log(hp + 1e-10)
        #HXentropy += (hp + xp)*log(hp + xp)
        Xentropy += xp*log(xp + 1e-10)
        Xsum += xp
        Hsum += hp
        #HXsum += xp + hp
    end

    @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        hp = H.nzval[i]/Hsum
        xp = H.x[i]/Xsum
        HXentropy += (hp + xp)*log(hp + xp + 1e-10)
        HXsum += xp + hp
    end
    Xentropy = log(Xsum) - Xentropy/Xsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    
    if (Xentropy < 3) | (Hentropy < 3)

        Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
        Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]^Hw#/Hsum
            xp = H.x[i]^Xw#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum 
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
        HXentropy = zero(Float32)
        HXsum = zero(Float32)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = (H.nzval[i]^Hw)/Hsum
            xp = (H.x[i]^Xw)/Xsum
            hxp = (hp + xp)^HXw
            HXentropy += (hxp)*log(hxp + 1e-10)
            HXsum += hxp
        end

        HXentropy = log(HXsum) - HXentropy/HXsum
    end

    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end
=#
function getEntropy(H::SparseArray{Ti, T}) where {Ti<:Integer,T<:AbstractFloat}
    entropy_sim = zeros(T, H.n)
    for col in range(1, H.n)

        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]#/Hsum
            xp = H.x[i]#/Xsum
            Hentropy += hp*log(hp + 1e-10)
            #HXentropy += (hp + xp)*log(hp + xp)
            Xentropy += xp*log(xp + 1e-10)
            Xsum += xp
            Hsum += hp
            #HXsum += xp + hp
        end

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            hp = H.nzval[i]/Hsum
            xp = H.x[i]/Xsum
            HXentropy += (hp + xp)*log(hp + xp + 1e-10)
            HXsum += xp + hp
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
        HXentropy = log(HXsum) - HXentropy/HXsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        
        if (Xentropy < 3) | (Hentropy < 3)

            Xw = Xentropy < 3 ? 0.25*(1 + Xentropy) : 1.0
            Hw = Hentropy < 3 ? 0.25*(1 + Hentropy) :  1.0
            HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0

            Hentropy = zero(T)
            Xentropy = zero(T)
            HXentropy = zero(T)
            Hsum = zero(T)
            Xsum = zero(T)
            HXsum = zero(T)
            @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
                #MASK is true for selected ions and false otherwise
                hp = H.nzval[i]^Hw#/Hsum
                xp = H.x[i]^Xw#/Xsum
                Hentropy += hp*log(hp + 1e-10)
                Xentropy += xp*log(xp + 1e-10)
                Xsum += xp
                Hsum += hp
            end

            @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
                #MASK is true for selected ions and false otherwise
                hp = (H.nzval[i]^Hw)/Hsum
                xp = (H.x[i]^Xw)/Xsum
                hxp = (hp + xp)^HXw
                HXentropy += (hxp)*log(hxp + 1e-10)
                HXsum += hxp
            end

            Xentropy = log(Xsum) - Xentropy/Xsum
            HXentropy = log(HXsum) - HXentropy/HXsum 
            Hentropy = log(Hsum) - Hentropy/Hsum
            HXw = HXentropy < 3 ? 0.25*(1 + HXentropy) : 1.0
            HXentropy = zero(Float32)
            HXsum = zero(Float32)
            @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
                #MASK is true for selected ions and false otherwise
                hp = (H.nzval[i]^Hw)/Hsum
                xp = (H.x[i]^Xw)/Xsum
                hxp = (hp + xp)^HXw
                HXentropy += (hxp)*log(hxp + 1e-10)
                HXsum += hxp
            end

            HXentropy = log(HXsum) - HXentropy/HXsum
        end

        entropy_sim[col] = (1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
    end
    return entropy_sim
end




function getEntropy(H::SparseArray{Ti, T}, col::Int64, included_indices) where {Ti<:Integer,T<:AbstractFloat}
    #println("col $col")
    Hsum = zero(T)
    Xsum = zero(T)
    HXsum = zero(T)

    Hentropy = zero(T)
    Xentropy = zero(T)
    HXentropy = zero(T)
    #@inbounds @fastmath 
    #println(" range(H.colptr[col], H.colptr[col + 1]-1) ", range(H.colptr[col], H.colptr[col + 1]-1))
    for i in range(H.colptr[col], H.colptr[col + 1]-1)
        if i ∉ included_indices
            continue
        end
    #for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #MASK is true for selected ions and false otherwise
        #iso = (one(UInt8)-H.isotope[i])
        hp = H.nzval[i]
        xp = H.x[i]
        Xsum += xp
        Hsum += hp
        HXentropy += (hp + xp)*log(hp + xp + Float32(1e-10))
        Xentropy += xp*log(xp + Float32(1e-10))
        Hentropy += hp*log(hp+Float32(1e-10))
    end
    HXsum = Hsum + Xsum
    Xentropy = log(Xsum) - Xentropy/Xsum
    Hentropy = log(Hsum) - Hentropy/Hsum
    HXentropy = log(HXsum) - HXentropy/HXsum
    if (Xentropy < 3) & ((Hentropy >= 3) & (HXentropy >= 3))
        Xw = Xentropy = T(0.25)*(1 + Xentropy)
        Xentropy = zero(T)
        Xsum = zero(T)
        for i in range(H.colptr[col], H.colptr[col + 1]-1)
            if i ∉ included_indices
                continue
            end
            #iso = (one(UInt8)-H.isotope[i])
            xp = (H.x[i]^Xw)#/Xsum
            Xentropy += xp*log(xp + Float32(1e-10))
            Xsum += xp
        end
        Xentropy = log(Xsum) - Xentropy/Xsum
    elseif (Xentropy < 3) | (Hentropy < 3) | (HXentropy < 3)

        Xw = Xentropy < 3 ? T(0.25)*(1 + Xentropy) : one(T)
        Hw = Hentropy < 3 ? T(0.25)*(1 + Hentropy) :  one(T)
        HXw = HXentropy < 3 ? T(0.25)*(1 + HXentropy) : one(T)

        Hentropy = zero(T)
        Xentropy = zero(T)
        HXentropy = zero(T)
        Hsum = zero(T)
        Xsum = zero(T)
        HXsum = zero(T)

        for i in range(H.colptr[col], H.colptr[col + 1]-1)
            if i ∉ included_indices
                continue
            end
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            #println("a H.nzval[i], ", H.nzval[i], "H.x[i], ",H.x[i])]        
            #iso = (one(UInt8)-H.isotope[i])
            hp = (H.nzval[i]^Hw) #/Hsum
            xp = (H.x[i]^Xw)#/Xsum
            hxp = ((H.nzval[i] + H.x[i])^HXw)
            Hentropy += hp*log(hp + Float32(1e-10))
            Xentropy += xp*log(xp + Float32(1e-10))
            HXentropy += hxp*log(hxp + Float32(1e-10))
            Xsum += xp
            Hsum += hp
            HXsum += hxp
        end
        #println("TRUE Hw $Hw Xw $Xw")
        #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
        Xentropy = log(Xsum) - Xentropy/Xsum
        Hentropy = log(Hsum) - Hentropy/Hsum
        HXentropy = log(HXsum) - HXentropy/HXsum
    end
    #println("XSum $Xsum Hsum $Hsum HXsum $HXsum")
    #println("Xentropy $Xentropy, Hentropy $Hentropy, HXentropy $HXentropy")
    return Float32(1 - (2*HXentropy - Xentropy - Hentropy)/(log(4)))
end