abstract type SpectralScores{T<:AbstractFloat} end

struct SpectralScoresComplex{T<:AbstractFloat} <: SpectralScores{T}
    scribe::T
    scribe_corrected::T
    scribe_fitted::T
    city_block::T
    city_block_fitted::T
    spectral_contrast::T
    spectral_contrast_corrected::T
    matched_ratio::T
    entropy_score::T
end

struct SpectralScoresSimple{T<:AbstractFloat} <: SpectralScores{T} 
    scribe::T
    city_block::T
    spectral_contrast::T
    matched_ratio::T
    entropy_score::T
end

function getDistanceMetrics(w::Vector{T}, H::SparseArray{Ti,T}, spectral_scores::Vector{SpectralScoresSimple{U}}) where {Ti<:Integer,T,U<:AbstractFloat}

    for col in range(1, H.n)
        H_sqrt_sum = zero(T)
        X_sqrt_sum = zero(T)

        X_sum = zero(T)
        H2_norm = zero(T)
        X2_norm = zero(T)

        dot_product = zero(T)

        scribe_score = zero(T)

        city_block_dist = zero(T)

        matched_sum = zero(T)
        unmatched_sum = zero(T)

        N = 0
        N_corrected = 0
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        #for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            X_sum += H.x[i]
            H_sqrt_sum += sqrt(H.nzval[i])
            #H_sqrt_sum_fitted += sqrt(w[col]*H.nzval[i])
            X_sqrt_sum += sqrt(H.x[i] + 1e-10)#/Xsum

            H2_norm += (H.nzval[i])^2 + 1e-10
            #H2_norm_fitted += (w[col]*H.nzval[i])^2 + 1e-10
            X2_norm += (H.x[i])^2 + 1e-10
            dot_product += H.nzval[i]*H.x[i]

            matched_sum += H.nzval[i]*H.matched[i]
            unmatched_sum += H.nzval[i]*(1 - H.matched[i])

            N += 1
            N_corrected += H.mask[i]
        end
          
        #Sqrt of sum of squares
        H2_norm = sqrt(H2_norm)
        #H2_norm_fitted = sqrt(H2_norm_fitted)
        X2_norm = sqrt(X2_norm)

        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            scribe_score +=  (
                                (sqrt(H.nzval[i])/H_sqrt_sum) - 
                                (sqrt(H.x[i])/X_sqrt_sum)
                                )^2  

            city_block_dist += abs(
                (H.nzval[i]/H2_norm) -
                (H.x[i]/X2_norm)
            )        

        end

        spectral_scores[col] = SpectralScoresSimple(
            Float16(-log((scribe_score)/N)), #scribe_score
            Float16(-log((city_block_dist)/N)), #city_block
            Float16(dot_product/(H2_norm*X2_norm)), #dot_p
            Float16(log2(matched_sum/unmatched_sum)), #matched_ratio
            Float16(getEntropy(H, col)) #entropy
        )
    end
end

function getDistanceMetrics(w::Vector{T}, H::SparseArray{Ti,T}, spectral_scores::Vector{SpectralScoresComplex{U}}) where {Ti<:Integer,T,U<:AbstractFloat}

    for col in range(1, H.n)
        
        H_sqrt_sum = zero(T)
        X_sqrt_sum = zero(T)

        H_sqrt_sum_corrected = zero(T)
        X_sqrt_sum_corrected = zero(T)

        H2_norm = zero(T)
        X2_norm = zero(T)
        X_sum = zero(T)
        #H2_norm_corrected = zero(T)
        #X2_norm_corrected = zero(T)

        dot_product = zero(T)
        #dot_product_corrected = zero(T)

        scribe_score = zero(T)
        scribe_score_corrected = zero(T)
        scribe_score_fitted = zero(T)

        city_block_dist = zero(T)
        city_block_fitted = zero(T)

        matched_sum = zero(T)
        unmatched_sum = zero(T)
    
        
        max_residual = zero(T)
        max_residual_int = 0

        for i in range(H.colptr[col], H.colptr[col + 1]-1)
            if abs(w[col]*H.nzval[i] - H.x[i]) > max_residual
                max_residual = abs(w[col]*H.nzval[i] - H.x[i])
                max_residual_int = i
            end
        end

        if !iszero(max_residual_int)
            setMask!(H, max_residual_int, false)
        end
    
        N = 0
        N_corrected = 0
        #@turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            X_sum += H.x[i]
            H_sqrt_sum += sqrt(H.nzval[i])
            #X_sqrt_sum += sqrt(H.x[i] + 1e-10)#/Xsum
            X_sqrt_sum += sqrt(H.x[i])#/Xsum
            scribe_score_fitted += (w[col]*H.nzval[i] - H.x[i])^2
            city_block_fitted += abs(w[col]*H.nzval[i] - H.x[i])
            H_sqrt_sum_corrected += sqrt(H.nzval[i]*H.mask[i])
            X_sqrt_sum_corrected += sqrt(H.x[i]*H.mask[i] + 1e-10)#/Xsum

            H2_norm += (H.nzval[i])^2# + 1e-10
            #H2_norm_fitted += (w[col]*H.nzval[i])^2 + 1e-10
            X2_norm += (H.x[i])^2# + 1e-10
            dot_product += H.nzval[i]*H.x[i]

            #H2_norm_corrected += (H.nzval[i]*H.mask[i])^2# + 1e-10
            #X2_norm_corrected += (H.x[i]*H.mask[i])^2# + 1e-10
            #dot_product_corrected += H.nzval[i]*H.x[i]*H.mask[i]

            matched_sum += H.nzval[i]*H.matched[i]
            unmatched_sum += H.nzval[i]*(1 - H.matched[i])

            N += 1
            N_corrected += H.mask[i]
        end
          
        #Sqrt of sum of squares
        H2_norm = sqrt(H2_norm)
        #H2_norm_fitted = sqrt(H2_norm_fitted)
        X2_norm = sqrt(X2_norm)
        #H2_norm_corrected = sqrt(H2_norm_corrected)
        #X2_norm_corrected = sqrt(X2_norm_corrected)

        #@turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            scribe_score +=  (
                                (sqrt(H.nzval[i])/H_sqrt_sum) - 
                                (sqrt(H.x[i])/X_sqrt_sum)
                                )^2  

            scribe_score_corrected +=  (
                (sqrt(H.nzval[i]*H.mask[i])/H_sqrt_sum_corrected) - 
                (sqrt(H.x[i]*H.mask[i])/X_sqrt_sum_corrected)
                )^2  

            city_block_dist += abs(
                (H.nzval[i]/H2_norm) -
                (H.x[i]/X2_norm)
            )          
        end

        spectral_scores[col] = SpectralScoresComplex(
            Float16(-log((scribe_score)/N)), #scribe_score
            Float16(-log((scribe_score_corrected)/N_corrected)), #scribe_score_corrected
            Float16(-log((scribe_score_fitted)/(X2_norm^2))), #scribe_score_corrected
            Float16(-log((city_block_dist)/N)), #city_block
            Float16(-log((city_block_fitted)/X_sum)), #city_block
            Float16(dot_product/(H2_norm*X2_norm)), #dot_p
            zero(Float16),#Float16(dot_product_corrected/(H2_norm_corrected*X2_norm_corrected)), #spectral_contrast_corrected
            Float16(log2(matched_sum/unmatched_sum)), #matched_ratio
            Float16(getEntropy(H, col)) #entropy
        )
    end
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



