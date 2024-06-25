abstract type SpectralScores{T<:AbstractFloat} end

struct SpectralScoresComplex{T<:AbstractFloat} <: SpectralScores{T}
    scribe::T
    scribe_corrected::T
    scribe_fitted::T
    gof::T
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

function getDistanceMetrics(H::SparseArray{Ti,T}, spectral_scores::Vector{SpectralScoresSimple{U}}) where {Ti<:Integer,T,U<:AbstractFloat}

    for col in range(1, H.n)
        H_sqrt_sum = zero(T)
        X_sqrt_sum = zero(T)

        H2_norm_m0 = zero(T)
        X2_norm_m0 = zero(T)

        dot_product_m0 = zero(T)

        scribe_score = zero(T)
        city_block_dist = zero(T)

        matched_sum = zero(T)
        unmatched_sum = zero(T)

        N_M0 = 0
        #@turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            M0 = one(UInt8) - H.isotope[i]

            H_sqrt_sum += M0*sqrt(H.nzval[i])
            X_sqrt_sum += M0*sqrt(H.x[i])#/Xsum

            H2_norm_m0 += M0*((H.nzval[i])^2)
            X2_norm_m0 += M0*((H.x[i])^2)
            dot_product_m0 += M0*(H.nzval[i]*H.x[i])

            matched_sum += M0*H.nzval[i]*H.matched[i]
            unmatched_sum += M0*H.nzval[i]*(1 - H.matched[i])

            N_M0 += M0
        end
          
        #Sqrt of sum of squares
        H2_norm_m0 = sqrt(H2_norm_m0)
        X2_norm_m0 = sqrt(X2_norm_m0)


        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            M0 = one(UInt8) - H.isotope[i]

            scribe_score +=  M0*((
                                (sqrt(H.nzval[i])/H_sqrt_sum) - 
                                (sqrt(H.x[i])/X_sqrt_sum)
                                )^2)

            city_block_dist += M0*abs(
                (H.nzval[i]/H2_norm_m0) -
                (H.x[i]/X2_norm_m0)
            )   
        end
        spectral_scores[col] = SpectralScoresSimple(
            Float16(-log((scribe_score)/N_M0)), #scribe_score
            Float16(-log((city_block_dist)/N_M0)), #city_block
            Float16(dot_product_m0/(H2_norm_m0*X2_norm_m0)), #dot_p
            Float16(log2(matched_sum/unmatched_sum)), #matched_ratio
            Float16(Float16(-1.0)*getEntropy(H, col)), #entropy
        )
    end
end

function getDistanceMetrics(w::Vector{T}, r::Vector{T}, H::SparseArray{Ti,T}, spectral_scores::Vector{SpectralScoresComplex{U}}) where {Ti<:Integer,T,U<:AbstractFloat}

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
        
        H_sqrt_sum = zero(T)
        X_sqrt_sum = zero(T)
        H2_norm = zero(T)
        X2_norm = zero(T)
        X_sum = zero(T)
        X2_sum = zero(T)
        dot_product = zero(T)
        scribe_score = zero(T)
        scribe_score_fitted = zero(T)
        city_block_dist = zero(T)
        city_block_fitted = zero(T)
        matched_sum = zero(T)
        unmatched_sum = zero(T)    
        sum_of_residuals = zero(T)
        wh_sum = zero(T)
        N = 0
        @turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
            #MASK is true for selected ions and false otherwise
            #=
            X_sum += H.x[i]
            X2_sum += H.x[i]^2
            =#

            X_sum += H.x[i]
            X2_sum += (H.x[i])^2
            H_sqrt_sum += sqrt(H.nzval[i])
            X_sqrt_sum += sqrt(H.x[i])
            scribe_score_fitted += (w[col]*H.nzval[i] - H.x[i])^2
            city_block_fitted += abs(w[col]*H.nzval[i] - H.x[i])
            wh_sum +=w[col]*H.nzval[i]
            sum_of_residuals += abs(r[H.rowval[i]])
            
            H2_norm += (H.nzval[i])^2
            X2_norm += (H.x[i])^2
            dot_product += H.nzval[i]*H.x[i]
            matched_sum += H.nzval[i]*H.matched[i]
            unmatched_sum += H.nzval[i]*(1 - H.matched[i])
            N += 1
        end
          
        #Sqrt of sum of squares
        H2_norm = sqrt(H2_norm)
        X2_norm = sqrt(X2_norm)
        #@turbo for i in range(H.colptr[col], H.colptr[col + 1]-1)
        for i in range(H.colptr[col], H.colptr[col + 1]-1)    
            #MASK is true for selected ions and false otherwise
            scribe_score +=  (
                                (sqrt(H.nzval[i])/H_sqrt_sum) - 
                                (sqrt(H.x[i])/X_sqrt_sum)
                                )^2  

            city_block_dist += abs(
                (H.nzval[i]/H2_norm) -
                (H.x[i]/X2_norm)
            )          
        end

        @turbo for i in range(1, H.m)
            r[i] = zero(T)
        end    

        for col in range(1, H.n)
            start = H.colptr[col]
            stop = H.colptr[col+1] - 1
            for n in start:stop
                r[H.rowval[n]] += w[col]*H.nzval[n]
            end
        end

        spectral_scores[col] = SpectralScoresComplex(
            Float16(-log((scribe_score)/N)), #scribe_score
            zero(Float16),#Float16(-log((scribe_score_corrected)/N_corrected)), #scribe_score_corrected
            Float16(-log((scribe_score_fitted)/(X2_sum))), #scribe_score_corrected
            Float16(-log((sum_of_residuals)/wh_sum)), #city_block
            Float16(-log((city_block_dist)/N)), #city_block
            Float16(-log((city_block_fitted)/X_sum)), #city_block
            Float16(dot_product/(H2_norm*X2_norm)), #dot_p
            zero(Float16),#Float16(dot_product_corrected/(H2_norm_corrected*X2_norm_corrected)), #spectral_contrast_corrected
            Float16(log2(matched_sum/unmatched_sum)), #matched_ratio
            Float16(Float16(-1.0)*getEntropy(H, r, col)) #entropy
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
        iso = (one(UInt8)-H.isotope[i])
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
#=
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



