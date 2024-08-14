function getBestPSMs!(psms::DataFrame,
                        prec_mz::Arrow.Primitive{T, Vector{T}}; 
                        max_q_val::Float32 = 0.01f0,
                        max_psms::Int64 = 250000) where {T<:AbstractFloat}

    #highest scoring psm for a given precursor
    psms[!,:best_psm] = zeros(Bool, size(psms, 1))
    #fwhm estimate of the precursor
    psms[!,:fwhm] = zeros(Union{Missing, Float32}, size(psms, 1))
    #number of scans below the q value threshold for hte precursor
    psms[!,:scan_count] = zeros(UInt16,size(psms, 1))

    #Get best psm for each precursor 
    #ASSUMES psms IS SORTED BY RT IN ASCENDING ORDER
    gpsms = groupby(psms,:precursor_idx)
    for (precursor_idx, prec_psms) in pairs(gpsms)

        #Get the best scoring psm, and the its row index. 
        #Get the maximum intensity psm (under the q_value threshold) and its row index. 
        max_irt, min_irt = missing, missing
        max_log2_intensity, max_idx = missing, one(Int64)
        best_psm_score, best_psm_idx = zero(Float32), one(Int64)
        scan_count = one(UInt8)
        for i in range(1, size(prec_psms, 1))
            if (prec_psms[i,:q_value].<=max_q_val) 
                if coalesce(max_log2_intensity, zero(Float32)) < prec_psms[i,:log2_summed_intensity]
                    max_log2_intensity = prec_psms[i,:log2_summed_intensity]
                    max_idx = i
                end
            end
            if prec_psms[i,:score]>best_psm_score
                best_psm_idx = i
                best_psm_score = prec_psms[i,:score]
            end
        end
        #Mark the best psm 
        prec_psms[best_psm_idx,:best_psm] = true

        #Try to estimate the fwhm. 
        i = max_idx - 1
        while i > 0
            #Is the i'th psm above half the maximum 
            if (prec_psms[i,:q_value].<=max_q_val) & (prec_psms[i,:log2_summed_intensity] > (max_log2_intensity - 1.0))
                scan_count += 1
                min_irt = prec_psms[i,:iRT]
            else
                break
            end
            i -= 1
        end
        
        i = max_idx + 1
        while i <= size(prec_psms, 1)
            #Is the i'th psm above half the maximum 
            if (prec_psms[i,:q_value].<=max_q_val) & (prec_psms[i,:log2_summed_intensity] > (max_log2_intensity - 1.0))
                scan_count += 1
                max_irt = prec_psms[i,:iRT]
            else
                break
            end
            i += 1
        end

        prec_psms[best_psm_idx,:fwhm] = max_irt - min_irt
        prec_psms[best_psm_idx,:scan_count] = scan_count
    end

    filter!(x->x.best_psm, psms);
    sort!(psms,:score, rev = true)
    n = size(psms, 1)
    select!(psms, [:precursor_idx,:RT,:iRT_predicted,:q_value,:score,:prob,:fwhm,:scan_count])
    delete!(psms, min(n, max_psms + 1):n)

    mz = zeros(T, size(psms, 1));
    precursor_idx = psms[!,:precursor_idx]::Vector{UInt32}
    Threads.@threads for i in range(1, size(psms, 1))
        mz[i] = prec_mz[precursor_idx[i]];
    end
    psms[!,:prec_mz] = mz

    return
end