
function refinePSMs!(PSMs::DataFrame, precursors::Vector{LibraryPrecursor}; loss::AbstractEstimator = TauEstimator{TukeyLoss}(), maxiter = 200, min_spectral_contrast::AbstractFloat = 0.8)
    transform!(PSMs, AsTable(:) => ByRow(psm -> isDecoy(precursors[psm[:precursor_idx]])) => :decoy)
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float64(getIRT(precursors[psm[:precursor_idx]]))) => :iRT)
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float64(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT)
    transform!(PSMs, AsTable(:) => ByRow(psm -> psm[:weight] < 10.0) => :nmf)
    function predictRTs!(PSMs::DataFrame; maxiter = 200, min_spectral_contrast::AbstractFloat = 0.95)
        best_PSMs = combine(sdf -> sdf[argmax(sdf.spectral_contrast_all), :], groupby(PSMs, [:scan_idx]))
        targets = best_PSMs[:,:decoy].==false
        
        spectral_contrast = best_PSMs[:,:spectral_contrast_all].>=min_spectral_contrast
        
        best_matches = targets .& spectral_contrast
    
        #Predicted iRTs
        iRTs = best_PSMs[:,:iRT][best_matches]#,# ones(Float32, sum(best_matches)))
        #Emperical retention times
        RTs = best_PSMs[:,:RT][best_matches]
        ns1 = Splines2.ns_(collect(range(-25.0, length=50, stop=160.0)),df=4,intercept=true);
        X = ns1(iRTs);
        fit1 = lm(X,RTs);
        PSMs[:,:RT_pred] =  GLM.predict(fit1, ns1(PSMs[:,:iRT]))
        PSMs[:,:RT_error] = abs.(PSMs[:,:RT_pred] .- PSMs[:,:RT])#abs.(PSMs[:,:RT] .- GLM.predict(fit1, ns1(PSMs[:,:iRT])))
        return # RTs, iRTs
    end

    predictRTs!(PSMs, loss = loss, maxiter = maxiter, min_spectral_contrast = min_spectral_contrast)
    sort!(PSMs, [:scan_idx, :total_ions]);

    # Group DataFrame by "day" column
    grouped_df = groupby(PSMs, :scan_idx);


    PSMs[:,:next_best] = (combine(grouped_df) do sub_df
        pushfirst!(diff(sub_df.total_ions), zero(UInt32))
        #next_scores = lead(sub_df.total_ions, default = missing)
        #coalesce(next_scores, 0)
    end)[:,:x1]

    PSMs[:,:diff_hyper] = (combine(grouped_df) do sub_df
        sort!(sub_df, :hyperscore)
        pushfirst!(diff(sub_df.hyperscore), zero(Float64))
        #next_scores = lead(sub_df.total_ions, default = missing)
        #coalesce(next_scores, 0)
    end)[:,:x1]

    PSMs[:,:diff_scribe] = (combine(grouped_df) do sub_df
        sort!(sub_df, :scribe_score)
        pushfirst!(diff(sub_df.scribe_score), zero(Float64))
        #next_scores = lead(sub_df.total_ions, default = missing)
        #coalesce(next_scores, 0)
    end)[:,:x1]

    PSMs[:,:median_ions] = (combine(grouped_df) do sub_df
        #sort!(sub_df, :hyperscore)
        #pushfirst!(diff(sub_df.hyperscore), zero(Float64))
        repeat([median(sub_df.total_ions)], size(sub_df)[1])
        #next_scores = lead(sub_df.total_ions, default = missing)
        #coalesce(next_scores, 0)
    end)[:,:x1]

    grouped_df = groupby(PSMs, :precursor_idx);

    PSMs[:,:n_obs] = (combine(grouped_df) do sub_df
        #sort!(sub_df, :hyperscore)
        #pushfirst!(diff(sub_df.hyperscore), zero(Float64))
        repeat([size(sub_df)[1]], size(sub_df)[1])
        #next_scores = lead(sub_df.total_ions, default = missing)
        #coalesce(next_scores, 0)
    end)[:,:x1]
    transform!(PSMs, AsTable(:) => ByRow(psm -> getCharge(prosit_precs[psm[:precursor_idx]])) => :charge)

end