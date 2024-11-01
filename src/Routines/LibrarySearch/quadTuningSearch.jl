function getNearestAdjacentScans(scan_idx::UInt32,
                            centerMz::AbstractArray{Union{Missing, T}},
                            isolationWidthMz::AbstractArray{Union{Missing, T}};
                            scans_to_search::Int64 = 500
        ) where {T<:AbstractFloat}
    upperBoundMz = centerMz[scan_idx] + isolationWidthMz[scan_idx]/T(2.0)
    min_diff, min_diff_idx = typemax(Float32), -1
    for near_scan_idx in range(scan_idx, min(scan_idx + scans_to_search, length(centerMz)))
        if ismissing(centerMz[near_scan_idx])
            continue
        end
        lowerBoundMz = centerMz[near_scan_idx] - isolationWidthMz[near_scan_idx]/T(2.0)
        if abs(upperBoundMz - lowerBoundMz) < min_diff
            min_diff_idx = near_scan_idx
            min_diff = abs(upperBoundMz - lowerBoundMz) 
        end
    end
    next_scan_idx = sign(min_diff_idx)==-1 ? scan_idx : min_diff_idx

    min_diff, min_diff_idx = typemax(Float32), -1
    lowerBoundMz = centerMz[scan_idx] - isolationWidthMz[scan_idx]/T(2.0)
    for near_scan_idx in range(scan_idx, max(scan_idx - scans_to_search, 1))
        if ismissing(centerMz[near_scan_idx])
            continue
        end
        upperBoundMz = centerMz[near_scan_idx] + isolationWidthMz[near_scan_idx]/T(2.0)
        if abs(upperBoundMz - lowerBoundMz) < min_diff
            min_diff_idx = near_scan_idx
            min_diff = abs(upperBoundMz - lowerBoundMz) 
        end
    end
    prev_scan_idx = sign(min_diff_idx)==-1 ? scan_idx : min_diff_idx



    return prev_scan_idx, next_scan_idx
end

function getScanToPrecIdx(
    scan_idxs::AbstractVector{UInt32},
    prec_idxs::AbstractVector{UInt32},
    centerMz::AbstractVector{Union{Missing, Float32}},
    isolationWidthMz::AbstractVector{Union{Missing, Float32}}
    )
    N = length(scan_idxs)
    scan_idx_to_prec_idx = Dictionary{UInt32, Vector{UInt32}}()
    for i in range(1, N)
        scan_idx = scan_idxs[i]
        prec_idx = prec_idxs[i]
        #Have encountered scan 
        if haskey(scan_idx_to_prec_idx, scan_idx)
            push!(scan_idx_to_prec_idx[scan_idx], prec_idx)#(zip(psms[!,:scan_idx], psms[!,:precursor_idx]))
            prev_scan_idx, next_scan_idx = getNearestAdjacentScans(
                scan_idx, centerMz, isolationWidthMz
            )
            #Have encountered nearest scan 
            if haskey(scan_idx_to_prec_idx, next_scan_idx)
                push!(scan_idx_to_prec_idx[next_scan_idx], prec_idx)#(zip(psms[!,:scan_idx], psms[!,:precursor_idx]))
            else
                insert!(scan_idx_to_prec_idx, next_scan_idx,[prec_idx])
            end
            if haskey(scan_idx_to_prec_idx, prev_scan_idx)
                push!(scan_idx_to_prec_idx[prev_scan_idx], prec_idx)#(zip(psms[!,:scan_idx], psms[!,:precursor_idx]))
            else
                insert!(scan_idx_to_prec_idx, prev_scan_idx,[prec_idx])
            end

        else
            insert!(scan_idx_to_prec_idx, scan_idx,[prec_idx])
            prev_scan_idx, next_scan_idx = getNearestAdjacentScans(
                scan_idx, centerMz, isolationWidthMz
            )
            if haskey(scan_idx_to_prec_idx, next_scan_idx)
                push!(scan_idx_to_prec_idx[next_scan_idx], prec_idx)#(zip(psms[!,:scan_idx], psms[!,:precursor_idx]))
            else
                insert!(scan_idx_to_prec_idx, next_scan_idx,[prec_idx])
            end
            if haskey(scan_idx_to_prec_idx, prev_scan_idx)
                push!(scan_idx_to_prec_idx[prev_scan_idx], prec_idx)#(zip(psms[!,:scan_idx], psms[!,:precursor_idx]))
            else
                insert!(scan_idx_to_prec_idx, prev_scan_idx,[prec_idx])
            end
        end
    end
    return scan_idx_to_prec_idx
end

function filterPSMs(
    iso_idx::AbstractVector{UInt8},
    n_matches::AbstractVector{UInt8},
    weight::AbstractVector{Float32})

    n = length(iso_idx)
    mask = Vector{Bool}(undef, n)
    @inbounds for i in 1:n
        mask[i] = (iso_idx[i] < 3) && 
                (n_matches[i] > 3) && 
                (weight[i]>0)
    end
    return mask
end

function addColumns(
    precursor_idx::AbstractVector{UInt32},
    lib_precursor_mz::AbstractVector{Float32},
    lib_prec_charge::AbstractVector{UInt8},
    lib_sulfur_count::AbstractVector{UInt8},
    iso_idx::AbstractVector{UInt8},
    center_mz::AbstractVector{Float32},
    iso_splines::IsotopeSplineModel)

    mono_mz = [lib_precursor_mz[pid] for pid in precursor_idx]
    prec_charge = [lib_prec_charge[pid] for pid in precursor_idx]
    sulfur_count = [lib_sulfur_count[pid] for pid in precursor_idx]
    iso_mz = Float32.(mono_mz .+ NEUTRON.*iso_idx./prec_charge)
    mz_offset = iso_mz .- center_mz
    δ = zeros(Float32, length(precursor_idx))
    for i in range(1, length(precursor_idx))
        s_count = min(Int64(sulfur_count[i]), 5)
        mono_mass = mono_mz[i]*prec_charge[i]
        δ[i] = iso_splines(s_count, 0, mono_mass)/iso_splines(s_count, 1, mono_mass)
    end
    return DataFrame((mono_mz = mono_mz, prec_charge = prec_charge, sulfur_count = sulfur_count, iso_mz = iso_mz, mz_offset = mz_offset, δ = δ))
end

function summarizePrecursor(
    iso_idx::AbstractVector{UInt8},
    center_mz::AbstractVector{Float32},
    iso_mz::AbstractVector{Float32},
    prec_charge::AbstractVector{UInt8},
    weight::AbstractVector{Float32},
    δ::AbstractVector{Float32})
    if length(iso_idx) == 2
        m0_idx, m1_idx = 0, 0
        if iso_idx[1] == 1
            m0_idx, m1_idx = 1, 2
        else
            m0_idx, m1_idx = 2, 1
        end
        return (center_mz = center_mz[m0_idx],
                δ = δ[m0_idx],
                yt = log(weight[m0_idx]/(weight[m1_idx]*δ[m0_idx])), 
                x0 = iso_mz[m0_idx]-center_mz[m0_idx], 
                x1 = iso_mz[m1_idx]-center_mz[m1_idx], 
                prec_charge = prec_charge[m0_idx])
    else
        return (center_mz = missing, 
                δ = missing, 
                yt = missing, 
                x0 = missing, 
                x1 = missing, 
                prec_charge = missing)
    end
end

function quadTuningSearch(    rt_to_irt_map_dict,
                                frag_err_dist_dict,
                                irt_errs,
                                MS_TABLE_PATHS,
                                params_,
                                spec_lib,
                                ionMatches,
                                ionMisses,
                                all_fmatches,
                                IDtoCOL,
                                ionTemplates,
                                iso_splines,
                                scored_psms,
                                unscored_psms,
                                spectral_scores,
                                precursor_weights,
                                precs)

    precursors = spec_lib["precursors"]
    model_fits = []
    quad_model_dict = Dict{Int64, }
    for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
        MS_TABLE = Arrow.Table(MS_TABLE_PATH)
        psms =  vcat(LibrarySearch(
                                MS_TABLE,
                                params_;
                                frag_index = spec_lib["presearch_f_index"],
                                precursors = spec_lib["precursors"],
                                fragment_lookup_table = spec_lib["f_det"],
                                rt_to_irt_spline =  rt_to_irt_map_dict[ms_file_idx],
                                ms_file_idx = UInt32(ms_file_idx),
                                irt_tol = irt_errs[ms_file_idx],
                                ion_matches = ionMatches,
                                ion_misses = ionMisses,
                                fmatches = all_fmatches,
                                id_to_col = IDtoCOL,
                                ion_templates = ionTemplates,
                                iso_splines = iso_splines,
                                scored_psms = scored_psms,
                                unscored_psms = unscored_psms,
                                spectral_scores = spectral_scores,
                                prec_to_score = precs,
                                mass_err_model = frag_err_dist_dict[ms_file_idx],
                                sample_rate = params_[:presearch_params]["sample_rate"],
                                params = params_[:presearch_params],
                                isotope_err_bounds = (0, 0),
                                mz_overhang = 1.0f0,
                                quad_transmission_model = SquareQuadModel(1.0f0)
                                                )...);

        psms[!,:best_psms ] .= false
        if iszero(size(psms, 1))
            continue
        end
        addPreSearchColumns!(psms, 
                                    MS_TABLE, 
                                    spec_lib["precursors"][:is_decoy],
                                    spec_lib["precursors"][:irt],
                                    spec_lib["precursors"][:prec_charge],
                                    MS_TABLE[:retentionTime],
                                    MS_TABLE[:TIC]
                                )
        #Score Precursors
        scorePresearch!(psms)
        getQvalues!(psms[!,:prob], psms[!,:target], psms[!,:q_value])
        filter!(:q_value => x -> x<=params_[:presearch_params]["max_qval"], psms)
        #Get best scan per passing precursor 
        psms[!,:best_psms] .= false
        grouped_psms = groupby(psms,:precursor_idx)
        for psms in grouped_psms
            best_idx = argmax(psms.prob)
            psms[best_idx,:best_psms] = true
        end
        filter!(x->x.best_psms, psms)
        filter!(x->x.target, psms)

        scan_idx_to_prec_idx = Dictionary{UInt32, Vector{UInt32}}()

        scan_idx_to_prec_idx = getScanToPrecIdx(
            psms[!,:scan_idx],
            psms[!,:precursor_idx],
            MS_TABLE[:centerMz], MS_TABLE[:isolationWidthMz]
        )

        scan_idxs = Set(keys(scan_idx_to_prec_idx))

        thread_tasks, total_peaks = partitionScansToThreads(MS_TABLE[:mz_array],
                                                            MS_TABLE[:retentionTime],
                                                            MS_TABLE[:centerMz],
                                                            MS_TABLE[:msOrder],
                                                            Threads.nthreads(),
                                                            1)

        tasks = map(thread_tasks) do thread_task
            Threads.@spawn begin 
                thread_id = first(thread_task)
                return QuadTransmissionSearch(
                                    MS_TABLE,
                                    last(thread_task), #getRange(thread_task),
                                    scan_idx_to_prec_idx,
                                    scan_idxs,
                                    precursors,
                                    spec_lib["f_det"], 
                                    UInt32(ms_file_idx),
                                    frag_err_dist_dict[ms_file_idx],
                                    Float32(100000.0f0), #Effectively use squared error 
                                    Int64(params_[:deconvolution_params]["max_iter_newton"]),
                                    Int64(params_[:deconvolution_params]["max_iter_bisection"]),
                                    Int64(params_[:deconvolution_params]["max_iter_outer"]),
                                    Float32(params_[:deconvolution_params]["accuracy_newton"]),
                                    Float32(params_[:deconvolution_params]["accuracy_bisection"]),
                                    Float32(params_[:deconvolution_params]["max_diff"]),
                                    ionMatches[thread_id],
                                    ionMisses[thread_id],
                                    IDtoCOL[thread_id],
                                    ionTemplates[thread_id],
                                    iso_splines,
                                    unscored_psms[thread_id],
                                    spectral_scores[thread_id],
                                    precursor_weights[thread_id],
                                    Set(2),
                                )
            end
        end
        psms = vcat(fetch.(tasks)...)


        psms = psms[filterPSMs(psms[!,:iso_idx], psms[!,:n_matches], psms[!,:weight]),:]
        sort!(psms, [:scan_idx,:precursor_idx,:iso_idx])
        psms = hcat(psms, addColumns(
            psms[!,:precursor_idx],
            precursors[:mz],
            precursors[:prec_charge],
            precursors[:sulfur_count],
            psms[!,:iso_idx],
            psms[!,:center_mz],
            iso_splines
        ))
        combined_psms = combine(groupby(psms, [:scan_idx,:precursor_idx])) do precursor_psms
            return summarizePrecursor(
                precursor_psms[!,:iso_idx],
                precursor_psms[!,:center_mz],
                precursor_psms[!,:iso_mz],
                precursor_psms[!,:prec_charge],
                precursor_psms[!,:weight],
                precursor_psms[!,:δ]
            )
        end
        filter!(x->!ismissing(x.yt), combined_psms)

        combined_psms[!,:prec_charge] = UInt8.(combined_psms[!,:prec_charge])
        combined_psms[!,:x0] = Float32.(combined_psms[!,:x0])
        combined_psms[!,:yt] = Float32.(combined_psms[!,:yt])
        combined_psms[!,:x1] = Float32.(combined_psms[!,:x1])

        combined_psms[!,:keep_charge] .= true
        for (prec_charge, sub_psms) in pairs(groupby(combined_psms, :prec_charge))
            median_yt = median(sub_psms[!,:yt]::AbstractVector{Float32})
            for i in range(1, size(sub_psms, 1))
                sub_psms[i,:yt] = sub_psms[i,:yt] .- median_yt
            end
            if size(sub_psms, 1) < 1000
                sub_psms[!,:keep_charge] .= false
            end
        end
        filter!(x->x.keep_charge::Bool, combined_psms)

        binned_psms = MergeBins(combined_psms,(-2.0, 2.0), min_bin_size = 10, min_bin_width = 0.1)
        fitted_rqm = fitRazoQuadModel(
            2.0,
            binned_psms[!,:median_x0],  binned_psms[!,:median_x1], binned_psms[!,:median_yt],
            λ0 = 1e-1, #Initial L-M damping coeficient 
            ϵ1 = 1e-5, #Convergence in the gradient
            ϵ2 = 1e-4, #Convergence in the coeficients
            ϵ3 = 1e-5 #Conergence in squared error
            )
        push!(model_fits, fitted_rqm)
    end
    return model_fits
end


