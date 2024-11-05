function addSecondSearchColumns!(psms::DataFrame, 
                        scan_retention_time::AbstractVector{Float32},
                        prec_charge::AbstractVector{UInt8},
                        prec_is_decoy::AbstractVector{Bool},
                        prec_id_to_cv_fold::Dictionary{UInt32, UInt8})
    
    ###########################
    #Correct Weights by base-peak intensity
    filter!(x->x.weight>0.0, psms);
    ###########################
    #Allocate new columns
   
    #Threads.@threads for i in ProgressBar(range(1, size(psms)[1]))
    N = size(psms, 1);
    decoys = zeros(Bool, N);
    rt = zeros(Float32, N);
    #TIC = zeros(Float16, N);
    total_ions = zeros(UInt16, N);
    err_norm = zeros(Float16, N);
    targets = zeros(Bool, N);
    prec_charges = zeros(UInt8, N);
    #prec_mzs = zeros(Float32, N);
    cv_fold = zeros(UInt8, N);
    scan_idxs::Vector{UInt32} = psms[!,:scan_idx]
    prec_idxs::Vector{UInt32} = psms[!,:precursor_idx]
    y_count::Vector{UInt8} = psms[!,:y_count]
    b_count::Vector{UInt8} = psms[!,:b_count]
    isotope_count::Vector{UInt8} = psms[!,:isotope_count]
    error::Vector{Float32} = psms[!,:error]
    #psms[!,:total_ions]
    #tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    #scan_retention_time = MS_TABLE[:retentionTime]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    matched_ratio::Vector{Float16} = psms[!,:matched_ratio]

    tasks_per_thread = 5
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                prec_idx = prec_idxs[i]
                scan_idx = scan_idxs[i]

                decoys[i] = prec_is_decoy[prec_idx];
                targets[i] = decoys[i] == false
                rt[i] = Float32(scan_retention_time[scan_idx]);
                prec_charges[i] = prec_charge[prec_idx];
                total_ions[i] = UInt16(y_count[i] + b_count[i] + isotope_count[i]);
                err_norm[i] = min(Float16((2^error[i])/(total_ions[i])), 6e4)
                if isinf(matched_ratio[i])
                    matched_ratio[i] = Float16(60000)*sign(matched_ratio[i])
                end
                cv_fold[i] = prec_id_to_cv_fold[prec_idx]
            end
        end
    end
    fetch.(tasks)
    psms[!,:matched_ratio] = matched_ratio
    psms[!,:decoy] = decoys
    psms[!,:rt] = rt
    #psms[!,:TIC] = TIC
    psms[!,:total_ions] = total_ions
    psms[!,:err_norm] = err_norm
    psms[!,:target] = targets
    psms[!,:charge] = prec_charges
    #psms[!,:prec_mz] = prec_mzs
    psms[!,:cv_fold] = cv_fold
    #######
    sort!(psms,:rt); #Sorting before grouping is critical. 
    return nothing
end

function addPostIntegrationFeatures!(psms::DataFrame, 
                                    MS_TABLE::Arrow.Table, 
                                    precursor_sequence::AbstractArray{String},
                                    structural_mods::AbstractArray{Union{String, Missing}},
                                    prec_mz::AbstractArray{T},
                                    prec_irt::AbstractArray{T},
                                    prec_charge::AbstractArray{UInt8},
                                    precursor_missed_cleavage::AbstractArray{UInt8},
                                    tic::AbstractVector{Float32},
                                    masses::AbstractArray,
                                    ms_file_idx::Integer,
                                    rt_to_irt_interp::UniformSpline,
                                    prec_id_to_irt::Dictionary{UInt32, @NamedTuple{irt::Float32, mz::Float32}}) where {T<:AbstractFloat}

    #filter!(x -> x.best_scan, psms);
    filter!(x->x.weight>0, psms);
    #filter!(x->x.data_points>0, psms)
    ###########################
    #Allocate new columns
    #println("TEST")
    N = size(psms, 1)
    irt_diff = zeros(Float32, N)
    irt_obs = zeros(Float32, N)
    irt_pred = zeros(Float32, N)
    irt_error = zeros(Float32, N)

    #psms[!,:missed_cleavage] .= zero(UInt8);
    #psms[!,:sequence] .= "";
    missed_cleavage = zeros(UInt8, N);
    #sequence = Vector{String}(undef, N);
    #stripped_sequence = Vector{String}(undef, N);
    adjusted_intensity_explained = zeros(Float16, N);
    #log2_base_peak_intensity = zeros(Float16, N);
    prec_charges = zeros(UInt8, N)
    sequence_length = zeros(UInt8, N);
    #b_y_overlap = zeros(Bool, N);
    spectrum_peak_count = zeros(Float16, N);
    prec_mzs = zeros(Float32, size(psms, 1));
    Mox = zeros(UInt8, N);
    TIC = zeros(Float16, N);

    #tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    precursor_idx::Vector{UInt32} = psms[!,:precursor_idx] 
    scan_idx::Vector{UInt32} = psms[!,:scan_idx]
    #masses = MS_TABLE[:mz_array]::Arrow.List{Union{Missing, SubArray{Union{Missing, Float32}, 1, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}, Tuple{UnitRange{Int64}}, true}}, Int64, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}}
    #longest_y::Vector{UInt8} = psms[!,:longest_y]
    #longest_b::Vector{UInt8} = psms[!,:longest_b]
    rt::Vector{Float32} = psms[!,:rt]
    #tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    log2_intensity_explained = psms[!,:log2_intensity_explained]::Vector{Float16}
    #precursor_idx = psms[!,:precursor_idx]::Vector{UInt32}
    function countMOX(seq::String)
        return UInt8(count("Unimod:35", seq))
    end

    function countMOX(seq::Missing)
        return zero(UInt8)
    end


    tasks_per_thread = 5
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin 
            for i in chunk
                prec_idx = precursor_idx[i]

                irt_obs[i] = rt_to_irt_interp(rt[i])
                irt_pred[i] = prec_irt[prec_idx]
                irt_diff[i] = abs(irt_obs[i] - first(prec_id_to_irt[prec_idx]))
                irt_error[i] = abs(irt_obs[i] - irt_pred[i])

                missed_cleavage[i] = precursor_missed_cleavage[prec_idx]
                #sequence[i] = precursor_sequence[prec_idx]
                sequence_length[i] = length(replace(precursor_sequence[prec_idx], r"\(.*?\)" => ""))#replace.(sequence[i], "M(ox)" => "M");
                Mox[i] = countMOX(structural_mods[prec_idx])::UInt8
                #sequence_length[i] = length(stripped_sequence[i])


                TIC[i] = Float16(log2(tic[scan_idx[i]]))
                adjusted_intensity_explained[i] = Float16(log2(TIC[i]) + log2_intensity_explained[i]);
                prec_charges[i] = prec_charge[prec_idx]
                #b_y_overlap[i] = ((sequence_length[i] - longest_y[i])>longest_b[i]) &  (longest_b[i] > 0) & (longest_y[i] > 0);

                spectrum_peak_count[i] = length(masses[scan_idx[i]])
         
                prec_mzs[i] = prec_mz[prec_idx];
            end
        end
    end
    fetch.(tasks)
    psms[!,:irt_obs] = irt_obs
    psms[!,:irt_pred] = irt_pred
    psms[!,:irt_diff] = irt_diff
    psms[!,:irt_error] = irt_error

    psms[!,:missed_cleavage] = missed_cleavage
    #psms[!,:sequence] = sequence
    #psms[!,:stripped_sequence] = stripped_sequence
    psms[!,:Mox] = Mox
    psms[!,:sequence_length] = sequence_length

    psms[!,:tic] = TIC
    psms[!,:adjusted_intensity_explained] = adjusted_intensity_explained
    psms[!,:charge] = prec_charges
    
    #psms[!,:b_y_overlap] = b_y_overlap
    psms[!,:spectrum_peak_count] = spectrum_peak_count
    psms[!,:prec_mz] = prec_mzs
    psms[!,:ms_file_idx] .= ms_file_idx
    return nothing
end

function getIsotopesCaptured!(chroms::DataFrame, 
                                quad_transmission_model::QuadTransmissionModel,
                                scan_idx::AbstractVector{UInt32},
                                prec_charge::AbstractArray{UInt8},
                                prec_mz::AbstractArray{Float32},
                                centerMz::AbstractVector{Union{Missing, Float32}},
                                isolationWidthMz::AbstractVector{Union{Missing, Float32}})
    #sum(MS2_CHROMS.weight.!=0.0)
    isotopes_captured = Vector{Tuple{Int8, Int8}}(undef, size(chroms, 1))
    
    tasks_per_thread = 5
    chunk_size = max(1, size(chroms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(chroms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                
                prec_id = chroms[i,:precursor_idx]
                mz = prec_mz[prec_id]
                charge = prec_charge[prec_id]

                scan_id = scan_idx[i]
                scan_mz = coalesce(centerMz[scan_id], zero(Float32))::Float32
                window_width = coalesce(isolationWidthMz[scan_id], zero(Float32))::Float32

                low_mz, high_mz = getQuadTransmissionBounds(quad_transmission_model, scan_mz, window_width)
                isotopes = getPrecursorIsotopeSet(mz, 
                                                    charge, 
                                                    low_mz, high_mz
                                                    )                
                isotopes_captured[i] = isotopes
            end
        end
        chroms[!,:isotopes_captured] = isotopes_captured
    end
    fetch.(tasks)
    return nothing
end

function initSummaryColumns!(
    psms::DataFrame,
    )

    new_cols = [
        (:max_entropy,              Float16)
        (:max_gof,         Float16)
        (:max_fitted_manhattan_distance,          Float16)
        (:max_fitted_spectral_contrast,         Float16)
        (:y_ions_sum,               UInt16)
        (:max_y_ions,               UInt16)
        (:max_matched_ratio,        Float16)
        ];

        N = size(psms, 1)
        for column in new_cols
            col_type = last(column);
            col_name = first(column)
            psms[!,col_name] = zeros(col_type, N)
        end
        return psms
end
