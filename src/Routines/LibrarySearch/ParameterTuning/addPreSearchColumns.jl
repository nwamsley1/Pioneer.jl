function addPreSearchColumns!(psms::DataFrame, 
                                MS_TABLE::Arrow.Table, 
                                prec_is_decoy::Arrow.BoolVector{Bool},
                                prec_irt::Arrow.Primitive{T, Vector{T}},
                                prec_charge::Arrow.Primitive{UInt8, Vector{UInt8}},
                                scan_retention_time::AbstractVector{Float32},
                                tic::AbstractVector{Float32}) where {T<:AbstractFloat}
    
    N = size(psms, 1)
    decoys = zeros(Bool, N);
    targets = zeros(Bool, N);
    TIC = zeros(Float16, N);
    charge = zeros(UInt8, N);
    spectrum_peak_count = zeros(UInt32, N);
    irt_pred = zeros(Float32, N);
    rt = zeros(Float32, N);
    scan_idx::Vector{UInt32} = psms[!,:scan_idx]
    precursor_idx::Vector{UInt32} = psms[!,:precursor_idx]
    matched_ratio::Vector{Float16} = psms[!,:matched_ratio]
    
    #tic = MS_TABLE[:TIC]::SentinelArrays.ChainedVector{Float32, Arrow.Primitive{Float32, Vector{Float32}}}
    #scan_retention_time = MS_TABLE[:retentionTime]::SentinelArrays.ChainedVector{Float32, Arrow.Primitive{Float32, Vector{Float32}}}

    tasks_per_thread = 10
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
            decoys[i] = prec_is_decoy[precursor_idx[i]];
            targets[i] = decoys[i] == false
            irt_pred[i] = Float32(prec_irt[precursor_idx[i]]);
            rt[i] = Float32(scan_retention_time[scan_idx[i]]);
            TIC[i] = Float16(log2(tic[scan_idx[i]]));
            charge[i] = UInt8(prec_charge[precursor_idx[i]]);
            matched_ratio[i] = Float16(min(matched_ratio[i], 6e4))
            end
        end
    end
    fetch.(tasks)
    psms[!,:matched_ratio] = matched_ratio
    psms[!,:decoy] = decoys
    psms[!,:irt_predicted] = irt_pred
    psms[!,:rt] = rt
    psms[!,:TIC] = TIC
    psms[!,:target] = targets
    psms[!,:charge] = charge
    psms[!,:spectrum_peak_count] = spectrum_peak_count
    psms[!,:intercept] = ones(Float16, size(psms, 1))
    psms[!,:q_value] = zeros(Float16, N);
    psms[!,:prob] = zeros(Float16, N);
end

function scorePresearch!(psms::DataFrame)
    
    features = [:entropy_score,:city_block,
                    :scribe,:spectral_contrast,
                    :y_count,:error,
                    :TIC,:intercept]
    psms[!,:prob] = zeros(Float32, size(psms, 1))

    M = size(psms, 1)
    tasks_per_thread = 10
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:M, chunk_size) # partition your data into chunks that
    β = zeros(Float64, length(features))
    β = ProbitRegression(β, psms[!,features], psms[!,:target], data_chunks, max_iter = 20)
    ModelPredictProbs!(psms[!,:prob], psms[!,features], β, data_chunks)
end