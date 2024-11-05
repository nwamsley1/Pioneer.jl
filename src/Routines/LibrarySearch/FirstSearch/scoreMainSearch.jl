function addMainSearchColumns!(psms::DataFrame, 
                                rt_irt::UniformSpline,
                                structural_mods::AbstractVector{Union{Missing, String}},
                                prec_missed_cleavages::Arrow.Primitive{UInt8, Vector{UInt8}},
                                prec_is_decoy::Arrow.BoolVector{Bool},
                                prec_irt::Arrow.Primitive{T, Vector{T}},
                                prec_charge::Arrow.Primitive{UInt8, Vector{UInt8}},
                                scan_retention_time::AbstractVector{Float32},
                                tic::AbstractVector{Float32},
                                masses::AbstractArray;
                                ) where {T<:AbstractFloat}
    
    ###########################
    #Allocate new columns
    N = size(psms, 1)
    missed_cleavage = zeros(UInt8, N);
    Mox = zeros(UInt8, N);
    irt_pred = zeros(Float32, N);
    rt = zeros(Float32, N);
    irt = zeros(Float32, N);
    TIC = zeros(Float16, N);
    charge = zeros(UInt8, N);
    err_norm = zeros(Float16, N);
    targets = zeros(Bool, N);
    spectrum_peak_count = zeros(UInt32, N);
    scan_idx::Vector{UInt32} = psms[!,:scan_idx]
    precursor_idx::Vector{UInt32} = psms[!,:precursor_idx] 
    error::Vector{Float32} = psms[!,:error]
    total_ions::Vector{UInt8} = psms[!,:y_count] .+ psms[!,:b_count]
    #tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    #scan_retention_time = MS_TABLE[:retentionTime]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    #masses = MS_TABLE[:mz_array]::Arrow.List{Union{Missing, SubArray{Union{Missing, Float32}, 1, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}, Tuple{UnitRange{Int64}}, true}}, Int64, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}}
    function countMOX(seq::String)
        return UInt8(count("Unimod:35", seq))
    end

    #Split data into chunk ranges
    tasks_per_thread = 10
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                prec_idx = precursor_idx[i]

                targets[i] = prec_is_decoy[prec_idx] == false;
                missed_cleavage[i] = prec_missed_cleavages[prec_idx]
                Mox[i] = countMOX(coalesce(structural_mods[prec_idx], ""))::UInt8 #UInt8(length(collect(eachmatch(r"ox",  precursors[precursor_idx[i]].sequence))))
                irt_pred[i] = Float32(prec_irt[prec_idx]);
                rt[i] = Float32(scan_retention_time[scan_idx[i]]);
                irt[i] = rt_irt(rt[i])
                TIC[i] = Float16(log2(tic[scan_idx[i]]));
                charge[i] = UInt8(prec_charge[prec_idx]);
                err_norm[i] = min(Float16.(abs(error[i]/total_ions[i])), Float16(6e4))#Float16(min(abs((error[i])/(total_ions[i])), 6e4))
                spectrum_peak_count[i] = UInt32(length(masses[scan_idx[i]]))
            end
        end
    end
    fetch.(tasks)
    psms[!,:irt_predicted] = irt_pred
    psms[!,:rt] = rt
    psms[!,:irt] = irt
    psms[!,:TIC] = TIC
    psms[!,:total_ions] = total_ions
    psms[!,:err_norm] = err_norm
    psms[!,:target] = targets
    psms[!,:missed_cleavage] = missed_cleavage
    psms[!,:Mox] = Mox
    psms[!,:charge] = charge
    psms[!,:spectrum_peak_count] = spectrum_peak_count
    psms[!,:score] = zeros(Float32, N);
    psms[!,:q_value] = zeros(Float16, N);
    psms[!,:intercept] = ones(Float16, N)
end

function scoreMainSearchpsms!(psms::DataFrame, column_names::Vector{Symbol};
                             n_train_rounds::Int64 = 2,
                             max_iter_per_round::Int64 = 50,
                             max_q_value::Float64 = 0.01)

    β = zeros(Float64, length(column_names));
    best_psms = nothing#ones(Bool, size(psms, 1))

    tasks_per_thread = 10
    M = size(psms, 1)
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:M, chunk_size) # partition your data into chunks that

    for i in range(1, n_train_rounds)
        if i < 2 #Train on all data during first round
            β = ProbitRegression(β, psms[!,column_names], psms[!,:target], data_chunks, max_iter = max_iter_per_round);
            ModelPredict!(psms[!,:score], psms[!,column_names], β, data_chunks); #Get Z-scores 
        else #Train using only decoys and high scoring targets from the previous round
            psms_targets = psms[best_psms,:target]
            M = size(psms_targets, 1)
            sub_chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
            sub_data_chunks = partition(1:M, sub_chunk_size) # partition your data into chunks that
            β = ProbitRegression(β, psms[best_psms,column_names], psms_targets, sub_data_chunks, max_iter = max_iter_per_round);
            ModelPredict!(psms[!,:score], psms[!,column_names], β, data_chunks); #Get Z-scores 
        end
        getQvalues!(psms[!,:score],psms[!,:target],psms[!,:q_value]);
        if i < n_train_rounds #Get Data to train on during subsequent round
            best_psms = ((psms[!,:q_value].<=max_q_value).&(psms[!,:target])) .| (psms[!,:target].==false);
        end
    end
    return 
end

function getProbs!(psms::DataFrame)
    psms[!,:prob] = zeros(Float32, size(psms, 1))
    tasks_per_thread = 10
    M = size(psms, 1)
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:M, chunk_size) # partition your data into chunks that
    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                psms[i,:prob] = (1 + SpecialFunctions.erf(psms[i,:score]/sqrt(2)))/2
            end
        end
    end
    fetch.(tasks)
    return nothing
end
