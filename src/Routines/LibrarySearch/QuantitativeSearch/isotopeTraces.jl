abstract type IsotopeTraceType end
struct CombineTraces <: IsotopeTraceType
    min_ft::Float32
end

seperateTraces(itt::CombineTraces) = false

function getPsmGroupbyCols(itt::CombineTraces)
    return [:precursor_idx]
end



function getIsotopesCaptured!(chroms::DataFrame, 
                                isotope_trace_type::CombineTraces,
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
                if first(isotopes) >= 2         
                    isotopes_captured[i] = isotopes
                else
                    isotopes_captured[i] = (Int8(-1), Int8(-1))
                end
            end
        end
    end
    fetch.(tasks)
    chroms[!,:isotopes_captured] = isotopes_captured
    return nothing
end



struct SeperateTraces <: IsotopeTraceType
end
 
seperateTraces(itt::SeperateTraces) = true

function getPsmGroupyCols(itt::SeperateTraces)
    return [:precursor_idx,:isotopes_captured]
end

function getIsotopesCaptured!(chroms::DataFrame, 
                                isotope_trace_type::SeperateTraces,
                                quad_transmission_model::QuadTransmissionModel,
                                scan_idx::AbstractVector{UInt32},
                                prec_charge::AbstractArray{UInt8},
                                prec_mz::AbstractArray{Float32},
                                centerMz::AbstractVector{Union{Missing, Float32}},
                                isolationWidthMz::AbstractVector{Union{Missing, Float32}})
    #sum(MS2_CHROMS.weight.!=0.0)
    chroms[!,:isotopes_captured] = Vector{Tuple{Int8, Int8}}(undef, size(chroms, 1))
    
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
                chroms[i,:isotopes_captured] = isotopes
            end
        end
    end
    fetch.(tasks)
    return nothing
end
