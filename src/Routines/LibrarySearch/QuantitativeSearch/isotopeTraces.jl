abstract type IsotopeTraceType end
struct CombineTraces <: IsotopeTraceType
    min_ft::Float32
end

seperateTraces(itt::CombineTraces) = false

function fractionTransmitted(
    isotope_trace_type::CombineTraces,
    mz::Float32,
    charge::UInt8,
    n_sulfur::UInt8,
    qtf::QuadTransmissionFunction
)
    ft = zero(Float32)
    iso_mz = mz
    prec_mass = Float32((mz - PROTON)*charge)
    iso_mass = prec_mass
    n_sulfur = min(Int64(n_sulfur) ,5)
    for iso_idx in range(0, 5)
        ft += iso_splines(n_sulfur, iso_idx, iso_mass)*qtf(iso_mz)
        iso_mass += Float32(NEUTRON)
        iso_mz += Float32(NEUTRON/charge)
    end
    return ft>=isotope_trace_type.min_ft 
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



struct SeperateTraces <: IsotopeTraceType
end
 
seperateTraces(itt::SeperateTraces) = true

function fractionTransmitted(
    isotope_trace_type::SeperateTraces,
    mz::Float32,
    charge::UInt8,
    n_sulfur::UInt8,
    qtf::QuadTransmissionFunction
)
    return true
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
