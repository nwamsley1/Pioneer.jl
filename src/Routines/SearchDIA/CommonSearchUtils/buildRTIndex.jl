# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

function buildRtIndex(rts::Vector{T}, prec_mzs::Vector{U}, prec_ids::Vector{I}, bin_rt_size::AbstractFloat) where {T,U<:AbstractFloat,I<:Integer}

    # Handle empty input gracefully
    if isempty(rts)
        return retentionTimeIndex(T, U)  # Return empty index
    end

    start_idx = 1
    start_rt =  rts[start_idx]
    rt_index = retentionTimeIndex(T, U) #Initialize retention time index
    i = 1
    while i < length(rts) + 1
        if ((rts[min(i + 1, length(rts))] - start_rt) > bin_rt_size) | (i == length(rts))
            push!(rt_index.rt_bins, 
                    rtIndexBin(rts[start_idx], #Retention time for first precursor in the bin
                          rts[i],     #Retention time for last precursor in the bin
                        [(zero(UInt32), zero(Float32)) for _ in 1:(i - start_idx + 1)] #Pre-allocate precursors 
                        )
                )

            n = 1 #n'th precursor 
            for idx in start_idx:(min(i, length(rts))) 
                rt_index.rt_bins[end].prec[n] = (prec_ids[idx], prec_mzs[idx]) #Add n'th precursor
                n += 1
            end

            sort!(rt_index.rt_bins[end].prec, by = x->last(x)) #Sort precursors by m/z
            i += 1
            start_idx = i
            start_rt = rts[min(start_idx, length(rts))]
            continue
        else
            i += 1
        end
    end


    function sortrtBins!(rt_index::retentionTimeIndex{T, U})
        for i in 1:length(rt_index.rt_bins)
            sort!(rt_index.rt_bins[i].prec, by = x->last(x));
        end
        return nothing
    end
    sortrtBins!(rt_index)
    return rt_index
end

buildRtIndex(PSMs::DataFrame; bin_rt_size::AbstractFloat = 0.1) = buildRtIndex(PSMs[:,:irt], PSMs[:,:prec_mz], PSMs[:,:precursor_idx], bin_rt_size)

buildRtIndex(PSMs::SubDataFrame; bin_rt_size::AbstractFloat = 0.1) = buildRtIndex(PSMs[:,:irt], PSMs[:,:prec_mz], PSMs[:,:precursor_idx], bin_rt_size)


function makeRTIndices(temp_folder::String,
                       psms_paths::Vector{String},
                       prec_to_irt::Dictionary{UInt32, @NamedTuple{irt::Float32, mz::Float32}},
                       rt_to_library_irt_splines::Any;
                       min_prob::AbstractFloat = 0.5)

    # Maps filepath to a retentionTimeIndex (see buildrtIndex.jl)
    rt_index_paths = Vector{String}(undef, length(psms_paths))
    # Fill retention time index for each file
    for (key, psms_path) in enumerate(psms_paths)
        psms = Arrow.Table(psms_path)
        rt_to_library_irt = rt_to_library_irt_splines[key]
        # Impute empirical library iRT value for psms with probability lower than the threshold
        irts = zeros(Float32, length(prec_to_irt))
        mzs = zeros(Float32, length(prec_to_irt))
        prec_ids = zeros(UInt32, length(prec_to_irt))
        # Map observed precursors to library iRT and probability score
        prec_set = Dict(zip(
            psms[:precursor_idx],
            map(x->(irt=first(x),prob=last(x)), zip(rt_to_library_irt.(psms[:rt]), psms[:prob]))
        ))

        Threads.@threads for (i, (prec_id, irt_mz)) in collect(enumerate(pairs(prec_to_irt)))
            prec_ids[i] = prec_id
            irt, mz = irt_mz::@NamedTuple{irt::Float32, mz::Float32}
            # Don't impute library iRT, use empirical
            if haskey(prec_set, prec_id)
                _irt_, prob = prec_set[prec_id]
                if (prob >= min_prob)
                    irts[i], mzs[i]  = _irt_, mz
                    continue
                end
            end
            # Impute library iRT from the best observed psm for the precursor across the experiment
            irts[i], mzs[i] = irt,mz
        end
        # Build RT index using library iRT values
        rt_df = DataFrame(Dict(:irt => irts,
                                :prec_mz => mzs,
                                :precursor_idx => prec_ids))
        sort!(rt_df, :irt)
        temp_path =joinpath(temp_folder, string(key)*"_rt_indices.arrow")
        Arrow.write(
            temp_path,
            rt_df,
            )
        rt_index_paths[key] = temp_path
    end
    return rt_index_paths
end
