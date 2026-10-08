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

# Streaming fragment stage: predict -> filter -> decode one batch of precursors at a time, appending the m/z-sorted
# fragments to detailed_fragments.bin and keeping only each precursor's fragment-index fragments (the
# IndexFragSelection the indexes are built from). Replaces predict_fragments + build_detailed_frags_from_raw +
# the sort/write of buildPionLib for spline models, without raw_fragments.arrow or the whole fragment table in memory.

"""
    stream_spline_fragments(precursors_path, lib_dir, model, filter_ctx, ion_dictionary, annotation_type,
                            immonium_path, mods_to_sulfur_diff, iso_mod_to_mass;
                            index_filters, batch_precs = 1_000_000, koina_batch = 1000, concurrency = 24)
        -> (selection, n_frags)

Writes `detailed_fragments.bin`, `spline_knots.jls`, `frag_name_to_idx.jls` and `ion_annotations.jls` to
`lib_dir`, with the same contents as the in-memory path, and returns the fragment-index selection
(`index_filters` = (y_start_index, b_start_index, include_p_index)) and the number of fragments.
"""
function stream_spline_fragments(precursors_path::String, lib_dir::String, model::SplineCoefficientModel,
                                 filter_ctx::SplineFragFilterCtx, ion_dictionary::Dict{Int32, String},
                                 annotation_type::FragAnnotation, immonium_path::String,
                                 mods_to_sulfur_diff::Dict{String, Int8}, iso_mod_to_mass::Dict{String, Float32};
                                 index_filters::Tuple{UInt8, UInt8, Bool}, batch_precs::Int = 1_000_000,
                                 koina_batch::Int = 1000, concurrency::Int = 24)
    prec = Arrow.Table(precursors_path)
    n_prec = length(prec.mz)
    # annotation lookup, as build_detailed_frags_from_raw builds it
    immonium_to_sulfur_count = get_immonium_sulfur_dict(immonium_path)
    annotations = Dict{Int32, PioneerFragAnnotation}()
    atype = typeof(annotation_type)
    for (ion_idx, ion_name) in ion_dictionary
        annotations[ion_idx] = parse_fragment_annotation(atype(ion_name); immonium_to_sulfur_count = immonium_to_sulfur_count)
    end
    prec_mzs = Vector{Float32}(prec.mz)
    prec_irts = Vector{Float32}(prec.irt)
    cols = (prec.sequence, prec.mods, prec.isotope_mods, prec.precursor_charge)
    sel_frags = SimpleFrag{Float32}[]
    sel_offsets = Int[1]
    ranges = Vector{UInt64}(undef, n_prec + 1)
    bin_path = joinpath(lib_dir, "detailed_fragments.bin")
    writer = nothing
    knots = nothing
    t_predict = 0.0; t_filter = 0.0; t_decode = 0.0; t_write = 0.0
    for lo in 1:batch_precs:n_prec
        hi = min(lo + batch_precs - 1, n_prec)
        t = time()
        input = DataFrame(koina_sequence = prec.koina_sequence[lo:hi], precursor_charge = prec.precursor_charge[lo:hi])
        results = koina_batch_results(model, input, KOINA_URLS[model.name]; batch_size = koina_batch,
                                      concurrency = concurrency)
        t_predict += time() - t; t = time()
        frags_df = _filter_spline_batches!(results, model, filter_ctx, lo, koina_batch)
        batch_knots = first(results).extra_data
        all(r -> r.extra_data == batch_knots, results) || error("Inconsistent knot vectors across batches")
        knots === nothing ? (knots = batch_knots) : (knots == batch_knots || error("Inconsistent knot vectors across batches"))
        t_filter += time() - t; t = time()
        detailed, local_ranges = _decode_spline_batch(frags_df, cols, annotations, mods_to_sulfur_diff, iso_mod_to_mass,
                                                      hi - lo + 1, lo - 1)
        _select_index_frags!(sel_frags, sel_offsets, detailed, local_ranges, lo - 1, prec_mzs, prec_irts, index_filters)
        sort_detailed_fragments_by_mz!(detailed, local_ranges)
        t_decode += time() - t; t = time()
        writer === nothing && (writer = DetailedFragsWriter{eltype(detailed)}(bin_path))
        base = UInt64(writer.n_frags)
        for k in 1:(hi - lo + 1)
            ranges[lo - 1 + k] = base + local_ranges[k]
        end
        append_frags!(writer, detailed)
        t_write += time() - t
    end
    writer === nothing && (writer = DetailedFragsWriter{SplineCompactFrag{4, Float32}}(bin_path))   # no precursors
    ranges[n_prec + 1] = UInt64(writer.n_frags + 1)
    finish_detailed_frags!(writer, ranges)
    serialize_to_jls(joinpath(lib_dir, "spline_knots.jls"), knots === nothing ? Float32[] : knots)
    serialize_to_jls(joinpath(lib_dir, "frag_name_to_idx.jls"), ion_dictionary)
    serialize_to_jls(joinpath(lib_dir, "ion_annotations.jls"), annotations)
    @user_info @sprintf("Streaming fragments: %d precursors, %d fragments (predict %.1f s, filter %.1f s, decode %.1f s, write %.1f s)",
                        n_prec, writer.n_frags, t_predict, t_filter, t_decode, t_write)
    return IndexFragSelection(sel_frags, sel_offsets, prec_mzs), writer.n_frags
end

"Number the predictions of each Koina batch (global precursor ids from `first_pid`) and filter them; one table."
function _filter_spline_batches!(results::Vector{KoinaBatchResult{Vector{Float32}}}, model::SplineCoefficientModel,
                                 filter_ctx::SplineFragFilterCtx, first_pid::Int, koina_batch::Int)
    dfs = Vector{DataFrame}(undef, length(results))
    Threads.@threads :dynamic for i in eachindex(results)
        dfs[i] = _filter_spline_batch!(results[i], model, filter_ctx, UInt32(first_pid + (i - 1) * koina_batch))
    end
    return vcat(dfs...)
end
function _filter_spline_batch!(r::KoinaBatchResult{Vector{Float32}}, model::SplineCoefficientModel,
                               filter_ctx::SplineFragFilterCtx, start::UInt32)
    df = r.fragments
    n = UInt32(fld(nrow(df), r.frags_per_precursor))
    df[!, :precursor_idx] = repeat(start:(start + n - one(UInt32)), inner = r.frags_per_precursor)
    filter_fragments!(df, model, filter_ctx)
    return df
end

"Decode one batch's filtered predictions (rank order) and their batch-local precursor -> fragment ranges."
function _decode_spline_batch(frags_df::DataFrame, cols::Tuple, annotations::Dict{Int32, PioneerFragAnnotation},
                              mods_to_sulfur_diff::Dict{String, Int8}, iso_mod_to_mass::Dict{String, Float32},
                              n_precs::Int, pid_offset::Int)
    coef = frags_df[!, :coefficients]
    return _decode_spline_batch(frags_df[!, :annotation], frags_df[!, :precursor_idx], frags_df[!, :mz], coef,
                                cols, annotations, mods_to_sulfur_diff, iso_mod_to_mass, n_precs, pid_offset)
end
function _decode_spline_batch(ann::AbstractVector, pids::AbstractVector, mzs::AbstractVector,
                              coef::AbstractVector{NTuple{N, T}}, cols::Tuple,
                              annotations::Dict{Int32, PioneerFragAnnotation}, mods_to_sulfur_diff::Dict{String, Int8},
                              iso_mod_to_mass::Dict{String, Float32}, n_precs::Int, pid_offset::Int) where {N, T}
    n_frags = length(mzs)
    detailed = Vector{SplineCompactFrag{N, T}}(undef, n_frags)
    local_ranges = zeros(UInt64, n_precs + 1)
    _fill_detailed_from_raw!(detailed, local_ranges, ann, pids, mzs, coef, cols[1], cols[2], cols[3], cols[4],
                             annotations, mods_to_sulfur_diff, iso_mod_to_mass, n_frags, n_precs, pid_offset)
    return detailed, local_ranges
end

"Append the batch's fragment-index fragments (select_index_fragments's rule) to the selection."
function _select_index_frags!(sel_frags::Vector{SimpleFrag{Float32}}, sel_offsets::Vector{Int},
                              detailed::Vector{F}, local_ranges::Vector{UInt64}, pid_offset::Int,
                              prec_mzs::Vector{Float32}, prec_irts::Vector{Float32},
                              filt::Tuple{UInt8, UInt8, Bool}) where {F}
    for k in 1:(length(local_ranges) - 1)
        pid = pid_offset + k
        push_frag = SelectionPush(sel_frags, UInt32(pid), prec_mzs[pid], prec_irts[pid])
        n = _visit_index_frags(push_frag, detailed, Int(local_ranges[k]):(Int(local_ranges[k + 1]) - 1), filt)
        push!(sel_offsets, sel_offsets[end] + n)
    end
    return nothing
end

"_visit_index_frags callback: one precursor's index fragment as select_index_fragments stores it."
struct SelectionPush
    sel_frags::Vector{SimpleFrag{Float32}}
    pid::UInt32
    prec_mz::Float32
    prec_irt::Float32
end
@inline (p::SelectionPush)(rank::Int, dfrag) =
    push!(p.sel_frags, SimpleFrag{Float32}(getMz(dfrag), p.pid, p.prec_mz, p.prec_irt, UInt8(0), UInt8(1) << UInt8(rank - 1)))
