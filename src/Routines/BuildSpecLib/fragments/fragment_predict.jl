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

# src/fragments/fragment_predict.jl

"""
    SplineFragFilterCtx

Per-batch filter context plumbed through the Altimeter prediction path so
filters that today fire at the END of the build (in `getDetailedFrags`) can
be applied at the START (in `filter_fragments!`), shrinking the
intermediate `raw_fragments.arrow` / `fragments_table.arrow` files.

Bundles per-fragment metadata decoder (`annotation_lookup`), per-precursor
lookups (`prec_len`, `prec_mz`), and the full set of filter knobs that
`getDetailedFrags` consults. The struct is plumbed through commit 2 of §9
without altering behavior; commit 3 puts the lookups to use inside
`filter_fragments!`.
"""
struct SplineFragFilterCtx
    annotation_lookup::Dict{Int32, PioneerFragAnnotation}
    prec_len::Vector{UInt8}
    prec_mz::Vector{Float32}
    y_start::UInt8
    b_start::UInt8
    include_p::Bool
    include_isotope::Bool
    include_immonium::Bool
    include_internal::Bool
    include_neutral_diff::Bool
    max_frag_charge::UInt8
    max_frag_rank::UInt8
    length_to_frag_count_multiple::Float32
    min_frag_intensity::Float32
    frag_bounds::FragBoundModel
end

"""
    AgnosticFragFilterCtx

Per-batch filter context for the Prosit (`InstrumentAgnosticModel`) path. Unlike
the Altimeter spline path, Prosit returns raw *monoisotopic* scalar intensities
with STRING annotations ("y7+1"), so this context carries what the predict-time
filter needs to (1) parse string annotations (`annotation_cache`, filled lazily
via the `GenericFragAnnotation` parser), (2) convert monoisotopic->total using the
SAME `iso_splines` the search uses (`getFragIsotopes!`), and (3) rank by the
converted total and keep the per-precursor topN.

The conversion needs per-fragment sulfur, computed on the fly from the precursor
`sequences`/`mods` (mirroring the decode kernel `_fill_detailed_from_raw!`), so
`raw_fragments.arrow` is written already-topN with `intensities` = total abundance.
"""
struct AgnosticFragFilterCtx
    annotation_cache::Dict{String, PioneerFragAnnotation}
    iso_splines::IsotopeSplineModel{40, Float32}
    sequences::Vector{String}
    mods::Vector{String}
    mods_to_sulfur_diff::Dict{String, Int8}
    prec_len::Vector{UInt8}
    prec_mz::Vector{Float32}
    y_start::UInt8
    b_start::UInt8
    include_p::Bool
    include_isotope::Bool
    include_immonium::Bool
    include_internal::Bool
    include_neutral_diff::Bool
    max_frag_charge::UInt8
    max_frag_rank::UInt8
    length_to_frag_count_multiple::Float32
    min_frag_intensity::Float32
    frag_bounds::FragBoundModel
end

"""
    predict_fragments(
        peptide_table_path::String,
        frags_out_path::String,
        model_type::KoinaModelType,
        filter_ctx::SplineFragFilterCtx,
        max_koina_batches::Int,
        batch_size::Int,
        model_name::String;
        intensity_threshold::Float32 = 0.001f0
    )

Fragment prediction dispatcher for the Altimeter (SplineCoefficientModel) Koina endpoint.
"""
# ── Env-gated fragment-prediction sub-step diagnostics (PIONEER_DIAG_BUILD=1) ──
# Accumulates @timed (seconds, bytes) per labelled sub-step; printed at the end of
# predict_fragments. Zero behavior change; only active when the env var is set.
const _BUILD_DIAG = Dict{String, Vector{Float64}}()
@inline _bdiag_on() = get(ENV, "PIONEER_DIAG_BUILD", "0") != "0"
@inline function _bdiag!(label::String, r)
    if _bdiag_on()
        v = get!(_BUILD_DIAG, label, [0.0, 0.0])
        v[1] += r.time; v[2] += r.bytes
    end
    return r.value
end
function _bdiag_report()
    (_bdiag_on() && !isempty(_BUILD_DIAG)) || return
    @user_info "[BUILD DIAG] fragment-prediction sub-steps (cumulative):"
    for (k, v) in sort(collect(_BUILD_DIAG), by = x -> -x[2][2])
        @user_info "  $(rpad(k, 26)) $(round(v[1], digits=2)) s   $(round(v[2]/1e9, digits=2)) GB"
    end
    empty!(_BUILD_DIAG)
end

function predict_fragments(
    peptide_table_path::String,
    frags_out_path::String,
    model_type::KoinaModelType,
    filter_ctx,  # SplineFragFilterCtx (Altimeter) or AgnosticFragFilterCtx (Prosit); dispatched on model_type in predict_fragments_batch
    max_koina_batches::Int,
    batch_size::Int,
    model_name::String;
    intensity_threshold::Float32 = 0.001f0
)
    # Verify model configuration
    if !haskey(KOINA_URLS, model_name)
        error("Invalid model name: $model_name. Valid options: $(join(keys(KOINA_URLS), ", "))")
    end

    # Load data
    peptides_df = DataFrame(Arrow.Table(peptide_table_path))

    # Process in batches
    koina_pool_size = max_koina_batches * 5
    nprecs = nrow(peptides_df)
    batch_size = min(batch_size, 1000)
    batch_start_idxs = collect(one(UInt32):UInt32(batch_size*koina_pool_size):UInt32(nprecs))

    rm(frags_out_path, force=true)

    # Track knot_vector across batches (only set by SplineCoefficientModel)
    knot_vector = nothing

    for start_idx in ProgressBar(batch_start_idxs)
        stop_idx = min(start_idx + batch_size*koina_pool_size - 1, nrow(peptides_df))
        batch_df = peptides_df[start_idx:stop_idx, :]

        # Generate predictions for batch - returns tuple (DataFrame, knot_vector or nothing)
        (frags_out, batch_knot_vector) = predict_fragments_batch(
            batch_df,
            model_type,
            filter_ctx,
            batch_size,
            max_koina_batches,
            start_idx,
        )

        # Capture knot_vector from first batch (if model provides one)
        if knot_vector === nothing && batch_knot_vector !== nothing
            knot_vector = batch_knot_vector
        end

        # Write or append results
        if start_idx == 1
            # Create file in stream format with metadata (knot_vector stored once, not per-row)
            metadata = if knot_vector !== nothing
                ["knot_vector" => JSON.json(knot_vector)]
            else
                String[]
            end
            open(frags_out_path, "w") do io
                Arrow.write(io, frags_out; file=false, metadata=metadata)
            end
        else
            _bdiag!("arrow_append", @timed Arrow.append(frags_out_path, frags_out))
        end
    end
    _bdiag_report()
end

"""
Fragment prediction for spline coefficient models (e.g., Altimeter).
Returns a tuple of (DataFrame, Vector{Float32}) where the second element is the knot_vector.
"""
function predict_fragments_batch(
    peptides_df::DataFrame,
    model::SplineCoefficientModel,
    filter_ctx::SplineFragFilterCtx,
    batch_size::Int,
    concurrent_koina_requests::Int,
    first_prec_idx::UInt32
)::Tuple{DataFrame, Vector{Float32}}
    json_batches = _bdiag!("prepare_koina_batch", @timed prepare_koina_batch(
        model,
        peptides_df,
        batch_size=batch_size
    ))

    responses = _bdiag!("koina_request", @timed make_koina_batch_requests(json_batches, KOINA_URLS[model.name]; concurrency=concurrent_koina_requests))

    batch_dfs = []
    knot_vectors = []

    for (i, response) in enumerate(responses)
        batch_result = _bdiag!("parse_koina_batch", @timed parse_koina_batch(model, response))
        # Keep precursor_idx UInt32 (not Int64) — it is a precursor row index; the
        # +1-1 cancels. Halves this per-fragment column on disk (8 B -> 4 B).
        start_idx = first_prec_idx + UInt32((i-1) * batch_size)

        batch_df = batch_result.fragments
        n_precursors_in_batch = UInt32(fld(size( batch_df , 1), batch_result.frags_per_precursor))
        batch_df[!, :precursor_idx] = repeat(start_idx:(start_idx + n_precursors_in_batch - one(UInt32)),
                                                inner=batch_result.frags_per_precursor)
        _bdiag!("filter_fragments", @timed filter_fragments!(batch_df, model, filter_ctx))
        push!(batch_dfs, batch_df)
        push!(knot_vectors, batch_result.extra_data)  # Store knot vectors
    end

    fragments_df = _bdiag!("vcat_batch_dfs", @timed vcat(batch_dfs...))

    # Verify knot vectors are consistent
    if !all(k == first(knot_vectors) for k in knot_vectors)
        error("Inconsistent knot vectors across batches")
    end

    # Return DataFrame and knot_vector separately (knot_vector stored as Arrow metadata, not per-row)
    return (fragments_df, first(knot_vectors))
end

"""
    build_spline_frag_filter_ctx(frag_bounds, precursors_df, ion_dict_path,
                                  immonium_data_path, frag_annotation_type;
                                  y_start, b_start, include_p, include_isotope,
                                  include_immonium, include_internal,
                                  include_neutral_diff, max_frag_charge,
                                  max_frag_rank, length_to_frag_count_multiple,
                                  min_frag_intensity) -> SplineFragFilterCtx

Construct the per-batch filter context used by `predict_fragments`.

The filter knobs are passed as keyword arguments rather than read from the
JSON config because BuildSpecLib hardcodes them at the call site (the same
constants flow into `buildPionLib` for the partitioned-index step). Accepting
them here as kwargs keeps the build-site as the single source of truth.
"""
function build_spline_frag_filter_ctx(
    frag_bounds::FragBoundModel,
    precursors_df::DataFrame,
    ion_dict_path::String,
    immonium_data_path::String,
    frag_annotation_type;
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8,
    max_frag_rank::UInt8,
    length_to_frag_count_multiple::Float32,
    min_frag_intensity::Float32,
)::SplineFragFilterCtx
    ion_dictionary = get_altimeter_ion_dict(ion_dict_path)
    immonium_to_sulfur_count = get_immonium_sulfur_dict(immonium_data_path)

    annotation_lookup = Dict{Int32, PioneerFragAnnotation}()
    annotation_type = typeof(frag_annotation_type)
    for (ion_idx, ion_name) in ion_dictionary
        instance = annotation_type(ion_name)
        annotation_lookup[ion_idx] = parse_fragment_annotation(
            instance;
            immonium_to_sulfur_count=immonium_to_sulfur_count,
        )
    end

    # prec_len + prec_mz are looked up by precursor_idx (1-based), so build
    # vectors sized to the table.
    n_precs = nrow(precursors_df)
    prec_len = Vector{UInt8}(undef, n_precs)
    prec_mz  = Vector{Float32}(undef, n_precs)
    seqs = precursors_df[!, :sequence]
    mzs  = precursors_df[!, :mz]
    @inbounds for i in 1:n_precs
        prec_len[i] = UInt8(min(255, length(seqs[i])))
        prec_mz[i]  = Float32(mzs[i])
    end

    return SplineFragFilterCtx(
        annotation_lookup, prec_len, prec_mz,
        y_start, b_start,
        include_p, include_isotope, include_immonium, include_internal,
        include_neutral_diff,
        max_frag_charge, max_frag_rank,
        length_to_frag_count_multiple, min_frag_intensity,
        frag_bounds,
    )
end

"""
    filter_fragments!(df, ::SplineCoefficientModel, ctx::SplineFragFilterCtx)

Apply per-fragment metadata filters and per-precursor rank caps in-place,
mirroring the (now-deleted) logic that previously ran inside
`getDetailedFrags()` / `fragFilter()` / `filterFrag()`. By filtering at
predict time rather than build time, the on-disk `raw_fragments.arrow` and
`fragments_table.arrow` contain only the fragments that will survive into
`detailed_fragments.jls` — typically ~12% of the predicted total.

Per-fragment filters (decided from `ctx.annotation_lookup[annotation]`):
  - mz > 0 (invalid Altimeter outputs sometimes carry zero/negative mz)
  - frag_bounds(prec_mz) m/z window
  - y-ion ⇒ frag_index ≥ y_start
  - b-ion ⇒ frag_index ≥ b_start
  - precursor ion ⇒ include_p
  - immonium ⇒ include_immonium
  - internal ⇒ include_internal
  - neutral-loss ⇒ include_neutral_diff
  - isotope > 0 ⇒ include_isotope
  - frag_charge ≤ max_frag_charge

Per-precursor rank cap (Altimeter returns fragments in importance order;
this preserves that ordering and keeps the first N per precursor):
  N = min(max_frag_rank, round(prec_len * length_to_frag_count_multiple) + 1)
"""
function filter_fragments!(df::DataFrame, model::SplineCoefficientModel,
                           ctx::SplineFragFilterCtx)
    n_rows = nrow(df)
    n_rows == 0 && return df

    # Concrete column types: df[!, :col] is inferred as AbstractVector, so col[i]
    # in the hot per-fragment loop below returns Any and boxes every read. The
    # ::Vector{T} assertions restore type stability (and act as a schema guard).
    annotations = df[!, :annotation]::Vector{Int32}
    mzs         = df[!, :mz]::Vector{Float32}
    pids        = df[!, :precursor_idx]::Vector{UInt32}

    keep = trues(n_rows)

    # --- Pass 1: per-fragment metadata filters ---
    @inbounds for i in 1:n_rows
        if mzs[i] <= 0
            keep[i] = false
            continue
        end

        pid = pids[i]
        prec_mz = ctx.prec_mz[pid]
        lo, hi = ctx.frag_bounds(prec_mz)
        if mzs[i] < lo || mzs[i] > hi
            keep[i] = false
            continue
        end

        ann = ctx.annotation_lookup[annotations[i]]
        bt = ann.base_type
        if bt == 'y' && ann.frag_index < ctx.y_start
            keep[i] = false; continue
        end
        if bt == 'b' && ann.frag_index < ctx.b_start
            keep[i] = false; continue
        end
        if bt == 'p' && !ctx.include_p
            keep[i] = false; continue
        end
        if ann.immonium && !ctx.include_immonium
            keep[i] = false; continue
        end
        if ann.internal && !ctx.include_internal
            keep[i] = false; continue
        end
        if ann.neutral_diff && !ctx.include_neutral_diff
            keep[i] = false; continue
        end
        if !ctx.include_isotope && ann.isotope > 0
            keep[i] = false; continue
        end
        if ann.charge > ctx.max_frag_charge
            keep[i] = false; continue
        end
    end

    # --- Pass 2: per-precursor rank cap (preserves Altimeter's input order) ---
    rank_per_pid = Dict{UInt32, Int}()
    @inbounds for i in 1:n_rows
        keep[i] || continue
        pid = pids[i]
        prec_len = ctx.prec_len[pid]
        max_n = min(Int(ctx.max_frag_rank),
                    round(Int, Float32(prec_len) * ctx.length_to_frag_count_multiple) + 1)
        cur = get(rank_per_pid, pid, 0) + 1
        rank_per_pid[pid] = cur
        if cur > max_n
            keep[i] = false
        end
    end

    deleteat!(df, findall(.!keep))
    return df
end


"""
Sort fragments by intensity within each precursor group.
"""
function sort_fragments!(df::DataFrame)
    sort!(df, [:precursor_idx, order(:intensities, rev=true)])
end

# ---------------------------------------------------------------------------
# Prosit (InstrumentAgnosticModel) prediction path
# ---------------------------------------------------------------------------

"""
Fragment prediction for instrument-agnostic scalar-intensity models (Prosit).
Returns `(DataFrame, nothing)` — no knot vector (contrast the spline path).
"""
function predict_fragments_batch(
    peptides_df::DataFrame,
    model::InstrumentAgnosticModel,
    filter_ctx::AgnosticFragFilterCtx,
    batch_size::Int,
    concurrent_koina_requests::Int,
    first_prec_idx::UInt32
)::Tuple{DataFrame, Nothing}
    json_batches = prepare_koina_batch(model, peptides_df, batch_size=batch_size)
    responses = make_koina_batch_requests(json_batches, KOINA_URLS[model.name];
                                          concurrency=concurrent_koina_requests)

    batch_dfs = DataFrame[]
    for (i, response) in enumerate(responses)
        batch_result = parse_koina_batch(model, response)
        start_idx = first_prec_idx + UInt32((i-1) * batch_size)
        batch_df = batch_result.fragments
        n_precursors_in_batch = UInt32(fld(nrow(batch_df), batch_result.frags_per_precursor))
        batch_df[!, :precursor_idx] = repeat(start_idx:(start_idx + n_precursors_in_batch - one(UInt32)),
                                             inner=batch_result.frags_per_precursor)
        batch_df = filter_fragments!(batch_df, model, filter_ctx)  # returns reordered topN table
        push!(batch_dfs, batch_df)
    end

    return (vcat(batch_dfs...), nothing)
end

"""
    build_agnostic_frag_filter_ctx(frag_bounds, precursors_df, immonium_data_path,
        mods_to_sulfur_diff; kwargs...) -> AgnosticFragFilterCtx

Construct the Prosit filter context. Loads the isotope splines (the SAME model
search uses, `assets/IsotopeSplines_10kDa_21isotopes.xml`) so the predict-time
monoisotopic->total conversion matches `getFragIsotopes!` exactly. Filter knobs
are passed as kwargs (BuildSpecLib hardcodes them at the call site, single source
of truth with `buildPionLib`).
"""
function build_agnostic_frag_filter_ctx(
    frag_bounds::FragBoundModel,
    precursors_df::DataFrame,
    mods_to_sulfur_diff::Dict{String, Int8};
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8,
    max_frag_rank::UInt8,
    length_to_frag_count_multiple::Float32,
    min_frag_intensity::Float32,
)::AgnosticFragFilterCtx
    iso_splines = parseIsoXML(isotope_spline_path())

    n_precs = nrow(precursors_df)
    prec_len = Vector{UInt8}(undef, n_precs)
    prec_mz  = Vector{Float32}(undef, n_precs)
    seqs = String.(precursors_df[!, :sequence])
    mods = String.(precursors_df[!, :mods])
    mzs  = precursors_df[!, :mz]
    @inbounds for i in 1:n_precs
        prec_len[i] = UInt8(min(255, length(seqs[i])))
        prec_mz[i]  = Float32(mzs[i])
    end

    return AgnosticFragFilterCtx(
        Dict{String, PioneerFragAnnotation}(),
        iso_splines, seqs, mods, mods_to_sulfur_diff,
        prec_len, prec_mz,
        y_start, b_start,
        include_p, include_isotope, include_immonium, include_internal,
        include_neutral_diff,
        max_frag_charge, max_frag_rank,
        length_to_frag_count_multiple, min_frag_intensity,
        frag_bounds,
    )
end

# Parse (and cache) a Prosit string annotation into a PioneerFragAnnotation.
# Accepts AbstractString so Arrow-backed columns (SubString/String) work at decode.
@inline function _agnostic_annotation!(cache::Dict{String, PioneerFragAnnotation}, s::AbstractString)
    pa = get(cache, s, nothing)
    pa === nothing || return pa
    pa = parse_fragment_annotation(GenericFragAnnotation(String(s)))
    cache[s] = pa
    return pa
end

"""
    filter_fragments!(df, ::InstrumentAgnosticModel, ctx::AgnosticFragFilterCtx)

Apply, per precursor and in place: (1) per-fragment metadata filters (mz>0,
frag_bounds, y/b start, charge, include_* flags, min intensity), (2) the
monoisotopic->total conversion `total = mono / iso_splines(min(sulfur,5), 0,
mz*charge)` — the SAME isotope model search uses — computing per-fragment sulfur
from the precursor sequence span, then (3) a per-precursor topN cap ranked by the
converted total. On return `df[:intensities]` holds the total abundance for the
surviving topN fragments (search reconstructs the monoisotopic peak exactly).

`df` is assumed grouped by ascending `precursor_idx` (parse_koina_batch emits
174 contiguous fragments per precursor), so per-precursor blocks are contiguous.
"""
function filter_fragments!(df::DataFrame, model::InstrumentAgnosticModel,
                           ctx::AgnosticFragFilterCtx)
    n_rows = nrow(df)
    n_rows == 0 && return df

    annotations = df[!, :annotation]::Vector{String}
    mzs         = df[!, :mz]::Vector{Float32}
    ints        = df[!, :intensities]::Vector{Float32}
    pids        = df[!, :precursor_idx]::Vector{UInt32}

    keep  = trues(n_rows)
    total = Vector{Float32}(undef, n_rows)
    # Captured for site-determining gap retention in Pass 2 (only kept rows are
    # read, via `block`). base_type/frag ordinal per fragment.
    base_types = Vector{Char}(undef, n_rows)
    frag_ords  = Vector{UInt8}(undef, n_rows)

    # --- Pass 1: metadata filters + sulfur + mono->total (per-precursor state) ---
    seq_idx_to_sulfur = zeros(UInt8, 255)
    last_pid = zero(UInt32)
    prec_sulfur = 0
    prec_length = zero(UInt8)
    @inbounds for i in 1:n_rows
        pid = pids[i]
        if pid != last_pid
            last_pid = pid
            fill!(seq_idx_to_sulfur, zero(UInt8))
            seq = ctx.sequences[pid]
            prec_sulfur = Int(count_sulfurs!(seq_idx_to_sulfur, seq,
                                             parseMods(ctx.mods[pid]),
                                             ctx.mods_to_sulfur_diff))
            prec_length = ctx.prec_len[pid]
        end

        if mzs[i] <= 0 || ints[i] <= ctx.min_frag_intensity
            keep[i] = false; continue
        end
        lo, hi = ctx.frag_bounds(ctx.prec_mz[pid])
        if mzs[i] < lo || mzs[i] > hi
            keep[i] = false; continue
        end

        pa = _agnostic_annotation!(ctx.annotation_cache, annotations[i])
        bt = pa.base_type
        base_types[i] = bt
        frag_ords[i]  = UInt8(min(255, pa.frag_index))
        if bt == 'y' && pa.frag_index < ctx.y_start; keep[i] = false; continue; end
        if bt == 'b' && pa.frag_index < ctx.b_start; keep[i] = false; continue; end
        if bt == 'p' && !ctx.include_p; keep[i] = false; continue; end
        if pa.immonium && !ctx.include_immonium; keep[i] = false; continue; end
        if pa.internal && !ctx.include_internal; keep[i] = false; continue; end
        if pa.neutral_diff && !ctx.include_neutral_diff; keep[i] = false; continue; end
        if !ctx.include_isotope && pa.isotope > 0; keep[i] = false; continue; end
        if pa.charge > ctx.max_frag_charge; keep[i] = false; continue; end

        # per-fragment sulfur over the covered sequence span (skip immonium)
        s_idx, e_idx = get_fragment_indices(bt, pa.frag_index, prec_length)
        sc = 0
        if !pa.immonium
            for k in s_idx:e_idx
                sc += seq_idx_to_sulfur[k]
            end
        end
        sc += pa.sulfur_diff
        sc = min(sc, prec_sulfur)

        f0 = ctx.iso_splines(min(sc, 5), 0, mzs[i] * pa.charge)
        total[i] = f0 > 0f0 ? ints[i] / f0 : ints[i]
    end

    # --- Pass 2: per-precursor topN, REORDERED by converted total so the decode's
    # positional rank equals the intensity rank (matches the spline path, where
    # Altimeter returns fragments already in importance order). Blocks are appended
    # in ascending precursor_idx order, so the result stays precursor-grouped. ---
    final_order = Int[]
    i = 1
    @inbounds while i <= n_rows
        pid = pids[i]
        j = i
        while j <= n_rows && pids[j] == pid
            j += 1
        end
        block = [k for k in i:(j-1) if keep[k]]
        sort!(block, by = k -> -total[k])
        max_n = min(Int(ctx.max_frag_rank),
                    round(Int, Float32(ctx.prec_len[pid]) * ctx.length_to_frag_count_multiple) + 1)
        for t in 1:min(length(block), max_n)
            push!(final_order, block[t])
        end
        i = j
    end

    # Emit a precursor-grouped, intensity-descending, topN table with `intensities`
    # holding total abundance (search reconstructs the monoisotopic peak). Returned
    # (not in-place) — predict_fragments_batch rebinds to the result.
    out = df[final_order, :]
    out[!, :intensities] = Float32[total[k] for k in final_order]
    return out
end

