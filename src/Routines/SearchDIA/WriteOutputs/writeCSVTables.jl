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



"""
Quantification columns whose zero values are a sentinel, not a measurement.

`integrate_chrom` writes `0.0f0` when a precursor could not be quantified in a run --
either integration found no peak, or `QUANT_MIN_AREA_SURVIVING_RATIO` withheld the area
because too little survived baseline subtraction. Downstream already reads zero as "not
quantified": `getS` maps any `abundance <= 0` to `missing` so MaxLFQ never ingests one.

Rendering that sentinel as `0.0` in the output tables is misleading -- it is
indistinguishable from a real measurement of zero. Blank it at the OUTPUT BOUNDARY only.
The chunk files keep plain `Float32` zeros, so MaxLFQ's input is byte-for-byte unchanged
and no `Union{Missing,Float32}` is introduced into any carried column.

`peak_area_normalized` carries the same sentinel when run-to-run normalization is enabled,
so it is blanked too. When normalization is *disabled* the column is never computed at all
(`ProteinQuantificationSearch` selects `:peak_area` as the quant column); it is then dropped from the output
by `drop_uncomputed_normalized` rather than emitted as an empty column.
"""
const OUTPUT_BLANKED_QUANT_COLUMNS = (:peak_area, :peak_area_normalized)

function _blank_unquantified(col::AbstractVector)
    out = Vector{Union{Missing, Float32}}(undef, length(col))
    @inbounds for i in eachindex(col)
        v = col[i]
        out[i] = (ismissing(v) || v <= 0) ? missing : Float32(v)
    end
    return out
end

"""Blank sentinel zero areas in a DataFrame, in place (TSV path)."""
function blank_unquantified_areas!(df)
    for c in OUTPUT_BLANKED_QUANT_COLUMNS
        hasproperty(df, c) && (df[!, c] = _blank_unquantified(df[!, c]))
    end
    return df
end

"""
    drop_uncomputed_normalized(cols, normalized::Bool)

`peak_area_normalized` is only populated when run-to-run normalization is enabled --
`ProteinQuantificationSearch` selects `:peak_area` as the quant column otherwise and never fills it. Emitting
it as an all-zero (or, after blanking, all-empty) column invites readers to treat uncomputed
values as measurements. Drop it from the output entirely when normalization is off.
"""
function drop_uncomputed_normalized(cols::NamedTuple, normalized::Bool)
    normalized && return cols
    haskey(cols, :peak_area_normalized) || return cols
    ks = filter(!=(:peak_area_normalized), keys(cols))
    return NamedTuple{ks}(map(k -> cols[k], ks))
end

"""Blank sentinel zero areas in a column table, returning a NamedTuple (Arrow path)."""
function blank_unquantified_areas(tbl)
    cols = Tables.columntable(tbl)
    ks = keys(cols)
    any(k -> k in OUTPUT_BLANKED_QUANT_COLUMNS, ks) || return cols
    return NamedTuple{ks}(map(k -> k in OUTPUT_BLANKED_QUANT_COLUMNS ?
                                   _blank_unquantified(cols[k]) : cols[k], ks))
end

function insert_at_indices(original::S, insertions::Vector{Tuple{String, UInt8}}) where {S<:AbstractString}
    # Convert the original string into an array of characters for easier manipulation
    char_array = collect(original)

    # Sort the insertions by index in ascending order
    sorted_insertions = sort(insertions, by = x -> x[2])

    # Adjust the index for each insertion
    offset = 0
    for (substr, idx) in sorted_insertions
        # Adjust the index with the current offset
        insertion_index = idx + offset
        # Insert each character of the substring at the specified index
        for (i, char) in enumerate(substr)
            insert!(char_array, insertion_index + i, char)
        end
        # Update the offset by the length of the inserted substring
        offset += length(substr)
    end

    # Join the array of characters back into a single string
    return join(char_array)
end


"""
    _map_accessions(accession_str::AbstractString, mapping::Dict{String,String})

Split a semicolon-delimited accession string and map each accession using
`mapping`. Missing accessions are ignored.
"""
function _map_accessions(accession_str::AbstractString, mapping::Dict{String,String})
    accs = split(accession_str, ';')
    names = String[]
    for acc in accs
        if acc != ""
            # Use mapped value if available and non-empty, otherwise fallback to accession
            mapped_value = get(mapping, acc, acc)
            # If mapped value is empty, use accession as fallback
            final_value = (mapped_value == "" || ismissing(mapped_value)) ? acc : mapped_value
            push!(names, final_value)
        end
    end
    return join(names, ';')
end
"""
    _map_accession_vector(column::AbstractVector{<:Union{Missing, AbstractString}}, mapping::Dict{String,String})

Return a vector where each element of `column` is mapped using `_map_accessions`.
Missing values are propagated. Results are cached for unique accession strings
to avoid redundant work.
"""
function _map_accession_vector(column::AbstractVector{T}, mapping::Dict{String,String}) where {T<:Union{Missing,AbstractString}}
    cache = Dict{String,String}()
    out = Vector{Union{Missing,String}}(undef, length(column))
    for i in eachindex(column)
        key = column[i]
        if ismissing(key)
            out[i] = missing
            continue
        end
        mapped = get!(cache, key) do
            _map_accessions(key, mapping)
        end
        out[i] = mapped
    end
    return out
end

function getModifiedSequence(
    sequence::S,
    isotope_mods::String,
    structural_mods::String,
    charge::UInt8) where {S<:AbstractString}

    mods = structural_mods*isotope_mods
    mods = [("("*getModName(mod.match)*")", getModIndex(mod.match)) for mod in parseMods(mods)]
    return "_"*insert_at_indices(sequence, mods)*"_."*string(charge)
end
function getModifiedSequence(
    sequence::S,
    isotope_mods::String,
    structural_mods::Missing,
    charge::UInt8) where {S<:AbstractString}

    mods = isotope_mods
    mods = [("("*getModName(mod.match)*")", getModIndex(mod.match)) for mod in parseMods(mods)]
    return "_"*insert_at_indices(sequence, mods)*"_."*string(charge)
end

function getModifiedSequence(
    sequence::S,
    isotope_mods::Missing,
    structural_mods::String,
    charge::UInt8) where {S<:AbstractString}

    mods = structural_mods
    mods = [("("*getModName(mod.match)*")", getModIndex(mod.match)) for mod in parseMods(mods)]
    return "_"*insert_at_indices(sequence, mods)*"_."*string(charge)
end
function getModifiedSequence(
    sequence::S,
    isotope_mods::Missing,
    structural_mods::Missing,
    charge::UInt8) where {S<:AbstractString}

    mods = ""
    mods = [("("*getModName(mod.match)*")", getModIndex(mod.match)) for mod in parseMods(mods)]
    return "_"*insert_at_indices(sequence, mods)*"_."*string(charge)
end

"""
    _n_distinct_precursors(peptide_lists) -> Int

Number of distinct precursor ids across the per-run `peptides` lists of one
protein group. Missing rows/entries are skipped.
"""
function _n_distinct_precursors(peptide_lists::AbstractVector)
    seen = Set{UInt32}()
    for pep in peptide_lists
        ismissing(pep) && continue
        for pid in pep
            ismissing(pid) || push!(seen, pid)
        end
    end
    return length(seen)
end

"""
    _extend_batch_to_group_end(keys, batch_end_idx, n_rows) -> Int

Walk `batch_end_idx` forward while `keys` keeps repeating, so a batch never splits one precursor
(or protein group) across two batches -- the wide format needs every row of a group together.

Takes the key column as a plain vector instead of indexing `df[i, :col]` inside the loop. That
form re-resolves the column by symbol on every iteration and infers `Any`, and this loop can run
once per row of the table.
"""
function _extend_batch_to_group_end(keys::AbstractVector, batch_end_idx::Int, n_rows::Int)
    @inbounds begin
        last_key = keys[batch_end_idx]
        while batch_end_idx < n_rows && isequal(keys[batch_end_idx + 1], last_key)
            batch_end_idx += 1
        end
    end
    return batch_end_idx
end

"""
    _modified_sequence_strings(peptide_lists, sequences, isotope_mods, structural_mods, charges)

Semicolon-joined modified sequences, one string per protein group.

Replaces a per-row loop that indexed `subdf[i, :peptides]` and assigned `subdf[i, :col] = ...`
(a symbol lookup per row, inferring `Any`), allocated a fresh `Vector{String}` per row, and then
allocated again in `filter(!isempty, ...)` before joining. Here the per-row scratch is one reused
buffer and the columns are typed vectors, so the only allocation left is the output string itself.
"""
function _modified_sequence_strings(
    peptide_lists::AbstractVector,
    sequences::AbstractVector,
    isotope_mods::AbstractVector,
    structural_mods::AbstractVector,
    charges::AbstractVector,
)
    out = Vector{String}(undef, length(peptide_lists))
    buf = String[]
    for i in eachindex(peptide_lists, out)
        empty!(buf)
        for pid in peptide_lists[i]
            ismissing(pid) && continue
            mod_seq = getModifiedSequence(
                sequences[pid], isotope_mods[pid], structural_mods[pid], charges[pid])
            isempty(mod_seq) || push!(buf, mod_seq)
        end
        out[i] = join(buf, ';')
    end
    return out
end

"""
    _ensure_typed_missing_file_columns!(df, file_names, T)

Ensure each file column in `file_names` exists in `df` and has element type
`Union{Missing,T}` when the column is fully missing for a batch.
"""
function _ensure_typed_missing_file_columns!(
    df::DataFrame,
    file_names::AbstractVector{<:AbstractString},
    ::Type{T}) where {T<:AbstractFloat}

    n = nrow(df)
    col_names = names(df)
    for fname in file_names
        if fname ∉ col_names || eltype(df[!, fname]) === Missing
            df[!, fname] = Vector{Union{Missing, T}}(missing, n)
        end
    end
    return df
end
# ---------------------------------------------------------------------------------------------
# DEAD CODE (commented out, not deleted).
#
# `writePrecursorCSV` is the unchunked precursor CSV writer. Nothing calls it: the live path is
# `writePrecursorCSV_chunked` below (ProteinQuantificationSearch.jl:228), and the only surviving mention of this
# name is the phrase "Chunked version of writePrecursorCSV" in that function's docstring.
#
# It also would not survive being called today -- it does the full-table materialisation the
# chunked version exists to avoid (`DataFrame(Arrow.Table(long_precursors_path))` over the whole
# precursors table) and it never received the batch-boundary and column-subset fixes applied to
# the chunked writer. If it is ever revived, port those first.
# ---------------------------------------------------------------------------------------------
# #Assume sorted by protein,peptide keys. Do this in batches and write a long and wide form .csv without
# #loading the entire table into memory. 
# function writePrecursorCSV(
#     long_precursors_path::String,
#     file_names::Vector{String},
#     normalized::Bool,
#     proteins::LibraryProteins;
#     output_schema_policy::OutputSchemaPolicy = OutputSchemaPolicy(),
#     write_csv::Bool = true,
#     batch_size::Int64 = 2000000)
#
#     function makeWideFormat(
#         longdf::DataFrame,
#         cols::AbstractVector{Symbol},
#         normalized::Bool)
#
#         value_col = normalized ? :peak_area_normalized : :peak_area
#
#         return unstack(longdf,
#                     cols,
#                     :file_name,
#                     value_col;
#                     combine = sum)      # or combine = maximum, first, etc.
#
#     end
#     # Replace empty strings with missing to avoid CSV.write indexing errors on empty Strings
#     function _sanitize_empty_strings!(df::DataFrame)
#         for nm in names(df)
#             col = df[!, nm]
#             if eltype(col) <: AbstractString
#                 if any(==(""), col)
#                     df[!, nm] = replace(col, "" => missing)
#                 end
#             elseif eltype(col) <: Union{Missing, AbstractString}
#                 has_empty = false
#                 @inbounds for v in col
#                     if !ismissing(v) && v == ""
#                         has_empty = true
#                         break
#                     end
#                 end
#                 if has_empty
#                     df[!, nm] = replace(col, "" => missing)
#                 end
#             end
#         end
#         return df
#     end
#     precursors_long = DataFrame(Arrow.Table(long_precursors_path))
#
#     unique_files_in_data = unique(precursors_long.file_name)
#
#     n_rows = size(precursors_long, 1)
#
#     out_dir, arrow_path = splitdir(long_precursors_path)
#     long_precursors_path = joinpath(out_dir,"precursors_long.tsv")
#     wide_precursors_path = joinpath(out_dir,"precursors_wide.tsv")
#     wide_precursors_arrow_path = joinpath(out_dir,"precursors_wide.arrow")
#     if isfile(wide_precursors_arrow_path)
#         safeRm(wide_precursors_arrow_path)
#     end
#     wide_columns = enabled_output_columns(output_schema_policy, :precursors, String[
#         "species",
#         "gene_names",
#         "inferred_protein_group",
#         "accession_numbers",
#         "sequence",
#         "charge",
#         "structural_mods",
#         "isotopic_mods",
#         "prec_mz",
#         "global_score",
#         "global_qval",
#         "use_for_protein_quant",
#         "precursor_idx",
#         "target",
#         "entrapment_group_id",
#     ])
#
#     long_columns_exclude = [:isotopes_captured, :scan_idx, :weight, :ms_file_idx]
#     select!(precursors_long, Not(long_columns_exclude))
#
#     accs = getAccession(proteins)
#     genes = getGeneName(proteins)
#     gene_map = Dict(accs[i] => genes[i] for i in eachindex(accs))
#     precursors_long[!, :gene_names] = _map_accession_vector(precursors_long.accession_numbers, gene_map)
#     # Build rename pairs dynamically to avoid conflicts
#     rename_pairs = Pair{Symbol,Symbol}[]
#     push!(rename_pairs, :new_best_scan => :apex_scan)
#     push!(rename_pairs, :prec_prob => :score)
#     push!(rename_pairs, :global_prob => :global_score)
#     push!(rename_pairs, :isotopes_captured_traces => :isotopes_captured)
#     push!(rename_pairs, :precursor_fraction_transmitted_traces => :precursor_fraction_transmitted)
#
#     # Apply all renames at once
#     if !isempty(rename_pairs)
#         rename!(precursors_long, rename_pairs)
#     end
#
#     requested_cols = enabled_output_columns(output_schema_policy, :precursors, Symbol[
#         :file_name,
#         :species,
#         :gene_names,
#         :inferred_protein_group,
#         :accession_numbers,
#         :sequence,
#         :charge,
#         :structural_mods,
#         :isotopic_mods,
#         :prec_mz,
#         :missed_cleavage,
#         :global_score,
#         :score,
#         :global_qval,
#         :qval,
#         :pep,
#         :peak_area,
#         :peak_area_normalized,
#         :points_integrated,
#         :precursor_fraction_transmitted,
#         :isotopes_captured,
#         :rt,
#         :apex_scan,
#         :global_pg_score,
#         :pg_score,
#         :use_for_protein_quant,
#         :precursor_idx,
#         :target,
#         :entrapment_group_id,
#     ])
#     available_cols = intersect(requested_cols, Symbol.(names(precursors_long)))
#     select!(precursors_long, available_cols)
#
#     sorted_columns = vcat(wide_columns, file_names)
#     open(long_precursors_path,"w") do io1
#         open(wide_precursors_path, "w") do io2
#             #Make file headers
#             if write_csv
#                 println(io1,join(names(precursors_long),"\t"))
#                 println(io2,join(sorted_columns,"\t"))
#             end
#             batch_start_idx, batch_end_idx = 1, min(batch_size+1,n_rows)
#             open(Arrow.Writer, wide_precursors_arrow_path; file=true) do arrow_writer
#                 while batch_start_idx <= n_rows
#                     #For the wide format, can't split a precursor between two batches.
#                     last_pid = precursors_long[batch_end_idx,:precursor_idx]
#                     while batch_end_idx < n_rows
#                         if precursors_long[batch_end_idx + 1,:precursor_idx] != last_pid
#                             break
#                         end
#                         batch_end_idx += 1
#                     end
#
#                     subdf =  precursors_long[range(batch_start_idx, batch_end_idx),:]
#                     batch_start_idx = batch_end_idx + 1
#                     batch_end_idx = min(batch_start_idx + batch_size, n_rows)
#                     if write_csv
#                         _sanitize_empty_strings!(subdf)
#                         CSV.write(io1, subdf, append=true, header=false, delim='\t')
#                     end
#                     subunstack = makeWideFormat(subdf, Symbol.(wide_columns), normalized)
#                     _ensure_typed_missing_file_columns!(subunstack, file_names, Float32)
#                     if write_csv
#                         _sanitize_empty_strings!(subunstack)
#                         CSV.write(io2, subunstack[!,sorted_columns], append=true,header=false,delim='\t')
#                     end
#                     # Normalize column types for consistent Arrow schema across batches
#                     allowmissing!(subunstack)
#                     Arrow.write(arrow_writer, subunstack[!,sorted_columns])
#                 end
#             end
#         end
#     end
#     if write_csv == false
#         safeRm(long_precursors_path; force = true)
#         safeRm(wide_precursors_path; force = true)
#     end
#     return wide_precursors_arrow_path
# end
"""
    _PRECURSOR_CSV_SOURCE_ALIASES

Output column name => the name it is read under in the merge chunk, for the columns the precursor
CSV writer renames. Needed so the writer can load only the columns it uses: `requested_cols` holds
post-rename names, while the chunk on disk still has the pre-rename ones.

Over-inclusion is safe here (a column that is read but unused costs only the read), whereas
under-inclusion would silently drop a column from the output, so the read set is the union of both
spellings intersected with what the chunk actually contains.
"""
const _PRECURSOR_CSV_SOURCE_ALIASES = Dict{Symbol, Symbol}(
    :global_score                   => :global_prob,
    :score                          => :prec_prob,
    :apex_scan                      => :new_best_scan,
    :isotopes_captured              => :isotopes_captured_traces,
    :precursor_fraction_transmitted => :precursor_fraction_transmitted_traces,
)

"""
    _precursor_csv_read_columns(tbl, requested_cols) -> Vector{Symbol}

The columns to pull out of a merge chunk: every requested output column, plus the pre-rename
spelling of each renamed one, plus the two the writer needs internally (`:precursor_idx` for batch
boundaries and `:accession_numbers` for the gene-name mapping), intersected with the chunk schema
and returned in the chunk's own column order.
"""
function _precursor_csv_read_columns(tbl, requested_cols)
    wanted = Set{Symbol}(requested_cols)
    for (out_name, src_name) in _PRECURSOR_CSV_SOURCE_ALIASES
        out_name in wanted && push!(wanted, src_name)
    end
    push!(wanted, :precursor_idx)
    push!(wanted, :accession_numbers)
    return Symbol[c for c in propertynames(tbl) if c in wanted]
end

"""
    writePrecursorCSV_chunked(chunk_refs, out_dir, file_names, normalized, proteins; ...)

Chunked version of writePrecursorCSV that processes merge chunks one at a time,
keeping memory bounded to ~1 chunk (~1 GB) instead of loading the full precursors table.
Produces identical output files: precursors_long.tsv, precursors_wide.tsv, precursors_wide.arrow.
"""
function writePrecursorCSV_chunked(
    chunk_refs::Vector{<:Any},
    out_dir::String,
    file_names::Vector{String},
    normalized::Bool,
    proteins::LibraryProteins;
    output_schema_policy::OutputSchemaPolicy = OutputSchemaPolicy(),
    write_csv::Bool = true,
    batch_size::Int64 = 2000000)

    function makeWideFormat(
        longdf::DataFrame,
        cols::AbstractVector{Symbol},
        normalized::Bool)

        value_col = normalized ? :peak_area_normalized : :peak_area

        return unstack(longdf,
                    cols,
                    :file_name,
                    value_col;
                    combine = sum)
    end

    function _sanitize_empty_strings!(df::DataFrame)
        for nm in names(df)
            col = df[!, nm]
            if eltype(col) <: AbstractString
                if any(==(""), col)
                    df[!, nm] = replace(col, "" => missing)
                end
            elseif eltype(col) <: Union{Missing, AbstractString}
                has_empty = false
                @inbounds for v in col
                    if !ismissing(v) && v == ""
                        has_empty = true
                        break
                    end
                end
                if has_empty
                    df[!, nm] = replace(col, "" => missing)
                end
            end
        end
        return df
    end

    # Build shared lookup maps once
    accs = getAccession(proteins)
    genes = getGeneName(proteins)
    gene_map = Dict(accs[i] => genes[i] for i in eachindex(accs))

    # Setup paths
    long_precursors_path = joinpath(out_dir, "precursors_long.tsv")
    wide_precursors_path = joinpath(out_dir, "precursors_wide.tsv")
    wide_precursors_arrow_path = joinpath(out_dir, "precursors_wide.arrow")
    isfile(wide_precursors_arrow_path) && safeRm(wide_precursors_arrow_path)

    wide_columns = enabled_output_columns(output_schema_policy, :precursors, String[
        "species",
        "gene_names",
        "inferred_protein_group",
        "accession_numbers",
        "sequence",
        "charge",
        "structural_mods",
        "isotopic_mods",
        "prec_mz",
        "peptide_start_positions",
        "missed_cleavage",
        "num_enzymatic_termini",
        "global_score",
        "global_qval",
        "use_for_protein_quant",
        "precursor_idx",
        "target",
        "entrapment_group_id",
    ])

    long_columns_exclude = [:isotopes_captured, :scan_idx, :weight, :ms_file_idx]

    requested_cols = enabled_output_columns(output_schema_policy, :precursors, Symbol[
        :file_name,
        :species,
        :gene_names,
        :inferred_protein_group,
        :accession_numbers,
        :sequence,
        :charge,
        :structural_mods,
        :isotopic_mods,
        :prec_mz,
        :peptide_start_positions,
        :missed_cleavage,
        :num_enzymatic_termini,
        :global_score,
        :score,
        :global_qval,
        :qval,
        :pep,
        :peak_area,
        :peak_area_normalized,
        :points_integrated,
        # True when QUANT_MIN_AREA_SURVIVING_RATIO withheld this row's area -- i.e. the
        # blank peak_area is a withheld measurement rather than a peak that was never found.
        # Recorded where the rule fires, so it classifies every blanked row exactly; the
        # previous Float32 peak_area_unsubtracted could not separate all cases.
        :quant_withheld,
        :precursor_fraction_transmitted,
        :isotopes_captured,
        :rt,
        :apex_scan,
        :global_pg_score,
        :pg_score,
        # The protein-level q-values. The table already carried the scores but not the q-values they
        # map to, so precursors_long.tsv gave no way to apply a protein FDR threshold.
        :global_pg_qval,
        :pg_qval,
        :use_for_protein_quant,
        :precursor_idx,
        :target,
        :entrapment_group_id,
    ])

    sorted_columns = vcat(wide_columns, file_names)
    n_chunks = length(chunk_refs)

    open(long_precursors_path, "w") do io1
        open(wide_precursors_path, "w") do io2
            open(Arrow.Writer, wide_precursors_arrow_path; file=true) do arrow_writer
                headers_written = false
                # Skip the progress bar when there's only one chunk —
                # ProgressBars displays "Inf:Inf, InfGs/it" on n=1 (rate divide-by-zero).
                pbar = n_chunks > 1 ? ProgressBar(total=n_chunks) : nothing
                pbar !== nothing && set_description(pbar, "Writing precursor CSV chunks:")

                for (ci, chunk_ref) in enumerate(chunk_refs)
                    # Load only the columns this writer uses. The chunk carries every column that
                    # survived to the final table; the TSV emits a curated subset, so materialising
                    # the whole chunk decoded columns that were dropped two lines later.
                    chunk_tbl = Arrow.Table(file_path(chunk_ref))
                    precursors_long = DataFrame()
                    for col in _precursor_csv_read_columns(chunk_tbl, requested_cols)
                        precursors_long[!, col] = chunk_tbl[col]
                    end
                    if :num_enzymatic_termini in requested_cols &&
                       !hasproperty(precursors_long, :num_enzymatic_termini)
                        # Resume/output compatibility for chunks written before
                        # semi-specific digestion existed. Those searches used
                        # fully specific libraries by construction.
                        precursors_long[!, :num_enzymatic_termini] = fill(
                            UInt8(2), length(chunk_tbl.precursor_idx)
                        )
                    end
                    n_rows = size(precursors_long, 1)
                    n_rows == 0 && continue

                    # Sentinel zero areas -> missing, so the TSV renders an empty cell
                    # rather than a value that looks like a measurement.
                    blank_unquantified_areas!(precursors_long)
                    if !normalized && hasproperty(precursors_long, :peak_area_normalized)
                        select!(precursors_long, Not(:peak_area_normalized))
                    end

                    # Drop excluded columns (only those that exist)
                    cols_to_drop = intersect(long_columns_exclude, Symbol.(names(precursors_long)))
                    if !isempty(cols_to_drop)
                        select!(precursors_long, Not(cols_to_drop))
                    end

                    # Add gene names from the protein table
                    precursors_long[!, :gene_names] = _map_accession_vector(precursors_long.accession_numbers, gene_map)

                    # Build rename pairs dynamically
                    rename_pairs = Pair{Symbol,Symbol}[]
                    col_names_set = Set(Symbol.(names(precursors_long)))
                    :new_best_scan ∈ col_names_set && push!(rename_pairs, :new_best_scan => :apex_scan)
                    :prec_prob ∈ col_names_set && push!(rename_pairs, :prec_prob => :score)
                    :global_prob ∈ col_names_set && push!(rename_pairs, :global_prob => :global_score)
                    :isotopes_captured_traces ∈ col_names_set && push!(rename_pairs, :isotopes_captured_traces => :isotopes_captured)
                    :precursor_fraction_transmitted_traces ∈ col_names_set && push!(rename_pairs, :precursor_fraction_transmitted_traces => :precursor_fraction_transmitted)
                    !isempty(rename_pairs) && rename!(precursors_long, rename_pairs)

                    # Reorder columns
                    available_cols = intersect(requested_cols, Symbol.(names(precursors_long)))
                    select!(precursors_long, available_cols)

                    # Write headers on first chunk
                    if !headers_written && write_csv
                        println(io1, join(names(precursors_long), "\t"))
                        println(io2, join(sorted_columns, "\t"))
                        headers_written = true
                    end

                    # Batch processing within chunk
                    pid_col = precursors_long[!, :precursor_idx]
                    batch_start_idx, batch_end_idx = 1, min(batch_size + 1, n_rows)
                    while batch_start_idx <= n_rows
                        batch_end_idx = _extend_batch_to_group_end(
                            pid_col, batch_end_idx, n_rows)

                        subdf = precursors_long[range(batch_start_idx, batch_end_idx), :]
                        batch_start_idx = batch_end_idx + 1
                        batch_end_idx = min(batch_start_idx + batch_size, n_rows)

                        if write_csv
                            _sanitize_empty_strings!(subdf)
                            CSV.write(io1, subdf, append=true, header=false, delim='\t')
                        end
                        subunstack = makeWideFormat(subdf, Symbol.(wide_columns), normalized)
                        _ensure_typed_missing_file_columns!(subunstack, file_names, Float32)
                        if write_csv
                            _sanitize_empty_strings!(subunstack)
                            CSV.write(io2, subunstack[!, sorted_columns], append=true, header=false, delim='\t')
                        end
                        # Normalize column types for consistent Arrow schema across batches
                        allowmissing!(subunstack)
                        Arrow.write(arrow_writer, subunstack[!, sorted_columns])
                    end

                    precursors_long = nothing  # Free chunk memory
                    pbar !== nothing && update(pbar)
                end
            end
        end
    end
    if !write_csv
        safeRm(long_precursors_path; force=true)
        safeRm(wide_precursors_path; force=true)
    end
    return wide_precursors_arrow_path
end

"""
    _protein_export_ranges(df, budget, row_limit; sequence_lengths=nothing)

Iterate row ranges sized for protein export buffers. One unusually large row
is indivisible; Arrow input record buffers are outside this workspace budget.
"""
struct ProteinExportRanges{T, L}
    data::T
    budget::Int
    row_limit::Int
    sequence_lengths::L
end

_protein_export_ranges(df, budget, row_limit; sequence_lengths=nothing) =
    ProteinExportRanges(df, budget, row_limit, sequence_lengths)

Base.IteratorSize(::Type{<:ProteinExportRanges}) = Base.SizeUnknown()
function Base.iterate(ranges::ProteinExportRanges, first_row::Int=1)
    df = ranges.data
    first_row > nrow(df) && return nothing
    bytes, i = 0, first_row
    while i <= nrow(df)
        row_bytes = 256 + 8length(df.peptides[i])
        if ranges.sequence_lengths !== nothing
            for id in df.peptides[i]
                ismissing(id) || (row_bytes += get(ranges.sequence_lengths, id, 0) + 1)
            end
        end
        if i > first_row && (bytes + row_bytes > ranges.budget || i - first_row >= ranges.row_limit)
            break
        end
        bytes += row_bytes
        i += 1
    end
    return first_row:i-1, i
end

function _spool_protein_export(long_path, dir, group_key, metadata_columns,
                               gene_map, protein_map, budget, row_limit, write_csv)
    catalog = DataFrame()
    seen = Set{Tuple}()
    precursor_ids = Set{UInt32}()
    current_key = nothing
    current_metadata = nothing
    current_path = ""
    writer = nothing
    column_types = Dict{Symbol, Type}(:gene_names => String, :protein_names => String)
    abundance_type = Float32
    function finish_group!()
        writer === nothing && return
        close(writer)
        push!(catalog, merge(current_metadata,
            (n_precursors_total=length(precursor_ids), export_path=current_path)); cols=:union)
    end
    try
        for batch in Arrow.Stream(long_path)
            df = DataFrame(batch; copycols=false)
            abundance_type = Base.nonmissingtype(eltype(df.abundance))
            for col in metadata_columns
                col in (:gene_names, :protein_names) && continue
                column_types[col] = Base.nonmissingtype(eltype(df[!, col]))
            end
            i = 1
            while i <= nrow(df)
                key = Tuple(df[i, col] for col in group_key)
                if !isequal(key, current_key)
                    finish_group!()
                    writer = nothing
                    key in seen && throw(ArgumentError("Protein export input must keep each protein group contiguous"))
                    push!(seen, key)
                    empty!(precursor_ids)
                    current_key = key
                    values = map(metadata_columns) do col
                        value = col == :gene_names ? _map_accessions(df.protein[i], gene_map) :
                                col == :protein_names ? _map_accessions(df.protein[i], protein_map) : df[i, col]
                        col in (:gene_names, :protein_names) && isempty(value) ? missing : value
                    end
                    current_metadata = NamedTuple{Tuple(metadata_columns)}(Tuple(values))
                    current_path = joinpath(dir, "group_$(length(seen)).arrow")
                    writer = open(Arrow.Writer, current_path; file=true, ntasks=0)
                end
                j = i
                while j <= nrow(df) && isequal(Tuple(df[j, col] for col in group_key), current_key)
                    for col in metadata_columns
                        col in (:gene_names, :protein_names) && continue
                        isequal(df[j, col], current_metadata[col]) || throw(ArgumentError(
                            "Inconsistent $col within protein group $(df.protein[j])"))
                    end
                    for id in df.peptides[j]
                        ismissing(id) || push!(precursor_ids, id)
                    end
                    j += 1
                end
                part = view(df, i:j-1, :)
                for rows in _protein_export_ranges(part, budget, row_limit)
                    Arrow.write(writer, view(part, rows, write_csv ? Colon() : [:file_name, :abundance]))
                end
                i = j
            end
        end
        finish_group!()
        writer = nothing
    finally
        writer === nothing || close(writer)
    end
    if !isempty(catalog)
        sort!(catalog, [order(:n_precursors_total, rev=true), order(:global_pg_score, rev=true), group_key...])
    end
    return catalog, column_types, abundance_type
end

"""
    writeProteinGroupsCSV(...; memory_budget_bytes=64*1024^2, batch_size=2000000)

Export protein-aligned Arrow batches through a temporary per-protein spool.
The group ordering index, per-protein precursor metadata, Arrow input records,
and bounded output buffers are retained. The original long Arrow table is not rewritten. TSV peptide strings
are constructed only when requested. Temporary files are removed on failure.
"""
function writeProteinGroupsCSV(
    long_pg_path::String,
    sequences::AbstractVector{<:AbstractString},
    isotope_mods::AbstractVector{Union{Missing, String}},
    structural_mods::AbstractVector{Union{String,Missing}},
    precursor_charge::AbstractVector{UInt8},
    file_names::Vector{String},
    proteins::LibraryProteins;
    output_schema_policy::OutputSchemaPolicy = OutputSchemaPolicy(),
    write_csv::Bool = true,
    batch_size::Int64 = 2000000,
    memory_budget_bytes::Int = 64 * 1024^2)

    batch_size > 0 || throw(ArgumentError("batch_size must be positive"))
    memory_budget_bytes > 0 || throw(ArgumentError("memory_budget_bytes must be positive"))
    length(unique(file_names)) == length(file_names) || throw(ArgumentError("Duplicate output file names"))
    out_dir = dirname(long_pg_path)
    isempty(out_dir) && (out_dir = ".")
    long_path = joinpath(out_dir, "protein_groups_long.tsv")
    wide_path = joinpath(out_dir, "protein_groups_wide.tsv")
    arrow_path = joinpath(out_dir, "protein_groups_wide.arrow")
    metadata_columns = enabled_output_columns(output_schema_policy, :protein_groups,
        Symbol[:species, :gene_names, :protein_names, :protein, :target, :entrap_id,
               :global_pg_score, :global_qval])
    group_key = enabled_output_columns(output_schema_policy, :protein_groups,
        Symbol[:species, :protein, :target, :entrap_id])
    long_columns = enabled_output_columns(output_schema_policy, :protein_groups, Symbol[
        :file_name, :species, :gene_names, :protein_names, :protein, :target, :entrap_id,
        :peptides, :n_precursors, :n_precursors_quantified, :n_modified_peptides, :n_peptides,
        :global_pg_score, :pg_score, :global_qval, :qval, :pg_pep, :total_peak_area, :abundance])
    file_indices = Dict(name => i for (i, name) in enumerate(file_names))
    accs, genes, names = getAccession(proteins), getGeneName(proteins), getProteinName(proteins)
    gene_map = Dict(accs[i] => genes[i] for i in eachindex(accs))
    protein_map = Dict(accs[i] => names[i] for i in eachindex(accs))
    budget = max(1, memory_budget_bytes ÷ 4)
    wide_capacity = max(1, min(batch_size, budget ÷ max(1, 16length(file_names) + 512)))

    mktempdir(out_dir; prefix=".protein_export_") do dir
        try
            catalog, column_types, abundance_type = _spool_protein_export(
                long_pg_path, dir, group_key, metadata_columns, gene_map, protein_map, budget, batch_size, write_csv)
            wide = DataFrame()
            for col in metadata_columns
                T = get(column_types, col, col in (:global_pg_score, :global_qval) ? Float32 :
                    col == :target ? Bool : col == :entrap_id ? UInt8 : String)
                wide[!, col] = Vector{Union{Missing, T}}()
            end
            for name in file_names
                Symbol(name) in metadata_columns && throw(ArgumentError("Run name conflicts with protein metadata: $name"))
                wide[!, name] = Vector{Union{Missing, abundance_type}}()
            end
            tmp_long, tmp_wide, tmp_arrow = joinpath.(dir, ("long.tsv", "wide.tsv", "wide.arrow"))
            long_io = write_csv ? open(tmp_long, "w") : nothing
            wide_io = write_csv ? open(tmp_wide, "w") : nothing
            try
                if write_csv
                    println(long_io, join(long_columns, '\t'))
                    println(wide_io, join(propertynames(wide), '\t'))
                end
                open(Arrow.Writer, tmp_arrow; file=true, ntasks=0) do arrow_writer
                    function flush_wide!()
                        write_csv && CSV.write(wide_io, wide; append=true, header=false, delim='\t')
                        Arrow.write(arrow_writer, wide)
                        wide = similar(wide, 0)
                    end
                    for group in eachrow(catalog)
                        areas = Vector{Union{Missing, abundance_type}}(missing, length(file_names))
                        occupied = falses(length(file_names))
                        sequence_cache = Dict{UInt32, String}()
                        sequence_lengths = Dict{UInt32, Int}()
                        for batch in Arrow.Stream(group.export_path)
                            df = DataFrame(batch; copycols=false)
                            for i in 1:nrow(df)
                                index = get(file_indices, df.file_name[i], 0)
                                index == 0 && throw(ArgumentError("Unknown run in protein export: $(df.file_name[i])"))
                                occupied[index] && throw(ArgumentError("Duplicate protein/run observation: $(group.protein), $(df.file_name[i])"))
                                occupied[index] = true
                                areas[index] = df.abundance[i]
                            end
                            if write_csv
                                for list in df.peptides, id in list
                                    ismissing(id) && continue
                                    if !haskey(sequence_cache, id)
                                        sequence_cache[id] = getModifiedSequence(
                                            sequences[id], isotope_mods[id], structural_mods[id], precursor_charge[id])
                                        sequence_lengths[id] = sizeof(sequence_cache[id])
                                    end
                                end
                                for rows in _protein_export_ranges(df, budget, batch_size; sequence_lengths)
                                    part = DataFrame(view(df, rows, :))
                                    part[!, :gene_names] = fill(group.gene_names, nrow(part))
                                    part[!, :protein_names] = fill(group.protein_names, nrow(part))
                                    part[!, :peptides] = map(part.peptides) do ids
                                        value = join((sequence_cache[id] for id in ids if !ismissing(id) && !isempty(sequence_cache[id])), ';')
                                        isempty(value) ? missing : value
                                    end
                                    CSV.write(long_io, select(part, long_columns); append=true, header=false, delim='\t')
                                end
                            end
                        end
                        for col in metadata_columns
                            push!(wide[!, col], group[col])
                        end
                        for (i, name) in enumerate(file_names)
                            push!(wide[!, name], areas[i])
                        end
                        nrow(wide) >= wide_capacity && flush_wide!()
                    end
                    (!isempty(wide) || isempty(catalog)) && flush_wide!()
                end
            finally
                long_io === nothing || close(long_io)
                wide_io === nothing || close(wide_io)
            end
            mv(tmp_arrow, arrow_path; force=true)
            if write_csv
                mv(tmp_long, long_path; force=true)
                mv(tmp_wide, wide_path; force=true)
            else
                safeRm(long_path; force=true)
                safeRm(wide_path; force=true)
            end
        finally
            # Release temporary Arrow mappings before directory cleanup on Windows.
            Sys.iswindows() && GC.gc()
        end
    end
    return arrow_path
end
