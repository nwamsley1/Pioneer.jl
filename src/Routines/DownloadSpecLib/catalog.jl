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

# Turning a repository listing into the catalog the CLI prints and the GUI
# renders.
#
# Two metadata sources, in priority order:
#   1. `libraries.json` at the repo root -- curation that a library cannot know
#      about itself: a human title, the prediction model, why you would pick it.
#   2. the library's own `config.json` -- everything derivable from the build.
#
# The manifest is optional and consulted first, so the usual case costs one
# request. A library missing from it still appears, described from its config,
# which is why uploading a library is a one-step publish.

"""
A library as presented to the user. `details` is display-ready
`label => value` pairs, already formatted.
"""
struct LibraryEntry
    name::String
    title::String
    model::String
    description::String
    recommended_for::String
    total_bytes::Int
    n_files::Int
    details::Vector{Pair{String, String}}
end

_str(x) = x === nothing ? "" : string(x)

"""
    _fmt_mods(mods) -> String

Format a `{"name": [...], "pattern": [...]}` modification group as
`"Unimod:4 (C), Unimod:35 (M)"`. Pioneer writes these three parallel arrays for
both fixed and variable mods.
"""
function _fmt_mods(mods)
    mods isa AbstractDict || return ""
    names = get(mods, "name", nothing)
    names isa AbstractVector && !isempty(names) || return "none"
    patterns = get(mods, "pattern", nothing)
    parts = String[]
    for (i, name) in enumerate(names)
        pattern = patterns isa AbstractVector && i <= length(patterns) ? _str(patterns[i]) : ""
        push!(parts, isempty(pattern) ? _str(name) : string(_str(name), " (", pattern, ")"))
    end
    return join(parts, ", ")
end

_fmt_range(lo, hi) =
    lo === nothing || hi === nothing ? "" : string(_str(lo), "-", _str(hi))

"""
    describe_config(config) -> Vector{Pair{String,String}}

Derive display fields from a library's `config.json`. Missing keys are skipped
rather than rendered blank, so an older library described by a thinner config
still produces a sensible table.

Pure, and the bulk of what the catalog shows.
"""
function describe_config(config::AbstractDict)
    details = Pair{String, String}[]
    push_if!(label, value) = isempty(_str(value)) || push!(details, label => _str(value))

    fastas = get(config, "fasta_names", nothing)
    fastas isa AbstractVector && push_if!("FASTA", join(string.(fastas), ", "))

    digest = get(config, "fasta_digest_params", nothing)
    if digest isa AbstractDict
        push_if!("Cleavage", get(digest, "cleavage_regex", nothing))
        push_if!("Missed cleavages", get(digest, "missed_cleavages", nothing))
        push_if!("Peptide length",
                 _fmt_range(get(digest, "min_length", nothing), get(digest, "max_length", nothing)))
        push_if!("Charge",
                 _fmt_range(get(digest, "min_charge", nothing), get(digest, "max_charge", nothing)))
        push_if!("Decoys", get(digest, "decoy_method", nothing))
        entrapment = get(digest, "entrapment_r", nothing)
        entrapment isa Number && entrapment > 0 && push_if!("Entrapment", entrapment)
    end

    push_if!("Fixed mods", _fmt_mods(get(config, "fixed_mods", nothing)))
    push_if!("Variable mods", _fmt_mods(get(config, "variable_mods", nothing)))

    nce = get(config, "nce_params", nothing)
    nce isa AbstractDict && push_if!("NCE", get(nce, "nce", nothing))

    library_params = get(config, "library_params", nothing)
    if library_params isa AbstractDict
        push_if!("Precursor m/z",
                 _fmt_range(get(library_params, "prec_mz_min", nothing),
                            get(library_params, "prec_mz_max", nothing)))
        push_if!("Fragment m/z",
                 _fmt_range(get(library_params, "frag_mz_min", nothing),
                            get(library_params, "frag_mz_max", nothing)))
    end

    return details
end

"""
    _manifest_details(entry) -> Union{Vector{Pair{String,String}}, Nothing}

Display fields carried by the manifest itself. When present these are used as-is
and no `config.json` is fetched, which is what keeps the common path to a single
request. `nothing` means "not curated -- go read the config".
"""
function _manifest_details(entry::AbstractDict)
    details = get(entry, "details", nothing)
    details isa AbstractDict || return nothing
    return Pair{String, String}[String(k) => _str(v) for (k, v) in details]
end

"""
    build_catalog(libraries, manifest; config_for) -> Vector{LibraryEntry}

Merge a tree listing with optional curation. `config_for(library)` is called
only for libraries the manifest does not already describe, and may return
`nothing` when the config cannot be read -- such a library still lists, just
without details.
"""
function build_catalog(libraries::Vector{HFLibrary},
                       manifest::Union{AbstractDict, Nothing};
                       config_for = _ -> nothing)
    curated = Dict{String, Any}()
    if manifest isa AbstractDict
        entries = get(manifest, "libraries", nothing)
        entries isa AbstractDict && (curated = Dict{String, Any}(entries))
    end

    catalog = LibraryEntry[]
    for library in libraries
        entry = get(curated, library.name, nothing)
        entry isa AbstractDict || (entry = Dict{String, Any}())

        details = _manifest_details(entry)
        if details === nothing
            config = config_for(library)
            details = config isa AbstractDict ? describe_config(config) : Pair{String, String}[]
        end

        push!(catalog, LibraryEntry(
            library.name,
            _str(get(entry, "title", library.name)),
            _str(get(entry, "model", "")),
            _str(get(entry, "description", "")),
            _str(get(entry, "recommended_for", "")),
            library.total_bytes,
            length(library.files),
            details,
        ))
    end
    return catalog
end

"""
    format_bytes(n) -> String

Human byte count, binary units, matching what the GUI shows next to a download.
"""
function format_bytes(n::Integer)
    n < 1024 && return string(n, " B")
    units = ("KiB", "MiB", "GiB", "TiB")
    value = float(n)
    for unit in units
        value /= 1024
        if value < 1024 || unit == last(units)
            return string(round(value; digits = 2), " ", unit)
        end
    end
    return string(n, " B")  # unreachable; keeps the return type concrete
end

"""
    catalog_to_json(catalog) -> String

The `--json` payload. The GUI parses exactly this, so the shape is the contract
between the two halves of the feature.
"""
function catalog_to_json(catalog::Vector{LibraryEntry})
    payload = [Dict{String, Any}(
        "name" => e.name,
        "title" => e.title,
        "model" => e.model,
        "description" => e.description,
        "recommended_for" => e.recommended_for,
        "total_bytes" => e.total_bytes,
        "size_human" => format_bytes(e.total_bytes),
        "n_files" => e.n_files,
        "details" => Dict{String, Any}(k => v for (k, v) in e.details),
    ) for e in catalog]
    return JSON.json(Dict{String, Any}("libraries" => payload))
end
