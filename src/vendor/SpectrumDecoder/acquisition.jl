# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# The acquisition mode an acquisition's embedded method states, as the label Pioneer reads from `.arrow` metadata.

"Version of the acquisition-mode metadata keys written by `write_pioneer_arrow`."
const ACQUISITION_LABEL_VERSION = "1"

const _MS1_ONLY_KINDS = ("MasterScan", "tSIM")

"""
    acquisition_mode(f) -> Union{Nothing, NamedTuple{(:mode, :evidence), Tuple{String, String}}}

The acquisition mode stated by the embedded method (`acquisition_mode_from_method` on `method_text(f)`), or `nothing`
when the method states none or the file has no method text.
"""
acquisition_mode(f::Acquisition) = acquisition_mode_from_method(join(method_text(f), "\n"))

"""
    acquisition_mode_from_method(text) -> Union{Nothing, NamedTuple{(:mode, :evidence)}}

Acquisition mode from a method summary. Three layouts are recognised:

- one `Scan <kind>` line per experiment: `ddMSnScan` → `dda` (`targeted` when the method also holds a
  `TargetedInclusionMassList` filter), `DIAScan` or `tMSn` → `dia`; a method whose only kinds are `MasterScan` and
  `tSIM` requests no MS2 scans → `not_searchable`;
- an `Experiment` / `Experiments` heading followed by experiment-type lines: `dd-MS²…` → `dda`, `DIA`, `vDIA`, `PRM` or
  `… / DIA` → `dia`;
- a `Data Dependent Settings:` section → `dda`.

Returns `nothing` when none applies (including unrecognised `Scan <kind>` values).
"""
function acquisition_mode_from_method(text::AbstractString)
    isempty(text) && return nothing
    lines = [strip(l) for l in split(text, r"\r?\n")]
    kinds = String[]
    for l in lines
        m = match(r"^Scan ([A-Za-z]\w*)$", l)
        m === nothing || m[1] in kinds || push!(kinds, m[1])
    end
    found(s) = (mode = s[1], evidence = "method: " * s[2])
    if "ddMSnScan" in kinds
        occursin("TargetedInclusionMassList", text) &&
            return found(("targeted", "Scan ddMSnScan with a targeted inclusion list"))
        return found(("dda", "Scan ddMSnScan"))
    end
    "DIAScan" in kinds && return found(("dia", "Scan DIAScan"))
    "tMSn" in kinds && return found(("dia", "Scan tMSn"))
    h = findfirst(l -> l == "Experiment" || l == "Experiments", lines)
    if h !== nothing
        after = lines[h+1:end]
        i = findfirst(l -> occursin(r"^(.*/\s*)?dd-MS"i, l), after)
        i === nothing || return found(("dda", "experiment \"$(after[i])\""))
        i = findfirst(l -> l in ("DIA", "vDIA", "PRM") || occursin(r"/\s*DIA$"i, l), after)
        i === nothing || return found(("dia", "experiment \"$(after[i])\""))
    end
    occursin("Data Dependent Settings:", text) && return found(("dda", "Data Dependent Settings"))
    !isempty(kinds) && all(in(_MS1_ONLY_KINDS), kinds) &&
        return found(("not_searchable", "no MS2 experiment (Scan $(join(kinds, ", ")))"))
    return nothing
end

"""
    acquisition_metadata(f) -> Union{Nothing, Vector{Pair{String, String}}}

`.arrow` schema metadata for the stated acquisition mode: `acquisition_mode`, `acquisition_mode_source` = "vendor",
`acquisition_mode_evidence` and `acquisition_mode_version`; `nothing` when the method states no mode.
"""
function acquisition_metadata(f::Acquisition)
    a = acquisition_mode(f)
    a === nothing && return nothing
    ["acquisition_mode" => a.mode, "acquisition_mode_source" => "vendor",
     "acquisition_mode_evidence" => a.evidence, "acquisition_mode_version" => ACQUISITION_LABEL_VERSION]
end
