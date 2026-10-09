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
    clamp_digest_length_to_model(model_name, min_length, max_length; rt_model = nothing) -> (min, max)

Reconcile the user's `fasta_digest_params` peptide-length bounds with what the
prediction models actually accept, returning the bounds the digest should use.

Both models a build sends every peptide to are checked: the fragment model
(`model_name`, from `MODEL_CONFIGS`) and, when given, the retention-time model
(`rt_model`, from `RT_MODEL_CONFIGS`). A model declares its accepted range as
`peptide_length = (min = m, max = M)`; a model that omits the field, sets it to
`nothing`, or is not listed is treated as unconstrained. The digest is narrowed
to the overlap of the declared ranges.

The problem this exists to prevent is silent loss. `digest_fasta` filters to the
requested `[min_length, max_length]` window, so nothing is dropped inside
Pioneer -- but a peptide that satisfies the user's bounds and *exceeds a
model's* is handed to Koina, which rejects or truncates it. The peptide then
disappears from the library with no entry in any log. Prosit's tokenizer, for
instance, hard-caps at 30 residues, so a build with `max_length = 40` loses
every 31-40mer without saying so.

Each override is reported with `@user_warn`, naming the model that sets the
bound and both values, so the narrowing appears in the build log rather than
being inferred later from a short library.

Returns a `Tuple{Int,Int}`. Throws if the models' range and the user's request
do not overlap at all, since that would otherwise produce an empty library.
"""
function clamp_digest_length_to_model(model_name::AbstractString,
                                      min_length::Integer,
                                      max_length::Integer;
                                      rt_model::Union{Nothing, AbstractString} = nothing)::Tuple{Int,Int}
    user_min, user_max = Int(min_length), Int(max_length)

    # (description, limits) for every model with a declared range. `get` with a
    # default rather than `cfg.peptide_length`: a model added without the field
    # stays unconstrained instead of erroring at build time.
    declared = Tuple{String, @NamedTuple{min::Int, max::Int}}[]
    for (what, name, configs) in (("prediction model", model_name, MODEL_CONFIGS),
                                  ("retention-time model", rt_model, RT_MODEL_CONFIGS))
        name === nothing && continue
        cfg = get(configs, String(name), nothing)
        cfg === nothing && continue
        limits = get(cfg, :peptide_length, nothing)
        limits === nothing || push!(declared, ("$what '$(name)'", limits))
    end
    isempty(declared) && return (user_min, user_max)

    # The binding bound at each end, and the model that sets it.
    lo_by, lo = declared[argmax([l.min for (_, l) in declared])]
    hi_by, hi = declared[argmin([l.max for (_, l) in declared])]
    model_min, model_max = lo.min, hi.max

    if user_min > model_max || user_max < model_min || model_min > model_max
        supported = join(["$(d) supports $(l.min)-$(l.max)" for (d, l) in declared], "; ")
        error("fasta_digest_params requests peptide lengths $(user_min)-$(user_max), " *
              "but $(supported). No peptide can satisfy all of them; " *
              "widen the digest range or choose different models.")
    end

    new_min, new_max = user_min, user_max

    if user_min < model_min
        new_min = model_min
        @user_warn "$(lo_by) does not support peptides shorter " *
                   "than $(model_min) residues; raising fasta_digest_params.min_length " *
                   "from $(user_min) to $(model_min)."
    end

    if user_max > model_max
        new_max = model_max
        @user_warn "$(hi_by) does not support peptides longer " *
                   "than $(model_max) residues; lowering fasta_digest_params.max_length " *
                   "from $(user_max) to $(model_max). Peptides longer than " *
                   "$(model_max) would otherwise be dropped silently."
    end

    return (new_min, new_max)
end
