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
    clamp_digest_length_to_model(model_name, min_length, max_length) -> (min, max)

Reconcile the user's `fasta_digest_params` peptide-length bounds with what the
fragment-prediction model actually accepts, returning the bounds the digest
should use.

A model declares its accepted range as `peptide_length = (min = m, max = M)` in
its `MODEL_CONFIGS` entry. Models that omit the field, or set it to `nothing`,
are treated as unconstrained and the user's bounds pass through untouched.

The problem this exists to prevent is silent loss. `digest_fasta` filters to the
requested `[min_length, max_length]` window, so nothing is dropped inside
Pioneer -- but a peptide that satisfies the user's bounds and *exceeds the
model's* is handed to Koina, which rejects or truncates it. The peptide then
disappears from the library with no entry in any log. Prosit's tokenizer, for
instance, hard-caps at 30 residues, so a build with `max_length = 40` loses
every 31-40mer without saying so.

Each override is reported with `@user_warn`, naming the model and both values,
so the narrowing appears in the build log rather than being inferred later from
a short library.

Returns a `Tuple{Int,Int}`. Throws if the model's range and the user's request
do not overlap at all, since that would otherwise produce an empty library.
"""
function clamp_digest_length_to_model(model_name::AbstractString,
                                      min_length::Integer,
                                      max_length::Integer)::Tuple{Int,Int}
    user_min, user_max = Int(min_length), Int(max_length)

    cfg = get(MODEL_CONFIGS, String(model_name), nothing)
    cfg === nothing && return (user_min, user_max)

    # `get` with a default rather than `cfg.peptide_length`: a model added
    # without the field stays unconstrained instead of erroring at build time.
    limits = get(cfg, :peptide_length, nothing)
    limits === nothing && return (user_min, user_max)

    if user_min > limits.max || user_max < limits.min
        error("fasta_digest_params requests peptide lengths $(user_min)-$(user_max), " *
              "but prediction model '$(model_name)' only supports " *
              "$(limits.min)-$(limits.max). No peptide can satisfy both; " *
              "widen the digest range or choose a different prediction_model.")
    end

    new_min, new_max = user_min, user_max

    if user_min < limits.min
        new_min = limits.min
        @user_warn "prediction model '$(model_name)' does not support peptides shorter " *
                   "than $(limits.min) residues; raising fasta_digest_params.min_length " *
                   "from $(user_min) to $(limits.min)."
    end

    if user_max > limits.max
        new_max = limits.max
        @user_warn "prediction model '$(model_name)' does not support peptides longer " *
                   "than $(limits.max) residues; lowering fasta_digest_params.max_length " *
                   "from $(user_max) to $(limits.max). Peptides longer than " *
                   "$(limits.max) would otherwise be dropped silently."
    end

    return (new_min, new_max)
end
