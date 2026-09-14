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
    parseFileNames(ms_table_paths::Vector{String})

Run names used for the `file_name` column and the wide-table headers: the
full file name with the extension removed. Names are never shortened, so
`20240101_lab_sampleA.arrow` stays `20240101_lab_sampleA`.
"""
function parseFileNames(
    ms_table_paths::Vector{String})
    file_names = first.(splitext.(basename.(ms_table_paths)))
    for (i, fname) in enumerate(file_names)
        if fname == ""
            file_names[i] = string(i)
        end
    end
    return file_names
end

"""
    distinguishingFileNames(file_names::Vector{String})

Short labels for plots: the `_`-delimited tokens that differ between the
runs, with every token shared by all of them dropped, so
`20240101_lab_rep1` / `20240101_lab_rep2` label as `rep1` / `rep2`. Falls
back to the full names when the runs do not split into the same number of
tokens, and to the full name when nothing distinguishes a run. Labels only:
the `file_name` column and the wide-table headers keep the full name
(`parseFileNames`).
"""
function distinguishingFileNames(file_names::Vector{String})
    split_names = split.(file_names, "_")
    length(unique(length.(split_names))) == 1 || return copy(file_names)
    n_tokens = first(length.(split_names))
    keep = [length(unique(s[i] for s in split_names)) > 1 for i in 1:n_tokens]
    labels = [join((s[i] for i in 1:n_tokens if keep[i]), "_") for s in split_names]
    for (i, label) in enumerate(labels)
        isempty(label) && (labels[i] = file_names[i])
    end
    return labels
end
