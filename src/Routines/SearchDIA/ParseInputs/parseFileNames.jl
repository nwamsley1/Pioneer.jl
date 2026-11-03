# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

"""
    parseFileNames(ms_table_paths::Vector{String})

Run names used for the `file_name` column and the wide-table headers: the
full file name with the extension removed. Names are never shortened, so
`20240101_lab_sampleA.arrow` stays `20240101_lab_sampleA`. The `.zt` that
convertSciex adds to ZT Scan DIA outputs is part of the extension:
`run1.zt.scxs` is run `run1`.
"""
function parseFileNames(
    ms_table_paths::Vector{String})
    file_names = [String(chopsuffix(first(splitext(basename(p))), ".zt")) for p in ms_table_paths]
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
