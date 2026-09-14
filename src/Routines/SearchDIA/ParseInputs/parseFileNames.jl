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
