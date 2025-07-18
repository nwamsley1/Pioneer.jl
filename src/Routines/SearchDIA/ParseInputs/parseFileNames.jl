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

function parseFileNames(
    ms_table_paths::Vector{String})
    file_names = first.(
                    split.(
                        basename.(ms_table_paths), '.'
                        ))
    split_file_names = split.(file_names, "_")

            #If file names have inconsistnat number of delimiters, give up on parsing and use the entire file name
            unique_split_file_name_lengths = unique(length.(split_file_names))
            parsed_file_names = copy(split_file_names)
            M = length(parsed_file_names)
    
            if length(unique_split_file_name_lengths) == 1
                N = first(unique_split_file_name_lengths)
                M = length(file_names)
                split_strings = Array{String}(undef, (M, N))
                for i in 1:N
                    for j in 1:M
                        split_strings[j, i] = split_file_names[j][i]
                    end
                end
                cols_to_keep = zeros(Bool, N)
                for (col_idx, col) in enumerate(eachcol(split_strings))
                    if length(unique(col))>1
                        cols_to_keep[col_idx] = true
                    end
                end
                parsed_file_names = Vector{String}(undef, M)
                split_strings = split_strings[:, cols_to_keep]
                for i in 1:M
                    parsed_file_names[i] = join(split_strings[i,:], "_")
                end
            else
                parsed_file_names = file_names
            end
    
            for (i, fname) in enumerate(parsed_file_names)
                if fname == ""
                    parsed_file_names[i] = string(i)
                end
            end
            
            return parsed_file_names#, parsed_fnames, file_path_to_parsed_name
end