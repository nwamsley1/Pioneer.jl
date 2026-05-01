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

#=
Training data selection strategies.
Uses AbstractPSMContainer for data access.
=#

"""
    select_training_data(psms::AbstractPSMContainer, strategy::TrainingDataStrategy, iteration::Int) -> AbstractPSMContainer

Select training data based on the strategy and current iteration.
Returns a copy of the training data.
"""
function select_training_data(psms::AbstractPSMContainer, ::AllDataSelection, ::Int)
    return copy_container(psms)
end

# DataFrame convenience wrapper
function select_training_data(psms::DataFrame, strategy::TrainingDataStrategy, iteration::Int)
    container = DataFramePSMContainer(psms, Val(:unsafe))
    result = select_training_data(container, strategy, iteration)
    return to_dataframe(result)
end
