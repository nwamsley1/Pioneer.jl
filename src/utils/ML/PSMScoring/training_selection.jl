# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

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
