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
Pipeline API for composable file operations.

Provides a framework for composing multiple file operations
that execute in a single pass for optimal performance and
memory efficiency.
"""

using DataFrames

#==========================================================
Pipeline Framework Types
==========================================================#

"""
    TransformPipeline

A composable pipeline of operations to apply to a DataFrame in a single pass.
Operations are executed in order, with optional post-actions for metadata updates.
"""
struct TransformPipeline
    operations::Vector{<:Pair{String}}  # description => operation
    post_actions::Vector{Function}  # Actions to run after transform (e.g., mark_sorted!)
end

# Constructor for empty pipeline
TransformPipeline() = TransformPipeline(Pair{String, Function}[], Function[])

"""
    PipelineOperation

Wrapper for operations that need post-processing actions (e.g., updating sort state).
"""
struct PipelineOperation
    operation::Pair{String, Function}
    post_action::Function
end

# Import Base operator to extend it
import Base: |>

# Override |> for pipeline composition
|>(pipeline::TransformPipeline, op::Pair{String, F}) where {F<:Function} = 
    TransformPipeline(
        vcat(pipeline.operations, op),
        pipeline.post_actions
    )

# Handle PipelineOperation specially
|>(pipeline::TransformPipeline, op::PipelineOperation) = 
    TransformPipeline(
        vcat(pipeline.operations, op.operation),
        vcat(pipeline.post_actions, op.post_action)
    )

