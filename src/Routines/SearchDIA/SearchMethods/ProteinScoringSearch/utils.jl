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

include(joinpath(@__DIR__, "feature_selection.jl"))
include(joinpath(@__DIR__, "protein_inference_pipeline.jl"))
include(joinpath(@__DIR__, "qc_plots.jl"))
include(joinpath(@__DIR__, "protein_groups.jl"))
include(joinpath(@__DIR__, "global_protein_scoring.jl"))
include(joinpath(@__DIR__, "qvalue_io.jl"))
include(joinpath(@__DIR__, "model_fit.jl"))
include(joinpath(@__DIR__, "run_protein_scoring.jl"))
