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
    MassErrSample

16-byte sample emitted by `run_fused_masserr!` (ParameterTuning's per-thread
mass-error collector) and consumed by `fit_mass_err_model`,
`fit_intensity_mass_error_model`, `fit_scout_calibrated_model`,
`extract_fragment_plot_data`, and `generate_wide_scout_plot`.

`observed_mz` is the **raw** (uncorrected) peak m/z — bias is learned downstream.
"""
struct MassErrSample
    theoretical_mz::Float32   # 4B — fragment library m/z
    observed_mz::Float32      # 4B — raw peak m/z
    intensity::Float32        # 4B — peak intensity
    rt::Float32               # 4B — scan retention time
end

MassErrSample() = MassErrSample(0f0, 0f0, 0f0, 0f0)
