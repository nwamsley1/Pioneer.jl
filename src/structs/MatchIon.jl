# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

"""
    MassErrSample

20-byte sample emitted by `run_fused_masserr!` (ParameterTuning's per-thread
mass-error collector) and consumed by `fit_mass_err_model`,
`fit_intensity_mass_error_model`, `fit_scout_calibrated_model`,
and `extract_fragment_plot_data`.

`observed_mz` is the **raw** (uncorrected) peak m/z — bias is learned downstream.
"""
struct MassErrSample
    theoretical_mz::Float32   # 4B — fragment library m/z
    observed_mz::Float32      # 4B — raw peak m/z
    intensity::Float32        # 4B — peak intensity
    rt::Float32               # 4B — scan retention time
    precursor_idx::UInt32     # peptide grouping for independent calibration QC
end

MassErrSample(mz::Float32, observed::Float32, intensity::Float32, rt::Float32) =
    MassErrSample(mz, observed, intensity, rt, UInt32(0))
MassErrSample() = MassErrSample(0f0, 0f0, 0f0, 0f0)
