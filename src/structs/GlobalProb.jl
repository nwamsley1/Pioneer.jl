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
    GlobalPrecFeatures{N}

Compact per-precursor feature summary for global probability model training.
Aggregates scores across runs into a fixed-size feature vector for ML-based
global precursor probability estimation.

# Fields
- `top_scores::NTuple{N, Float32}`: Top N :prec_prob values across runs (descending), padded with -1.0f0
- `n_scores_available::Int32`: Actual number of runs with scores (≤ N)
- `stats::NTuple{5, Float32}`: (mean, max, min, std, skewness); -1.0f0 for undefined
- `counts::NTuple{3, Int32}`: (n_runs_with_score, n_runs_total, n_above_thresh)
- `deltas::NTuple{2, Float32}`: (top1 - top2, top2 - top3); -1.0f0 if undefined
- `logodds_baseline::Float32`: Baseline probability from logodds(p, sqrt_n_runs)

# Notes
- Sentinel value -1.0f0 indicates missing/undefined features (not 0.0 or NaN)
- Statistics computed using streaming algorithms (no full array materialization)
- `logodds_baseline` serves dual purpose: feature input and fallback score
"""
struct GlobalPrecFeatures{N}
    top_scores::NTuple{N, Float32}
    n_scores_available::Int32
    stats::NTuple{5, Float32}
    counts::NTuple{3, Int32}
    deltas::NTuple{2, Float32}
    logodds_baseline::Float32
end

# Convenience constructor with validation
function GlobalPrecFeatures{N}(
    top_scores::NTuple{N, Float32},
    n_scores_available::Int32,
    meanp::Float32, maxp::Float32, minp::Float32, stdp::Float32, skewp::Float32,
    n_runs_with_score::Int32, n_runs_total::Int32, n_above_thresh::Int32,
    delta_12::Float32, delta_23::Float32,
    logodds_baseline::Float32
) where {N}
    GlobalPrecFeatures{N}(
        top_scores,
        n_scores_available,
        (meanp, maxp, minp, stdp, skewp),
        (n_runs_with_score, n_runs_total, n_above_thresh),
        (delta_12, delta_23),
        logodds_baseline
    )
end
