# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

"""Build score calibration from Arrow record batches using bounded workspace."""
function _score_calibration_from_arrow(path, score_col; kwargs...)
    build_score_calibration(; kwargs...) do emit
        # Scores are emitted one by one, so the file is unmapped afterwards.
        with_arrow_stream(path) do stream
            for table in stream
                _emit_score_arrays(emit, Tables.getcolumn(table, score_col), Tables.getcolumn(table, :target))
            end
        end
    end
end

"""Return a callable grouped PAVA fit, or `nothing` for empty input."""
function get_pep_interpolation(path::String, score_col::Symbol;
    fdr_scale_factor::Float32=1.0f0, memory_budget_bytes::Int=SCORE_WORKSPACE_BYTES)
    fit = _score_calibration_from_arrow(path, score_col; fdr_scale_factor, memory_budget_bytes)
    return fit === nothing ? nothing : fit.pep_interp
end
