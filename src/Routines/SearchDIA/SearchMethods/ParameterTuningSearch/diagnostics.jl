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

# Type definitions moved to types.jl

"""
    record_tuning_status!(diagnostics, status)

Record the parameter tuning status for a file.
"""
function record_tuning_status!(diagnostics::ParameterTuningDiagnostics, status::ParameterTuningStatus)
    diagnostics.file_statuses[status.file_idx] = status
    
    if status.converged && !status.used_fallback
        diagnostics.n_successful += 1
    elseif status.used_fallback
        diagnostics.n_fallback += 1
    else
        diagnostics.n_failed += 1
    end
end

"""
    generate_parameter_tuning_report(diagnostics, output_path)

Generate a summary report of parameter tuning outcomes.
"""
function generate_parameter_tuning_report(diagnostics::ParameterTuningDiagnostics, output_path::String)
    open(output_path, "w") do io
        println(io, "# Parameter Tuning Summary Report")
        println(io)
        println(io, "## Overall Statistics")
        println(io, "- Total files: $(length(diagnostics.file_statuses))")
        println(io, "- Successful: $(diagnostics.n_successful)")
        println(io, "- Used fallback: $(diagnostics.n_fallback)")
        println(io, "- Failed: $(diagnostics.n_failed)")
        println(io)
        
        if diagnostics.n_fallback > 0
            println(io, "## Files Using Fallback Parameters")
            println(io)
            println(io, "| File | Reason | Iterations | PSMs |")
            println(io, "|------|--------|------------|------|")
            
            for (idx, status) in diagnostics.file_statuses
                if status.used_fallback
                    println(io, "| $(status.file_name) | $(status.fallback_reason) | $(status.n_iterations) | $(status.final_psm_count) |")
                end
            end
            println(io)
        end
        
        # Mass calibration statistics
        successful_offsets = [s.final_mass_offset for s in values(diagnostics.file_statuses) if s.converged && !s.used_fallback]
        if !isempty(successful_offsets)
            println(io, "## Mass Calibration Statistics (Successful Files)")
            println(io, "- Median offset: $(round(median(successful_offsets), digits=2)) ppm")
            println(io, "- MAD: $(round(mad(successful_offsets, normalize=false), digits=2)) ppm")
            println(io, "- Range: [$(round(minimum(successful_offsets), digits=2)), $(round(maximum(successful_offsets), digits=2))] ppm")
        end
    end
end

"""
    has_parameter_tuning_warnings(diagnostics)

Check if any files used fallback parameters.
"""
function has_parameter_tuning_warnings(diagnostics::ParameterTuningDiagnostics)
    return diagnostics.n_fallback > 0
end