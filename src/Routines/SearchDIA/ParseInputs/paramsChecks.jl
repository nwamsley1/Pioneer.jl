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

using JSON
#see src/Routines/BuildspecKib/utils
#struct InvalidParametersError <: Exception
#    message::String
#    params::Dict{String, Any}
#end

# Parameter default functions are loaded from paramDefaults.jl

function checkParams(json_path::String)
    # Read user params
    user_params = JSON.parsefile(json_path, dicttype=Dict{String,Any})
    
    # Apply defaults before validation
    defaults = get_default_parameters()
    params = merge_with_defaults(user_params, defaults)
    # Helper function to check if a key exists and has the correct type
    function check_param(dict, key, expected_type)
        if !haskey(dict, key)
            throw(InvalidParametersError("Missing parameter: $key", dict))
        elseif !(dict[key] isa expected_type)
            throw(InvalidParametersError("Invalid type for parameter $key: expected $(expected_type), got $(typeof(dict[key]))", dict))
        end
    end

    # Check top-level sections
    required_sections = [
        "global", "search",
        "acquisition", "optimization", "proteinScoring", "maxLFQ", "output", "paths"
    ]

    for section in required_sections
        check_param(params, section, Dict)
    end

    # Validate global parameters. Old configs had scoring.q_value_threshold
    # nested; the flat global.q_value_threshold form is now canonical, with a
    # backwards-compat fallback in the consumer (_resolve_q_value_threshold).
    global_params = params["global"]
    if !haskey(global_params, "q_value_threshold") &&
       !(haskey(global_params, "scoring") &&
         haskey(global_params["scoring"], "q_value_threshold"))
        throw(InvalidParametersError(
            "Missing parameter: global.q_value_threshold", global_params))
    end

    # Validate search parameters
    search_section = params["search"]
    # n_isotopes lives at search.n_isotopes (flattened from search.fragment_settings.n_isotopes).
    # Old configs nested under fragment_settings still parse — see backwards-compat
    # fallback in MainSearchParameters / IntegrateChromatogramSearchParameters.
    if haskey(search_section, "n_isotopes")
        check_param(search_section, "n_isotopes", Integer)
    elseif haskey(search_section, "fragment_settings") &&
           haskey(search_section["fragment_settings"], "n_isotopes")
        check_param(search_section["fragment_settings"], "n_isotopes", Integer)
    else
        throw(InvalidParametersError(
            "Missing parameter: search.n_isotopes", search_section))
    end

    # Validate acquisition parameters
    acq_params = params["acquisition"]
    check_param(acq_params, "nce", Integer)

    # Validate optimization parameters
    opt_params = params["optimization"]
    check_param(opt_params, "machine_learning", Dict)

    ml_params = opt_params["machine_learning"]
    check_param(ml_params, "max_psm_memory_mb", Real)
    check_param(ml_params, "pep_bin_size", Integer)

    check_param(opt_params, "chromatogram_integration", Dict)
    chrom_params = opt_params["chromatogram_integration"]
    check_param(chrom_params, "trace_mode", String)
    if !(chrom_params["trace_mode"] in ("combined", "separate"))
        throw(InvalidParametersError(
            "Invalid parameter optimization.chromatogram_integration.trace_mode: expected \"combined\" or \"separate\"",
            chrom_params,
        ))
    end
    check_param(chrom_params, "deconvolution_solver", String)
    if !(chrom_params["deconvolution_solver"] in ("huber", "pmm"))
        throw(InvalidParametersError(
            "Invalid parameter optimization.chromatogram_integration.deconvolution_solver: expected \"huber\" or \"pmm\"",
            chrom_params,
        ))
    end

    # Validate protein scoring parameters
    protein_scoring = params["proteinScoring"]
    check_param(protein_scoring, "min_peptides", Integer)
    check_param(protein_scoring, "global_protein_inference", Bool)
    check_param(protein_scoring, "write_qc_plots", Bool)

    # Validate MaxLFQ parameters
    output = params["maxLFQ"]
    check_param(output, "run_to_run_normalization", Bool)

    # Validate output parameters
    output = params["output"]
    check_param(output, "write_csv", Bool)
    check_param(output, "delete_temp", Bool)
    check_param(output, "write_decoys", Bool)

    # Validate path parameters
    paths = params["paths"]
    check_param(paths, "library", String)
    check_param(paths, "ms_data", String)
    check_param(paths, "results", String)

    return params
end
