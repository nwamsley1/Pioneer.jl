function load_predict_sources()
    files_loaded = Set{String}()

    function get_julia_files(dir::String)
        julia_files = String[]
        for (root, _, files) in walkdir(dir)
            for file in files
                if endswith(file, ".jl")
                    push!(julia_files, joinpath(root, file))
                end
            end
        end
        return julia_files
    end

    function safe_include!(file_path::String)::Bool
        if file_path in files_loaded
            return false
        end

        if !isfile(file_path)
            @user_warn "File not found: $file_path"
            return false
        end

        try
            include(file_path)
            push!(files_loaded, file_path)
            return true
        catch e
            @user_warn "Failed to include $file_path: $e"
            return false
        end
    end

    function include_files!(file_dir::String, file_names::Vector{String})
        for fname in file_names
            safe_include!(joinpath(file_dir, fname))
        end
        return nothing
    end

    function safe_include_directory!(dir::String)
        for jfile in get_julia_files(dir)
            safe_include!(jfile)
        end
        return nothing
    end

    include_files!(
        joinpath(REPO_ROOT, "src", "structs"),
        [
            "MassSpecData.jl",
            "FilteredMassSpecData.jl",
            "ChromObject.jl",
            "ArrayDict.jl",
            "Counter.jl",
            "Ion.jl",
            "LibraryIon.jl",
            "LibraryProteins.jl",
            "LibraryFragmentIndex.jl",
            "IsotopeTraceType.jl",
            "MatchIon.jl",
            "SparseArray.jl",
            "FragBoundModel.jl",
            "RetentionTimeIndex.jl",
            "MassErrorModel.jl",
            "RetentionTimeConversionModel.jl",
            "protein_inference_types.jl",
        ],
    )

    include_files!(
        joinpath(REPO_ROOT, "src", "utils"),
        [
            "serialization.jl",
            "isotopes.jl",
            "writeArrow.jl",
            "safeFileOps.jl",
            "profile.jl",
            "pdfUtils.jl",
        ],
    )

    safe_include!(joinpath(REPO_ROOT, "src", "utils", "FileOperations", "FileOperations.jl"))
    safe_include_directory!(joinpath(REPO_ROOT, "src", "structs", "KoinaStructs"))

    root_path = joinpath(REPO_ROOT, "src", "Routines", "BuildSpecLib")

    safe_include!(joinpath(root_path, "structs", "mods.jl"))
    safe_include!(joinpath(root_path, "fasta", "fasta_parser.jl"))
    safe_include!(joinpath(root_path, "fasta", "fasta_digest.jl"))
    safe_include!(joinpath(root_path, "fasta", "fasta_utils.jl"))
    safe_include!(joinpath(root_path, "fasta", "fasta_protein_table.jl"))

    safe_include!(joinpath(root_path, "fragments", "get_frag_bounds.jl"))
    safe_include!(joinpath(root_path, "fragments", "fragment_parse.jl"))
    safe_include!(joinpath(root_path, "fragments", "fragment_index.jl"))
    safe_include!(joinpath(root_path, "fragments", "fragment_annotation.jl"))
    safe_include!(joinpath(root_path, "fragments", "fragment_predict.jl"))

    safe_include!(joinpath(root_path, "koina", "koina_api.jl"))
    safe_include!(joinpath(root_path, "koina", "koina_batch_prep.jl"))
    safe_include!(joinpath(root_path, "koina", "koina_batch_parse.jl"))

    safe_include!(joinpath(root_path, "utils", "io.jl"))
    safe_include!(joinpath(root_path, "utils", "estimate_collision_ev.jl"))
    safe_include!(joinpath(root_path, "utils", "math.jl"))
    safe_include!(joinpath(root_path, "utils", "get_mz.jl"))
    safe_include!(joinpath(root_path, "utils", "parse_isotope_mods.jl"))
    safe_include!(joinpath(root_path, "utils", "check_params.jl"))
    safe_include!(joinpath(root_path, "utils", "essential_mods.jl"))
    safe_include!(joinpath(root_path, "utils", "parse_mods.jl"))

    safe_include!(joinpath(root_path, "build", "build_poin_lib.jl"))

    safe_include!(joinpath(root_path, "chronologer", "pair_decoys.jl"))
    safe_include!(joinpath(root_path, "chronologer", "chronologer_prep.jl"))
    safe_include!(joinpath(root_path, "chronologer", "chronologer_predict.jl"))
    safe_include!(joinpath(root_path, "chronologer", "chronologer_parse.jl"))

    return nothing
end

load_predict_sources()
