function importScripts()
    package_root = dirname(dirname(dirname(@__DIR__)))
    #package_root = dirname(@__DIR__)
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
    
    root_path = joinpath(package_root, "src", "Routines", "BuildSpecLib")
    [include(fname) for fname in get_julia_files(joinpath(package_root, "src", "structs", "KoinaStructs"))]
    # FASTA processing
    include(joinpath(root_path, "structs", "mods.jl"))
    include(joinpath(root_path, "fasta", "fasta_parser.jl"))
    include(joinpath(root_path, "fasta", "fasta_digest.jl"))
    include(joinpath(root_path, "fasta", "fasta_utils.jl"))
    # Fragment handling
    include(joinpath(root_path, "fragments", "get_frag_bounds.jl"))
    include(joinpath(root_path, "fragments", "fragment_parse.jl"))
    include(joinpath(root_path, "fragments", "fragment_index.jl"))
    include(joinpath(root_path, "fragments", "fragment_annotation.jl"))
    include(joinpath(root_path, "fragments", "fragment_predict.jl"))
    # Koina integration
    include(joinpath(root_path, "koina", "koina_api.jl"))
    include(joinpath(root_path, "koina", "koina_batch_prep.jl"))
    include(joinpath(root_path, "koina", "koina_batch_parse.jl"))
    # Utilities
    include(joinpath(root_path, "utils", "io.jl"))
    include(joinpath(root_path, "utils", "estimate_collision_ev.jl"))
    include(joinpath(root_path, "utils", "math.jl"))
    include(joinpath(root_path, "utils", "get_mz.jl"))
    include(joinpath(root_path, "utils", "parse_isotope_mods.jl"))
    include(joinpath(root_path, "utils", "check_params.jl"))
    #structs 
    include(joinpath(root_path, "structs", "EmpiricalLibrary.jl"))
    include(joinpath(root_path, "utils", "parse_mods.jl"))
    # Library building
    include(joinpath(root_path, "build", "build_poin_lib.jl"))
    #Chronologer Methods 
    include(joinpath(root_path, "chronologer", "chronologer_prep.jl"))
    include(joinpath(root_path, "chronologer", "chronologer_predict.jl"))
    include(joinpath(root_path, "chronologer", "chronologer_parse.jl"))
    # Profiling
    include(joinpath(package_root, "src", "utils", "profile.jl"))
end