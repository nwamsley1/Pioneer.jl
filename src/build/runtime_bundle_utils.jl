using Pkg
using TOML

const RUNTIME_BUNDLE_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function rewrite_local_sources(project_data::Dict{String, Any})
    sources = get(project_data, "sources", nothing)
    if !(sources isa Dict)
        return project_data
    end

    rewritten = Dict{String, Any}()
    for (name, spec) in sources
        if spec isa Dict && haskey(spec, "path")
            rewritten[name] = Dict("path" => joinpath("dev", name))
        else
            rewritten[name] = spec
        end
    end

    project_data["sources"] = rewritten
    return project_data
end

function vendor_local_packages!(bundle_share_dir::String; repo_root::String=RUNTIME_BUNDLE_REPO_ROOT)
    source_root = joinpath(repo_root, "packages")
    vendored_root = joinpath(bundle_share_dir, "dev")
    mkpath(vendored_root)

    for package_name in sort(readdir(source_root))
        package_src = joinpath(source_root, package_name)
        isdir(package_src) || continue
        package_dest = joinpath(vendored_root, package_name)
        isdir(package_dest) && rm(package_dest; force=true, recursive=true)
        cp(package_src, package_dest; force=true)
    end

    return vendored_root
end

function vendor_repo_source_tree!(bundle_share_dir::String; repo_root::String=RUNTIME_BUNDLE_REPO_ROOT)
    source_tree = joinpath(repo_root, "src")
    vendored_source_tree = joinpath(bundle_share_dir, "src")
    isdir(vendored_source_tree) && rm(vendored_source_tree; force=true, recursive=true)
    cp(source_tree, vendored_source_tree; force=true)
    return vendored_source_tree
end

function write_runtime_bundle_project!(app_project_dir::String, bundle_share_dir::String)
    project_path = joinpath(app_project_dir, "Project.toml")
    project_data = TOML.parsefile(project_path)
    rewrite_local_sources(project_data)
    delete!(project_data, "name")
    delete!(project_data, "uuid")
    delete!(project_data, "version")
    delete!(project_data, "authors")
    mkpath(bundle_share_dir)
    open(joinpath(bundle_share_dir, "Project.toml"), "w") do io
        TOML.print(io, project_data; sorted=true)
    end
    return joinpath(bundle_share_dir, "Project.toml")
end

function instantiate_runtime_bundle!(bundle_share_dir::String)
    path_sep = Sys.iswindows() ? ";" : ":"
    julia_bin = joinpath(Sys.BINDIR, Base.julia_exename())
    load_path = string(bundle_share_dir, path_sep, "@stdlib")
    cmd = addenv(
        Cmd([
            julia_bin,
            "--startup-file=no",
            "--project=$bundle_share_dir",
            "-e",
            "using Pkg; Pkg.instantiate(); Pkg.precompile()",
        ]),
        "JULIA_LOAD_PATH" => load_path,
        "JULIA_DEPOT_PATH" => bundle_share_dir,
    )
    run(cmd)
    return bundle_share_dir
end

function bundle_runtime_project!(app_project_dir::String, bundled_app_dir::String; repo_root::String=RUNTIME_BUNDLE_REPO_ROOT)
    bundle_share_dir = joinpath(bundled_app_dir, "share", "julia")
    vendor_local_packages!(bundle_share_dir; repo_root)
    vendor_repo_source_tree!(bundle_share_dir; repo_root)
    write_runtime_bundle_project!(app_project_dir, bundle_share_dir)
    return bundle_share_dir
end

function bundle_and_instantiate_runtime_project!(app_project_dir::String, bundled_app_dir::String; repo_root::String=RUNTIME_BUNDLE_REPO_ROOT)
    bundle_share_dir = bundle_runtime_project!(app_project_dir, bundled_app_dir; repo_root)
    instantiate_runtime_bundle!(bundle_share_dir)
    return bundle_share_dir
end
