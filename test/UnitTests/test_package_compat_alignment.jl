using Test
using TOML

@testset "Package compat alignment" begin
    repo_root = normpath(joinpath(@__DIR__, "../.."))
    root_project_path = joinpath(repo_root, "Project.toml")
    root_project = TOML.parsefile(root_project_path)
    root_compat = get(root_project, "compat", Dict{String, Any}())

    project_paths = vcat(
        [root_project_path],
        sort(readdir(joinpath(repo_root, "packages"); join=true)) .|> x -> joinpath(x, "Project.toml"),
        sort(readdir(joinpath(repo_root, "apps"); join=true)) .|> x -> joinpath(x, "Project.toml"),
    )
    filter!(isfile, project_paths)

    for project_path in project_paths
        project = TOML.parsefile(project_path)
        compat = get(project, "compat", Dict{String, Any}())
        deps = get(project, "deps", Dict{String, Any}())

        @test haskey(compat, "julia")
        @test compat["julia"] == root_compat["julia"]

        for dep_name in keys(deps)
            if haskey(root_compat, dep_name)
                @test haskey(compat, dep_name)
                @test compat[dep_name] == root_compat[dep_name]
            end
        end
    end
end

@testset "Public Pioneer registration surface has no internal package dependencies" begin
    repo_root = normpath(joinpath(@__DIR__, "../.."))
    root_project = TOML.parsefile(joinpath(repo_root, "Project.toml"))
    root_deps = get(root_project, "deps", Dict{String, Any}())
    root_sources = get(root_project, "sources", Dict{String, Any}())

    internal_package_names = Set(
        basename(pkg_dir) for pkg_dir in readdir(joinpath(repo_root, "packages"); join=true)
        if isfile(joinpath(pkg_dir, "Project.toml"))
    )

    @test isempty(intersect(keys(root_deps), internal_package_names))
    @test isempty(intersect(keys(root_sources), internal_package_names))
end
