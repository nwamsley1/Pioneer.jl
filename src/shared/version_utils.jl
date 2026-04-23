"""
    get_pioneer_version()

Return the Pioneer version string used for user-facing logs and CLI output.
Lookup order:
1. `ENV["PIONEER_VERSION"]` (if set)
2. Bundled `VERSION` file near the executable/app root
3. `VERSION` or `Project.toml` near the active package root
4. Repository-root `VERSION` or `Project.toml`
"""
function get_pioneer_version()
    env_version = strip(get(ENV, "PIONEER_VERSION", ""))
    if !isempty(env_version)
        return env_version
    end

    version_candidates = String[]

    if !isempty(PROGRAM_FILE)
        exe = PROGRAM_FILE
        if !isabspath(exe)
            exe_full = Sys.which(exe)
            exe = exe_full === nothing ? exe : exe_full
        end
        if ispath(exe)
            exe_real = realpath(exe)
            exe_dir = dirname(exe_real)
            push!(version_candidates, joinpath(exe_dir, "VERSION"))
            push!(version_candidates, joinpath(exe_dir, "..", "VERSION"))
        end
    end

    package_root = try
        pkgdir(@__MODULE__)
    catch
        nothing
    end

    if package_root !== nothing
        push!(version_candidates, joinpath(package_root, "VERSION"))
        push!(version_candidates, joinpath(package_root, "Project.toml"))
    end

    repo_root = normpath(joinpath(@__DIR__, "..", ".."))
    push!(version_candidates, joinpath(repo_root, "VERSION"))
    push!(version_candidates, joinpath(repo_root, "Project.toml"))

    for path in unique(version_candidates)
        if !isfile(path)
            continue
        end

        try
            if endswith(path, "Project.toml")
                content = read(path, String)
                version_match = match(r"^version[ \t]*=[ \t]*\"([^\"]+)\""m, content)
                if version_match !== nothing
                    return strip(version_match[1])
                end
            else
                v = strip(read(path, String))
                if !isempty(v)
                    return first(split(v, '\n'))
                end
            end
        catch
            # Ignore malformed or unreadable candidates and continue fallback chain.
        end
    end

    return "unknown"
end
