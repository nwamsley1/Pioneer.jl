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

# src/Routines/DownloadSpecLib.jl

# Entry point for PackageCompiler
function main_DownloadSpecLib(argv=ARGS)::Cint

    settings = ArgParseSettings(; autofix_names = true)
    @add_arg_table! settings begin
        "library"
            help = "Name of the library directory to download (see --list)"
            arg_type = String
            required = false
        "--dest"
            help = "Local directory the library is written into"
            arg_type = String
        "--repo"
            help = "Hugging Face dataset repository"
            arg_type = String
            default = HF_DEFAULT_REPO
        "--list"
            help = "List available libraries instead of downloading"
            action = :store_true
        "--json"
            help = "With --list, emit JSON instead of a table"
            action = :store_true
        "--force"
            help = "Replace an existing library directory at the destination"
            action = :store_true
    end
    parsed_args = parse_args(argv, settings; as_symbols = true)

    try
        if parsed_args[:list]
            ListSpecLibs(parsed_args[:repo]; as_json = parsed_args[:json])
        else
            library = parsed_args[:library]
            library === nothing && throw(ArgumentError(
                "no library named. Run with --list to see what is available."))
            dest = parsed_args[:dest]
            dest === nothing && throw(ArgumentError(
                "--dest is required: name the directory the library is written into."))
            DownloadSpecLib(library, dest;
                            repo = parsed_args[:repo], force = parsed_args[:force])
        end
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

"""
    ListSpecLibs(repo=HF_DEFAULT_REPO; as_json=false)

Print the libraries available in `repo`.

`as_json=true` emits the machine-readable catalog the GUI parses; the default is
a table for terminal use. Returns the catalog either way.
"""
function ListSpecLibs(repo::AbstractString = HF_DEFAULT_REPO; as_json::Bool = false)
    libraries = list_libraries(repo)
    manifest = fetch_json(repo, HF_MANIFEST_FILE)
    catalog = build_catalog(libraries, manifest;
                            config_for = lib -> fetch_json(repo, lib.name * "/config.json"))

    if as_json
        # Deliberately `print`, not a logging macro: this is a data channel the
        # GUI parses, so nothing else may share stdout.
        print(catalog_to_json(catalog))
        return catalog
    end

    if isempty(catalog)
        @user_info "No libraries found in $(repo)."
        return catalog
    end

    for entry in catalog
        @user_print "$(entry.name)"
        entry.title == entry.name || @user_print "    $(entry.title)"
        isempty(entry.model) || @user_print "    model: $(entry.model)"
        @user_print "    size: $(format_bytes(entry.total_bytes)) in $(entry.n_files) files"
        isempty(entry.description) || @user_print "    $(entry.description)"
        isempty(entry.recommended_for) || @user_print "    use when: $(entry.recommended_for)"
        for (label, value) in entry.details
            @user_print "    $(label): $(value)"
        end
        @user_print ""
    end
    return catalog
end

"""
    DownloadSpecLib(library, dest; repo=HF_DEFAULT_REPO, force=false)

Download the spectral library named `library` from `repo` into the directory
`dest`, returning the path of the downloaded library.

The transfer stages into `<name>.partial` and is renamed only after every file
has been fetched and checked against the repository's `SHA256SUMS`, so an
interrupted download never leaves a directory that looks like a usable library.
Re-running starts over: any leftover staging directory is discarded first,
rather than resumed into.
"""
function DownloadSpecLib(library::AbstractString, dest::AbstractString;
                         repo::AbstractString = HF_DEFAULT_REPO,
                         force::Bool = false)
    libraries = list_libraries(repo)
    idx = findfirst(l -> l.name == library, libraries)
    if idx === nothing
        available = isempty(libraries) ? "none" : join((l.name for l in libraries), ", ")
        throw(ArgumentError("no library named $(library) in $(repo). Available: $(available)"))
    end
    entry = libraries[idx]

    @user_info "Downloading $(entry.name) ($(format_bytes(entry.total_bytes)), " *
               "$(length(entry.files)) files) from $(repo)"

    last_report = Ref(0)
    path = download_library(entry, dest, repo;
        force = force,
        on_file = (name, i, n) -> @user_info("[$(i)/$(n)] $(name)"),
        on_progress = function (done, total)
            # Report per whole percent: a byte-level callback would emit
            # millions of lines into the GUI's log drawer.
            percent = total > 0 ? floor(Int, 100 * done / total) : 0
            if percent > last_report[]
                last_report[] = percent
                @user_info "  $(percent)% ($(format_bytes(done)) of $(format_bytes(total)))"
            end
        end)

    @user_info "Library ready at $(path)"
    return path
end
