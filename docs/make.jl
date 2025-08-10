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

using Pioneer 
using Documenter 

makedocs(
    modules = [Pioneer],
    authors = "N.T. Wamsley",
    sitename = "Pioneer.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold = 300_000,  # Increase threshold to 300KB
        size_threshold_warn = 150_000  # Warning threshold to 150KB
    ),
    pages = [
        "Home" => "index.md",
        "User Guide" => [
            "Installation Guide" => "user_guide/installation.md",
            "Quick Start Tutorial" => "user_guide/quickstart.md",
            "Parameter Configuration" => "user_guide/parameters.md",
        ]
    ],
    doctest = true,
    clean = true,
    checkdocs = :exports
)

deploydocs(
    repo = "github.com/nwamsley1/Pioneer.jl",
    devbranch = "develop",
    devurl    = "dev",          # shows up at /dev
    versions  = [
        "stable" => "v^",       # latest release tag (vX.Y.Z) becomes /stable
        "v#.#",                 # publish each minor series: /v0.1, /v0.2, ...
        "dev"    => "develop",  # explicitly map the dev label to develop
    ]
)