# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

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
            "Search Results" => "user_guide/outputs.md",
        ],
        "Advanced" => [
            "Parameter Configuration" => "user_guide/parameters.md",
        ],
        "Function Reference" => "reference.md",
    ],
    doctest = true,
    clean = true,
    checkdocs = :none
)

deploydocs(
    repo = "github.com/nwamsley1/Pioneer.jl",
    devbranch = "develop",
    devurl    = "dev",          # shows up at /dev
    dirname   = "docs",
    versions  = [
        "stable" => "v^",       # latest release tag (vX.Y.Z) becomes /stable
        "v#.#",                 # publish each minor series: /v0.1, /v0.2, ...
        "dev"    => "develop",  # explicitly map the dev label to develop
    ]
)
