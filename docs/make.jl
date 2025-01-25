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
    devbranch = "main"
)