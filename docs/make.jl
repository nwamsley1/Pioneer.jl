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
        "Getting Started" => [
            "Installation" => "getting_started/installation.md",
            "Quick Start" => "getting_started/quickstart.md",
            "Basic Usage" => "getting_started/basic_usage.md"
        ],
        "User Guide" => [
            "Overview" => "user_guide/overview.md",
            "Parameter Configuration" => "user_guide/parameters.md",
            "Search Modes" => "user_guide/search_modes.md",
            "Output Files" => "user_guide/outputs.md"
        ],
        "API Reference" => [
            "Core Functions" => "api/core.md",
            "Search Types" => "api/search_types.md",
            "Parameter Types" => "api/parameter_types.md",
            "Utilities" => "api/utilities.md"
        ],
        "Advanced Topics" => [
            "Algorithm Details" => "advanced/algorithms.md",
            "Performance Tuning" => "advanced/performance.md",
            "Custom Extensions" => "advanced/extensions.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/nwamsley1/Pioneer.jl",
    devbranch = "main"
)