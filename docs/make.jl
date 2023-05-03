using Titus
using Documenter
#using DemoCards
DocMeta.setdocmeta!(Titus, :DocTestSetup, :(using Titus); recursive=true)

#op_templates, op_theme = cardtheme("grid")
#operations, operations_cb = makedemos("operations", op_templates)

About = "Introduction" => "../README.md"

GettingStarted = "gettingstarted.md"

UserGuide = "User's guide" => [
        "interface.md"
    ]

DevGuide = "Developer's guide" => [
        "wrappers.md"
    ]

Examples = "Examples" => [
        "examples/test.md"
    ]

Index = "Index" => "index.md"

License = "License" => "license.md"

PAGES = [
    About,
    GettingStarted,
    UserGuide,
    DevGuide,
    Examples,
    Index,
    License
    ]

makedocs(;
    modules=[Titus],
    authors="Nathan-Wamsley",
    repo="https://github.com/nwamsley1/Titus.jl/blob/{commit}{path}#{line}",
    sitename="Titus.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://nwamsley1.github.io/Titus.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=
        PAGES
    ,
)

deploydocs(
        devbranch="main",
    repo="github.com/nwamsley1/Titus.jl"
    
)
