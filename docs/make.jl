using Titus
using Documenter

DocMeta.setdocmeta!(Titus, :DocTestSetup, :(using Titus); recursive=true)

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
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nwamsley1/Titus.jl",
    devbranch="main",
)
