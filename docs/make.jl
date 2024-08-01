using Pioneer
using Documenter

DocMeta.setdocmeta!(Pioneer, :DocTestSetup, :(using Pioneer); recursive=true)

makedocs(;
    modules=[Pioneer],
    authors="Natahan Wamsley",
    sitename="Pioneer.jl",
    format=Documenter.HTML(;
        canonical="https://nwamsley1.github.io/Pioneer.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nwamsley1/Pioneer.jl",
    devbranch="main",
)
