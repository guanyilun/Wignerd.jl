using Wignerd
using Documenter

DocMeta.setdocmeta!(Wignerd, :DocTestSetup, :(using Wignerd); recursive=true)

makedocs(;
    modules=[Wignerd],
    authors="Yilun Guan",
    repo="https://github.com/guanyilun/Wignerd.jl/blob/{commit}{path}#{line}",
    sitename="Wignerd.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://guanyilun.github.io/Wignerd.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/guanyilun/Wignerd.jl",
    devbranch="main",
)
