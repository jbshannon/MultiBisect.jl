using MultiBisect
using Documenter

DocMeta.setdocmeta!(MultiBisect, :DocTestSetup, :(using MultiBisect); recursive=true)

makedocs(;
    modules=[MultiBisect],
    authors="Jack Shannon <j.b.shannon@lse.ac.uk> and contributors",
    repo="https://github.com/jbshannon/MultiBisect.jl/blob/{commit}{path}#{line}",
    sitename="MultiBisect.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
