using MultiBisect
using Documenter

DocTestSetup = quote
    using MultiBisect
    using MultiBisect: forwardinds, expandzeros, edgetuple, edgebounds, edgedim, domainindex, revert
    import MultiBisect.Roots as Roots
end
DocMeta.setdocmeta!(MultiBisect, :DocTestSetup, DocTestSetup; recursive=true)

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
        "API Reference" => "api.md"
    ],
)

deploydocs(
    repo="github.com/jbshannon/MultiBisect.jl.git",
    push_preview=true,
)
