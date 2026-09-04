using Documenter
import GridapTrilinos

DocMeta.setdocmeta!(
    GridapTrilinos,
    :DocTestSetup,
    :(import GridapTrilinos);
    recursive=true,
)

makedocs(;
    modules=[GridapTrilinos],
    sitename="GridapTrilinos.jl",
    format=Documenter.HTML(; edit_link="main"),
    pages=[
        "Home" => "index.md",
    ],
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/shreyas02/GridapTrilinos.jl.git",
    devbranch="main",
)
