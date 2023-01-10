using CrystallographyRecipes
using Documenter

DocMeta.setdocmeta!(CrystallographyRecipes, :DocTestSetup, :(using CrystallographyRecipes); recursive=true)

makedocs(;
    modules=[CrystallographyRecipes],
    authors="singularitti <singularitti@outlook.com> and contributors",
    repo="https://github.com/MineralsCloud/CrystallographyRecipes.jl/blob/{commit}{path}#{line}",
    sitename="CrystallographyRecipes.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/CrystallographyRecipes.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/CrystallographyRecipes.jl",
    devbranch="main",
)
