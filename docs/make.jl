using ParameterEstimation
using Documenter

DocMeta.setdocmeta!(ParameterEstimation, :DocTestSetup, :(using ParameterEstimation); recursive=true)

makedocs(;
    modules=[ParameterEstimation],
    authors="Daniel Pals <Daniel.Pals@tum.de>",
    repo="https://github.com/DanielJonathanPals/ParameterEstimation.jl/blob/{commit}{path}#{line}",
    sitename="ParameterEstimation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
