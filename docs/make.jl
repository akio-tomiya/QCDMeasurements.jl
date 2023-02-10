using QCDMeasurements
using Documenter

DocMeta.setdocmeta!(QCDMeasurements, :DocTestSetup, :(using QCDMeasurements); recursive=true)

makedocs(;
    modules=[QCDMeasurements],
    authors="Yuki Nagai <cometscome@gmail.com> and contributors",
    repo="https://github.com/cometscome/QCDMeasurements.jl/blob/{commit}{path}#{line}",
    sitename="QCDMeasurements.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cometscome.github.io/QCDMeasurements.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cometscome/QCDMeasurements.jl",
    devbranch="main",
)
