using Documenter, SurvivalAnalysis, DataFrames, Distributions


DocMeta.setdocmeta!(SurvivalAnalysis, :DocTestSetup, :(using SurvivalAnalysis); recursive=true)

makedocs(
    sitename = "SurvivalAnalysis",
    modules  = [SurvivalAnalysis],
    pages=[
        "Home" => "index.md"
        "Examples" => "examples.md"
        "API" => "api.md"
    ],
    doctest = true,
    strict = :doctest,
)

deploydocs(
    repo = "github.com/RaphaelS1/SurvivalAnalysis.jl.git"
)
