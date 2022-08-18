using Documenter, SurvivalAnalysis


DocMeta.setdocmeta!(SurvivalAnalysis, :DocTestSetup, :(using SurvivalAnalysis); recursive=true)

makedocs(
    sitename = "SurvivalAnalysis.jl",
    modules  = [SurvivalAnalysis, Distributions, DataFrames, Random],
    pages=[
        "Home" => "index.md"
        "Examples" => "examples.md"
        "API" => "api.md"
        ]
)

deploydocs(
    repo = "github.com/RaphaelS1/SurvivalAnalysis.jl.git"
)
