using Documenter, SurvivalAnalysis, DataFrames, Distributions


DocMeta.setdocmeta!(SurvivalAnalysis, :DocTestSetup, :(using SurvivalAnalysis, Distributions); recursive=true)

makedocs(
    sitename = "SurvivalAnalysis",
    modules  = [SurvivalAnalysis],
    pages=[
        "Home" => "index.md"
        "Tutorials" => "tutorials.md"
        "How-to guides" => "howto.md"
        "Explanations" => "explanations.md"
        "Reference (API)" => "api.md"
        "Changelog" => "news.md"
    ],
    doctest = true,
    strict = :doctest,
)

deploydocs(
    repo = "github.com/RaphaelS1/SurvivalAnalysis.jl.git"
)
