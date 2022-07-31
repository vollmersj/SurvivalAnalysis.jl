using Documenter, SurvivalAnalysis

makedocs(
    sitename = "SurvivalAnalysis.jl",
    modules  = [SurvivalAnalysis],
    pages=[
        "Home" => "index.md"
        "Examples" => "examples.md"
        "API" => "api.md"
        ]
)

deploydocs(
    repo = "github.com/RaphaelS1/SurvivalAnalysis.jl.git"
)
