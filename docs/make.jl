using Documenter, SurvivalAnalysis

makedocs(
    sitename = "SurvivalAnalysis.jl",
    modules  = [SurvivalAnalysis],
    pages=[
        "Home" => "index.md"
        ]
)

deploydocs(
    repo = "github.com/RaphaelS1/SurvivalAnalysis.jl.git"
)
