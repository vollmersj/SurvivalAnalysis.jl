using SurvivalAnalysis
using Test
using Random
using Random: seed!
using StableRNGs
using FiniteDifferences
using RCall
using Distributions
using DataFrames
using StatsModels
using Optim

const files = setdiff(readdir(), ["runtests.jl", "Project.toml", "Manifest.toml"])

for f in files
    @testset "$(f)" begin
        include(f)
    end
end
