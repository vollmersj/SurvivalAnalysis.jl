using SurvivalAnalysis
using Test
using Random
using Random: seed!
using RCall
using Distributions
using DataFrames
using StatsModels

const files = setdiff(readdir(), ["runtests.jl"])

for f in files
    @testset "$(f)" begin
        @test include(f)
    end
end
