using SurvivalAnalysis
using Test
using Random
using Random: seed!
using RCall
using Distributions
using DataFrames
using StatsModels

const files = filter(x -> occursin("test_", x), readdir())

for f in files
    @testset "$(titlecase(replace(f, r"test_|.jl" => "", "_" => " ")))" begin
        @test include(f)
    end
end
