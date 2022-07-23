using SurvivalAnalysis
using Test
using Random
using Random: seed!
using RCall
using Distributions
using DataFrames
using StatsModels

@testset "Parametric models" begin
    @test include("test_parametricmodels.jl")
end
