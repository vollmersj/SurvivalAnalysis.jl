using SurvivalAnalysis
using Test
using Random
using RCall
using Distributions
using DataFrames

@testset "Parametric models" begin
    @test include("test_parametricmodels.jl")
end
