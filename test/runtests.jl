using SurvivalAnalysis
using Test
using Random
using Random: seed!
using RCall
using Distributions
using DataFrames
using StatsModels

@testset "Nonparametric estimators" begin
    @test include("test_survivalestimators.jl")
end

@testset "Continuous parametric distributions" begin
    @test include("test_continuousparametric.jl")
end

@testset "Parametric models" begin
    @test include("test_parametricmodels.jl")
end

@testset "Utils" begin
    @test include("test_utils.jl")
end

@testset "Surv" begin
    @test include("test_surv.jl")
end
