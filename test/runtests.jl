using Survival
using Test
using Random
using RCall
using Distributions

@testset "Survival.jl" begin
    @test typeof(Surv(1, 1)) == Survival.intSurv
    @test typeof(Surv(1, 1, "right")) == Survival.rcSurv
    @test typeof(Surv(1, 1, "left")) == Survival.lcSurv
end

@testset "Parametric models" begin
    @test include("test_parametricmodels.jl")
end
