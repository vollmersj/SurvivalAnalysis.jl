using Survival
using Test

@testset "Survival.jl" begin
    @test typeof(Surv(1, 1)) == Survival.intSurv
    @test typeof(Surv(1, 1, "right")) == Survival.rcSurv
    @test typeof(Surv(1, 1, "left")) == Survival.lcSurv
end
