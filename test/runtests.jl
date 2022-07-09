using Survival
using Test

@testset "Survival.jl" begin
    @test Surv(1, 1).times = 1
end
