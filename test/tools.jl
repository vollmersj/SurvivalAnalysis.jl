@testset "Test operators" begin
    @test SurvivalAnalysis.∑(1, 2, 3, 4) == 10
    @test SurvivalAnalysis.∑([1, 2, 3, 4]) == 10
    @test SurvivalAnalysis.∏(1, 2, 3, 4) == 24
    @test SurvivalAnalysis.∏([1, 2, 3, 4]) == 24
end

@testset "Survival distribution methods" begin
    d = Weibull()
    @test fₜ(d, π) == pₜ(d, π) == pdf(d, π)
    @test Fₜ(d, π) == cdf(d, π)
    @test survival(d, π) == Sₜ(d, π) == ccdf(d, π)
    @test hazard(d, π) == hₜ(d, π) == pdf(d, π) / ccdf(d, π)
    @test cum_hazard(d, π) == Hₜ(d, π) == -log(ccdf(d, π))
end

true
