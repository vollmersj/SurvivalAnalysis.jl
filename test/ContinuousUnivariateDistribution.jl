@testset "Test PH" begin
    d = SurvivalAnalysis.ContinuousPHDistribution(Exponential(1/5), 0.1)
    @test params(d) == (d.distr, d.lp)
    @test partype(d) == Float64
    @test hazard(d, 5) == 5 * exp(0.1)
    @test ccdf(d, 5) == exp(-25)^exp(0.1)
    @test cdf(d, 5) == 1 - (exp(-25)^exp(0.1))
    @test pdf(d, 5) == (5 * exp(0.1)) * (exp(-25)^exp(0.1))
    @test cdf.(d, quantile.(d, [0.2, 0.5, 0.8])) ≈ [0.2, 0.5, 0.8]
    @test show(d) === nothing ## TODO - HOW TO TEST OUTPUT?
end

@testset "Test AFT" begin
    d = SurvivalAnalysis.ContinuousAFTDistribution(Exponential(1/5), 0.1)
    @test params(d) == (d.distr, d.lp)
    @test partype(d) == Float64
    @test hazard(d, 5) == 5 * exp(-0.1)
    @test ccdf(d, 5) == exp(-5*(5/exp(0.1)))
    @test cdf(d, 5) == 1 - exp(-5*(5/exp(0.1)))
    @test pdf(d, 5) == (5 * exp(-0.1)) * exp(-5*(5/exp(0.1)))
    @test cdf.(d, quantile.(d, [0.2, 0.5, 0.8])) ≈ [0.2, 0.5, 0.8]
    @test show(d) === nothing ## TODO - HOW TO TEST OUTPUT?
end
