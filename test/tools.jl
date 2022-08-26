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

@testset "pstring" begin
    @test SurvivalAnalysis.pstring("G", 0) === "1"
    @test SurvivalAnalysis.pstring("G", 1) === "G"
    @test SurvivalAnalysis.pstring("G", -1) === "/G"
    @test SurvivalAnalysis.pstring("G", -1234567890) === "/G¹²³⁴⁵⁶⁷⁸⁹⁰"
end

@testset "c_pstring" begin
    @test SurvivalAnalysis.c_pstring("1", "1") === "1"
    @test SurvivalAnalysis.c_pstring("A", "/B") === "A/B"
    @test SurvivalAnalysis.c_pstring("/A", "B") === "B/A"
    @test SurvivalAnalysis.c_pstring("/A", "/B") === "1/A/B"
end

true
