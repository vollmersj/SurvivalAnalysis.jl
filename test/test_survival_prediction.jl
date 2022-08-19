@testset "Continuous distr only" begin
    sp = SurvivalPrediction(distr = [Exponential()])
    @test sp isa SurvivalAnalysis.ContinuousSurvivalPrediction
    @test all(isnan.(sp.lp))
    @test all(isnan.(sp.crank))
    @test all(isnan.(sp.time))
    @test sp.distr == [Exponential()]
    @test show(sp) === nothing
end

@testset "Discrete distr only" begin
    @test_throws ArgumentError SurvivalPrediction(
        distr = [DiscreteNonParametric([1.0], [1.0]), DiscreteNonParametric([2.0], [1.0])])
    distr = [DiscreteNonParametric([1.0, 2.0], [0.5, 0.5]),
            DiscreteNonParametric([1.0, 2.0], [0.3, 0.7])]
    sp = SurvivalPrediction(distr = distr)
    @test sp isa SurvivalAnalysis.DiscreteSurvivalPrediction
    @test all(isnan.(sp.lp))
    @test all(isnan.(sp.crank))
    @test all(isnan.(sp.time))
    @test sp.distr == distr
    @test sp.survival_matrix.time == [1.0, 2.0] # TODO - MAKE TIME/TIMES CONSISTENT
    @test sp.survival_matrix.surv == [[0.5,0.7] [0.0,0.0]]
    @test show(sp) === nothing
end

@testset "lp only" begin
    sp = SurvivalPrediction(lp = [1.0])
    sp isa SurvivalAnalysis.DeterministicSurvivalPrediction
    @test sp.lp == [1]
    @test all(isnan.(sp.crank))
    @test all(isnan.(sp.time))
    @test show(sp) === nothing
end

@testset "crank only" begin
    sp = SurvivalPrediction(crank = [1.0])
    sp isa SurvivalAnalysis.DeterministicSurvivalPrediction
    @test sp.crank == [1]
    @test all(isnan.(sp.lp))
    @test all(isnan.(sp.time))
end

@testset "time only" begin
    sp = SurvivalPrediction(time = [1.0])
    sp isa SurvivalAnalysis.DeterministicSurvivalPrediction
    @test sp.time == [1]
    @test all(isnan.(sp.lp))
    @test all(isnan.(sp.crank))
end

@testset "fit_times and survival_matrix only" begin
    @test_throws ArgumentError SurvivalPrediction(fit_times = [1.0])
    @test_throws ArgumentError SurvivalPrediction(survival_matrix = ones((1, 5)))
    @test_throws ArgumentError SurvivalPrediction(fit_times = [1.0], survival_matrix = ones((1, 5)))
    M = randn((5, 2))
    sp = SurvivalPrediction(fit_times = [1.0, 5.0], survival_matrix = M)
    sp isa SurvivalAnalysis.DiscreteSurvivalPrediction
    @test all(isnan.(sp.time))
    @test all(isnan.(sp.lp))
    @test all(isnan.(sp.crank))
    @test sp.survival_matrix.time == [1.0, 5.0]
    @test sp.survival_matrix.surv == M ## FIXME - SURV/SURVIVAL CONSISTENCY
    @test unique(support.(sp.distr)) == [[0.0, 1.0, 5.0]]
    @test probs.(sp.distr) == map(x -> [0, abs.(diff(x))...], eachrow(hcat(ones(5), M)))
end

@testset "SurvivalPrediction errors as expected" begin
    @test_throws ArgumentError SurvivalPrediction(lp = [1], crank = [1, 2])
end

true
