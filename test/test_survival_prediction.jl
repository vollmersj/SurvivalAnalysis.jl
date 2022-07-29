@testset "Continuous ζ only" begin
    sp = SurvivalPrediction(ζ = [Exponential()])
    @test sp isa SurvivalAnalysis.ContinuousSurvivalPrediction
    @test all(isnan.(sp.lp))
    @test all(isnan.(sp.crank))
    @test all(isnan.(sp.time))
    @test sp.distr == [Exponential()]
    @test show(sp) === nothing
end

@testset "Discrete ζ only" begin
    @test_throws ArgumentError SurvivalPrediction(
        ζ = [DiscreteNonParametric([1.0], [1.0]), DiscreteNonParametric([2.0], [1.0])])
    ζ = [DiscreteNonParametric([1.0, 2.0], [0.5, 0.5]),
            DiscreteNonParametric([1.0, 2.0], [0.3, 0.7])]
    sp = SurvivalPrediction(ζ = ζ)
    @test sp isa SurvivalAnalysis.DiscreteSurvivalPrediction
    @test all(isnan.(sp.lp))
    @test all(isnan.(sp.crank))
    @test all(isnan.(sp.time))
    @test sp.distr == ζ
    @test sp.survival_matrix.time == [1.0, 2.0] # TODO - MAKE TIME/TIMES CONSISTENT
    @test sp.survival_matrix.surv == [[0.5,0.7] [0.0,0.0]]
    @test show(sp) === nothing
end

@testset "η only" begin
    sp = SurvivalPrediction(η = [1.0])
    sp isa SurvivalAnalysis.DeterministicSurvivalPrediction
    @test sp.lp == [1]
    @test all(isnan.(sp.crank))
    @test all(isnan.(sp.time))
    @test show(sp) === nothing
end

@testset "ϕ only" begin
    sp = SurvivalPrediction(ϕ = [1.0])
    sp isa SurvivalAnalysis.DeterministicSurvivalPrediction
    @test sp.crank == [1]
    @test all(isnan.(sp.lp))
    @test all(isnan.(sp.time))
end

@testset "T̂ only" begin
    sp = SurvivalPrediction(T̂ = [1.0])
    sp isa SurvivalAnalysis.DeterministicSurvivalPrediction
    @test sp.time == [1]
    @test all(isnan.(sp.lp))
    @test all(isnan.(sp.crank))
end

@testset "Ts and Ŝ only" begin
    @test_throws ArgumentError SurvivalPrediction(Ts = [1.0])
    @test_throws ArgumentError SurvivalPrediction(Ŝ = ones((1, 5)))
    @test_throws ArgumentError SurvivalPrediction(Ts = [1.0], Ŝ = ones((1, 5)))
    M = randn((5, 2))
    sp = SurvivalPrediction(Ts = [1.0, 5.0], Ŝ = M)
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
    @test_throws ArgumentError SurvivalPrediction(η = [1], ϕ = [1, 2])
end

true
